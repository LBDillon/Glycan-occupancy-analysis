"""Context features v2 — the environment of a whole sequon, not one residue.

`features.py` describes a single asparagine, and exists to *match* sites. This
describes the three-residue sequon and exists to *characterise* it, so that a
distribution of natural occupied environments can be built and an arbitrary site
located within it.

It is a separate module rather than an extension of `features.py` because that
file feeds the frozen matching. Changing its outputs would move the pairs every
comparison rests on. This shares its parsing and SASA cache and adds nothing to
it.

## What is measured, and why

**At all three sequon positions, not just the Asn.** Oligosaccharyltransferase
recognises the motif through a catalytic geometry that constrains the +2 hydroxyl
as much as the asparagine, so describing only the Asn discards the position most
likely to carry signal.

**Secondary structure from DSSP, RSA from Shrake-Rupley.** Deliberately mixed
sources. DSSP's own RSA runs about 1.22x ours (different probe, radii and max-ASA
reference; Spearman 0.988), and the frozen matching used the Shrake-Rupley scale
— so taking RSA from DSSP would put this atlas on a different scale from the
benchmark it is meant to explain. SS has no such conflict.

**DSSP on a single-chain extract.** Whole-assembly DSSP fails on large
structures — 38 of 50 failures in the survey were `Empty file`, and failing
structures were twice the median size. Extracting the chain rescues about a
third of them and, verified on 39 residues where both work, changes no SS call.
Coverage is still not equal between arms (~95% occupied, ~89% control), so
`dssp_ok` and `dssp_reason` are recorded per site and any analysis using SS must
report coverage per arm rather than dropping quietly.

**Backbone dihedrals** as a continuous read on local conformation, coarser than
SS but available whenever the neighbouring residues are.

**Side-chain contacts, not just CA neighbours.** A CA within 8 A says little
about whether a glycan could be accommodated; heavy-atom contacts around the Asn
side chain say more.
"""
from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

from .features import MAX_ASA, _clean_icode, _model_with_sasa, _rsa_bin

MMCIF_SUFFIXES = {".cif", ".mmcif"}

# Chemistry classes, matching the ortholog module's vocabulary so the two
# analyses describe neighbourhoods in the same words.
AA_CLASS = {**{a: "acidic" for a in "DE"}, **{a: "basic" for a in "KRH"},
            **{a: "polar" for a in "STNQ"}, **{a: "hydrophobic" for a in "AVLIM"},
            **{a: "aromatic" for a in "FWY"},
            "G": "glycine", "P": "proline", "C": "cysteine"}
CLASS_ORDER = ("acidic", "basic", "polar", "hydrophobic", "aromatic",
               "glycine", "proline", "cysteine", "other")
CHARGE = {"D": -1, "E": -1, "K": +1, "R": +1}

# DSSP's eight states, collapsed. Both are kept: the coarse one for stratifying,
# the raw code because collapsing loses the G/I distinction.
SS_COARSE = {"H": "helix", "G": "helix", "I": "helix",
             "E": "sheet", "B": "sheet",
             "T": "loop", "S": "loop", "-": "loop", " ": "loop"}

NEIGHBOUR_RADIUS = 8.0        # CA-CA, matching features.py
CONTACT_RADIUS = 5.0          # heavy atom, around the Asn side chain
AROMATIC = set("FWY")

_DSSP_CACHE: dict = {}
_DSSP_LIMIT = 32


class _OneChain(Select):
    """Protein residues of one chain: what DSSP is asked to look at."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " "


def _parse_model(path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = (MMCIFParser(QUIET=True)
                  if Path(path).suffix.lower() in MMCIF_SUFFIXES else PDBParser(QUIET=True))
        return next(iter(parser.get_structure("s", str(path))), None)


def dssp_for_chain(path: Path, chain_id: str) -> tuple[dict, str]:
    """`{resseq: (ss_code, dssp_rsa)}` for one chain, and a failure reason.

    Runs on a single-chain extract rather than the whole assembly: whole-assembly
    DSSP fails on large structures, and the extract produces identical SS calls
    where both work. Only SS is taken from it; RSA is returned for diagnostics
    and should not be mixed with the Shrake-Rupley values.
    """
    key = (str(path), str(chain_id))
    if key in _DSSP_CACHE:
        return _DSSP_CACHE[key]

    from tempfile import NamedTemporaryFile

    from Bio.PDB.DSSP import DSSP

    result, reason = {}, ""
    try:
        model = _parse_model(path)
        if model is None:
            raise ValueError("structure unreadable")
        io = PDBIO()
        io.set_structure(model)
        with NamedTemporaryFile(suffix=".pdb", delete=True) as handle:
            io.save(handle.name, _OneChain(str(chain_id)))
            extracted = _parse_model(Path(handle.name))
            if extracted is None:
                raise ValueError("chain extract unreadable")
            dssp = DSSP(extracted, handle.name, dssp="mkdssp")
            for (chain, residue_id), entry in dssp.property_dict.items():
                if chain == str(chain_id) and residue_id[0] == " ":
                    result[int(residue_id[1])] = (entry[2], entry[3])
    except Exception as exc:
        reason = f"{type(exc).__name__}: {str(exc)[:60]}"

    _DSSP_CACHE[key] = (result, reason)
    while len(_DSSP_CACHE) > _DSSP_LIMIT:
        _DSSP_CACHE.pop(next(iter(_DSSP_CACHE)))
    return result, reason


def _dihedral(p0, p1, p2, p3) -> float | None:
    """Torsion angle in degrees, or None if any atom is absent.

    Note the sign of `b0`: it is p0 - p1, NOT p1 - p0. Getting that backwards
    negates every angle, which is easy to miss because the values stay in range
    and look plausible — phi simply comes out positive, where L-amino acids are
    almost always negative. Verified against Bio.PDB's `get_phi_psi_list`.
    """
    if any(p is None for p in (p0, p1, p2, p3)):
        return None
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    norm = np.linalg.norm(b1)
    if norm == 0:
        return None
    b1 = b1 / norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return round(math.degrees(math.atan2(y, x)), 1)


def _atom(residue, name):
    return residue[name].coord if residue is not None and name in residue else None


def _phi_psi(previous, residue, following):
    phi = _dihedral(_atom(previous, "C"), _atom(residue, "N"),
                    _atom(residue, "CA"), _atom(residue, "C"))
    psi = _dihedral(_atom(residue, "N"), _atom(residue, "CA"),
                    _atom(residue, "C"), _atom(following, "N"))
    return phi, psi


def _one_letter(residue) -> str:
    try:
        return seq1(residue.get_resname(), undef_code="X")
    except Exception:
        return "X"


def _residue_at(chain, resseq: int, icode: str = ""):
    icode = _clean_icode(icode)
    for residue in chain:
        if not is_aa(residue, standard=False):
            continue
        _, seq, code = residue.id
        if int(seq) == int(resseq) and str(code).strip() == icode:
            return residue
    return None


def _position_features(residue, dssp: dict, prefix: str) -> dict:
    """RSA, secondary structure and dihedrals at one sequon position."""
    out = {f"{prefix}_residue": None, f"{prefix}_rsa": None, f"{prefix}_rsa_bin": None,
           f"{prefix}_ss": None, f"{prefix}_ss_coarse": None,
           f"{prefix}_phi": None, f"{prefix}_psi": None}
    if residue is None:
        return out
    code = _one_letter(residue)
    sasa = float(getattr(residue, "sasa", 0.0) or 0.0)
    max_asa = MAX_ASA.get(code)
    rsa = min(sasa / max_asa, 1.5) if max_asa else None
    entry = dssp.get(int(residue.id[1]))
    out[f"{prefix}_residue"] = code
    out[f"{prefix}_rsa"] = round(rsa, 4) if rsa is not None else None
    out[f"{prefix}_rsa_bin"] = _rsa_bin(rsa)
    if entry:
        out[f"{prefix}_ss"] = entry[0]
        out[f"{prefix}_ss_coarse"] = SS_COARSE.get(entry[0], "unknown")
    return out


def sequon_context(path: Path, chain_id: str, resseq: int, icode: str = "",
                   pdb_id: "str | None" = None) -> dict | None:
    """Everything measurable about one sequon's environment, or None if absent.

    `resseq` is the asparagine's author residue number; +1 and +2 are located by
    walking the chain's resolved residues, not by adding 1 and 2, because
    insertion codes and numbering gaps make arithmetic wrong.
    """
    model = _model_with_sasa(Path(path))
    if model is None:
        return None
    try:
        chain = model[str(chain_id)]
    except KeyError:
        return None

    resolved = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    if not resolved:
        return None
    order = {id(r): i for i, r in enumerate(resolved)}

    asn = _residue_at(chain, resseq, icode)
    if asn is None or id(asn) not in order:
        return None
    i = order[id(asn)]
    plus1 = resolved[i + 1] if i + 1 < len(resolved) else None
    plus2 = resolved[i + 2] if i + 2 < len(resolved) else None

    dssp, dssp_reason = dssp_for_chain(Path(path), chain_id)
    out: dict = {"dssp_ok": bool(dssp), "dssp_reason": dssp_reason}

    for prefix, residue in (("n", asn), ("plus1", plus1), ("plus2", plus2)):
        out.update(_position_features(residue, dssp, prefix))
        if residue is not None:
            j = order[id(residue)]
            phi, psi = _phi_psi(resolved[j - 1] if j > 0 else None, residue,
                                resolved[j + 1] if j + 1 < len(resolved) else None)
            out[f"{prefix}_phi"], out[f"{prefix}_psi"] = phi, psi

    out["triplet_observed"] = "".join(
        out.get(f"{p}_residue") or "?" for p in ("n", "plus1", "plus2"))
    out["triplet_resolved"] = all(out.get(f"{p}_residue") for p in ("n", "plus1", "plus2"))

    # --- neighbourhood, by CA within 8 A ------------------------------------
    origin = asn["CA"].coord
    neighbours = []
    for other in model:
        for residue in other:
            if residue is asn or not is_aa(residue, standard=False) or "CA" not in residue:
                continue
            if float(np.linalg.norm(residue["CA"].coord - origin)) <= NEIGHBOUR_RADIUS:
                neighbours.append((_one_letter(residue), residue))
    codes = [c for c, _ in neighbours]
    n = len(codes)
    out["n_neighbours_8a"] = n
    counts = {k: 0 for k in CLASS_ORDER}
    for c in codes:
        counts[AA_CLASS.get(c, "other")] += 1
    for k in CLASS_ORDER:
        out[f"neighbour_{k}_count_8a"] = counts[k]
        out[f"neighbour_{k}_fraction_8a"] = round(counts[k] / n, 4) if n else None
    out["neighbour_net_charge_8a"] = sum(CHARGE.get(c, 0) for c in codes)

    # --- heavy-atom contacts around the Asn side chain ----------------------
    # A CA within 8 A says little about whether a glycan fits; side-chain
    # contacts say more.
    side_chain = [a for a in asn if a.element != "H"
                  and a.get_id() not in ("N", "CA", "C", "O")]
    contacts = 0
    if side_chain:
        for _, residue in neighbours:
            for atom in residue:
                if atom.element == "H":
                    continue
                if any(float(np.linalg.norm(atom.coord - s.coord)) <= CONTACT_RADIUS
                       for s in side_chain):
                    contacts += 1
                    break
    out["sidechain_contacts_5a"] = contacts if side_chain else None
    out["has_sidechain_atoms"] = bool(side_chain)

    # --- nearest aromatic and nearest disulfide -----------------------------
    def nearest(predicate):
        best = None
        for other in model:
            for residue in other:
                if residue is asn or not is_aa(residue, standard=False) or "CA" not in residue:
                    continue
                if not predicate(residue):
                    continue
                d = float(np.linalg.norm(residue["CA"].coord - origin))
                best = d if best is None else min(best, d)
        return round(best, 2) if best is not None else None

    out["nearest_aromatic_ca"] = nearest(lambda r: _one_letter(r) in AROMATIC)
    out["aromatic_within_8a"] = counts["aromatic"] > 0

    def is_cystine(residue):
        if _one_letter(residue) != "C" or "SG" not in residue:
            return False
        for other in model:
            for partner in other:
                if partner is residue or _one_letter(partner) != "C" or "SG" not in partner:
                    continue
                if float(np.linalg.norm(partner["SG"].coord - residue["SG"].coord)) <= 2.5:
                    return True
        return False

    out["nearest_disulfide_ca"] = nearest(is_cystine)

    # --- position in the chain, and QC -------------------------------------
    numbers = [int(r.id[1]) for r in resolved]
    out["chain_length_resolved"] = len(resolved)
    out["distance_to_n_terminus_resolved"] = int(resseq) - min(numbers)
    out["distance_to_c_terminus_resolved"] = max(numbers) - int(resseq)
    bfactors = [float(a.get_bfactor()) for a in asn if a.get_bfactor() is not None]
    out["mean_bfactor_asn"] = round(float(np.mean(bfactors)), 2) if bfactors else None
    out["structure_pdb_id"] = str(pdb_id or Path(path).stem).upper()
    out["structure_chain_id"] = str(chain_id)
    out["structure_format"] = Path(path).suffix.lstrip(".")
    return out


FEATURE_COLUMNS = (
    "dssp_ok", "dssp_reason", "triplet_observed", "triplet_resolved",
    *[f"{p}_{f}" for p in ("n", "plus1", "plus2")
      for f in ("residue", "rsa", "rsa_bin", "ss", "ss_coarse", "phi", "psi")],
    "n_neighbours_8a",
    *[f"neighbour_{k}_{s}_8a" for k in CLASS_ORDER for s in ("count", "fraction")],
    "neighbour_net_charge_8a", "sidechain_contacts_5a", "has_sidechain_atoms",
    "nearest_aromatic_ca", "aromatic_within_8a", "nearest_disulfide_ca",
    "chain_length_resolved", "distance_to_n_terminus_resolved",
    "distance_to_c_terminus_resolved", "mean_bfactor_asn",
    "structure_pdb_id", "structure_chain_id", "structure_format",
)
