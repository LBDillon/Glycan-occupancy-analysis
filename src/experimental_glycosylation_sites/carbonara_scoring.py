"""CARBonAra scoring — a fourth opinion on the sequon, from a fourth parser.

CARBonAra (Krapp et al., Nat Commun 2024) is a geometric transformer over atomic
coordinates and element names, built on PeSTo. It has no amino-acid-specific
parameterisation, which is why it can read arbitrary molecular context — and why
this module has to be careful about what it is handed.

This integration is deliberately narrow. It asks only whether CARBonAra assigns a
different occupancy-associated statistical preference to experimentally occupied
sequons than to the frozen matched controls, on **exactly the input the other
structure models get**: one protein chain, no glycans, no ligands, no water. The
context-aware scoring CARBonAra is actually known for is not done here, and
nothing below should be read as evidence about glycosylation biology.

## What the conditional is

CARBonAra is not autoregressive and has no decoding order. Its residue-identity
input is an "imprint": a per-residue one-hot `yt`, zero wherever the identity is
withheld. `apply_model` sees element types plus `Mr @ yt` and nothing else — side
chains are stripped and an idealised C-beta is placed at every residue, including
glycine — so a position with `yt = 0` genuinely leaks no identity.

Scoring position *i* therefore means imprinting every other mappable canonical
residue natively and zeroing row *i*:

    P(residue at i | backbone, all other native residues)

which is one deterministic forward pass per position. Recorded as
`conditional_all_other_native` with `n_orders = 1`, so it can never be pooled
with ProteinMPNN's eight-order average or ESM-IF's causal prefix. As always, raw
magnitudes are not comparable between models; the SD-standardised matched-pair
contrast within each model is.

## The two things that will bite

**The alphabet is sorted by abundance, not alphabetically.** Upstream
`src/data_encoding.py`:

    std_aminoacids = np.array([
        'LEU', 'GLU', 'ARG', 'LYS', 'VAL', 'ILE', 'PHE', 'ASP', 'TYR',
        'ALA', 'THR', 'SER', 'GLN', 'ASN', 'PRO', 'GLY', 'HIS', 'TRP',
        'MET', 'CYS'])

Asparagine is column 13, not column 11. Assuming otherwise is the 2026-08-20
defect exactly — that one read P(Asp) as P(Asn) for a month. `verify_alphabet`
runs against the checkout at load time so a change upstream stops the run.

**MSE has no special case.** ProteinMPNN maps `HETATM MSE` to methionine.
CARBonAra's `clean_structure` does not: dropped as a heteroatom the residue
vanishes and every index after it shifts, kept as MSE it is not in
`std_aminoacids` and becomes a ligand. Preparation converts it here.

The residue enumeration itself is friendlier than ProteinMPNN's. `clean_structure`
renumbers observed residues consecutively from 1 in file order, so there is no
gap-filling and no placeholder token — a numbering gap shifts nothing. That makes
the mapping an identity, but it is still *verified* against what CARBonAra parsed
rather than assumed, because an unverified identity mapping is exactly what the
25.3% misindexing looked like before anyone checked.
"""
from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import numpy as np

# Shared with ProteinMPNN so the score means the same thing under every model:
# the same clamp before the logit, and the same tolerance on a distribution.
from .mpnn_scoring import EPSILON, PROBABILITY_SUM_TOLERANCE, logit  # noqa: F401

# Upstream `src/data_encoding.py`, quoted in the module docstring above. Order is
# load-bearing: it is the column order of every probability vector this module
# returns. Verified against the checkout by `verify_alphabet` before any scoring.
CARBONARA_RESNAMES = (
    "LEU", "GLU", "ARG", "LYS", "VAL", "ILE", "PHE", "ASP", "TYR", "ALA",
    "THR", "SER", "GLN", "ASN", "PRO", "GLY", "HIS", "TRP", "MET", "CYS",
)

# Upstream `src/structure.py: res3to1`, restricted to the twenty and applied in
# CARBonAra's order rather than a conventional one.
_THREE_TO_ONE = {
    "CYS": "C", "ASP": "D", "SER": "S", "GLN": "Q", "LYS": "K",
    "ILE": "I", "PRO": "P", "THR": "T", "PHE": "F", "ASN": "N",
    "GLY": "G", "HIS": "H", "LEU": "L", "ARG": "R", "TRP": "W",
    "ALA": "A", "VAL": "V", "GLU": "E", "TYR": "Y", "MET": "M",
}
RES3TO1 = {name: _THREE_TO_ONE[name] for name in CARBONARA_RESNAMES}
ALPHABET = "".join(RES3TO1[name] for name in CARBONARA_RESNAMES)
AA_INDEX = {letter: index for index, letter in enumerate(ALPHABET)}

# `add_virtual_cb` needs N, CA and C to place the C-beta; O completes the
# backbone encoding. A residue missing any of the four is an exclusion.
BACKBONE_ATOMS = ("N", "CA", "C", "O")

DEFAULT_MODEL = "s_v6_4_2022-09-16_11-51"

# Named for what it is, so it can never be pooled with the other models'.
CONDITIONING = "conditional_all_other_native"
N_ORDERS = 1
SCORE_SD = 0.0


class ChainUnreadableError(ValueError):
    """CARBonAra could not be given a chain whose indices we can justify."""


class IncompleteBackboneError(ValueError):
    """A requested residue has no usable CARBonAra counterpart."""


class SequonMismatchError(ValueError):
    """CARBonAra read different residues than the manifest recorded."""


class InvalidProbabilityVector(ValueError):
    """A row that is not a probability distribution over the twenty."""


class AlphabetMismatchError(RuntimeError):
    """The checkout's amino-acid order is not the one these columns assume."""


def verify_alphabet(std_aminoacids) -> None:
    """Confirm the checkout's `std_aminoacids` is the order this module assumes.

    Called at model load. A silent reordering upstream would not raise anywhere
    else — every vector would still be twenty finite numbers summing to one, and
    every score would be wrong.
    """
    observed = tuple(str(name) for name in std_aminoacids)
    if observed != CARBONARA_RESNAMES:
        raise AlphabetMismatchError(
            "CARBonAra's std_aminoacids is not the order this module was "
            f"written against.\n  checkout: {observed}\n  expected: "
            f"{CARBONARA_RESNAMES}\nScoring would read the wrong columns; "
            "update carbonara_scoring.CARBONARA_RESNAMES and its tests.")


# --------------------------------------------------------------------------
# Locating and importing the checkout.
# --------------------------------------------------------------------------

def _import_carbonara(carbonara_dir=None):
    """Import the upstream package, lazily and from its checkout.

    Upstream is a script layout rather than an installable package: `carbonara.py`
    sits at the checkout root and falls back to `from src...` when its relative
    import fails, so putting the checkout on `sys.path` is what makes it work.

    Import is deferred to call time because `src/structure_io.py` needs gemmi and
    `carbonara.py` needs blosum, and neither belongs in this package's core
    dependencies. Every non-CARBonAra command must keep working without them.
    """
    import sys

    if carbonara_dir is None:
        from .runner_support import carbonara_dir as locate

        carbonara_dir = locate()
    path = str(Path(carbonara_dir).expanduser().resolve())
    if path not in sys.path:
        sys.path.insert(0, path)
    import carbonara

    return carbonara


def load_model(carbonara_dir=None, model_name: str = DEFAULT_MODEL,
               device: str = "cpu"):
    """Load a checkpoint, refusing one whose alphabet is not the assumed one."""
    upstream = _import_carbonara(carbonara_dir)
    verify_alphabet(upstream.std_aminoacids)
    return upstream.CARBonAra(model_name=model_name, device_name=device)


def _load_structure(pdb_text: str, carbonara_dir=None):
    """Parse prepared PDB text through CARBonAra's own reader.

    Through *its* parser rather than ours, because the mapping has to be checked
    against what CARBonAra actually read. `rm_hetatm` and `rm_wat` are set even
    though preparation already removed both: the baseline arm must be
    protein-only whatever reaches it, and a belt-and-braces flag costs nothing.

    The temporary file exists because the upstream reader takes a path. It is
    removed before returning, so no converted structure outlives the chain.
    """
    upstream = _import_carbonara(carbonara_dir)
    handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    try:
        handle.write(pdb_text)
        handle.close()
        return upstream.load_structure(handle.name, rm_hetatm=True, rm_wat=True)
    finally:
        Path(handle.name).unlink(missing_ok=True)


# --------------------------------------------------------------------------
# Protein-only structure preparation.
# --------------------------------------------------------------------------

def _atom_line(serial: int, name: str, resname: str, chain: str, resseq: int,
               xyz, element: str) -> str:
    """One fixed-column ATOM record.

    Atom names are placed by the PDB convention — four-character names start in
    column 13, shorter ones in column 14 — so a reader that slices by column
    sees `CD1` and not `D1`.
    """
    padded = name if len(name) >= 4 else f" {name:<3s}"
    x, y, z = (float(v) for v in xyz)
    return (f"ATOM  {serial:5d} {padded:4s} {resname:>3s} {chain:1s}{resseq:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"          {element:>2s}")


def protein_only_pdb(structure_path, chain_id: str, pdb_id: "str | None" = None):
    """One protein chain as PDB text, plus the metadata the mapping needs.

    Returns `(pdb_text, residue_ids, native_sequence)`, where `residue_ids` are
    the author `(resseq, icode)` pairs in manifest order and `native_sequence` is
    the manifest's own sequence for the chain.

    The residue filter is `_parse_chains`' filter, and the result is checked
    against `_parse_chains` before being returned — the manifest's `model_index`
    is an ordinal into that list, so a preparation that kept a different set of
    residues would renumber the sequon before CARBonAra ever saw it.

    What is removed, and why each matters:

    - every other chain, and every glycan, ligand, ion and water, so the input
      matches what ProteinMPNN is given rather than CARBonAra's context-aware
      default;
    - hydrogens and deuterium, which `clean_structure` would drop anyway;
    - the selenium of a selenomethionine, which is not a methionine atom name.

    What is converted: `MSE` becomes `MET`, emitted as `ATOM`. Everything is
    emitted as `ATOM`, because anything left as `HETATM` would be removed by
    `rm_hetatm` and take its backbone position with it.

    Side chains are kept. CARBonAra strips them itself for canonical residues and
    places its own idealised C-beta; doing it here instead would make this arm
    differ from every other use of the model.
    """
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Polypeptide import is_aa

    from .structures import MMCIF_SUFFIXES, _one_letter, _parse_chains

    path = Path(structure_path)
    identifier = str(pdb_id or path.stem)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        warnings.simplefilter("ignore")
        parser = (MMCIFParser(QUIET=True)
                  if path.suffix.lower() in MMCIF_SUFFIXES else PDBParser(QUIET=True))
        try:
            structure = parser.get_structure(identifier, str(path))
        except Exception as exc:
            raise ChainUnreadableError(
                f"{type(exc).__name__}: {str(exc)[:120]}") from exc

    chain = next((c for c in next(iter(structure)) if c.id == str(chain_id)), None)
    if chain is None:
        raise ChainUnreadableError(
            f"chain {chain_id!r} absent from {path.name}")

    lines, residue_ids, sequence, serial, index = [], [], [], 1, 0
    for residue in chain:
        if not is_aa(residue, standard=False) or "CA" not in residue:
            continue
        _, resseq, icode = residue.id
        resname = residue.get_resname().strip().upper()
        selenomethionine = resname == "MSE"
        if selenomethionine:
            resname = "MET"

        index += 1
        for atom in residue:
            name = atom.get_name().strip()
            element = (atom.element or name[:1]).strip().upper()
            if element in ("H", "D"):
                continue
            # By name as well as by element: depositions are unreliable about the
            # element column, and an atom called SE in an MSE is a selenium
            # whatever columns 77-78 claim.
            if selenomethionine and (element == "SE" or name.upper() == "SE"):
                continue
            lines.append(_atom_line(serial, name, resname,
                                    str(chain_id)[:1], index, atom.get_coord(),
                                    element.title() if len(element) > 1 else element))
            serial += 1

        residue_ids.append((int(resseq), str(icode).strip()))
        sequence.append("M" if selenomethionine else _one_letter(residue.get_resname()))

    if not lines:
        raise ChainUnreadableError(
            f"chain {chain_id!r} of {path.name} has no usable protein residues")

    native_sequence = "".join(sequence)
    chains = _parse_chains(path, identifier)
    native = next((c for c in chains if c.chain_id == str(chain_id)), None)
    if native is not None and native.sequence != native_sequence:
        raise ChainUnreadableError(
            f"preparation kept {len(native_sequence)} residues where the "
            f"manifest's parse lists {len(native.sequence)}; the manifest's "
            "indices would not address the prepared chain")

    return "\n".join(lines) + "\nTER\nEND\n", residue_ids, native_sequence


# --------------------------------------------------------------------------
# Index mapping.
# --------------------------------------------------------------------------

def residue_view(structure) -> "tuple[str, list[set[str]]]":
    """CARBonAra's residues as a sequence and a set of atom names each.

    Reproduces how `process_structure` reads identity — `yr` takes each residue's
    resname from its first atom — so what this returns is what the model will be
    imprinted with and read at. A residue outside `std_aminoacids` is `X`: it
    keeps its slot, so nothing shifts, but it has no identity to imprint.
    """
    resids = np.asarray(structure["resid"])
    resnames = np.asarray(structure["resname"])
    names = np.asarray(structure["name"])

    letters, atoms = [], []
    for resid in np.unique(resids):
        selected = resids == resid
        letters.append(RES3TO1.get(str(resnames[selected][0]).strip().upper(), "X"))
        atoms.append({str(name).strip() for name in names[selected]})
    return "".join(letters), atoms


def build_index_map(native_sequence: str, model_sequence: str) -> "dict[int, int]":
    """Manifest index -> CARBonAra index, or `{}` if the two cannot be reconciled.

    `clean_structure` renumbers observed residues consecutively in file order and
    the prepared chain contains exactly the manifest's residues in that order, so
    the correspondence is the identity. That is a claim about the other parser,
    not a fact, so it is checked residue by residue and any disagreement returns
    an empty mapping for the caller to drop the chain.

    A non-canonical residue is expected to read as `X` on CARBonAra's side. One
    that comes back as a canonical letter means the two parses disagree about
    what the residue *is*, which is not something to score around.
    """
    if not native_sequence or not model_sequence:
        return {}
    if len(native_sequence) != len(model_sequence):
        return {}
    for index, letter in enumerate(native_sequence):
        expected = letter if letter in AA_INDEX else "X"
        if model_sequence[index] != expected:
            return {}
    return {index: index for index in range(len(native_sequence))}


class CarbonaraMapping:
    """The correspondence between manifest indices and CARBonAra's residues.

    `pdb_text` is the prepared protein-only chain, carried so that scoring parses
    exactly what the mapping was verified against.
    """

    __slots__ = ("to_model", "model_seq", "backbone_ok", "residue_ids",
                 "manifest_length", "pdb_text")

    def __init__(self, to_model, model_seq, backbone_ok, residue_ids, pdb_text):
        self.to_model = to_model
        self.model_seq = model_seq
        self.backbone_ok = backbone_ok
        self.residue_ids = residue_ids
        self.manifest_length = len(residue_ids)
        self.pdb_text = pdb_text

    def map_indices(self, indices) -> "list[int]":
        """Manifest indices -> CARBonAra indices, or raise."""
        mapped = [self.to_model.get(int(i)) for i in indices]
        missing = [int(i) for i, m in zip(indices, mapped) if m is None]
        if missing:
            raise IncompleteBackboneError(
                f"manifest {'index' if len(missing) == 1 else 'indices'} "
                f"{missing} have no CARBonAra counterpart")
        rejected = [m for m in mapped if not bool(self.backbone_ok[m])]
        if rejected:
            raise IncompleteBackboneError(
                "CARBonAra backbone incomplete (missing one of "
                f"{', '.join(BACKBONE_ATOMS)}) at "
                f"{'index' if len(rejected) == 1 else 'indices'} {rejected}")
        return mapped

    def check_triplet(self, mapped, expected: str) -> None:
        observed = "".join(self.model_seq[m] for m in mapped)
        if observed != expected:
            raise SequonMismatchError(
                f"CARBonAra reads {observed!r} where the manifest records "
                f"{expected!r}")


def chain_mapping(structure_path, chain_id: str, pdb_id: "str | None" = None,
                  carbonara_dir=None) -> CarbonaraMapping:
    """Build and verify the mapping for one chain. Needs no model pass."""
    pdb_text, residue_ids, native_sequence = protein_only_pdb(
        structure_path, chain_id, pdb_id)
    structure = _load_structure(pdb_text, carbonara_dir)
    model_sequence, atoms = residue_view(structure)

    to_model = build_index_map(native_sequence, model_sequence)
    if not to_model:
        raise ChainUnreadableError(
            f"chain {chain_id} of {Path(structure_path).name} cannot be mapped "
            f"onto CARBonAra's parse: it reads {len(model_sequence)} residues "
            f"where the manifest lists {len(native_sequence)}")

    required = set(BACKBONE_ATOMS)
    backbone_ok = np.array([required <= present for present in atoms], dtype=bool)
    return CarbonaraMapping(to_model, model_sequence, backbone_ok, residue_ids,
                            pdb_text)


def decodable_positions(structure_path, chain_id: str,
                        pdb_id: "str | None" = None,
                        carbonara_dir=None) -> np.ndarray:
    """Which manifest indices CARBonAra can evaluate, in MANIFEST index space.

    Manifest space because that is what `05_scoreability.py` indexes with, and
    because scoreability has to be comparable across models. A chain CARBonAra
    cannot read is entirely unscoreable rather than an exception, so one bad
    structure does not stop a sweep of thousands.
    """
    try:
        mapping = chain_mapping(structure_path, chain_id, pdb_id, carbonara_dir)
    except ChainUnreadableError:
        from .structures import _parse_chains

        try:
            chains = _parse_chains(Path(structure_path),
                                   str(pdb_id or Path(structure_path).stem))
            native = next((c for c in chains if c.chain_id == str(chain_id)), None)
            length = len(native.sequence) if native else 0
        except Exception:
            length = 0
        return np.zeros(length, dtype=bool)

    decodable = np.zeros(mapping.manifest_length, dtype=bool)
    for manifest_index, model_index in mapping.to_model.items():
        if (0 <= manifest_index < mapping.manifest_length
                and bool(mapping.backbone_ok[model_index])):
            decodable[manifest_index] = True
    return decodable


# --------------------------------------------------------------------------
# The conditional.
# --------------------------------------------------------------------------

def conditional_probabilities(mapping: CarbonaraMapping, model, positions,
                              device: str = "cpu", carbonara_dir=None
                              ) -> "dict[int, np.ndarray]":
    """P(residue at i | backbone, every other native residue), for each position.

    One forward pass per position, because the imprint differs for each: the
    scored row is zeroed and every other mappable canonical residue carries its
    native identity. Batching them is not safe — the passes differ only in `yt`,
    and sharing one would leak the scored residue into its own conditional.

    Returns calibrated, normalised twenty-vectors keyed by CARBonAra index.
    `apply_model` emits independent sigmoids that do not sum to one, so
    `model.conf` — the checkpoint's own empirical confidence map — runs before
    anything is scored.
    """
    structure = _load_structure(mapping.pdb_text, carbonara_dir)
    X, qe, _, _, Mr, _, y, mr_aa = model.process_structure(structure)

    # Canonical residues only: a residue outside `std_aminoacids` has an all-zero
    # row in `y`, so there is no identity to imprint and claiming one would be a
    # statement about the parser rather than the protein.
    imprintable = mr_aa.float()

    probabilities = {}
    for position in sorted({int(p) for p in positions}):
        known = imprintable.clone()
        known[position] = 0.0
        raw = model.apply_model(X, qe, Mr, yt=known.unsqueeze(1) * y)
        calibrated = model.conf(
            np.asarray(raw.detach().cpu().numpy(), dtype=float))
        probabilities[position] = np.asarray(calibrated[position], dtype=float)
    return probabilities


def check_scoreable(probabilities: "dict[int, np.ndarray]", indices) -> None:
    """Raise unless all three rows are genuine distributions over the twenty.

    Four separate checks, because each catches something the others pass. A
    missing key is a position the model never evaluated; a wrong length is an
    alphabet disagreement; a non-finite entry is a failed calibration; a row that
    does not sum to one is not a probability vector, and putting a logit on it
    would produce a plausible number for nothing at all.
    """
    for index in indices:
        if index not in probabilities:
            raise IncompleteBackboneError(
                f"CARBonAra index {index} was not evaluated")
        row = np.asarray(probabilities[index], dtype=float)
        if row.shape != (len(ALPHABET),):
            raise InvalidProbabilityVector(
                f"row at CARBonAra index {index} has shape {row.shape}, "
                f"not ({len(ALPHABET)},)")
        if not np.all(np.isfinite(row)):
            raise InvalidProbabilityVector(
                f"row at CARBonAra index {index} has non-finite entries")
        if float(row.min()) < 0.0 or float(row.max()) > 1.0 + PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at CARBonAra index {index} has values outside [0, 1]: "
                f"min {float(row.min()):.6g}, max {float(row.max()):.6g}")
        total = float(row.sum())
        if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at CARBonAra index {index} sums to {total:.6g}, not 1")


def sequon_score(probabilities: "dict[int, np.ndarray]", n_index: int,
                 plus1_index: int, plus2_index: int) -> dict:
    """Score one sequon, with the column names every other adapter uses.

    The same statistic: the mean of the log odds of asparagine at the first
    position and of serine-or-threonine at the third. The middle residue is
    excluded because any residue except proline permits a sequon, so a preference
    there is not a preference for the motif; proline is kept as a diagnostic.

    `conditional_sequon_score_sd` is structurally zero and `n_decoding_orders` is
    one. Those columns exist so the models share a schema — they are not evidence
    that CARBonAra is the more precise model.
    """
    check_scoreable(probabilities, (n_index, plus1_index, plus2_index))

    n_row = np.asarray(probabilities[n_index], dtype=float)
    plus1_row = np.asarray(probabilities[plus1_index], dtype=float)
    plus2_row = np.asarray(probabilities[plus2_index], dtype=float)

    p_n = float(n_row[AA_INDEX["N"]])
    p_s = float(plus2_row[AA_INDEX["S"]])
    p_t = float(plus2_row[AA_INDEX["T"]])
    p_pro = float(plus1_row[AA_INDEX["P"]])

    return {
        "conditional_sequon_score": 0.5 * (logit(p_n) + logit(p_s + p_t)),
        "conditional_sequon_score_sd": SCORE_SD,
        "n_decoding_orders": N_ORDERS,
        "p_asn_at_n": p_n,
        "p_ser_at_plus2": p_s,
        "p_thr_at_plus2": p_t,
        "p_ser_or_thr_at_plus2": p_s + p_t,
        "p_pro_at_plus1": p_pro,
        "logit_p_asn": logit(p_n),
        "logit_p_ser_or_thr": logit(p_s + p_t),
        "probs_n": n_row.tolist(),
        "probs_plus1": plus1_row.tolist(),
        "probs_plus2": plus2_row.tolist(),
    }
