"""Structural context features for scoring experiments.

These describe the local environment of a site — how exposed it is, what
surrounds it, where it sits in the chain. They are matching variables, not
evidence: occupied sites sit disproportionately in exposed loops, so any
comparison of model scores between occupied and observed-unmodified sites has to
control for context before a residual difference can mean anything.

Features come only from experimental coordinates. Predicted models are excluded
deliberately: a predicted structure contains no glycans, and its local
conformation near a real glycosylation site may be systematically wrong in a way
that correlates with the label being tested.
"""
from __future__ import annotations

import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

# Theoretical maximum solvent accessibility, Tien et al. 2013 (PLoS ONE),
# empirical Gly-X-Gly values. Used to turn absolute SASA into RSA.
MAX_ASA = {
    "A": 121.0, "R": 265.0, "N": 187.0, "D": 187.0, "C": 148.0,
    "E": 214.0, "Q": 214.0, "G": 97.0, "H": 216.0, "I": 195.0,
    "L": 191.0, "K": 230.0, "M": 203.0, "F": 228.0, "P": 154.0,
    "S": 143.0, "T": 163.0, "W": 264.0, "Y": 255.0, "V": 165.0,
}

HYDROPHOBIC = set("AVLIMFWCY")
CHARGED = set("DEKR")
NEIGHBOUR_RADIUS = 8.0


def _rsa_bin(rsa: float | None) -> str:
    if rsa is None:
        return ""
    if rsa < 0.09:
        return "buried"
    if rsa < 0.36:
        return "intermediate"
    return "exposed"


# Bounded deliberately. Sites are processed in accession order, so consecutive
# lookups hit the same structure and a small cache captures nearly all the reuse.
# An unbounded one would hold every parsed model in memory at once — thousands of
# structures across the control sets — and exhaust RAM long before it helped.
_MODEL_CACHE: "OrderedDict[str, object]" = OrderedDict()
_MODEL_CACHE_LIMIT = 4


# Solvent accessibility is computed over the whole structure, so cost scales with
# the entire assembly rather than the one residue of interest. Cytosolic proteins
# are frequently components of ribosomes and proteasomes, whose entries run to
# hundreds of megabytes and take many minutes each. The residue environment is the
# same in a smaller structure of the same protein, so oversized entries are skipped
# rather than allowed to dominate a run.
MAX_STRUCTURE_BYTES = 20_000_000


def _model_with_sasa(path: Path):
    """Parsed first model with per-residue SASA attached, cached per file.

    SASA is a whole-structure calculation; doing it per site would repeat the
    same work for every residue of the same protein.
    """
    key = str(path)
    if key in _MODEL_CACHE:
        return _MODEL_CACHE[key]

    if Path(path).exists() and Path(path).stat().st_size > MAX_STRUCTURE_BYTES:
        _MODEL_CACHE[key] = None
        return None

    model = None
    with warnings.catch_warnings():
        simplefilter = warnings.simplefilter
        simplefilter("ignore", PDBConstructionWarning)
        simplefilter("ignore")
        try:
            parser = (
                MMCIFParser(QUIET=True)
                if Path(path).suffix.lower() in {".cif", ".mmcif"}
                else PDBParser(QUIET=True)
            )
            structure = parser.get_structure("s", str(path))
            model = next(iter(structure), None)
            if model is not None:
                # Strip non-amino-acid residues BEFORE computing accessibility.
                # An occupied site carries its glycan in the same structure, and
                # that sugar occludes the residue it is attached to — so leaving
                # it in makes occupied sites look systematically more buried than
                # unmodified ones. The burial would then be an artefact of the
                # label, not a property of the protein, and would confound the
                # very comparison these features exist to support.
                for chain in list(model):
                    for residue in list(chain):
                        if not is_aa(residue, standard=False):
                            chain.detach_child(residue.id)
                ShrakeRupley().compute(model, level="R")
        except Exception:
            model = None

    _MODEL_CACHE[key] = model
    while len(_MODEL_CACHE) > _MODEL_CACHE_LIMIT:
        _MODEL_CACHE.popitem(last=False)
    return model


def _clean_icode(value) -> str:
    """Insertion codes arrive as NaN from CSV; str(nan) would never match."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def residue_features(
    path: Path, chain_id: str, resseq: int, icode: str = ""
) -> dict | None:
    """Local environment of one residue, or None if it cannot be located."""
    model = _model_with_sasa(path)
    if model is None:
        return None

    icode = _clean_icode(icode)
    try:
        chain = model[chain_id]
    except KeyError:
        return None

    target = None
    for residue in chain:
        if not is_aa(residue, standard=False):
            continue
        _, seq, code = residue.id
        if int(seq) == int(resseq) and str(code).strip() == icode:
            target = residue
            break
    if target is None or "CA" not in target:
        return None

    from Bio.SeqUtils import seq1

    try:
        code1 = seq1(target.get_resname(), undef_code="X")
    except Exception:
        code1 = "X"

    sasa = float(getattr(target, "sasa", 0.0) or 0.0)
    max_asa = MAX_ASA.get(code1)
    rsa = min(sasa / max_asa, 1.5) if max_asa else None

    # Neighbourhood composition, a coarse read on what the site sits against.
    origin = target["CA"].coord
    neighbours = []
    for other_chain in model:
        for residue in other_chain:
            if residue is target or not is_aa(residue, standard=False):
                continue
            if "CA" not in residue:
                continue
            if float(np.linalg.norm(residue["CA"].coord - origin)) <= NEIGHBOUR_RADIUS:
                try:
                    neighbours.append(seq1(residue.get_resname(), undef_code="X"))
                except Exception:
                    neighbours.append("X")

    resolved = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    order = [int(r.id[1]) for r in resolved]
    n_neighbours = len(neighbours)

    return {
        "sasa": round(sasa, 2),
        "rsa": round(rsa, 4) if rsa is not None else None,
        "rsa_bin": _rsa_bin(rsa),
        "observed_residue": code1,
        "n_neighbours_8a": n_neighbours,
        "hydrophobic_fraction_8a": (
            round(sum(1 for r in neighbours if r in HYDROPHOBIC) / n_neighbours, 4)
            if n_neighbours else None
        ),
        "charged_fraction_8a": (
            round(sum(1 for r in neighbours if r in CHARGED) / n_neighbours, 4)
            if n_neighbours else None
        ),
        "chain_length_resolved": len(resolved),
        "distance_to_chain_terminus": (
            min(abs(int(resseq) - min(order)), abs(max(order) - int(resseq)))
            if order else None
        ),
    }


FEATURE_COLUMNS = [
    "sasa", "rsa", "rsa_bin", "observed_residue", "n_neighbours_8a",
    "hydrophobic_fraction_8a", "charged_fraction_8a",
    "chain_length_resolved", "distance_to_chain_terminus",
]


def build_features(sites: pd.DataFrame, structure_paths: dict[str, Path]) -> pd.DataFrame:
    """One row per site, with features where the residue could be located.

    `sites` needs accession, position, structure_pdb_id, structure_chain_id,
    structure_resseq and structure_icode. Sites without mapped coordinates are
    kept with empty features so the output stays a complete account of the input
    rather than a silently filtered subset.
    """
    cache: dict[tuple, dict | None] = {}
    rows = []

    for row in sites.itertuples(index=False):
        record = {
            "accession": row.accession,
            "position": int(row.position),
            "occupancy_status": getattr(row, "occupancy_status", ""),
            "structure_pdb_id": getattr(row, "structure_pdb_id", ""),
            "structure_chain_id": getattr(row, "structure_chain_id", ""),
            "structure_resseq": getattr(row, "structure_resseq", None),
        }
        pdb_id = str(record["structure_pdb_id"] or "").upper()
        resseq = record["structure_resseq"]
        path = structure_paths.get(pdb_id)

        features = None
        if path is not None and pd.notna(resseq):
            key = (pdb_id, record["structure_chain_id"], int(resseq))
            if key not in cache:
                cache[key] = residue_features(
                    path, str(record["structure_chain_id"]), int(resseq),
                    _clean_icode(getattr(row, "structure_icode", "")),
                )
            features = cache[key]

        record["features_available"] = features is not None
        for column in FEATURE_COLUMNS:
            record[column] = features.get(column) if features else None
        rows.append(record)

    return (
        pd.DataFrame(rows)
        .sort_values(["accession", "position"])
        .reset_index(drop=True)
    )
