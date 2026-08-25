"""Holding positions fixed during fixed-backbone design.

ProteinMPNN reads `chain_M_pos` to decide what it may change: 1 designs the
position, 0 keeps the native residue. Fixing a sequon is therefore a mask
applied before sampling, not a repair of the output afterwards. The distinction
matters — repairing afterwards lets the model condition on residues it was
going to overwrite, which measures something else.

The policy names follow SugarFix, so the two codebases mean the same thing:

    full_sequon           fix N, X and S/T
    functional_preserve   fix N and S/T, leaving X free to change
"""
from __future__ import annotations

import numpy as np

POLICIES = ("full_sequon", "functional_preserve")


def sequon_positions(n_index: int, policy: str = "full_sequon") -> "list[int]":
    """Chain indices to hold fixed for one sequon, given the Asn's index."""
    if policy == "full_sequon":
        return [n_index, n_index + 1, n_index + 2]
    if policy == "functional_preserve":
        return [n_index, n_index + 2]
    raise ValueError(f"unknown policy {policy!r}; expected one of {POLICIES}")


from experimental_glycosylation_sites.retention import design_mask  # re-exported


def native_sequence(pdb_path, chain_id: str) -> str:
    """The chain's deposited sequence, indexed exactly as designs are.

    Read through ProteinMPNN's own parser rather than BioPython, because every
    index in this experiment refers to that parse. Taking it from anywhere else
    would reintroduce the alignment question the shell gate exists to close.
    """
    from experimental_glycosylation_sites.mpnn_scoring import (_ensure_importable,
                                                               _prepare_environment)
    from experimental_glycosylation_sites.runner_support import proteinmpnn_dir
    _prepare_environment()
    _ensure_importable(proteinmpnn_dir())
    from protein_mpnn_utils import StructureDatasetPDB, parse_PDB

    parsed = parse_PDB(str(pdb_path), input_chain_list=[str(chain_id)])
    protein = StructureDatasetPDB(parsed, truncate=None, max_length=20000)[0]
    return protein[f"seq_chain_{chain_id}"]


def verify_sequon_index(sequence: str, n_index: int, expected_triplet: str) -> bool:
    """Whether the claimed index really points at the expected sequon.

    Checked directly against the decoded sequence, independently of any
    structure parse. This is the assertion that catches an off-by-one before it
    protects the wrong three residues and reports a clean run -- the failure
    that produced a day of rework when it went unchecked.
    """
    if not (0 <= n_index and n_index + 3 <= len(sequence)):
        return False
    observed = sequence[n_index:n_index + 3]
    if len(expected_triplet) == 3 and observed != expected_triplet:
        return False
    return observed[0] == "N" and observed[2] in ("S", "T") and observed[1] != "P"
