"""Invariants the context table must satisfy, checked rather than assumed.

Each of these corresponds to a defect that reached a finished table once. A QC
report that only counts rows would have passed all of them.
"""
from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.context_qc import check_invariants

CLEAN = {
    "accession": "P1", "position": 1,
    "distance_to_n_terminus_resolved": 10, "distance_to_c_terminus_resolved": 89,
    "chain_length_resolved": 100,
    "triplet_observed": "NAS", "has_nd2": True,
    "nd2_atoms_8a_same_chain": 20, "nd2_residues_8a_same_chain": 5,
    "loop_run_length": 4, "n_dssp_ok": True, "n_ss": "H",
}


def _names(violations):
    return {v["invariant"] for v in violations}


def test_clean_row_has_no_violations():
    assert check_invariants(pd.DataFrame([CLEAN])) == []


def test_terminal_distances_must_sum_to_the_chain_length():
    """The defect that left 1,218 rows inconsistent: author-number arithmetic
    against a residue count."""
    row = {**CLEAN, "distance_to_c_terminus_resolved": 500}
    assert "terminal_distances_sum_to_chain_length" in _names(
        check_invariants(pd.DataFrame([row])))


def test_nd2_counts_must_be_absent_when_there_is_no_nd2():
    """A missing side chain must not read as an uncrowded site."""
    row = {**CLEAN, "has_nd2": False}
    assert "nd2_counts_present_iff_nd2" in _names(
        check_invariants(pd.DataFrame([row])))


def test_secondary_structure_requires_dssp_availability():
    row = {**CLEAN, "n_dssp_ok": False, "n_ss": "H"}
    assert "ss_requires_dssp_ok" in _names(check_invariants(pd.DataFrame([row])))


def test_triplet_must_have_three_positions():
    row = {**CLEAN, "triplet_observed": "NA"}
    assert "triplet_has_three_positions" in _names(
        check_invariants(pd.DataFrame([row])))


def test_violations_report_how_many_rows_and_an_example():
    row = {**CLEAN, "distance_to_c_terminus_resolved": 500}
    violation = check_invariants(pd.DataFrame([row, row]))[0]
    assert violation["rows"] == 2
    assert violation["example"]["accession"] == "P1"
