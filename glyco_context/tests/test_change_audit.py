"""Attributing every difference between the old and corrected tables.

The point is not to summarise the diff but to prove it: each changed row must
trace to a named fix. A row that changed for no reason we can name is the
signal that something else moved, and it is reported rather than absorbed.
"""
from __future__ import annotations

import pandas as pd

from glyco_context.change_audit import attribute_changes

OLD = {"accession": "P1", "position": 1, "population": "occupied",
       "triplet_observed": "NAS", "triplet_matches": True, "dssp_ok": True,
       "n_phi": -60.0, "distance_to_n_terminus_resolved": 10,
       "distance_to_c_terminus_resolved": 89}
NEW = {**OLD, "mapping_continuous": True, "n_icode": ""}


def _audit(old_row, new_row, **kwargs):
    return attribute_changes(pd.DataFrame([old_row]), pd.DataFrame([new_row]), **kwargs)


def test_identical_rows_produce_no_findings():
    assert len(_audit(OLD, NEW)) == 0


def test_insertion_code_change_is_named():
    old = {**OLD, "triplet_observed": "LKN", "triplet_matches": False}
    new = {**NEW, "triplet_observed": "NAS", "n_icode": "B"}
    audit = _audit(old, new)
    assert audit.iloc[0].category == "insertion_code_propagation"


def test_invisible_gap_jump_is_named_separately():
    """Matched before, discontinuous now: the nine silent rows."""
    new = {**NEW, "triplet_observed": "N??", "triplet_matches": False,
           "mapping_continuous": False}
    assert _audit(OLD, new).iloc[0].category == "invisible_gap_jump"


def test_gap_aware_mapping_change_is_named():
    old = {**OLD, "triplet_observed": "NEY", "triplet_matches": False}
    new = {**NEW, "triplet_observed": "N??", "triplet_matches": False,
           "mapping_continuous": False}
    assert _audit(old, new).iloc[0].category == "gap_aware_mapping"


def test_dssp_gained_is_named():
    old = {**OLD, "dssp_ok": False}
    assert _audit(old, NEW).iloc[0].category == "dssp_multichain_repair"


def test_dihedral_change_is_named():
    new = {**NEW, "n_phi": None}
    assert _audit(OLD, new).iloc[0].category == "dihedral_no_longer_crosses_gap"


def test_terminal_distance_change_is_named():
    new = {**NEW, "distance_to_c_terminus_resolved": 89 + 40}
    assert _audit(OLD, new).iloc[0].category == "terminal_distance_units"


def test_unexplained_change_is_flagged_not_absorbed():
    new = {**NEW, "triplet_observed": "NAS", "n_rsa": 0.9}
    old = {**OLD, "n_rsa": 0.1}
    audit = _audit(old, new)
    assert audit.iloc[0].category == "UNEXPLAINED"
    assert "n_rsa" in audit.iloc[0].changed_columns
