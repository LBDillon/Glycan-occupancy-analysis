"""Row-level checks on the sequence context behind each site.

A population-level check ("does this population have any sequences?") passes a
table that lost half a cache. These are per row, so partial loss is detectable,
and every failure carries a reason a person can act on.
"""
from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.sequence_qc import sequence_context_failures

CLEAN = {"accession": "P1", "position": 10, "population": "occupied",
         "sequon_triplet": "NAS", "uniprot_length": 100}


def _reasons(rows):
    failures = sequence_context_failures(pd.DataFrame(rows))
    return dict(zip(failures.accession, failures.reason))


def test_clean_row_passes():
    assert len(sequence_context_failures(pd.DataFrame([CLEAN]))) == 0


def test_missing_sequence_is_reported():
    """Partial cache loss shows up as some rows without a length, not all."""
    rows = [CLEAN, {**CLEAN, "accession": "P2", "uniprot_length": None}]
    assert _reasons(rows) == {"P2": "no_uniprot_sequence"}


def test_incomplete_triplet_is_reported():
    rows = [{**CLEAN, "accession": "P2", "sequon_triplet": "NA"}]
    assert _reasons(rows) == {"P2": "incomplete_expected_triplet"}


def test_triplet_must_start_with_asparagine():
    rows = [{**CLEAN, "accession": "P2", "sequon_triplet": "QAS"}]
    assert _reasons(rows) == {"P2": "expected_triplet_not_a_sequon"}


def test_triplet_must_end_in_serine_or_threonine():
    rows = [{**CLEAN, "accession": "P2", "sequon_triplet": "NAA"}]
    assert _reasons(rows) == {"P2": "expected_triplet_not_a_sequon"}


def test_proline_at_the_x_position_is_not_a_sequon():
    """N-P-S/T is not glycosylated; if one is in the set the mapping is wrong."""
    rows = [{**CLEAN, "accession": "P2", "sequon_triplet": "NPS"}]
    assert _reasons(rows) == {"P2": "expected_triplet_not_a_sequon"}


def test_sequon_running_past_the_sequence_end_is_reported():
    rows = [{**CLEAN, "accession": "P2", "position": 99, "uniprot_length": 100}]
    assert _reasons(rows) == {"P2": "sequon_outside_sequence"}


def test_non_positive_position_is_reported():
    rows = [{**CLEAN, "accession": "P2", "position": 0}]
    assert _reasons(rows) == {"P2": "sequon_outside_sequence"}


def test_failures_carry_the_population_for_reporting():
    rows = [{**CLEAN, "accession": "P2", "population": "secretory_unannotated",
             "uniprot_length": None}]
    failures = sequence_context_failures(pd.DataFrame(rows))
    assert failures.iloc[0].population == "secretory_unannotated"


def test_sequence_distances_are_derived_from_the_manifest():
    from experimental_glycosylation_sites.sequence_qc import add_sequence_distances
    frame = add_sequence_distances(pd.DataFrame([CLEAN]))
    assert frame.iloc[0].uniprot_residues_after_asn == 90     # 100 - 10
    assert frame.iloc[0].uniprot_residues_after_sequon == 88  # 100 - 12


def test_sequence_distances_reject_a_negative_result():
    """Negative means the sequon runs past the sequence end, which the
    row-level checks should already have excluded -- so it is a bug, not data."""
    import pytest
    from experimental_glycosylation_sites.sequence_qc import add_sequence_distances
    bad = pd.DataFrame([{**CLEAN, "position": 120, "uniprot_length": 100}])
    with pytest.raises(ValueError, match="negative"):
        add_sequence_distances(bad)


def test_sequence_distances_are_absent_without_a_length():
    from experimental_glycosylation_sites.sequence_qc import add_sequence_distances
    frame = add_sequence_distances(pd.DataFrame([{**CLEAN, "uniprot_length": None}]))
    assert pd.isna(frame.iloc[0].uniprot_residues_after_asn)
