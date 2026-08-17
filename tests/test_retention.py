from __future__ import annotations

import pytest

from experimental_glycosylation_sites.retention import (
    PREPRINT_CONDITION,
    RETENTION_CATEGORIES,
    STANDARD_CONDITION,
    classify_retention,
)

# index 0 = N position, 1 = X, 2 = S/T position
IDX = (0, 1, 2)


def test_intact_sequon_counts_as_full_retention():
    r = classify_retention(["NKT", "NAS"], *IDX)
    assert r["frac_full_sequon_retained"] == 1.0
    assert r["frac_complete_motif_loss"] == 0.0


def test_proline_at_x_abolishes_the_sequon_even_with_n_and_t():
    """N-P-T is not a sequon: proline at +1 blocks glycosylation."""
    r = classify_retention(["NPT"], *IDX)
    assert r["frac_full_sequon_retained"] == 0.0
    assert r["frac_proline_introduced_at_x"] == 1.0
    # the individual residues did survive, and that is reported separately
    assert r["frac_asn_retained_motif_lost"] == 1.0
    assert r["frac_ser_thr_retained_motif_lost"] == 1.0


def test_losing_the_asparagine_is_not_complete_loss_if_hydroxyl_remains():
    r = classify_retention(["QKT"], *IDX)
    assert r["frac_full_sequon_retained"] == 0.0
    assert r["frac_ser_thr_retained_motif_lost"] == 1.0
    assert r["frac_complete_motif_loss"] == 0.0


def test_losing_both_is_complete_loss():
    r = classify_retention(["QKA"], *IDX)
    assert r["frac_complete_motif_loss"] == 1.0
    assert r["frac_asn_retained_motif_lost"] == 0.0
    assert r["frac_ser_thr_retained_motif_lost"] == 0.0


def test_full_retention_and_complete_loss_are_mutually_exclusive():
    r = classify_retention(["NKT", "QKA", "NKA", "QKT"], *IDX)
    assert r["frac_full_sequon_retained"] + r["frac_complete_motif_loss"] <= 1.0
    assert r["frac_full_sequon_retained"] == 0.25
    assert r["frac_complete_motif_loss"] == 0.25


def test_fractions_use_the_number_actually_scored():
    r = classify_retention(["NKT", "NKT", "QKA", "QKA"], *IDX)
    assert r["n_designs_scored"] == 4
    assert r["frac_full_sequon_retained"] == 0.5


def test_designs_too_short_are_skipped_not_miscounted():
    r = classify_retention(["NKT", "NK"], *IDX)
    assert r["n_designs_scored"] == 1
    assert r["frac_full_sequon_retained"] == 1.0


def test_no_scoreable_designs_returns_nulls_not_zeros():
    """Zero would read as 'never retained'; absent must read as absent."""
    r = classify_retention(["NK"], *IDX)
    assert r["n_designs_scored"] == 0
    assert all(r[f"frac_{k}"] is None for k in RETENTION_CATEGORIES)


def test_conditions_match_the_preprint_and_raise_only_the_sample_count():
    assert PREPRINT_CONDITION["temperature"] == STANDARD_CONDITION["temperature"] == 0.1
    assert PREPRINT_CONDITION["n_designs"] == 8
    assert STANDARD_CONDITION["n_designs"] > PREPRINT_CONDITION["n_designs"]
