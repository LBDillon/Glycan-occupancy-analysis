from __future__ import annotations

import numpy as np
import pandas as pd

from experimental_glycosylation_sites.matching import (
    MATCH_FEATURES,
    balance_report,
    match_controls,
    standardised_mean_difference,
)


def frame(n, rsa, neighbours, hydro=0.4, prefix="C", seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "accession": [f"{prefix}{i}" for i in range(n)],
        "position": range(1, n + 1),
        "rsa": rng.normal(rsa, 0.05, n),
        "n_neighbours_8a": rng.normal(neighbours, 1.0, n),
        "hydrophobic_fraction_8a": rng.normal(hydro, 0.05, n),
    })


def test_smd_is_zero_for_identical_groups():
    a = pd.Series([1.0, 2.0, 3.0, 4.0])
    assert standardised_mean_difference(a, a.copy()) == 0.0


def test_smd_sign_says_which_group_is_higher():
    low, high = pd.Series([1.0, 2.0, 3.0]), pd.Series([4.0, 5.0, 6.0])
    assert standardised_mean_difference(high, low) > 0
    assert standardised_mean_difference(low, high) < 0


def test_matching_reduces_imbalance():
    """The whole point: cases and controls differ before, and should not after."""
    cases = frame(60, rsa=0.45, neighbours=8, prefix="A", seed=1)
    # a control pool spanning a wider range, so well-matched controls exist inside it
    controls = frame(600, rsa=0.32, neighbours=10, prefix="B", seed=2)
    pairs = match_controls(cases, controls, k=3)
    report = balance_report(cases, controls, pairs)
    before = abs(report["features"]["rsa"]["smd_before"])
    after = abs(report["features"]["rsa"]["smd_after"])
    assert before > 0.5, "fixture should start imbalanced"
    assert after < before, "matching must reduce imbalance"


def test_controls_are_not_reused():
    """Reusing a control would inflate apparent power."""
    cases = frame(20, rsa=0.4, neighbours=8, prefix="A", seed=3)
    controls = frame(200, rsa=0.4, neighbours=8, prefix="B", seed=4)
    pairs = match_controls(cases, controls, k=4)
    used = pairs[["control_accession", "control_position"]]
    assert not used.duplicated().any()


def test_caliper_rejects_distant_controls():
    """An unmatched case is visible; a badly matched one is not."""
    cases = frame(20, rsa=0.9, neighbours=3, prefix="A", seed=5)
    controls = frame(200, rsa=0.1, neighbours=25, prefix="B", seed=6)
    pairs = match_controls(cases, controls, k=3, caliper=0.05)
    assert len(pairs) == 0


def test_unmatched_cases_are_counted_not_hidden():
    cases = frame(20, rsa=0.9, neighbours=3, prefix="A", seed=7)
    controls = frame(200, rsa=0.1, neighbours=25, prefix="B", seed=8)
    report = balance_report(cases, controls, match_controls(cases, controls, caliper=0.05))
    assert report["cases_unmatched"] == 20
    assert report["cases_matched"] == 0


def test_matching_is_reproducible():
    cases = frame(30, rsa=0.4, neighbours=8, prefix="A", seed=9)
    controls = frame(300, rsa=0.35, neighbours=9, prefix="B", seed=10)
    first = match_controls(cases, controls, seed=42)
    second = match_controls(cases, controls, seed=42)
    pd.testing.assert_frame_equal(first, second)


def test_rows_with_missing_features_are_excluded():
    cases = frame(10, rsa=0.4, neighbours=8, prefix="A", seed=11)
    cases.loc[0, "rsa"] = np.nan
    controls = frame(100, rsa=0.4, neighbours=8, prefix="B", seed=12)
    pairs = match_controls(cases, controls)
    assert "A0" not in set(pairs.case_accession)
