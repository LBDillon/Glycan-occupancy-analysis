from __future__ import annotations

import numpy as np
import pandas as pd

from experimental_glycosylation_sites.matching import (
    MATCH_FEATURES,
    balance_report,
    match_controls,
    match_controls_optimal,
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


def test_exact_column_blocks_a_closer_but_mismatched_control():
    """Exact matching must override distance, not merely be one term within it."""
    cases = pd.DataFrame([{"accession": "P1", "position": 10, "rsa": 0.50,
                           "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.5,
                           "subtype": "NXT"}])
    controls = pd.DataFrame([
        # nearer, wrong subtype — must not be chosen
        {"accession": "C1", "position": 1, "rsa": 0.50, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.5, "subtype": "NXS"},
        # further, right subtype
        {"accession": "C2", "position": 2, "rsa": 0.52, "n_neighbours_8a": 11.0,
         "hydrophobic_fraction_8a": 0.52, "subtype": "NXT"},
    ])
    pairs = match_controls(cases, controls, k=1, caliper=10.0, exact=("subtype",))
    assert list(pairs.control_accession) == ["C2"]


def test_no_pair_is_emitted_when_no_control_shares_the_subtype():
    """A case with no same-subtype control is dropped, not matched across."""
    cases = pd.DataFrame([{"accession": "P1", "position": 10, "rsa": 0.5,
                           "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.5,
                           "subtype": "NXT"}])
    controls = pd.DataFrame([{"accession": "C1", "position": 1, "rsa": 0.5,
                              "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.5,
                              "subtype": "NXS"}])
    assert match_controls(cases, controls, k=1, caliper=10.0, exact=("subtype",)).empty


def test_every_matched_pair_shares_the_exact_column():
    """The invariant the primary comparison relies on, over a larger draw."""
    rng = np.random.default_rng(0)
    cases = pd.DataFrame({
        "accession": [f"P{i}" for i in range(40)], "position": range(40),
        "rsa": rng.uniform(0, 1, 40), "n_neighbours_8a": rng.uniform(5, 20, 40),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 40),
        "subtype": rng.choice(["NXS", "NXT"], 40)})
    controls = pd.DataFrame({
        "accession": [f"C{i}" for i in range(60)], "position": range(60),
        "rsa": rng.uniform(0, 1, 60), "n_neighbours_8a": rng.uniform(5, 20, 60),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 60),
        "subtype": rng.choice(["NXS", "NXT"], 60)})
    pairs = match_controls(cases, controls, k=3, caliper=1.0, exact=("subtype",))
    assert len(pairs) > 0
    merged = (pairs.merge(cases[["accession", "position", "subtype"]],
                          left_on=["case_accession", "case_position"],
                          right_on=["accession", "position"])
                   .merge(controls[["accession", "position", "subtype"]],
                          left_on=["control_accession", "control_position"],
                          right_on=["accession", "position"], suffixes=("_case", "_ctrl")))
    assert (merged.subtype_case == merged.subtype_ctrl).all()


# --------------------------------------------------------------------------
# Deterministic optimal matching
# --------------------------------------------------------------------------

def test_optimal_matching_is_independent_of_input_order():
    """The property greedy matching lacks: no dependence on ordering or seed."""
    rng = np.random.default_rng(1)
    cases = pd.DataFrame({
        "accession": [f"P{i}" for i in range(25)], "position": range(25),
        "rsa": rng.uniform(0, 1, 25), "n_neighbours_8a": rng.uniform(5, 20, 25),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 25),
        "subtype": rng.choice(["NXS", "NXT"], 25)})
    controls = pd.DataFrame({
        "accession": [f"C{i}" for i in range(12)], "position": range(12),
        "rsa": rng.uniform(0, 1, 12), "n_neighbours_8a": rng.uniform(5, 20, 12),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 12),
        "subtype": rng.choice(["NXS", "NXT"], 12)})

    a = match_controls_optimal(cases, controls, caliper=0.6, exact=("subtype",))
    b = match_controls_optimal(cases.iloc[::-1], controls.iloc[::-1],
                               caliper=0.6, exact=("subtype",))
    pd.testing.assert_frame_equal(a, b)


def test_optimal_matching_never_beats_greedy_on_pair_count():
    """Maximum cardinality: greedy can tie it, never exceed it."""
    rng = np.random.default_rng(2)
    cases = pd.DataFrame({
        "accession": [f"P{i}" for i in range(30)], "position": range(30),
        "rsa": rng.uniform(0, 1, 30), "n_neighbours_8a": rng.uniform(5, 20, 30),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 30),
        "subtype": rng.choice(["NXS", "NXT"], 30)})
    controls = pd.DataFrame({
        "accession": [f"C{i}" for i in range(14)], "position": range(14),
        "rsa": rng.uniform(0, 1, 14), "n_neighbours_8a": rng.uniform(5, 20, 14),
        "hydrophobic_fraction_8a": rng.uniform(0, 1, 14),
        "subtype": rng.choice(["NXS", "NXT"], 14)})

    best = match_controls_optimal(cases, controls, caliper=0.5, exact=("subtype",))
    for seed in range(40):
        greedy = match_controls(cases, controls, k=1, caliper=0.5,
                                seed=seed, exact=("subtype",))
        assert len(greedy) <= len(best)


def test_optimal_matching_respects_caliper_and_exact_columns():
    best = match_controls_optimal(
        pd.DataFrame([{"accession": "P1", "position": 1, "rsa": 0.5,
                       "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.5,
                       "subtype": "NXT"}]),
        pd.DataFrame([{"accession": "C1", "position": 1, "rsa": 0.5,
                       "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.5,
                       "subtype": "NXS"}]),
        caliper=10.0, exact=("subtype",))
    assert best.empty


def test_optimal_matching_uses_each_control_once():
    rng = np.random.default_rng(3)
    cases = pd.DataFrame({
        "accession": [f"P{i}" for i in range(20)], "position": range(20),
        "rsa": rng.uniform(0.4, 0.6, 20), "n_neighbours_8a": rng.uniform(9, 11, 20),
        "hydrophobic_fraction_8a": rng.uniform(0.4, 0.6, 20), "subtype": "NXT"})
    controls = pd.DataFrame({
        "accession": [f"C{i}" for i in range(8)], "position": range(8),
        "rsa": rng.uniform(0.4, 0.6, 8), "n_neighbours_8a": rng.uniform(9, 11, 8),
        "hydrophobic_fraction_8a": rng.uniform(0.4, 0.6, 8), "subtype": "NXT"})
    best = match_controls_optimal(cases, controls, caliper=5.0, exact=("subtype",))
    assert len(best) == 8
    assert best.control_accession.nunique() == 8
    assert best.case_accession.nunique() == 8


def test_optimal_matching_breaks_ties_canonically():
    """Two controls equidistant from one case: the choice must not depend on
    the caller's row order, so it falls to accession then position."""
    case = pd.DataFrame([{"accession": "P1", "position": 10, "rsa": 0.50,
                          "n_neighbours_8a": 10.0, "hydrophobic_fraction_8a": 0.50,
                          "subtype": "NXT"}])
    # symmetric about the case in the first feature, identical in the others
    controls = pd.DataFrame([
        {"accession": "C_bbb", "position": 5, "rsa": 0.60, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
        {"accession": "C_aaa", "position": 5, "rsa": 0.40, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
    ])
    chosen = {
        tuple(match_controls_optimal(case, controls.iloc[order], caliper=5.0,
                                     exact=("subtype",)).control_accession)
        for order in ([0, 1], [1, 0])
    }
    assert len(chosen) == 1, "tie broken differently for different input orders"
    assert chosen.pop() == ("C_aaa",)


def test_optimal_matching_recovers_a_pair_greedy_loses():
    """The concrete cost of greedy matching on a scarce control pool.

    Both cases sit close to C_shared, and only P_first is within the caliper of
    C_far. A greedy pass that reaches P_first first takes C_shared — its nearest
    — and strands P_second, which has no admissible partner left: one pair. The
    optimal assignment spends P_first on the further-but-admissible C_far and
    keeps C_shared for P_second: two pairs.

    Only `rsa` varies, so distances are proportional to the gaps below.
    """
    cases = pd.DataFrame([
        {"accession": "P_first", "position": 1, "rsa": 0.33, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
        {"accession": "P_second", "position": 2, "rsa": 0.28, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
    ])
    controls = pd.DataFrame([
        {"accession": "C_shared", "position": 1, "rsa": 0.30, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
        {"accession": "C_far", "position": 2, "rsa": 0.62, "n_neighbours_8a": 10.0,
         "hydrophobic_fraction_8a": 0.50, "subtype": "NXT"},
    ])
    # scaled: P_first-C_far 1.82, P_second-C_far 2.13, so the caliper admits the
    # first and excludes the second
    caliper = 2.0

    best = match_controls_optimal(cases, controls, caliper=caliper, exact=("subtype",))
    assert len(best) == 2
    assert dict(zip(best.case_accession, best.control_accession)) == {
        "P_first": "C_far", "P_second": "C_shared"}

    greedy_counts = {len(match_controls(cases, controls, k=1, caliper=caliper,
                                        seed=seed, exact=("subtype",)))
                     for seed in range(30)}
    assert 1 in greedy_counts, "expected some seed ordering to strand a case"
    assert max(greedy_counts) <= len(best)
