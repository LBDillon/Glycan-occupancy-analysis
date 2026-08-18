"""Contrast construction, resampling units and the bootstrap.

These are shared by the primary analysis and the matching-sensitivity sweep, so
a defect here would move a result and its own sensitivity check together and
look like agreement.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from experimental_glycosylation_sites.contrasts import (
    assign_resample_units,
    build_contrasts,
    classify,
    cluster_bootstrap,
)


def site_frame(rows):
    """rows: (accession, position, score, subtype, ortholog_clusters)"""
    return pd.DataFrame(
        [{"accession": a, "position": p, "conditional_sequon_score": s,
          "subtype": t, "ortholog_clusters": c} for a, p, s, t, c in rows]
    )


def pair_frame(rows):
    """rows: (case_accession, case_position, control_accession, control_position)"""
    return pd.DataFrame(
        [{"case_accession": ca, "case_position": cp,
          "control_accession": na, "control_position": np_, "distance": 0.1,
          "match_rank": 1} for ca, cp, na, np_ in rows]
    )


# --------------------------------------------------------------------------
# build_contrasts
# --------------------------------------------------------------------------

def test_contrast_is_case_minus_control():
    site = site_frame([("P1", 10, 1.5, "NXT", "cl1"), ("C1", 5, 0.5, "NXT", None)])
    result = build_contrasts(pair_frame([("P1", 10, "C1", 5)]), site)
    assert len(result) == 1
    assert result.contrast.iloc[0] == pytest.approx(1.0)
    assert result.case_score.iloc[0] == pytest.approx(1.5)
    assert result.control_mean_score.iloc[0] == pytest.approx(0.5)


def test_several_controls_are_averaged_into_one_contrast():
    """One contrast per occupied site, whatever its matched-set size."""
    site = site_frame([("P1", 10, 2.0, "NXT", "cl1"),
                       ("C1", 1, 0.0, "NXT", None), ("C2", 2, 1.0, "NXT", None)])
    result = build_contrasts(pair_frame([("P1", 10, "C1", 1), ("P1", 10, "C2", 2)]), site)
    assert len(result) == 1
    assert result.n_controls.iloc[0] == 2
    assert result.control_mean_score.iloc[0] == pytest.approx(0.5)
    assert result.contrast.iloc[0] == pytest.approx(1.5)
    assert result.control_proteins.iloc[0] == "C1;C2"


def test_case_without_a_scored_control_is_dropped_not_zeroed():
    """A missing score must remove the contrast, never contribute as zero."""
    site = site_frame([("P1", 10, 2.0, "NXT", "cl1"), ("P2", 20, 3.0, "NXT", "cl2"),
                       ("C1", 1, 0.5, "NXT", None)])
    pairs = pair_frame([("P1", 10, "C1", 1), ("P2", 20, "C_absent", 9)])
    result = build_contrasts(pairs, site)
    assert list(result.case_accession) == ["P1"]


def test_unscored_case_is_dropped():
    site = site_frame([("C1", 1, 0.5, "NXT", None)])
    result = build_contrasts(pair_frame([("P_absent", 10, "C1", 1)]), site)
    assert result.empty


def test_missing_ortholog_cluster_falls_back_to_the_site_itself():
    """A site with no cluster must not be pooled with every other such site."""
    site = site_frame([("P1", 10, 1.0, "NXT", None), ("P2", 20, 1.0, "NXT", None),
                       ("C1", 1, 0.0, "NXT", None), ("C2", 2, 0.0, "NXT", None)])
    result = build_contrasts(pair_frame([("P1", 10, "C1", 1), ("P2", 20, "C2", 2)]), site)
    assert result.ortholog_cluster.nunique() == 2
    assert result.resample_unit.nunique() == 2


# --------------------------------------------------------------------------
# resampling units
# --------------------------------------------------------------------------

def base_contrasts(rows):
    """rows: (case, cluster, control_proteins)"""
    return pd.DataFrame(
        [{"case_accession": c, "case_position": i, "subtype": "NXT",
          "case_score": 1.0, "control_mean_score": 0.0, "n_controls": 1,
          "control_proteins": p, "ortholog_cluster": cl, "contrast": 1.0}
         for i, (c, cl, p) in enumerate(rows)]
    )


def test_shared_control_protein_joins_two_units():
    units = assign_resample_units(
        base_contrasts([("P1", "clA", "Cx"), ("P2", "clB", "Cx")]))
    assert units.resample_unit.nunique() == 1


def test_shared_ortholog_cluster_joins_two_units():
    units = assign_resample_units(
        base_contrasts([("P1", "clA", "Cx"), ("P2", "clA", "Cy")]))
    assert units.resample_unit.nunique() == 1


def test_dependence_chains_transitively():
    """A and B share a control; B and C share a cluster; all three move together."""
    units = assign_resample_units(base_contrasts([
        ("P1", "clA", "Cx"), ("P2", "clB", "Cx"), ("P3", "clB", "Cz")]))
    assert units.resample_unit.nunique() == 1


def test_independent_contrasts_stay_separate():
    units = assign_resample_units(base_contrasts([
        ("P1", "clA", "Cx"), ("P2", "clB", "Cy"), ("P3", "clC", "Cz")]))
    assert units.resample_unit.nunique() == 3


def test_a_case_matched_to_several_controls_links_all_of_them():
    units = assign_resample_units(base_contrasts([
        ("P1", "clA", "Cx;Cy"), ("P2", "clB", "Cy")]))
    assert units.resample_unit.nunique() == 1


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------

def spread_contrasts(values, units):
    return pd.DataFrame({
        "contrast": values,
        "resample_unit": units,
    })


def test_bootstrap_is_deterministic_for_a_given_seed():
    frame = spread_contrasts([1.0, 2.0, -1.0, 0.5], ["u1", "u2", "u3", "u4"])
    a = cluster_bootstrap(frame, 500, seed=7)
    b = cluster_bootstrap(frame, 500, seed=7)
    np.testing.assert_array_equal(a, b)


def test_bootstrap_differs_across_seeds():
    frame = spread_contrasts([1.0, 2.0, -1.0, 0.5], ["u1", "u2", "u3", "u4"])
    assert not np.array_equal(cluster_bootstrap(frame, 500, seed=7),
                              cluster_bootstrap(frame, 500, seed=8))


def test_bootstrap_resamples_units_not_rows():
    """With every contrast in one unit there is nothing to vary, so a
    row-level bootstrap would produce spread and a unit-level one cannot."""
    frame = spread_contrasts([1.0, 2.0, -1.0, 0.5], ["u1"] * 4)
    draws = cluster_bootstrap(frame, 200, seed=1)
    assert np.allclose(draws, frame.contrast.mean())


def test_bootstrap_mean_is_centred_near_the_observed_mean():
    rng = np.random.default_rng(0)
    values = rng.normal(0.4, 1.0, 40)
    frame = spread_contrasts(values, [f"u{i}" for i in range(40)])
    draws = cluster_bootstrap(frame, 2000, seed=3)
    assert draws.mean() == pytest.approx(values.mean(), abs=0.05)


# --------------------------------------------------------------------------
# verdicts
# --------------------------------------------------------------------------

@pytest.mark.parametrize("low,high,expected", [
    (-0.10, 0.10, "equivalent within the margin"),
    (0.50, 1.00, "difference beyond the margin"),
    (-1.00, -0.50, "difference beyond the margin"),
    (0.05, 1.20, "directional, magnitude undetermined"),
    (-0.30, 1.10, "inconclusive"),
])
def test_classify(low, high, expected):
    assert classify(low, high, margin=0.2663) == expected
