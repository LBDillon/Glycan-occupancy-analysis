"""Effect sizes and uncertainty for the context comparison.

The property being bought here is the cluster bootstrap. One protein contributes
up to 19 sequons, and treating those as 19 independent observations narrows
every interval by a factor that has nothing to do with the biology. The decisive
test is that perfectly correlated sites within a protein widen the interval
rather than narrowing it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from glyco_context.context_stats import (benjamini_hochberg,
                                                            cluster_bootstrap_smd,
                                                            smd)


def test_smd_is_zero_for_identical_groups():
    x = pd.Series([1.0, 2.0, 3.0, 4.0])
    assert smd(x, x) == pytest.approx(0.0)


def test_smd_scales_by_pooled_spread():
    """Means differ by 2, each group has unit variance -> SMD of 2."""
    a = pd.Series([1.0, 3.0])      # mean 2, var 2
    b = pd.Series([-1.0, 1.0])     # mean 0, var 2
    assert smd(a, b) == pytest.approx(2.0 / np.sqrt(2.0))


def test_smd_is_nan_without_spread():
    a = pd.Series([1.0, 1.0])
    assert np.isnan(smd(a, pd.Series([2.0, 2.0])))


def test_smd_ignores_missing_values():
    a = pd.Series([1.0, 3.0, np.nan])
    b = pd.Series([-1.0, 1.0])
    assert smd(a, b) == pytest.approx(2.0 / np.sqrt(2.0))


def _clustered(sites_per_protein, n_proteins, seed=0):
    """Every site in a protein carries that protein's value exactly."""
    rng = np.random.default_rng(seed)
    rows = []
    for group in ("occupied", "control"):
        shift = 1.0 if group == "occupied" else 0.0
        for p in range(n_proteins):
            value = rng.normal(shift, 1.0)
            for _ in range(sites_per_protein):
                rows.append({"accession": f"{group}{p}", "group": group, "x": value})
    return pd.DataFrame(rows)


def test_cluster_bootstrap_does_not_gain_precision_from_repeated_sites():
    """Twenty copies of one protein's value is still one observation.

    Row resampling would treat it as twenty and shrink the interval; resampling
    proteins must not.
    """
    one = cluster_bootstrap_smd(_clustered(1, 30), "x", "group", "occupied",
                                "control", "accession", n_boot=400, seed=1)
    twenty = cluster_bootstrap_smd(_clustered(20, 30), "x", "group", "occupied",
                                   "control", "accession", n_boot=400, seed=1)
    one_width = one["ci_high"] - one["ci_low"]
    twenty_width = twenty["ci_high"] - twenty["ci_low"]
    assert twenty_width > one_width * 0.6, (
        "duplicating sites within a protein must not sharply narrow the interval")


def test_cluster_bootstrap_is_reproducible_under_a_seed():
    frame = _clustered(3, 20)
    a = cluster_bootstrap_smd(frame, "x", "group", "occupied", "control",
                              "accession", n_boot=200, seed=7)
    b = cluster_bootstrap_smd(frame, "x", "group", "occupied", "control",
                              "accession", n_boot=200, seed=7)
    assert a == b


def test_cluster_bootstrap_reports_group_sizes_and_coverage():
    frame = _clustered(2, 10)
    frame.loc[frame.index[:4], "x"] = np.nan
    result = cluster_bootstrap_smd(frame, "x", "group", "occupied", "control",
                                   "accession", n_boot=100, seed=3)
    assert result["n_occupied"] + result["n_comparison"] == int(frame.x.notna().sum())
    assert result["proteins_occupied"] <= 10


def test_benjamini_hochberg_matches_a_worked_example():
    p = [0.01, 0.02, 0.03, 0.04, 0.05]
    q = benjamini_hochberg(p)
    assert q[0] == pytest.approx(0.05)
    assert all(q[i] <= q[i + 1] + 1e-12 for i in range(len(q) - 1))


def test_benjamini_hochberg_leaves_a_single_p_value_alone():
    assert benjamini_hochberg([0.031])[0] == pytest.approx(0.031)


def test_benjamini_hochberg_handles_missing_values():
    q = benjamini_hochberg([0.01, float("nan"), 0.05])
    assert np.isnan(q[1])
    assert not np.isnan(q[0])


def test_ramachandran_regions():
    """Coarse regions only -- enough to say 'helical' from 'extended'."""
    from glyco_context.context_stats import ramachandran_region
    assert ramachandran_region(-63, -43) == "alpha_R"      # right-handed helix
    assert ramachandran_region(-120, 130) == "beta"        # extended
    assert ramachandran_region(57, 40) == "alpha_L"        # left-handed, often Gly/Asn
    assert ramachandran_region(None, -43) is None
    assert ramachandran_region(-63, None) is None


def test_ramachandran_region_series_handles_a_column():
    from glyco_context.context_stats import ramachandran_region_series
    out = ramachandran_region_series(pd.Series([-63.0, -120.0, np.nan]),
                                     pd.Series([-43.0, 130.0, 20.0]))
    assert list(out[:2]) == ["alpha_R", "beta"]
    assert pd.isna(out.iloc[2])


def test_proportion_difference_on_a_binary_indicator():
    """For a 0/1 column the mean difference IS the difference in proportion."""
    from glyco_context.context_stats import cluster_bootstrap_difference
    rows = ([{"accession": f"o{i}", "group": "occupied", "x": 1.0} for i in range(20)]
            + [{"accession": f"o{i}", "group": "occupied", "x": 0.0} for i in range(20, 30)]
            + [{"accession": f"c{i}", "group": "control", "x": 1.0} for i in range(10)]
            + [{"accession": f"c{i}", "group": "control", "x": 0.0} for i in range(10, 40)])
    out = cluster_bootstrap_difference(pd.DataFrame(rows), "x", "group", "occupied",
                                       "control", "accession", n_boot=200, seed=2,
                                       statistic="mean_difference")
    # 20/30 occupied vs 10/40 control -> 0.6667 - 0.25
    assert out["estimate"] == pytest.approx(20/30 - 0.25, abs=1e-9)
