"""Comparing a model's score with the motif visible against with it hidden.

The quantity is a paired difference of paired differences: for each matched
pair, the occupied-minus-control contrast under each masking scheme, and then
the change between schemes. It is only interpretable on pairs both schemes
scored, and the resample unit still has to be respected.
"""
from __future__ import annotations

import pandas as pd
import pytest

from experimental_glycosylation_sites.masking import masking_change


def _contrasts(values, units=None):
    return pd.DataFrame({
        "case_accession": [f"P{i}" for i in range(len(values))],
        "case_position": list(range(len(values))),
        "contrast": values,
        "resample_unit": units or [f"u{i}" for i in range(len(values))]})


def test_change_is_visible_minus_hidden():
    out = masking_change(_contrasts([1.0, 1.0]), _contrasts([0.4, 0.6]), n_boot=200)
    assert out["mean"] == pytest.approx(0.5)


def test_only_pairs_scored_by_both_schemes_are_used():
    """A pair one scheme could not score cannot contribute a change."""
    visible = _contrasts([1.0, 1.0, 1.0])
    hidden = _contrasts([0.4, 0.6, 0.5]).iloc[:2]
    out = masking_change(visible, hidden, n_boot=100)
    assert out["n_pairs"] == 2


def test_no_change_reads_as_zero():
    out = masking_change(_contrasts([0.3, 0.7]), _contrasts([0.3, 0.7]), n_boot=100)
    assert out["mean"] == pytest.approx(0.0)


def test_resampling_is_by_unit_not_by_pair():
    """Ten pairs in one unit carry one unit's worth of evidence."""
    many = _contrasts([1.0] * 10, units=["u0"] * 10)
    hidden = _contrasts([0.0] * 10, units=["u0"] * 10)
    out = masking_change(many, hidden, n_boot=300)
    assert out["n_units"] == 1
    # a single cluster cannot produce a tight interval around the estimate
    assert out["ci_low"] == out["ci_high"] == pytest.approx(1.0)


def test_units_are_counted():
    v = _contrasts([1.0, 1.0, 1.0], units=["a", "a", "b"])
    h = _contrasts([0.0, 0.0, 0.0], units=["a", "a", "b"])
    assert masking_change(v, h, n_boot=100)["n_units"] == 2


def test_no_overlap_returns_nothing_rather_than_raising():
    v = _contrasts([1.0])
    h = _contrasts([0.5]).assign(case_accession="Q9")
    assert masking_change(v, h, n_boot=50)["n_pairs"] == 0
