"""Distance from the natural occupied distribution.

The reference must never contain the site being scored, or a wild-type site is
partly compared against itself and looks more natural than it is. Held-out means
held out by protein, not by row: a protein with several sequons would otherwise
leak into its own reference.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from glyco_context.context_distance import context_distance, reference_moments

FEATURES = ["a", "b"]
REF = pd.DataFrame([
    {"accession": "P1", "a": 1.0, "b": 10.0},
    {"accession": "P1", "a": 1.0, "b": 10.0},
    {"accession": "P2", "a": 3.0, "b": 20.0},
    {"accession": "P3", "a": 5.0, "b": 30.0},
])


def test_reference_excludes_the_protein_being_scored():
    mu, sigma = reference_moments(REF, FEATURES, exclude_accession="P1")
    assert mu["a"] == pytest.approx(4.0)      # only P2 and P3 remain


def test_reference_uses_everything_when_nothing_is_excluded():
    mu, _ = reference_moments(REF, FEATURES, exclude_accession=None)
    assert mu["a"] == pytest.approx(2.5)


def test_distance_is_zero_at_the_reference_mean():
    mu, sigma = reference_moments(REF, FEATURES, exclude_accession=None)
    row = {"a": mu["a"], "b": mu["b"]}
    assert context_distance(row, mu, sigma, FEATURES) == pytest.approx(0.0)


def test_distance_grows_with_departure_from_the_reference():
    mu, sigma = reference_moments(REF, FEATURES, exclude_accession=None)
    near = context_distance({"a": mu["a"], "b": mu["b"] + sigma["b"]}, mu, sigma, FEATURES)
    far = context_distance({"a": mu["a"], "b": mu["b"] + 4 * sigma["b"]}, mu, sigma, FEATURES)
    assert far > near


def test_features_without_spread_are_skipped_not_infinite():
    """A constant reference feature cannot standardise anything."""
    flat = pd.DataFrame([{"accession": "P1", "a": 1.0, "b": 7.0},
                         {"accession": "P2", "a": 3.0, "b": 7.0}])
    mu, sigma = reference_moments(flat, FEATURES, exclude_accession=None)
    d = context_distance({"a": 3.0, "b": 99.0}, mu, sigma, FEATURES)
    assert np.isfinite(d)


def test_missing_values_are_skipped():
    mu, sigma = reference_moments(REF, FEATURES, exclude_accession=None)
    d = context_distance({"a": mu["a"], "b": None}, mu, sigma, FEATURES)
    assert d == pytest.approx(0.0)


def test_distance_is_nan_when_nothing_is_measurable():
    mu, sigma = reference_moments(REF, FEATURES, exclude_accession=None)
    assert np.isnan(context_distance({"a": None, "b": None}, mu, sigma, FEATURES))
