"""Aggregating design-site rows into the reported retention statistics.

The quantities are simple; the ways of getting them wrong are not. Thirty-two
designs of one chain are replicates of one draw, so a mean over rows weights
chains by how many sites they carry. Intervals must resample proteins, because
sites in one protein are not independent observations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from glyco_context.retention_rate import summarise


def _rows(spec):
    """spec: (accession, position, sequon_exact, control, background) per design."""
    return pd.DataFrame([
        {"accession": a, "position": p, "structure_pdb_id": a,
         "sequon_exact": se, "sequon_pattern": se,
         "control_triplet_exact": c, "background_mutation_rate": b}
        for a, p, se, c, b in spec])


def test_designs_are_averaged_within_site_before_anything_else():
    """Two designs at one site, one retaining: the site's rate is 0.5, not two
    observations of 1 and 0."""
    table = _rows([("P1", 1, True, 0.1, 0.5), ("P1", 1, False, 0.1, 0.5)])
    out = summarise(table, n_boot=50)
    assert out["sites"] == 1
    assert out["sequon_exact"]["mean"] == pytest.approx(0.5)


def test_a_protein_with_many_sites_does_not_dominate():
    """Ten sites in one protein and one in another is two proteins' evidence."""
    many = [("P1", i, True, 0.1, 0.5) for i in range(10)]
    one = [("P2", 1, False, 0.1, 0.5)]
    out = summarise(_rows(many + one), n_boot=200)
    # site-level mean would be 10/11 = 0.91; protein-level resampling must not
    # report near-certainty from what is really two clusters
    assert out["sequon_exact"]["ci_low"] < 0.9


def test_proteins_and_sites_are_counted_separately():
    out = summarise(_rows([("P1", 1, True, 0.1, 0.5), ("P1", 2, False, 0.1, 0.5),
                           ("P2", 1, True, 0.1, 0.5)]), n_boot=50)
    assert out["sites"] == 3
    assert out["proteins"] == 2


def test_excess_is_control_minus_sequon():
    out = summarise(_rows([("P1", 1, False, 0.2, 0.5), ("P2", 1, False, 0.2, 0.5)]),
                    n_boot=100)
    assert out["control_minus_sequon"]["mean"] == pytest.approx(0.2)


def test_interval_brackets_the_estimate():
    rng = np.random.default_rng(0)
    spec = [(f"P{i}", 1, bool(rng.integers(0, 2)), 0.1, 0.5) for i in range(40)]
    out = summarise(_rows(spec), n_boot=500)
    e = out["sequon_exact"]
    assert e["ci_low"] <= e["mean"] <= e["ci_high"]


def test_empty_input_does_not_raise():
    out = summarise(_rows([]), n_boot=10)
    assert out["sites"] == 0
