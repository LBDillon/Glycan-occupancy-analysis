"""Aggregating design-site rows into the reported retention statistics.

Kept apart from the stage that generates the designs so the numbers can be
recomputed, and tested, without a ProteinMPNN run.

Two things this exists to get right. Thirty-two designs of one chain are
replicates of a single draw, so they are averaged within site before anything
else; a mean over raw rows would weight a chain by how many sequons it carries.
And sites within a protein are not independent, so intervals resample proteins.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

QUANTITIES = {
    "sequon_exact": "sequon retained, exact triplet",
    "sequon_pattern": "sequon retained, N-X-S/T pattern",
    "control_triplet_exact": "control triplet retained, exact",
    "background_mutation_rate": "background mutation rate",
}


def _bootstrap(values: pd.Series, clusters: pd.Series, n_boot: int, seed: int = 0):
    values = values.dropna()
    if len(values) < 2:
        return {"mean": float(values.mean()) if len(values) else float("nan"),
                "ci_low": float("nan"), "ci_high": float("nan"), "n": int(len(values))}
    groups = [g.to_numpy(float) for _, g in values.groupby(clusters.loc[values.index])]
    rng = np.random.default_rng(seed)
    draws = np.array([np.concatenate([groups[i] for i in
                                      rng.integers(0, len(groups), len(groups))]).mean()
                      for _ in range(n_boot)])
    return {"mean": float(values.mean()),
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)),
            "n": int(len(values))}


def summarise(rows: pd.DataFrame, n_boot: int = 4000, seed: int = 0) -> dict:
    """Per-site means, then protein-resampled intervals, for each quantity."""
    columns = [c for c in QUANTITIES if c in rows.columns]
    if not len(rows) or not columns:
        return {"sites": 0, "proteins": 0, "design_rows": int(len(rows))}

    site = rows.groupby(["accession", "position"])[columns].mean().reset_index()
    proteins = site.accession
    out = {"sites": int(len(site)), "proteins": int(site.accession.nunique()),
           "design_rows": int(len(rows))}
    for column in columns:
        out[column] = _bootstrap(site[column], proteins, n_boot, seed)

    if {"control_triplet_exact", "sequon_exact"} <= set(columns):
        excess = site.control_triplet_exact - site.sequon_exact
        out["control_minus_sequon"] = _bootstrap(excess, proteins, n_boot, seed)
    return out
