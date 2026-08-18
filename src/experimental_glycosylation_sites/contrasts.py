"""Turning matched pairs into contrasts, and contrasts into an interval.

Extracted so the primary analysis and the matching-sensitivity sweep cannot
drift apart. A sensitivity analysis that quietly computes its interval a
different way from the result it is testing is worse than no sensitivity
analysis, because it looks like corroboration.

The unit of analysis is one occupied site:

    site_contrast = occupied_score - mean(scores of its matched controls)

The resampling unit is coarser. Two dependencies run through these contrasts:
occupied sites in the same ortholog cluster are near copies of one another, and
one control protein can be matched to several occupied cases. Resampling on
either alone leaves the other unaccounted for, so contrasts are grouped into the
connected components of the bipartite graph joining occupied clusters to shared
control proteins, and components are resampled whole.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def build_contrasts(pairs: pd.DataFrame, site: pd.DataFrame) -> pd.DataFrame:
    """One row per occupied case, with its control mean and provenance."""
    lookup = site.set_index(["accession", "position"])

    def score_of(accession, position):
        try:
            return float(lookup.loc[(str(accession), int(position))].conditional_sequon_score)
        except KeyError:
            return np.nan

    rows = []
    for (accession, position), group in pairs.groupby(["case_accession", "case_position"]):
        case_score = score_of(accession, position)
        controls = [(str(r.control_accession), int(r.control_position),
                     score_of(r.control_accession, r.control_position))
                    for r in group.itertuples(index=False)]
        controls = [c for c in controls if np.isfinite(c[2])]
        if not np.isfinite(case_score) or not controls:
            continue
        record = lookup.loc[(str(accession), int(position))]
        clusters = record.get("ortholog_clusters", pd.NA)
        rows.append({
            "case_accession": str(accession), "case_position": int(position),
            "subtype": record.subtype,
            "case_score": case_score,
            "control_mean_score": float(np.mean([c[2] for c in controls])),
            "n_controls": len(controls),
            "control_proteins": ";".join(sorted({c[0] for c in controls})),
            "ortholog_cluster": (str(clusters).split(";")[0]
                                 if pd.notna(clusters) and str(clusters)
                                 else f"solo:{accession}"),
        })

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["contrast"] = frame.case_score - frame.control_mean_score
    return assign_resample_units(frame)


def assign_resample_units(contrasts: pd.DataFrame) -> pd.DataFrame:
    """Group contrasts that share an ortholog cluster or a control protein."""
    parent: dict[str, str] = {}

    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a, b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_a] = root_b

    for row in contrasts.itertuples(index=False):
        anchor = f"cluster:{row.ortholog_cluster}"
        find(anchor)
        for protein in row.control_proteins.split(";"):
            union(anchor, f"control:{protein}")

    contrasts = contrasts.copy()
    contrasts["resample_unit"] = [find(f"cluster:{r.ortholog_cluster}")
                                  for r in contrasts.itertuples(index=False)]
    return contrasts


def cluster_bootstrap(contrasts: pd.DataFrame, n_boot: int, seed: int) -> np.ndarray:
    """Bootstrap means, resampling whole components with replacement."""
    units = contrasts.resample_unit.unique()
    by_unit = [contrasts.contrast[contrasts.resample_unit == u].to_numpy() for u in units]
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        picked = rng.integers(0, len(by_unit), len(by_unit))
        draws[i] = np.concatenate([by_unit[j] for j in picked]).mean()
    return draws


def classify(low: float, high: float, margin: float) -> str:
    """Four outcomes, not two.

    A confidence interval that excludes zero but reaches past the margin has
    established a direction and not a magnitude. Calling that "inconclusive"
    discards the informative half; calling it "a difference" overstates the
    other.
    """
    if low >= -margin and high <= margin:
        return "equivalent within the margin"
    if low > margin or high < -margin:
        return "difference beyond the margin"
    if not (low <= 0 <= high):
        return "directional, magnitude undetermined"
    return "inconclusive"
