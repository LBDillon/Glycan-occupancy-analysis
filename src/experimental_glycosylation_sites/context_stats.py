"""Effect sizes and uncertainty for the context comparison.

Two decisions are load-bearing here and both are pre-specified.

**Sites are not independent.** One protein contributes up to 7 occupied sequons
and up to 19 secretory ones. Treating those as separate observations narrows
every interval by a factor that reflects how many sequons a protein happens to
have, not how certain the comparison is. Every interval is therefore a cluster
bootstrap over proteins.

**The standardised mean difference describes separation, not probability.** It
says how far apart two distributions sit in pooled standard deviations. It is
not a probability that a site is occupied, and must never be reported as one.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def smd(occupied: pd.Series, comparison: pd.Series) -> float:
    """(mean_occ - mean_cmp) / sqrt((var_occ + var_cmp) / 2), NaN without spread."""
    a = pd.Series(occupied).dropna().astype(float)
    b = pd.Series(comparison).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = (a.var(ddof=1) + b.var(ddof=1)) / 2.0
    if not np.isfinite(pooled) or pooled <= 0:
        return float("nan")
    return float((a.mean() - b.mean()) / np.sqrt(pooled))


def _mean_difference(a: pd.Series, b: pd.Series) -> float:
    """Difference in means. On a 0/1 indicator this is the proportion difference."""
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if not len(a) or not len(b):
        return float("nan")
    return float(a.mean() - b.mean())


STATISTICS = {"smd": smd, "mean_difference": _mean_difference}


def cluster_bootstrap_difference(frame: pd.DataFrame, feature: str, group_column: str,
                                 occupied_label: str, comparison_label: str,
                                 cluster_column: str = "accession",
                                 n_boot: int = 2000, seed: int = 0,
                                 alpha: float = 0.05,
                                 statistic: str = "smd") -> dict:
    """SMD with a percentile interval from resampling clusters, not rows.

    Clusters are resampled with replacement within each arm separately, so the
    arm sizes stay comparable across replicates. The p-value is the two-sided
    bootstrap value 2 * min(P(SMD* <= 0), P(SMD* >= 0)), which is the fraction
    of replicates falling on the wrong side of no effect.
    """
    data = frame[[feature, group_column, cluster_column]].copy()
    data = data[data[feature].notna()]
    occ = data[data[group_column] == occupied_label]
    cmp_ = data[data[group_column] == comparison_label]

    estimator = STATISTICS[statistic]
    result = {
        "feature": feature,
        "statistic": statistic,
        "comparison": comparison_label,
        "n_occupied": int(len(occ)),
        "n_comparison": int(len(cmp_)),
        "proteins_occupied": int(occ[cluster_column].nunique()),
        "proteins_comparison": int(cmp_[cluster_column].nunique()),
        "mean_occupied": float(occ[feature].mean()) if len(occ) else float("nan"),
        "mean_comparison": float(cmp_[feature].mean()) if len(cmp_) else float("nan"),
        "estimate": estimator(occ[feature], cmp_[feature]),
        "ci_low": float("nan"), "ci_high": float("nan"), "p": float("nan"),
        "n_boot": int(n_boot),
    }
    if not np.isfinite(result["estimate"]):
        return result

    rng = np.random.default_rng(seed)
    # Values grouped by cluster once, so a replicate is an index draw rather
    # than a repeated groupby over the frame.
    def by_cluster(part):
        groups = part.groupby(cluster_column)[feature]
        return [g.to_numpy(dtype=float) for _, g in groups]

    occ_clusters, cmp_clusters = by_cluster(occ), by_cluster(cmp_)
    if len(occ_clusters) < 2 or len(cmp_clusters) < 2:
        return result

    draws = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        a = np.concatenate([occ_clusters[j] for j in
                            rng.integers(0, len(occ_clusters), len(occ_clusters))])
        b = np.concatenate([cmp_clusters[j] for j in
                            rng.integers(0, len(cmp_clusters), len(cmp_clusters))])
        draws[i] = estimator(pd.Series(a), pd.Series(b))

    finite = draws[np.isfinite(draws)]
    if len(finite) < 2:
        return result
    result["ci_low"] = float(np.percentile(finite, 100 * alpha / 2))
    result["ci_high"] = float(np.percentile(finite, 100 * (1 - alpha / 2)))
    below = float(np.mean(finite <= 0))
    above = float(np.mean(finite >= 0))
    # Bounded away from zero: with B replicates nothing is resolvable past 1/B.
    result["p"] = float(min(1.0, max(2 * min(below, above), 1.0 / len(finite))))
    return result


def benjamini_hochberg(pvalues) -> "list[float]":
    """BH-adjusted values, preserving input order; NaN passes through as NaN."""
    values = np.asarray(pvalues, dtype=float)
    out = np.full(values.shape, np.nan)
    present = np.where(np.isfinite(values))[0]
    if len(present) == 0:
        return out.tolist()

    ordered = present[np.argsort(values[present])]
    m = len(ordered)
    running = 1.0
    for rank in range(m - 1, -1, -1):
        candidate = values[ordered[rank]] * m / (rank + 1)
        running = min(running, candidate)
        out[ordered[rank]] = running
    return out.tolist()


def ramachandran_region(phi, psi) -> "str | None":
    """Coarse backbone region, or None if either angle is missing.

    Deliberately coarse -- enough to separate helical from extended from
    left-handed, which is what the panel needs. Dihedrals are circular, so they
    cannot enter the comparison as means; this is how backbone geometry is
    tested instead.

    Left-handed alpha matters here specifically: asparagine and glycine occupy
    it far more often than other residues, so it is not a rounding category.
    """
    if phi is None or psi is None:
        return None
    try:
        phi, psi = float(phi), float(psi)
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(phi) and np.isfinite(psi)):
        return None
    if phi > 0:
        return "alpha_L" if -60.0 <= psi <= 90.0 else "other"
    if -120.0 <= psi <= 50.0:
        return "alpha_R"
    return "beta"


def ramachandran_region_series(phi: pd.Series, psi: pd.Series) -> pd.Series:
    """`ramachandran_region` over two columns, preserving the index."""
    return pd.Series([ramachandran_region(a, b) for a, b in zip(phi, psi)],
                     index=phi.index, dtype="object")


def cluster_bootstrap_smd(frame: pd.DataFrame, feature: str, group_column: str,
                          occupied_label: str, comparison_label: str,
                          cluster_column: str = "accession",
                          n_boot: int = 2000, seed: int = 0,
                          alpha: float = 0.05) -> dict:
    """`cluster_bootstrap_difference` with the standardised mean difference."""
    out = cluster_bootstrap_difference(frame, feature, group_column, occupied_label,
                                       comparison_label, cluster_column, n_boot,
                                       seed, alpha, statistic="smd")
    out["smd"] = out["estimate"]
    return out
