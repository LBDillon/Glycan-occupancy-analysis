"""Significance testing for both outcomes, across all four control sets.

Three deliberate choices.

**A cluster-level permutation test, not Wilcoxon.** Wilcoxon and the sign test
assume the pairs are independent. They are not: occupied sites in one ortholog
cluster are near-copies, and one control protein can serve several occupied
cases. The test here flips the sign of every contrast within a whole resampling
unit at once — the same connected components used by the bootstrap — so the null
respects the dependency. Under the null of no effect, the sign of a unit's
contrasts is exchangeable; under clustering, individual contrasts are not.

**Correction across the whole family.** Eight tests are reported (four control
sets, two outcomes). Holm controls the family-wise error rate and Benjamini–
Hochberg the false discovery rate; both are shown, because the honest reading of
eight related tests is not the same as eight independent ones.

**p-values are reported beside effect sizes and intervals, never instead.** For
the conditional score the pre-specified inference is an equivalence assessment
against a ±0.2 SD margin, which a p-value cannot deliver. A small p there means
"not identical", not "meaningfully different".
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

N_PERM, SEED = 20000, 20260819
MARGIN_SD = 0.2

SCORE = [("internal control", "results/analysis/contrasts_optimal.csv"),
         ("eukaryotic secretory", "results/analysis/contrasts_secretory.csv"),
         ("bacterial", "results/analysis/contrasts_bacterial.csv"),
         ("cytosolic", "results/analysis/contrasts_cytosolic.csv")]
RETENTION = [("internal control", "results/analysis/retention_paired_internal.csv"),
             ("eukaryotic secretory", "results/analysis/retention_paired_eukaryotic.csv"),
             ("bacterial", "results/analysis/retention_paired_bacterial.csv"),
             ("cytosolic", "results/analysis/retention_paired_cytosolic.csv")]

dataset = pd.read_csv("results/scores/scores_dataset.csv", low_memory=False)
REFERENCE_SD = float(dataset.conditional_sequon_score.std(ddof=1))


def cluster_permutation(frame, n_perm=N_PERM, seed=SEED):
    """Two-sided p by flipping the sign of whole resampling units."""
    rng = np.random.default_rng(seed)
    units = [g.contrast.to_numpy() for _, g in frame.groupby("resample_unit")]
    sizes = np.array([len(u) for u in units])
    values = np.concatenate(units)
    observed = values.mean()
    unit_index = np.repeat(np.arange(len(units)), sizes)
    extreme = 0
    for _ in range(n_perm):
        flips = rng.choice([-1.0, 1.0], size=len(units))
        if abs((values * flips[unit_index]).mean()) >= abs(observed) - 1e-15:
            extreme += 1
    # add-one correction: a permutation p is never exactly zero
    return observed, (extreme + 1) / (n_perm + 1)


def holm(pvalues):
    order = np.argsort(pvalues)
    n = len(pvalues)
    adjusted = np.empty(n)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (n - rank) * pvalues[idx])
        adjusted[idx] = min(1.0, running)
    return adjusted


def benjamini_hochberg(pvalues):
    order = np.argsort(pvalues)
    n = len(pvalues)
    adjusted = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        idx = order[rank]
        running = min(running, n * pvalues[idx] / (rank + 1))
        adjusted[idx] = min(1.0, running)
    return adjusted


rows = []
for outcome, table in (("conditional score", SCORE), ("design retention", RETENTION)):
    for label, path in table:
        if not Path(path).exists():
            continue
        frame = pd.read_csv(path)
        observed, p_perm = cluster_permutation(frame)
        informative = int((frame.contrast != 0).sum())
        rows.append({
            "outcome": outcome, "comparison": label,
            "n_pairs": int(len(frame)),
            "n_units": int(frame.resample_unit.nunique()),
            "n_informative": informative,
            "effect": observed,
            "effect_sd_units": observed / REFERENCE_SD if outcome == "conditional score" else np.nan,
            "p_permutation": p_perm,
            "p_wilcoxon": (float(st.wilcoxon(frame.contrast).pvalue)
                           if informative else np.nan),
        })

results = pd.DataFrame(rows)
results["p_holm"] = holm(results.p_permutation.to_numpy())
results["p_bh"] = benjamini_hochberg(results.p_permutation.to_numpy())
results.to_csv("results/analysis/significance.csv", index=False)

print(f"cluster-level sign-flip permutation, {N_PERM:,} draws, {len(results)} tests\n")
for outcome in results.outcome.unique():
    part = results[results.outcome == outcome]
    print(f"=== {outcome} ===")
    print(f"{'comparison':22s} {'pairs':>6s} {'units':>6s} {'effect':>9s} "
          f"{'p (perm)':>9s} {'p (Wilcox)':>11s} {'p Holm':>8s} {'p BH':>7s}")
    for r in part.itertuples(index=False):
        effect = (f"{r.effect_sd_units:+.3f} SD" if outcome == "conditional score"
                  else f"{r.effect:+.4f}")
        print(f"{r.comparison:22s} {r.n_pairs:>6d} {r.n_units:>6d} {effect:>9s} "
              f"{r.p_permutation:>9.4f} {r.p_wilcoxon:>11.4f} {r.p_holm:>8.4f} {r.p_bh:>7.4f}")
    print()

print("equivalence assessment (conditional score only; the pre-specified inference):")
for label, path in SCORE:
    key = {"internal control": "optimal", "eukaryotic secretory": "secretory",
           "bacterial": "bacterial", "cytosolic": "cytosolic"}[label]
    js = Path(f"results/analysis/analysis_{key}.json")
    if js.exists():
        d = json.loads(js.read_text())
        lo, hi = d["ci95_sd"]
        print(f"  {label:22s} [{lo:+.3f}, {hi:+.3f}] SD vs margin +/-{MARGIN_SD}  ->  {d['verdict']}")

Path("results/analysis/significance.json").write_text(
    results.to_json(orient="records", indent=2))
print("\nwrote results/analysis/significance.csv and results/analysis/significance.json")
