"""How much of the primary result is the matching algorithm's arbitrary choices?

Greedy nearest-neighbour matching walks the cases in a seeded random order. With
28 controls for 314 occupied sites that order matters: an early case can take the
only admissible partner for a later one. The reported seed-0 result was one draw
from this distribution, presented as a value.

This runs 200 of them, and reports the distribution of the mean contrast and of
whether the interval excludes zero. The deterministic optimal matching is shown
alongside as the primary reference.

Uses the same contrast construction and bootstrap as the primary analysis, from
`experimental_glycosylation_sites.contrasts`, so the comparison is like for like.
Bootstrap draws are reduced to 2,000 per seed because 200 intervals are being
summarised rather than any one of them quoted.
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.contrasts import (
    build_contrasts, classify, cluster_bootstrap)
from experimental_glycosylation_sites.matching import (
    DEFAULT_CALIPER, DEFAULT_EXACT, MATCH_FEATURES,
    match_controls, match_controls_optimal)

N_SEEDS, N_BOOT, BOOT_SEED, MARGIN_SD = 200, 2000, 20260818, 0.2
KEY = ["accession", "position"]

manifest = pd.read_csv("results/manifests/candidate_manifest_dataset.csv", low_memory=False)
manifest = manifest[manifest.scoreable.astype(bool)].copy()
scores = pd.read_csv("results/scores/scores_dataset.csv", low_memory=False)
for frame in (manifest, scores):
    frame["accession"] = frame.accession.astype(str)
    frame["position"] = frame.position.astype(int)
site = manifest.merge(scores[KEY + ["conditional_sequon_score"]], on=KEY, how="inner")

reference_sd = float(site.conditional_sequon_score.std(ddof=1))
margin_raw = MARGIN_SD * reference_sd
cases = site[site.occupancy_status == "occupied_supported"].copy()
controls = site[site.occupancy_status == "observed_unmodified"].copy()
print(f"reference SD {reference_sd:.4f}   margin +/-{margin_raw:.4f}")
print(f"{len(cases)} occupied, {len(controls)} internal controls, {N_SEEDS} seeds\n")


def evaluate(pairs, boot_seed):
    contrasts = build_contrasts(pairs, site)
    if contrasts.empty:
        return None
    draws = cluster_bootstrap(contrasts, N_BOOT, boot_seed)
    low, high = np.percentile(draws, [2.5, 97.5])
    return {"n_pairs": int(len(contrasts)),
            "mean": float(contrasts.contrast.mean()),
            "ci_low": float(low), "ci_high": float(high),
            "excludes_zero": bool(not (low <= 0 <= high)),
            "verdict": classify(low, high, margin_raw)}


rows = []
for seed in range(N_SEEDS):
    pairs = match_controls(cases, controls, features=MATCH_FEATURES, k=1,
                           caliper=DEFAULT_CALIPER, seed=seed, exact=DEFAULT_EXACT)
    result = evaluate(pairs, BOOT_SEED)
    if result:
        rows.append({"seed": seed, **result})
    if (seed + 1) % 50 == 0:
        print(f"  {seed + 1}/{N_SEEDS} seeds", flush=True)

sweep = pd.DataFrame(rows)
sweep.to_csv("results/analysis/matching_seed_sweep.csv", index=False)

best = evaluate(match_controls_optimal(cases, controls, features=MATCH_FEATURES,
                                       caliper=DEFAULT_CALIPER, exact=DEFAULT_EXACT),
                BOOT_SEED)

print(f"\n=== 200 greedy seeds ===")
print(f"  pairs            {sweep.n_pairs.min()}-{sweep.n_pairs.max()} "
      f"(median {int(sweep.n_pairs.median())})")
print(f"  mean contrast    {sweep['mean'].min():+.4f} to {sweep['mean'].max():+.4f} "
      f"(median {sweep['mean'].median():+.4f})")
print(f"  in SD units      {sweep['mean'].min()/reference_sd:+.3f} to "
      f"{sweep['mean'].max()/reference_sd:+.3f}")
print(f"  ALL point estimates positive: {bool((sweep['mean'] > 0).all())}")
print(f"  CI excludes zero in {int(sweep.excludes_zero.sum())} of {len(sweep)} seeds "
      f"({100 * sweep.excludes_zero.mean():.0f}%)")
print("\n  verdicts across seeds:")
for verdict, count in sweep.verdict.value_counts().items():
    print(f"    {verdict:38s} {count}")

print(f"\n=== deterministic optimal matching (primary) ===")
print(f"  pairs {best['n_pairs']}   mean {best['mean']:+.4f} "
      f"({best['mean']/reference_sd:+.3f} SD)")
print(f"  95% CI [{best['ci_low']:+.4f}, {best['ci_high']:+.4f}]   {best['verdict']}")

Path("results/analysis/matching_sensitivity.json").write_text(json.dumps({
    "purpose": "how much the primary result depends on the matching algorithm",
    "n_seeds": N_SEEDS, "n_boot_per_seed": N_BOOT,
    "reference_sd": round(reference_sd, 6), "margin_raw": round(margin_raw, 6),
    "greedy": {
        "pairs_min": int(sweep.n_pairs.min()), "pairs_max": int(sweep.n_pairs.max()),
        "mean_min": round(float(sweep["mean"].min()), 6),
        "mean_max": round(float(sweep["mean"].max()), 6),
        "mean_median": round(float(sweep["mean"].median()), 6),
        "all_point_estimates_positive": bool((sweep["mean"] > 0).all()),
        "fraction_ci_excludes_zero": round(float(sweep.excludes_zero.mean()), 4),
        "verdicts": sweep.verdict.value_counts().to_dict()},
    "optimal": best,
    "interpretation": "every point estimate is positive; whether the interval "
                      "excludes zero depends on which matching is selected",
}, indent=2))
print("\nwrote results/analysis/matching_seed_sweep.csv and results/analysis/matching_sensitivity.json")
