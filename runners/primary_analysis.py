"""The primary comparison, end to end and reproducible from committed inputs.

    Given the native backbone and the surrounding native sequence, does
    ProteinMPNN assign a higher conditional sequon score to experimentally
    occupied sites than to structurally matched sites with no modelled glycan?

Usage:  primary_analysis.py [optimal|greedy_seed0|k5|bacterial|cytosolic]

`optimal` is the primary comparison. The others are sensitivity analyses and
diagnostics; see docs/primary_result.md for which conclusions each supports.

Contrast construction, the resampling unit and the bootstrap live in
`experimental_glycosylation_sites.contrasts`, shared with the matching
sensitivity sweep so the two cannot compute their intervals differently.
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, "src")
from experimental_glycosylation_sites.contrasts import (
    build_contrasts, classify, cluster_bootstrap)

MARGIN_SD = 0.2                  # frozen, exploratory
N_BOOT, BOOT_SEED = 10000, 20260818
KEY = ["accession", "position"]

LABEL = sys.argv[1] if len(sys.argv) > 1 else "optimal"
PAIRS = Path(f"results/matched_pairs_{LABEL}.csv")

# Only the diagnostics draw on the external control pool. Loading it for the
# primary comparison would make an analysis of 16 dataset pairs depend on files
# it never reads, so a stale control file could silently move a result that does
# not involve controls at all.
NEEDS_CONTROL_POOL = LABEL in ("bacterial", "cytosolic")

ROLE = {"optimal": "PRIMARY", "greedy_seed0": "sensitivity (greedy, seed 0)",
        "k5": "sensitivity (greedy, k=5)",
        "bacterial": "DIAGNOSTIC (bacterial extracytoplasmic)",
        "cytosolic": "DIAGNOSTIC (cytosolic eukaryotic)"}.get(LABEL, LABEL)


def load(manifest_path, score_path):
    if not Path(manifest_path).exists() or not Path(score_path).exists():
        return pd.DataFrame()
    manifest = pd.read_csv(manifest_path, low_memory=False)
    manifest = manifest[manifest.scoreable.astype(bool)].copy()
    scores = pd.read_csv(score_path, low_memory=False)
    for frame in (manifest, scores):
        frame["accession"] = frame.accession.astype(str)
        frame["position"] = frame.position.astype(int)
    return manifest.merge(
        scores[KEY + ["conditional_sequon_score", "conditional_sequon_score_sd"]],
        on=KEY, how="inner")


dataset = load("results/candidate_manifest_dataset.csv", "results/scores_dataset.csv")
controls = (load("results/manifest_matched_controls.csv", "results/scores_controls.csv")
            if NEEDS_CONTROL_POOL else pd.DataFrame())
site = pd.concat([dataset, controls], ignore_index=True) if len(controls) else dataset.copy()
if "ortholog_clusters" not in site.columns:
    site["ortholog_clusters"] = pd.NA

# The frozen rule: pooled across all structurally scoreable DATASET sites,
# without reference to their labels. Controls never enter it, so every
# comparison sits on one scale that cannot move when a control pool is rebuilt.
reference_sd = float(dataset.conditional_sequon_score.std(ddof=1))
margin_raw = MARGIN_SD * reference_sd
reference_population = (
    f"{len(dataset)} scoreable dataset sites "
    f"({int((dataset.occupancy_status == 'occupied_supported').sum())} occupied + "
    f"{int((dataset.occupancy_status == 'observed_unmodified').sum())} internal control); "
    "labels not consulted; control sites excluded by construction")

print(f"comparison                     : {LABEL}  [{ROLE}]")
print(f"reference population           : {reference_population}")
if len(controls):
    print(f"matched control sites loaded   : {len(controls)}")
print(f"reference SD                   : {reference_sd:.4f} log-odds")
print(f"equivalence margin (+/-0.2 SD) : +/-{margin_raw:.4f} log-odds\n")

pairs = pd.read_csv(PAIRS, low_memory=False)
contrasts = build_contrasts(pairs, site)
contrasts.to_csv(f"results/contrasts_{LABEL}.csv", index=False)

reuse = contrasts.control_proteins.str.split(";").explode().value_counts()
print(f"contrasts {len(contrasts)}   occupied proteins {contrasts.case_accession.nunique()}   "
      f"ortholog clusters {contrasts.ortholog_cluster.nunique()}   "
      f"resampling units {contrasts.resample_unit.nunique()}")
print(f"control proteins used {len(reuse)}   most reused appears in {int(reuse.max())} contrasts\n")

observed = float(contrasts.contrast.mean())
draws = cluster_bootstrap(contrasts, N_BOOT, BOOT_SEED)
low, high = np.percentile(draws, [2.5, 97.5])
verdict = classify(low, high, margin_raw)

rng = np.random.default_rng(BOOT_SEED + 1)
values = contrasts.contrast.to_numpy()
site_draws = np.array([rng.choice(values, len(values), replace=True).mean()
                       for _ in range(N_BOOT)])
site_low, site_high = np.percentile(site_draws, [2.5, 97.5])

def sd(x):
    return x / reference_sd

n_positive = int((values > 0).sum())
loo = np.array([np.delete(values, i).mean() for i in range(len(values))])
robustness = {
    "median": round(float(np.median(values)), 4),
    "trimmed_mean_10pct": round(float(st.trim_mean(values, 0.1)), 4),
    "n_positive": n_positive, "n_contrasts": int(len(values)),
    "sign_test_p": round(float(st.binomtest(n_positive, len(values)).pvalue), 4),
    "wilcoxon_p": round(float(st.wilcoxon(values).pvalue), 4),
    "leave_one_out_mean_min": round(float(loo.min()), 4),
    "leave_one_out_mean_max": round(float(loo.max()), 4),
}

print(f"=== occupied vs matched control [{ROLE}] ===")
print(f"  mean paired difference   {observed:+.4f} log-odds   ({sd(observed):+.3f} SD)")
print(f"  95% CI, cluster-aware    [{low:+.4f}, {high:+.4f}]   "
      f"([{sd(low):+.3f}, {sd(high):+.3f}] SD)")
print(f"  95% CI, site-level       [{site_low:+.4f}, {site_high:+.4f}]   "
      f"([{sd(site_low):+.3f}, {sd(site_high):+.3f}] SD)   <- ignores clustering, not used")
print(f"  equivalence margin       +/-{margin_raw:.4f} ({MARGIN_SD} SD)")
print(f"  VERDICT                  {verdict}")
print(f"  occupied scores higher in {n_positive} of {len(values)} contrasts")
print("\nrobustness — does the mean survive its own skew?")
print(f"  median {robustness['median']:+.4f}   10% trimmed mean "
      f"{robustness['trimmed_mean_10pct']:+.4f}")
print(f"  leave-one-out mean stays in [{robustness['leave_one_out_mean_min']:+.4f}, "
      f"{robustness['leave_one_out_mean_max']:+.4f}]")
print(f"  sign test p={robustness['sign_test_p']:.3f} (direction alone)   "
      f"Wilcoxon signed-rank p={robustness['wilcoxon_p']:.3f}")

by_subtype = {name: {"n": int(len(group)),
                     "mean": round(float(group.contrast.mean()), 4),
                     "mean_sd_units": round(sd(float(group.contrast.mean())), 4)}
              for name, group in contrasts.groupby("subtype")}
print("\nby subtype (matched exactly, so these are like-for-like):")
for name, stats in by_subtype.items():
    print(f"  {name}  n={stats['n']:2d}  {stats['mean']:+.4f} log-odds "
          f"({stats['mean_sd_units']:+.3f} SD)")

Path(f"results/analysis_{LABEL}.json").write_text(json.dumps({
    "comparison": LABEL, "role": ROLE,
    "question": "Does ProteinMPNN score occupied sequons higher than matched "
                "sequons with no modelled glycan, given the native backbone and sequence?",
    "model": "ProteinMPNN v_48_020", "conditioning": "conditional", "n_decoding_orders": 8,
    "reference_sd": round(reference_sd, 6),
    "reference_sd_population": reference_population,
    "margin_standardised": MARGIN_SD, "margin_raw": round(margin_raw, 6),
    "n_contrasts": int(len(contrasts)),
    "n_occupied_proteins": int(contrasts.case_accession.nunique()),
    "n_ortholog_clusters": int(contrasts.ortholog_cluster.nunique()),
    "n_resampling_units": int(contrasts.resample_unit.nunique()),
    "max_control_protein_reuse": int(reuse.max()),
    "mean_difference_raw": round(observed, 6),
    "mean_difference_sd": round(sd(observed), 6),
    "ci95_raw": [round(low, 6), round(high, 6)],
    "ci95_sd": [round(sd(low), 6), round(sd(high), 6)],
    "ci95_site_level_not_used": [round(site_low, 6), round(site_high, 6)],
    "verdict": verdict, "robustness": robustness, "by_subtype": by_subtype,
    "bootstrap": {"n": N_BOOT, "seed": BOOT_SEED,
                  "unit": "connected component of occupied ortholog clusters "
                          "and shared control proteins"},
}, indent=2))
print(f"\nwrote results/analysis_{LABEL}.json and results/contrasts_{LABEL}.csv")
