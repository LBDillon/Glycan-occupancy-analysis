import json, sys
import numpy as np, pandas as pd
from pathlib import Path

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
SCORE = "conditional_sequon_score"
MARGIN_SD = 0.2
RNG = np.random.default_rng(20260817)
N_BOOT = 5000

scores = pd.concat([
    pd.read_csv("results/scores/mpnn_conditional_scores.csv", low_memory=False),
    pd.read_csv("results/scores/mpnn_conditional_scores_unmatched.csv", low_memory=False),
], ignore_index=True).drop_duplicates(KEY)

feats = pd.read_csv("results/datasets/site_structural_features.csv", low_memory=False)
manifest = pd.read_csv("results/manifests/scoring_manifest.csv", low_memory=False)

# ---- reference SD: all scoreable dataset sites, pooled, labels not consulted ----
dataset_keys = set(zip(feats[feats.features_available &
                             feats.occupancy_status.isin(["occupied_supported",
                                                          "observed_unmodified"])].accession,
                       feats[feats.features_available &
                             feats.occupancy_status.isin(["occupied_supported",
                                                          "observed_unmodified"])].position))
dataset_scores = scores[[(a, p) in dataset_keys for a, p in zip(scores.accession, scores.position)]]
REF_SD = float(dataset_scores[SCORE].std(ddof=1))
MARGIN_RAW = MARGIN_SD * REF_SD

print("=" * 74)
print("REFERENCE SCALE (labels not used)")
print(f"  dataset sites scored          {len(dataset_scores)}")
print(f"  reference SD                  {REF_SD:.4f} log-odds")
print(f"  equivalence margin +/-0.2 SD  {MARGIN_RAW:.4f} log-odds")
print("=" * 74)

lookup = scores.set_index(["accession", "position"])[SCORE].to_dict()
clusters = manifest.drop_duplicates(["accession", "position"]).set_index(
    ["accession", "position"]).ortholog_clusters.to_dict()

pairs = pd.read_csv("results/matching/matched_pairs.csv", low_memory=False)
report = {"reference_sd": round(REF_SD, 4), "margin_raw": round(MARGIN_RAW, 4),
          "margin_sd_units": MARGIN_SD, "n_dataset_sites_for_sd": len(dataset_scores),
          "comparisons": {}}

for comparison, group in pairs.groupby("comparison"):
    contrasts = []
    for (acc, pos), g in group.groupby(["case_accession", "case_position"]):
        case = lookup.get((acc, pos))
        ctl = [lookup.get((a, p)) for a, p in zip(g.control_accession, g.control_position)]
        ctl = [c for c in ctl if c is not None and np.isfinite(c)]
        if case is None or not np.isfinite(case) or not ctl:
            continue
        contrasts.append({"accession": acc, "position": pos,
                          "contrast": case - float(np.mean(ctl)), "n_controls": len(ctl),
                          "cluster": clusters.get((acc, pos), "")})
    c = pd.DataFrame(contrasts)
    if c.empty:
        continue

    # cluster bootstrap: resample whole ortholog clusters (proteins where absent)
    units = c.cluster.where(c.cluster.astype(str) != "", c.accession)
    groups_list = [c[units == u].contrast.values for u in units.unique()]
    # Resample cluster INDICES: the clusters hold different numbers of sites, so
    # sampling the arrays directly fails and sampling sites would ignore the
    # non-independence the cluster bootstrap exists to handle.
    n_units = len(groups_list)
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        picks = RNG.integers(0, n_units, size=n_units)
        boot[b] = np.concatenate([groups_list[i] for i in picks]).mean()
    lo, hi = np.percentile(boot, [2.5, 97.5])

    mean = float(c.contrast.mean())
    # Four distinct states, because "not equivalent" and "different" are not the
    # same claim and neither is "we cannot tell".
    inside = (lo > -MARGIN_RAW) and (hi < MARGIN_RAW)
    excludes_zero = (lo > 0) or (hi < 0)
    beyond_margin = (lo > MARGIN_RAW) or (hi < -MARGIN_RAW)
    if inside:
        verdict = "equivalent within the margin"
    elif beyond_margin:
        verdict = "difference established beyond the margin"
    elif excludes_zero:
        verdict = ("difference excludes zero but the interval straddles the margin "
                   "- too imprecise to call equivalence either way")
    else:
        verdict = ("INCONCLUSIVE - interval includes zero and extends beyond the "
                   "margin; not evidence of no effect")

    report["comparisons"][comparison] = {
        "n_sites": len(c), "n_proteins": int(c.accession.nunique()),
        "n_clusters": int(units.nunique()),
        "mean_contrast_logodds": round(mean, 4),
        "median_contrast_logodds": round(float(c.contrast.median()), 4),
        "ci95_logodds": [round(float(lo), 4), round(float(hi), 4)],
        "standardised_effect": round(mean / REF_SD, 4),
        "ci95_standardised": [round(float(lo / REF_SD), 4), round(float(hi / REF_SD), 4)],
        "prob_occupied_higher": round(float((c.contrast > 0).mean()), 4),
        "verdict": verdict,
        "ci_width_sd_units": round(float((hi - lo) / REF_SD), 4),
        "ci_excludes_zero": bool(excludes_zero),
    }
    # stratify by sequon subtype, since NXS and NXT are chemically distinct
    subtype = manifest.drop_duplicates(["accession", "position"]).set_index(
        ["accession", "position"]).subtype.to_dict()
    c["subtype"] = [subtype.get((a, p), "") for a, p in zip(c.accession, c.position)]
    report["comparisons"][comparison]["by_subtype"] = {
        st: {"n": int(len(g)), "mean_logodds": round(float(g.contrast.mean()), 4),
             "standardised": round(float(g.contrast.mean() / REF_SD), 4)}
        for st, g in c.groupby("subtype") if st
    }
    c.to_csv(f"results/analysis/contrasts_{comparison}.csv", index=False)

for name, r in report["comparisons"].items():
    print(f"\n{name}")
    print(f"  sites {r['n_sites']}, proteins {r['n_proteins']}, clusters {r['n_clusters']}")
    print(f"  mean contrast   {r['mean_contrast_logodds']:+.4f} log-odds  "
          f"95% CI [{r['ci95_logodds'][0]:+.4f}, {r['ci95_logodds'][1]:+.4f}]")
    print(f"  standardised    {r['standardised_effect']:+.4f} SD  "
          f"95% CI [{r['ci95_standardised'][0]:+.4f}, {r['ci95_standardised'][1]:+.4f}]")
    print(f"  P(occupied higher) {r['prob_occupied_higher']:.3f}")
    print(f"  -> {r['verdict']}")

Path("results/phase4_primary_analysis.json").write_text(json.dumps(report, indent=2))
