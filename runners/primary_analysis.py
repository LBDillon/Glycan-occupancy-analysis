"""The primary comparison, end to end and reproducible from committed inputs.

    Given the native backbone and the surrounding native sequence, does
    ProteinMPNN assign a higher conditional sequon score to experimentally
    occupied sites than to structurally matched sites with no modelled glycan?

One contrast per occupied site, as the frozen configuration defines it:

    site_contrast = occupied_score - mean(scores of its matched controls)

Uncertainty is the part most easily got wrong. Two distinct dependencies run
through these contrasts. Occupied sites in the same ortholog cluster are near
copies of one another, and a single control protein can be matched to several
different occupied cases. Resampling on either alone leaves the other
unaccounted for, so the resampling unit here is the connected component of the
bipartite graph joining occupied clusters to the control proteins they share.
Contrasts that are linked, however indirectly, move together or not at all.
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

MARGIN_SD = 0.2                  # frozen, exploratory
N_BOOT, BOOT_SEED = 10000, 20260818
KEY = ["accession", "position"]

LABEL = sys.argv[1] if len(sys.argv) > 1 else "primary"
PAIRS = Path(f"results/matched_pairs_{LABEL}.csv")

def load(manifest_path, score_path):
    if not Path(manifest_path).exists() or not Path(score_path).exists():
        return pd.DataFrame()
    m = pd.read_csv(manifest_path, low_memory=False)
    m = m[m.scoreable.astype(bool)].copy()
    sc = pd.read_csv(score_path, low_memory=False)
    for d in (m, sc):
        d["accession"] = d.accession.astype(str)
        d["position"] = d.position.astype(int)
    return m.merge(sc[KEY + ["conditional_sequon_score", "conditional_sequon_score_sd"]],
                   on=KEY, how="inner")

dataset = load("results/candidate_manifest_dataset.csv", "results/scores_dataset.csv")
controls = load("results/manifest_matched_controls.csv", "results/scores_controls.csv")
site = pd.concat([dataset, controls], ignore_index=True)
if "ortholog_clusters" not in site.columns:
    site["ortholog_clusters"] = pd.NA

# ---------------------------------------------------------------- reference SD
# The frozen rule: pooled across all structurally scoreable DATASET sites,
# computed without reference to their labels. Controls never enter it, so every
# comparison is expressed on one common scale.
reference_sd = float(dataset.conditional_sequon_score.std(ddof=1))
margin_raw = MARGIN_SD * reference_sd
print(f"scoreable dataset sites scored : {len(dataset)} "
      f"({int((dataset.occupancy_status=='occupied_supported').sum())} occupied, "
      f"{int((dataset.occupancy_status=='observed_unmodified').sum())} internal control)")
if len(controls):
    print(f"matched control sites scored   : {len(controls)}")
print(f"reference SD                   : {reference_sd:.4f} log-odds")
print(f"equivalence margin (+/-0.2 SD) : +/-{margin_raw:.4f} log-odds\n")

# ------------------------------------------------------------------- contrasts
pairs = pd.read_csv(PAIRS, low_memory=False)
pairs["case_accession"] = pairs.case_accession.astype(str)
pairs["control_accession"] = pairs.control_accession.astype(str)

lookup = site.set_index(["accession", "position"])
def score_of(acc, pos):
    try:
        return float(lookup.loc[(acc, int(pos))].conditional_sequon_score)
    except KeyError:
        return np.nan

rows = []
for (acc, pos), grp in pairs.groupby(["case_accession", "case_position"]):
    case_score = score_of(acc, pos)
    ctrl = [(c.control_accession, int(c.control_position),
             score_of(c.control_accession, c.control_position))
            for c in grp.itertuples(index=False)]
    ctrl = [c for c in ctrl if np.isfinite(c[2])]
    if not np.isfinite(case_score) or not ctrl:
        continue
    clusters = lookup.loc[(acc, int(pos))].get("ortholog_clusters", pd.NA)
    rows.append({
        "case_accession": acc, "case_position": int(pos),
        "subtype": lookup.loc[(acc, int(pos))].subtype,
        "case_score": case_score,
        "control_mean_score": float(np.mean([c[2] for c in ctrl])),
        "n_controls": len(ctrl),
        "control_proteins": ";".join(sorted({c[0] for c in ctrl})),
        "ortholog_cluster": (str(clusters).split(";")[0]
                             if pd.notna(clusters) and str(clusters) else f"solo:{acc}"),
    })
contrasts = pd.DataFrame(rows)
contrasts["contrast"] = contrasts.case_score - contrasts.control_mean_score
contrasts.to_csv(f"results/contrasts_{LABEL}.csv", index=False)

# ------------------------------------------------- resampling unit: components
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

for r in contrasts.itertuples(index=False):
    anchor = f"cluster:{r.ortholog_cluster}"
    find(anchor)
    for protein in r.control_proteins.split(";"):
        union(anchor, f"control:{protein}")
contrasts["resample_unit"] = [find(f"cluster:{r.ortholog_cluster}")
                              for r in contrasts.itertuples(index=False)]

n_units = contrasts.resample_unit.nunique()
print(f"contrasts {len(contrasts)}   occupied proteins {contrasts.case_accession.nunique()}   "
      f"ortholog clusters {contrasts.ortholog_cluster.nunique()}   resampling units {n_units}")
reuse = contrasts.control_proteins.str.split(";").explode().value_counts()
print(f"control proteins used {len(reuse)}   most reused appears in {int(reuse.max())} contrasts\n")

# ------------------------------------------------------------------- bootstrap
rng = np.random.default_rng(BOOT_SEED)
units = contrasts.resample_unit.unique()
by_unit = {u: contrasts.contrast[contrasts.resample_unit == u].to_numpy() for u in units}

def boot(values_by_unit, unit_list):
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        picked = rng.integers(0, len(unit_list), len(unit_list))
        pooled = np.concatenate([values_by_unit[unit_list[j]] for j in picked])
        draws[i] = pooled.mean()
    return draws

observed = float(contrasts.contrast.mean())
draws = boot(by_unit, list(units))
lo, hi = np.percentile(draws, [2.5, 97.5])

# a site-level bootstrap, shown only to expose how much the clustering matters
site_draws = np.array([rng.choice(contrasts.contrast.to_numpy(),
                                  len(contrasts), replace=True).mean()
                       for _ in range(N_BOOT)])
site_lo, site_hi = np.percentile(site_draws, [2.5, 97.5])

def sd(x): return x / reference_sd

# Four outcomes, not two. A confidence interval that excludes zero but reaches
# past the margin has established a direction and not a magnitude, and calling
# that "inconclusive" discards the half of it that is informative — just as
# calling it "a difference" would overstate the half that is not.
within_margin = (lo >= -margin_raw) and (hi <= margin_raw)
excludes_zero = not (lo <= 0 <= hi)
beyond_margin = (lo > margin_raw) or (hi < -margin_raw)
if within_margin:
    verdict = "equivalent within the margin"
elif beyond_margin:
    verdict = "difference beyond the margin"
elif excludes_zero:
    verdict = "directional, magnitude undetermined"
else:
    verdict = "inconclusive"

# The mean is the pre-specified estimand, but with sixteen contrasts it can be
# carried by a skewed few. These say whether it is.
x = contrasts.contrast.to_numpy()
n_pos = int((x > 0).sum())
loo = np.array([np.delete(x, i).mean() for i in range(len(x))])
robustness = {
    "median": round(float(np.median(x)), 4),
    "trimmed_mean_10pct": round(float(st.trim_mean(x, 0.1)), 4),
    "n_positive": n_pos, "n_contrasts": int(len(x)),
    "sign_test_p": round(float(st.binomtest(n_pos, len(x)).pvalue), 4),
    "wilcoxon_p": round(float(st.wilcoxon(x).pvalue), 4),
    "leave_one_out_mean_min": round(float(loo.min()), 4),
    "leave_one_out_mean_max": round(float(loo.max()), 4),
}

CONTROL_NAME = {"primary": "internal control (no modelled glycan)",
                "k5": "internal control (no modelled glycan), k=5 sensitivity",
                "bacterial": "bacterial extracytoplasmic control [DIAGNOSTIC]",
                "cytosolic": "cytosolic eukaryotic control [DIAGNOSTIC]"}
print(f"=== occupied vs {CONTROL_NAME.get(LABEL, LABEL)} ===")
print(f"  mean paired difference   {observed:+.4f} log-odds   "
      f"({sd(observed):+.3f} SD)")
print(f"  95% CI, cluster-aware    [{lo:+.4f}, {hi:+.4f}]   "
      f"([{sd(lo):+.3f}, {sd(hi):+.3f}] SD)")
print(f"  95% CI, site-level       [{site_lo:+.4f}, {site_hi:+.4f}]   "
      f"([{sd(site_lo):+.3f}, {sd(site_hi):+.3f}] SD)   <- ignores clustering, not used")
print(f"  equivalence margin       +/-{margin_raw:.4f} ({MARGIN_SD} SD)")
print(f"  VERDICT                  {verdict}")
print(f"  occupied scores higher in {n_pos} of {len(contrasts)} contrasts")
print("\nrobustness — does the mean survive its own skew?")
print(f"  median {robustness['median']:+.4f}   10% trimmed mean "
      f"{robustness['trimmed_mean_10pct']:+.4f}")
print(f"  leave-one-out mean stays in [{robustness['leave_one_out_mean_min']:+.4f}, "
      f"{robustness['leave_one_out_mean_max']:+.4f}]")
print(f"  sign test p={robustness['sign_test_p']:.3f}  "
      f"(direction alone)   Wilcoxon signed-rank p={robustness['wilcoxon_p']:.3f}")

by_subtype = {st: {"n": int(len(g)), "mean": round(float(g.contrast.mean()), 4),
                   "mean_sd_units": round(sd(float(g.contrast.mean())), 4)}
              for st, g in contrasts.groupby("subtype")}
print("\nby subtype (matched exactly, so these are like-for-like):")
for st, s in by_subtype.items():
    print(f"  {st}  n={s['n']:2d}  {s['mean']:+.4f} log-odds ({s['mean_sd_units']:+.3f} SD)")

Path(f"results/analysis_{LABEL}.json").write_text(json.dumps({
    "question": "Does ProteinMPNN score occupied sequons higher than matched "
                "sequons with no modelled glycan, given the native backbone and sequence?",
    "model": "ProteinMPNN v_48_020", "conditioning": "conditional", "n_decoding_orders": 8,
    "reference_sd": round(reference_sd, 6),
    "reference_sd_population": f"{len(site)} scoreable dataset sites, labels not consulted",
    "margin_standardised": MARGIN_SD, "margin_raw": round(margin_raw, 6),
    "n_contrasts": int(len(contrasts)),
    "n_occupied_proteins": int(contrasts.case_accession.nunique()),
    "n_ortholog_clusters": int(contrasts.ortholog_cluster.nunique()),
    "n_resampling_units": int(n_units),
    "max_control_protein_reuse": int(reuse.max()),
    "mean_difference_raw": round(observed, 6),
    "mean_difference_sd": round(sd(observed), 6),
    "ci95_raw": [round(lo, 6), round(hi, 6)],
    "ci95_sd": [round(sd(lo), 6), round(sd(hi), 6)],
    "ci95_site_level_not_used": [round(site_lo, 6), round(site_hi, 6)],
    "verdict": verdict,
    "robustness": robustness,
    "by_subtype": by_subtype,
    "bootstrap": {"n": N_BOOT, "seed": BOOT_SEED,
                  "unit": "connected component of occupied ortholog clusters and shared control proteins"},
}, indent=2))
print(f"\nwrote results/analysis_{LABEL}.json and results/contrasts_{LABEL}.csv")
