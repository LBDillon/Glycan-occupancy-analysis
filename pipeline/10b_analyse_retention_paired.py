"""Design retention as a PAIRED contrast, not a class average.

The class-average view (pipeline/10_analyse_retention_by_class.py) compares whole
populations, so any difference between them in fold, size or composition rides
along. Every control set here was matched site-by-site to occupied sites on
local structure and sequon subtype, so the paired contrast is available and is
the better test: it asks whether THIS occupied sequon is retained more often
than the structurally matched control chosen for it.

Retention remains a secondary, descriptive outcome. It has no pre-specified
equivalence margin, and this analysis was not pre-registered — it was run after
the eukaryotic secretory set was added. Read it as a lead, not a conclusion.
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, "src")
from experimental_glycosylation_sites.contrasts import assign_resample_units, cluster_bootstrap

FULL = "std_frac_full_sequon_retained"
N_BOOT, BOOT_SEED = 10000, 20260818

COMPARISONS = [
    ("internal control", "results/matching/matched_pairs_optimal.csv"),
    ("eukaryotic secretory", "results/matching/matched_pairs_secretory.csv"),
    ("bacterial", "results/matching/matched_pairs_bacterial.csv"),
    ("cytosolic", "results/matching/matched_pairs_cytosolic.csv"),
]

ret = pd.read_csv("results/designs/retention_all_classes.csv", low_memory=False)
ret["accession"] = ret.accession.astype(str)
ret["position"] = ret.position.astype(int)
# A site can appear under more than one structure; keep one retention value each.
ret = ret.drop_duplicates(["accession", "position"])
retention = dict(zip(zip(ret.accession, ret.position), ret[FULL]))

man = pd.read_csv("results/manifests/candidate_manifest_dataset.csv", low_memory=False)
man["accession"] = man.accession.astype(str)
man["position"] = man.position.astype(int)
man = man.drop_duplicates(["accession", "position"])
clusters = dict(zip(zip(man.accession, man.position),
                    man.get("ortholog_clusters", pd.Series(dtype=object))))

summary = {}
print(f"{'comparison':24s} {'pairs':>6s} {'occ':>7s} {'ctl':>7s} {'diff':>8s} "
      f"{'95% CI':>20s} {'informative':>12s} {'wilcoxon':>9s}")
print("-" * 104)
for label, path in COMPARISONS:
    if not Path(path).exists():
        continue
    pairs = pd.read_csv(path)
    rows = []
    for r in pairs.itertuples(index=False):
        case = (str(r.case_accession), int(r.case_position))
        control = (str(r.control_accession), int(r.control_position))
        if case in retention and control in retention:
            cluster = clusters.get(case)
            rows.append({
                "case_accession": case[0], "case_position": case[1],
                "occ": float(retention[case]), "ctl": float(retention[control]),
                "control_proteins": control[0],
                "ortholog_cluster": (str(cluster).split(";")[0]
                                     if pd.notna(cluster) and str(cluster) != "nan"
                                     else f"solo:{case[0]}"),
            })
    if len(rows) < 5:
        print(f"{label:24s} {len(rows):>6d}  too few pairs designed on both sides")
        continue
    frame = pd.DataFrame(rows)
    frame["contrast"] = frame.occ - frame.ctl
    frame = assign_resample_units(frame)
    draws = cluster_bootstrap(frame, N_BOOT, BOOT_SEED)
    low, high = np.percentile(draws, [2.5, 97.5])
    informative = int((frame.contrast != 0).sum())
    wilcoxon = float(st.wilcoxon(frame.contrast).pvalue) if informative else float("nan")
    print(f"{label:24s} {len(frame):>6d} {frame.occ.mean():>7.4f} {frame.ctl.mean():>7.4f} "
          f"{frame.contrast.mean():>+8.4f} [{low:>+7.4f},{high:>+7.4f}] "
          f"{informative:>12d} {wilcoxon:>9.4f}")
    summary[label] = {
        "n_pairs": int(len(frame)), "n_informative": informative,
        "n_tied": int((frame.contrast == 0).sum()),
        "occupied_mean": round(float(frame.occ.mean()), 4),
        "control_mean": round(float(frame.ctl.mean()), 4),
        "paired_difference": round(float(frame.contrast.mean()), 4),
        "ci95": [round(low, 4), round(high, 4)],
        "excludes_zero": bool(not (low <= 0 <= high)),
        "wilcoxon_p": round(wilcoxon, 4),
        "n_resampling_units": int(frame.resample_unit.nunique()),
    }
    frame.to_csv(f"results/analysis/retention_paired_{label.split()[0]}.csv", index=False)

Path("results/analysis/retention_paired.json").write_text(json.dumps({
    "status": "secondary and descriptive; no pre-specified margin; not pre-registered",
    "unit": "one occupied site, against the control matched to it",
    "bootstrap": "connected components of ortholog clusters and shared control proteins",
    "comparisons": summary,
}, indent=2))
print("\nwrote results/analysis/retention_paired.json")
