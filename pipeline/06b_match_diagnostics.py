"""Bacterial and cytosolic matched sets — diagnostics, not the primary result.

Each of these control sets is confounded by construction. The cytosolic set
matches taxonomy and differs in subcellular compartment; the bacterial set
matches compartment and differs in taxonomy. Neither can substitute for the
internal controls, and they cannot compensate for how few of those there are:
a large, well-matched comparison against the wrong population is still the
wrong comparison. They are read against one another, in an appendix.

Same pool discipline as the primary: scoreable sites only, exact subtype, k=1.
"""
import json, sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.matching import (
    DEFAULT_CALIPER, DEFAULT_EXACT, MATCH_FEATURES,
    balance_report, match_controls_optimal, weighted_balance_report)

SEED = 0
data = pd.read_csv("results/manifests/candidate_manifest_dataset.csv", low_memory=False)
cases = data[data.scoreable.astype(bool)
             & (data.occupancy_status == "occupied_supported")].copy()

POOLS = {
    "bacterial": ("results/manifests/candidate_manifest_controls.csv",
                  "control_bacterial_extracytoplasmic"),
    "cytosolic": ("results/manifests/candidate_manifest_controls.csv",
                  "control_cytosolic_eukaryotic"),
    "secretory": ("results/manifests/candidate_manifest_secretory.csv",
                  "control_secretory_eukaryotic_unannotated"),
}

for label, (source, status) in POOLS.items():
    if not Path(source).exists():
        print(f"\n--- {label}: {source} not built yet, skipping ---")
        continue
    ctrl = pd.read_csv(source, low_memory=False)
    pool = ctrl[ctrl.scoreable.astype(bool) & (ctrl.occupancy_status == status)].copy()
    if pool.empty:
        print(f"\n--- {label}: no scoreable sites with status {status}, skipping ---")
        continue
    # Deterministic, matching the primary comparison (Amendment 2).
    pairs = match_controls_optimal(cases, pool, features=MATCH_FEATURES,
                                   caliper=DEFAULT_CALIPER, exact=DEFAULT_EXACT)
    pairs["comparison"] = f"vs_{label}_control"
    pairs.to_csv(f"results/matching/matched_pairs_{label}.csv", index=False)

    report = {"comparison": f"vs_{label}_control", "status": "diagnostic only",
              "matching": {"features": list(MATCH_FEATURES), "exact": list(DEFAULT_EXACT),
                           "algorithm": "maximum-cardinality, minimum-total-distance",
                           "deterministic": True, "k": 1,
                           "caliper": DEFAULT_CALIPER},
              "unweighted": balance_report(cases, pool, pairs),
              "weighted": weighted_balance_report(cases, pool, pairs)}
    Path(f"results/matching/matching_balance_{label}.json").write_text(json.dumps(report, indent=2))

    print(f"\n--- {label}: {len(pool)} scoreable controls ---")
    print(f"pairs {len(pairs)}   occupied cases "
          f"{pairs.drop_duplicates(['case_accession','case_position']).shape[0]}")
    for feature in MATCH_FEATURES:
        u = report["unweighted"]["features"][feature]
        print(f"  {feature:26s} before {u['smd_before']:+.3f}   after {u['smd_after']:+.3f}")
