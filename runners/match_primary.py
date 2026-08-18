"""Match occupied sites to internal controls, among scoreable sites only.

The primary comparison. Internal controls are sequons in the same solved
glycoprotein structures that carry no modelled glycan under conditions where a
glycan would have been visible had it been there. They are the only controls
that hold organism, subcellular compartment and experimental context constant,
and there are very few of them, which is the binding constraint on this study.

Two things are fixed here that were wrong before:

  * only sites ProteinMPNN can actually decode enter the pool, so matched sets
    cannot lose members afterwards;
  * subtype is matched exactly, so an occupied NXS is never paired against an
    unoccupied NXT.

Structural features and caliper are unchanged. The neighbour count is set to
one, and the reason is pre-specified here rather than chosen afterwards: there
are 28 internal controls against 314 occupied sites, so greedy matching without
replacement at k=5 lets the first few cases absorb the whole control pool. That
yielded 14 occupied cases from 314. At k=1 every control supports a distinct
occupied case, which is what the analysis unit — one contrast per occupied site
— is counting. k=5 is retained as a sensitivity analysis and both are written.

This was fixed before any score for these sites existed under the corrected
scorer. It is recorded as a deviation because the earlier, defective run's
primary estimate had already been seen, so this choice cannot claim to have been
made blind.
"""
import json, sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.matching import (
    DEFAULT_CALIPER, DEFAULT_EXACT, MATCH_FEATURES,
    balance_report, match_controls, weighted_balance_report)

PRIMARY_K, SENSITIVITY_K, SEED = 1, 5, 0
man = pd.read_csv("results/candidate_manifest_dataset.csv", low_memory=False)

total = len(man)
man = man[man.scoreable.astype(bool)].copy()
print(f"candidate sites {total} -> {len(man)} scoreable "
      f"({total - len(man)} dropped before matching, not after)")

cases = man[man.occupancy_status == "occupied_supported"].copy()
controls = man[man.occupancy_status == "observed_unmodified"].copy()
print(f"occupied {len(cases)}   internal controls {len(controls)}")
print("control subtypes: " + controls.subtype.value_counts().to_dict().__str__())

for k, label in ((PRIMARY_K, "primary"), (SENSITIVITY_K, "k5")):
    pairs = match_controls(cases, controls, features=MATCH_FEATURES, k=k,
                           caliper=DEFAULT_CALIPER, seed=SEED, exact=DEFAULT_EXACT)
    pairs["comparison"] = "vs_internal_control"
    pairs.to_csv(f"results/matched_pairs_{label}.csv", index=False)

    report = {
        "comparison": "vs_internal_control", "label": label,
        "matching": {"features": list(MATCH_FEATURES), "exact": list(DEFAULT_EXACT),
                     "k": k, "caliper": DEFAULT_CALIPER, "seed": SEED,
                     "replacement": False,
                     "pool": "structurally scoreable sites only, determined before matching"},
        "unweighted": balance_report(cases, controls, pairs),
        "weighted": weighted_balance_report(cases, controls, pairs),
    }
    Path(f"results/matching_balance_{label}.json").write_text(json.dumps(report, indent=2))

    n_cases = pairs.drop_duplicates(["case_accession", "case_position"]).shape[0]
    n_ctrl = pairs.drop_duplicates(["control_accession", "control_position"]).shape[0]
    print(f"\n--- k={k} ({label}) ---")
    print(f"pairs {len(pairs)}   distinct occupied cases {n_cases}   distinct controls {n_ctrl}")

    merged = (pairs.merge(cases[["accession", "position", "subtype"]],
                          left_on=["case_accession", "case_position"],
                          right_on=["accession", "position"])
                   .merge(controls[["accession", "position", "subtype"]],
                          left_on=["control_accession", "control_position"],
                          right_on=["accession", "position"], suffixes=("_case", "_ctrl")))
    assert (merged.subtype_case == merged.subtype_ctrl).all(), "subtype mismatch survived matching"
    print("subtype identical in every pair: yes   "
          + merged.subtype_case.value_counts().to_dict().__str__())

    print("balance (standardised mean difference, occupied minus control):")
    for feature in MATCH_FEATURES:
        unw = report["unweighted"]["features"][feature]
        wt = report["weighted"]["features"][feature]
        print(f"  {feature:26s} before {unw['smd_before']:+.3f}   "
              f"after {unw['smd_after']:+.3f}   weighted {wt['smd_weighted']:+.3f}")
