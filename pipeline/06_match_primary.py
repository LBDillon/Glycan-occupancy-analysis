"""Matched sets for the primary comparison.

Internal controls are sequons in solved glycoprotein structures carrying no
modelled glycan under conditions where a glycan would have been visible. They
hold organism, compartment and experimental context roughly constant, there are
28 of them that can be scored, and that scarcity is the binding constraint on
the study.

Matching is now deterministic. Greedy nearest-neighbour matching walks the cases
in a seeded random order, so with a pool this small an early case can consume the
only admissible partner for a later one — both the number of pairs and their
identity depend on the seed. The primary set is instead the assignment that
maximises the number of admissible pairs and, among those, minimises total
distance. It has no seed and no ordering.

Three sets are written:

  optimal        the primary set, deterministic
  greedy_seed0   the previously reported greedy result, retained as sensitivity
  k5             greedy allowing up to five controls per case, retained likewise

Features, caliper and the exact NXS/NXT requirement are identical across all
three. Only the algorithm that assigns pairs differs.
"""
import json, sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.matching import (
    DEFAULT_CALIPER, DEFAULT_EXACT, MATCH_FEATURES,
    balance_report, match_controls, match_controls_optimal, weighted_balance_report)

SEED = 0
man = pd.read_csv("results/candidate_manifest_dataset.csv", low_memory=False)

total = len(man)
man = man[man.scoreable.astype(bool)].copy()
print(f"candidate sites {total} -> {len(man)} scoreable "
      f"({total - len(man)} dropped before matching, not after)")

cases = man[man.occupancy_status == "occupied_supported"].copy()
controls = man[man.occupancy_status == "observed_unmodified"].copy()
print(f"occupied {len(cases)}   internal controls {len(controls)} "
      f"{controls.subtype.value_counts().to_dict()}")

recipes = {
    "optimal": lambda: match_controls_optimal(
        cases, controls, features=MATCH_FEATURES,
        caliper=DEFAULT_CALIPER, exact=DEFAULT_EXACT),
    "greedy_seed0": lambda: match_controls(
        cases, controls, features=MATCH_FEATURES, k=1,
        caliper=DEFAULT_CALIPER, seed=SEED, exact=DEFAULT_EXACT),
    "k5": lambda: match_controls(
        cases, controls, features=MATCH_FEATURES, k=5,
        caliper=DEFAULT_CALIPER, seed=SEED, exact=DEFAULT_EXACT),
}

for label, build in recipes.items():
    pairs = build()
    pairs["comparison"] = "vs_internal_control"
    pairs.to_csv(f"results/matched_pairs_{label}.csv", index=False)

    report = {
        "comparison": "vs_internal_control", "label": label,
        "role": "primary" if label == "optimal" else "sensitivity",
        "matching": {
            "algorithm": ("maximum-cardinality, minimum-total-distance assignment"
                          if label == "optimal" else "greedy nearest neighbour"),
            "deterministic": label == "optimal",
            "features": list(MATCH_FEATURES), "exact": list(DEFAULT_EXACT),
            "caliper": DEFAULT_CALIPER, "replacement": False,
            "k": 5 if label == "k5" else 1,
            "seed": None if label == "optimal" else SEED,
            "pool": "structurally scoreable sites only, determined before matching"},
        "unweighted": balance_report(cases, controls, pairs),
        "weighted": weighted_balance_report(cases, controls, pairs),
    }
    Path(f"results/matching_balance_{label}.json").write_text(json.dumps(report, indent=2))

    n_cases = pairs.drop_duplicates(["case_accession", "case_position"]).shape[0]
    print(f"\n--- {label} ({report['role']}) ---")
    print(f"pairs {len(pairs)}   occupied cases {n_cases}   "
          f"controls used {pairs.drop_duplicates(['control_accession','control_position']).shape[0]}"
          f" of {len(controls)}")

    merged = (pairs.merge(cases[["accession", "position", "subtype"]],
                          left_on=["case_accession", "case_position"],
                          right_on=["accession", "position"])
                   .merge(controls[["accession", "position", "subtype"]],
                          left_on=["control_accession", "control_position"],
                          right_on=["accession", "position"], suffixes=("_case", "_ctrl")))
    assert (merged.subtype_case == merged.subtype_ctrl).all(), "subtype mismatch survived matching"
    print(f"subtype identical in every pair: yes   {merged.subtype_case.value_counts().to_dict()}")
    print(f"mean matching distance {pairs.distance.mean():.4f}")
    for feature in MATCH_FEATURES:
        u = report["unweighted"]["features"][feature]
        print(f"  {feature:26s} before {u['smd_before']:+.3f}   after {u['smd_after']:+.3f}")
