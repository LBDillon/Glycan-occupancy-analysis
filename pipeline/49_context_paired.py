"""Amendment 1 — occupied versus its matched control, within pair.

The population-level Step 3 compares sets of sites that share no proteins and no
chains, so occupancy is confounded with protein identity. This compares each
occupied site to the control it was matched to, using the matching the occupancy
benchmark already rests on: RSA, neighbour count and hydrophobic fraction, never
model output.

Those three are the matching variables. They are balanced by construction, so a
null on them means the matching worked -- not that occupancy is unrelated to
exposure. They are reported as a balance check and excluded from the family.

Uncertainty is a cluster bootstrap over `resample_unit` -- 72 units behind 262
pairs -- because ortholog-connected sites are not independent.

Usage:
    49_context_paired.py [--frozen benchmark_frozen/2026-08-23] [--boot 2000]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_stats import (benjamini_hochberg,
                                                            ramachandran_region_series)
from experimental_glycosylation_sites.provenance import hash_file, _git_state

MATCHING_VARIABLES = ["n_rsa", "n_neighbours_8a", "neighbour_hydrophobic_fraction_8a"]
TESTED = ["plus1_rsa", "plus2_rsa", "loop_run_length", "nd2_atoms_8a_same_chain",
          "nd2_residues_8a_same_chain", "nd2_atoms_8a_other_chain",
          "sidechain_neighbour_residues_5a", "neighbour_net_charge_8a",
          "neighbour_aromatic_fraction_8a", "nearest_aromatic_sidechain_nd2",
          "uniprot_residues_after_asn", "uniprot_residues_after_sequon",
          "distance_to_n_terminus_resolved", "distance_to_c_terminus_resolved"]
CATEGORICAL = {"n_ss_coarse": ("loop", "helix", "sheet"),
               "plus1_ss_coarse": ("loop", "helix", "sheet"),
               "plus2_ss_coarse": ("loop", "helix", "sheet"),
               "aromatic_within_8a": (True,)}
COMPARISONS = {"secretory": "matched_pairs_secretory.csv",
               "internal": "matched_pairs_optimal.csv"}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--frozen", default="benchmark_frozen/2026-08-23")
parser.add_argument("--views", default="results/datasets")
parser.add_argument("--boot", type=int, default=2000)
parser.add_argument("--out", default="results/analysis")
args = parser.parse_args()

core = pd.read_csv(Path(args.views) / "context_triplet_core.csv", low_memory=False)
core["n_rama_region"] = ramachandran_region_series(core.n_phi, core.n_psi)
for column, levels in CATEGORICAL.items():
    for level in levels:
        indicator = f"{column}=={level}"
        core[indicator] = (core[column] == level).astype(float)
        core.loc[core[column].isna(), indicator] = np.nan
FEATURES = TESTED + [f"{c}=={l}" for c, levels in CATEGORICAL.items() for l in levels]

by_site = core.drop_duplicates(["accession", "position"]).set_index(["accession", "position"])


def paired_frame(pairs: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    """One row per pair: the case value, the control value and their difference."""
    rows = []
    for pair in pairs.itertuples():
        case_key = (pair.case_accession, pair.case_position)
        ctl_key = (pair.control_accession, pair.control_position)
        if case_key not in by_site.index or ctl_key not in by_site.index:
            continue
        case, control = by_site.loc[case_key], by_site.loc[ctl_key]
        row = {"case_accession": pair.case_accession, "case_position": pair.case_position}
        for feature in FEATURES + MATCHING_VARIABLES:
            if feature in by_site.columns:
                row[feature] = float(case[feature]) - float(control[feature]) \
                    if pd.notna(case[feature]) and pd.notna(control[feature]) else np.nan
        rows.append(row)
    frame = pd.DataFrame(rows)
    return frame.merge(units, on=["case_accession", "case_position"], how="left")


def bootstrap_mean(differences: pd.Series, unit: pd.Series, n_boot: int,
                   seed: int) -> dict:
    values = differences.dropna()
    if len(values) < 3:
        return {"n_pairs": int(len(values)), "mean_difference": np.nan,
                "ci_low": np.nan, "ci_high": np.nan, "p": np.nan, "units": 0}
    keys = unit.loc[values.index].fillna("_none").astype(str)
    clusters = [g.to_numpy(dtype=float) for _, g in values.groupby(keys)]
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(clusters), len(clusters))
        draws[i] = np.concatenate([clusters[j] for j in pick]).mean()
    below, above = float(np.mean(draws <= 0)), float(np.mean(draws >= 0))
    # Standardised by the spread of the differences, so it is comparable across
    # features -- a paired effect size, not a raw unit difference.
    sd = values.std(ddof=1)
    return {"n_pairs": int(len(values)), "units": len(clusters),
            "mean_difference": float(values.mean()),
            "standardised": float(values.mean() / sd) if sd > 0 else np.nan,
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)),
            "p": float(min(1.0, max(2 * min(below, above), 1.0 / n_boot)))}


results, balance_rows = [], []
for name, filename in COMPARISONS.items():
    pairs = pd.read_csv(Path(args.frozen) / "matching" / filename, low_memory=False)
    variant = "optimal" if name == "internal" else name
    contrasts = pd.read_csv(
        Path(args.frozen) / "analysis" / f"contrasts_{variant}_alphabet_corrected.csv",
        low_memory=False)
    units = contrasts[["case_accession", "case_position", "resample_unit"]]
    frame = paired_frame(pairs, units)
    print(f"\n{name}: {len(frame)} pairs joined of {len(pairs)}, "
          f"{frame.resample_unit.nunique()} resample units")

    for feature in MATCHING_VARIABLES:
        if feature in frame.columns:
            stats = bootstrap_mean(frame[feature], frame.resample_unit, args.boot, 5)
            balance_rows.append({"comparison": name, "feature": feature, **stats})
    for feature in [f for f in FEATURES if f in frame.columns]:
        stats = bootstrap_mean(frame[feature], frame.resample_unit, args.boot, 5)
        results.append({"comparison": name, "feature": feature, **stats})

table = pd.DataFrame(results)
table["q"] = np.nan
for name, group in table.groupby("comparison"):
    table.loc[group.index, "q"] = benjamini_hochberg(group.p.tolist())
balance = pd.DataFrame(balance_rows)

print("\n=== balance check on the matching variables (should be near zero) ===")
print(f"{'comparison':12}{'feature':36}{'mean diff':>11}{'95% CI':>20}")
for row in balance.itertuples():
    print(f"{row.comparison:12}{row.feature:36}{row.mean_difference:11.3f}"
          f"  [{row.ci_low:7.3f},{row.ci_high:7.3f}]")

for name in COMPARISONS:
    subset = table[table.comparison == name].sort_values("q")
    if not len(subset):
        continue
    print(f"\n=== occupied minus matched control, {name} "
          f"({int(subset.n_pairs.max())} pairs, {int(subset.units.max())} units) ===")
    print(f"{'feature':34}{'mean diff':>11}{'std':>7}{'95% CI':>20}{'q':>9}")
    for row in subset.itertuples():
        flag = "*" if row.q < 0.05 else " "
        print(f"{row.feature:34}{row.mean_difference:11.3f}{row.standardised:7.2f}"
              f"  [{row.ci_low:7.3f},{row.ci_high:7.3f}]{row.q:9.4f}{flag}")

outdir = Path(args.out)
outdir.mkdir(parents=True, exist_ok=True)
table.to_csv(outdir / "context_paired_differences.csv", index=False)
balance.to_csv(outdir / "context_paired_balance.csv", index=False)
(outdir / "context_paired_summary.json").write_text(json.dumps({
    "amendment": "docs/prespecification_amendment_1_matched_pairs.md",
    "n_boot": args.boot, "git": _git_state()}, indent=2, default=str))
print(f"\ntable -> {(outdir/'context_paired_differences.csv').resolve()}")
