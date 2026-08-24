"""Step 3 — occupancy-associated context differences, as pre-specified.

Compares occupied sites against the two comparison sets *separately*, following
docs/prespecification_2026-08-24_context_differences.md. The two are never
pooled: the internal controls carry the better label and almost no power, the
secretory set the reverse, and each is informative exactly where the other is
weak.

Bacterial and cytosolic are reported as diagnostics and excluded from the
confirmatory family. They differ from occupied sites for compartment and
composition reasons unrelated to occupancy, and including them would
mis-attribute a confound to glycosylation.

Every interval resamples proteins rather than rows, because one protein
contributes up to 7 occupied sequons and up to 19 secretory ones.

Usage:
    48_context_differences.py [--views results/datasets] [--boot 2000]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_stats import (benjamini_hochberg,
                                                            cluster_bootstrap_difference,
                                                            ramachandran_region_series)
from experimental_glycosylation_sites.provenance import hash_file, _git_state

CONTINUOUS = ["n_rsa", "plus1_rsa", "plus2_rsa", "loop_run_length",
              "n_neighbours_8a", "nd2_atoms_8a_same_chain",
              "nd2_residues_8a_same_chain", "nd2_atoms_8a_other_chain",
              "sidechain_neighbour_residues_5a", "neighbour_net_charge_8a",
              "neighbour_hydrophobic_fraction_8a", "neighbour_aromatic_fraction_8a",
              "nearest_aromatic_sidechain_nd2", "uniprot_residues_after_asn",
              "uniprot_residues_after_sequon", "distance_to_n_terminus_resolved",
              "distance_to_c_terminus_resolved"]
CATEGORICAL = {"n_ss_coarse": ("loop", "helix", "sheet"),
               "plus1_ss_coarse": ("loop", "helix", "sheet"),
               "plus2_ss_coarse": ("loop", "helix", "sheet"),
               "aromatic_within_8a": (True,)}
CONFIRMATORY = ("internal_control", "secretory_unannotated")
DIAGNOSTIC = ("bacterial", "cytosolic")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--views", default="results/datasets")
parser.add_argument("--diagnostic",
                    default="results/datasets/context_features_diagnostic.csv")
parser.add_argument("--boot", type=int, default=2000)
parser.add_argument("--out", default="results/analysis")
args = parser.parse_args()

core = pd.read_csv(Path(args.views) / "context_triplet_core.csv", low_memory=False)
frame = core
diagnostic_path = Path(args.diagnostic)
if diagnostic_path.exists():
    diag = pd.read_csv(diagnostic_path, low_memory=False)
    keep = (diag.triplet_matches.fillna(False).astype(bool)
            & diag.mapping_continuous.fillna(False).astype(bool))
    frame = pd.concat([core, diag[keep]], ignore_index=True)
    print(f"diagnostic populations added: {int(keep.sum())} of {len(diag)} rows "
          "meeting triplet_core criteria")
frame["n_rama_region"] = ramachandran_region_series(frame.n_phi, frame.n_psi)

rows = []
for comparison in CONFIRMATORY + DIAGNOSTIC:
    if comparison not in set(frame.population):
        continue
    for feature in [c for c in CONTINUOUS if c in frame.columns]:
        result = cluster_bootstrap_difference(
            frame, feature, "population", "occupied", comparison,
            n_boot=args.boot, seed=11, statistic="smd")
        result["family"] = ("confirmatory" if comparison in CONFIRMATORY
                            else "diagnostic")
        rows.append(result)
    for column, levels in CATEGORICAL.items():
        if column not in frame.columns:
            continue
        for level in levels:
            indicator = f"{column}=={level}"
            work = frame.assign(**{indicator: (frame[column] == level).astype(float)})
            work.loc[frame[column].isna(), indicator] = float("nan")
            result = cluster_bootstrap_difference(
                work, indicator, "population", "occupied", comparison,
                n_boot=args.boot, seed=11, statistic="mean_difference")
            result["family"] = ("confirmatory" if comparison in CONFIRMATORY
                                else "diagnostic")
            rows.append(result)

table = pd.DataFrame(rows)
# BH within each comparison, across that comparison's whole family.
table["q"] = float("nan")
for comparison, group in table.groupby("comparison"):
    table.loc[group.index, "q"] = benjamini_hochberg(group.p.tolist())

outdir = Path(args.out)
outdir.mkdir(parents=True, exist_ok=True)
table.to_csv(outdir / "context_differences.csv", index=False)

def show(subset, title):
    print(f"\n=== {title} ===")
    print(f"{'feature':34}{'est':>7}{'95% CI':>18}{'q':>9}{'n_occ':>7}{'n_cmp':>7}")
    for row in subset.sort_values("q").itertuples():
        flag = "*" if row.q < 0.05 else " "
        print(f"{row.feature:34}{row.estimate:7.2f}"
              f"  [{row.ci_low:6.2f},{row.ci_high:6.2f}]{row.q:9.4f}{flag}"
              f"{row.n_occupied:6d}{row.n_comparison:7d}")

for comparison in CONFIRMATORY:
    subset = table[table.comparison == comparison]
    if len(subset):
        show(subset, f"occupied vs {comparison} (confirmatory)")

# Direction agreement is the pre-specified reading, so it is computed rather
# than left to the eye.
print("\n=== direction agreement between the two confirmatory comparisons ===")
pivot = table[table.family == "confirmatory"].pivot_table(
    index="feature", columns="comparison", values=["estimate", "q"], aggfunc="first")
agree = []
for feature in pivot.index:
    try:
        e_int = pivot.loc[feature, ("estimate", "internal_control")]
        e_sec = pivot.loc[feature, ("estimate", "secretory_unannotated")]
        q_int = pivot.loc[feature, ("q", "internal_control")]
        q_sec = pivot.loc[feature, ("q", "secretory_unannotated")]
    except KeyError:
        continue
    if pd.isna(e_int) or pd.isna(e_sec):
        continue
    same = (e_int > 0) == (e_sec > 0)
    reading = ("both significant, same direction" if same and q_int < 0.05 and q_sec < 0.05
               else "secretory only" if q_sec < 0.05 and q_int >= 0.05
               else "internal only" if q_int < 0.05 and q_sec >= 0.05
               else "opposite directions" if not same
               else "neither")
    agree.append({"feature": feature, "internal": round(float(e_int), 3),
                  "secretory": round(float(e_sec), 3), "reading": reading})
agreement = pd.DataFrame(agree).sort_values("reading")
print(agreement.to_string(index=False))

for comparison in DIAGNOSTIC:
    subset = table[table.comparison == comparison]
    if len(subset):
        strongest = subset.reindex(subset.estimate.abs().sort_values(ascending=False).index).head(5)
        print(f"\n=== {comparison} (diagnostic, NOT a test of occupancy) ===")
        for row in strongest.itertuples():
            print(f"  {row.feature:34}{row.estimate:7.2f}  q={row.q:.4f}")

record = {
    "prespecification": "docs/prespecification_2026-08-24_context_differences.md",
    "n_boot": args.boot,
    "direction_agreement": agree,
    "git": _git_state(),
}
(outdir / "context_differences_summary.json").write_text(json.dumps(record, indent=2, default=str))
agreement.to_csv(outdir / "context_direction_agreement.csv", index=False)
print(f"\ntable -> {(outdir/'context_differences.csv').resolve()}")
