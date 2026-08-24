"""Step 2 — describe the natural occupied context. No comparison, no test.

The reference picture. Before asking whether occupied sites differ from anything,
say what they actually look like: how exposed the asparagine is, what secondary
structure it sits in, how crowded its attachment point is, what is nearby.

Deliberately descriptive. No classifier, no composite score, and no comparison
group -- those are steps 3 and 5. Coverage is printed beside every distribution,
because missingness here is not random: flexible and glycosylated regions are
preferentially unresolved.

Usage:
    47_context_describe.py [--views results/datasets] [--out results/analysis]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_stats import ramachandran_region_series
from experimental_glycosylation_sites.provenance import hash_file, _git_state

CONTINUOUS = ["n_rsa", "plus1_rsa", "plus2_rsa", "loop_run_length",
              "n_neighbours_8a", "nd2_atoms_8a_same_chain",
              "nd2_residues_8a_same_chain", "nd2_atoms_8a_other_chain",
              "sidechain_neighbour_residues_5a", "neighbour_net_charge_8a",
              "neighbour_hydrophobic_fraction_8a", "neighbour_aromatic_fraction_8a",
              "nearest_aromatic_sidechain_nd2", "uniprot_residues_after_asn",
              "uniprot_residues_after_sequon", "distance_to_n_terminus_resolved",
              "distance_to_c_terminus_resolved"]
CATEGORICAL = ["n_ss_coarse", "plus1_ss_coarse", "plus2_ss_coarse",
               "n_rama_region", "aromatic_within_8a", "subtype"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--views", default="results/datasets")
parser.add_argument("--out", default="results/analysis")
args = parser.parse_args()

core = pd.read_csv(Path(args.views) / "context_triplet_core.csv", low_memory=False)
core["n_rama_region"] = ramachandran_region_series(core.n_phi, core.n_psi)
occ = core[core.population == "occupied"].copy()

report = {"sites": int(len(occ)), "proteins": int(occ.accession.nunique())}
print(f"occupied sites in triplet_core: {len(occ)} across "
      f"{occ.accession.nunique()} proteins")
print(f"sequon subtype: {occ.subtype.value_counts().to_dict()}")
print(f"most sequons contributed by one protein: {occ.accession.value_counts().max()}")

print("\n=== continuous features (occupied sites) ===")
print(f"{'feature':34}{'cov%':>6}{'median':>9}{'IQR':>19}{'mean':>9}")
continuous = {}
for column in [c for c in CONTINUOUS if c in occ.columns]:
    values = occ[column].dropna().astype(float)
    coverage = 100 * len(values) / len(occ)
    if not len(values):
        continue
    q1, q2, q3 = values.quantile([0.25, 0.5, 0.75])
    continuous[column] = {"coverage_percent": round(coverage, 1),
                          "median": round(float(q2), 3),
                          "q1": round(float(q1), 3), "q3": round(float(q3), 3),
                          "mean": round(float(values.mean()), 3),
                          "sd": round(float(values.std(ddof=1)), 3), "n": int(len(values))}
    print(f"{column:34}{coverage:6.1f}{q2:9.2f}   [{q1:7.2f},{q3:7.2f}]{values.mean():9.2f}")
report["continuous"] = continuous

print("\n=== categorical features (occupied sites, % of non-missing) ===")
categorical = {}
for column in [c for c in CATEGORICAL if c in occ.columns]:
    counts = occ[column].value_counts(dropna=True)
    coverage = 100 * counts.sum() / len(occ)
    categorical[column] = {"coverage_percent": round(coverage, 1),
                           "counts": {str(k): int(v) for k, v in counts.items()}}
    shares = ", ".join(f"{k}={100*v/counts.sum():.1f}%" for k, v in counts.head(5).items())
    print(f"  {column:26} cov {coverage:5.1f}%   {shares}")
report["categorical"] = categorical

# Loop membership and length are separate facts: how often the Asn is in a loop,
# and how long that loop is *given* it is in one. Censored runs are a lower
# bound, so they are reported apart rather than averaged in.
print("\n=== loop membership and length ===")
in_loop = occ.loop_run_length.notna()
censored = occ.loop_run_censored.fillna(False).astype(bool)
lengths = occ.loc[in_loop & ~censored, "loop_run_length"].astype(float)
report["loop"] = {
    "asn_in_loop": int(in_loop.sum()),
    "asn_in_loop_percent": round(100 * float(in_loop.mean()), 1),
    "censored_runs": int((in_loop & censored).sum()),
    "uncensored_median_length": float(lengths.median()) if len(lengths) else None,
}
print(f"  Asn in a loop: {in_loop.sum()}/{len(occ)} ({100*in_loop.mean():.1f}%)")
print(f"  of those, censored at an unresolved boundary: {(in_loop & censored).sum()}")
print(f"  median length among uncensored runs: {lengths.median() if len(lengths) else 'n/a'}")

print("\n=== stratified by sequon subtype ===")
strata = {}
for subtype, group in occ.groupby("subtype"):
    strata[str(subtype)] = {c: round(float(group[c].median()), 3)
                            for c in ("n_rsa", "plus1_rsa", "plus2_rsa",
                                      "nd2_atoms_8a_same_chain")
                            if c in group.columns and group[c].notna().any()}
    print(f"  {subtype}: n={len(group)}  " +
          "  ".join(f"{k} median {v}" for k, v in strata[str(subtype)].items()))
report["by_subtype"] = strata

print("\n=== evidence tiers (sensitivity, described not tested) ===")
tiers = {}
for column in ("support_count", "glycan_modelled_at_site", "structure_choice"):
    if column in occ.columns:
        counts = occ[column].value_counts(dropna=False)
        tiers[column] = {str(k): int(v) for k, v in counts.items()}
        print(f"  {column:26}{tiers[column]}")
report["evidence_tiers"] = tiers

print("\n=== disulfide proximity: described within occupied only ===")
if "nearest_disulfide_sg_nd2" in occ.columns:
    values = occ.nearest_disulfide_sg_nd2.dropna()
    report["disulfide_within_occupied"] = {
        "coverage_percent": round(100 * len(values) / len(occ), 1),
        "median": round(float(values.median()), 2) if len(values) else None}
    print(f"  coverage {100*len(values)/len(occ):.1f}%, median "
          f"{values.median() if len(values) else 'n/a'} A "
          "(excluded from cross-arm tests by pre-specification)")

outdir = Path(args.out)
outdir.mkdir(parents=True, exist_ok=True)
source = Path(args.views) / "context_triplet_core.csv"
report["source"] = {str(source.resolve()): hash_file(source)}
report["git"] = _git_state()
path = outdir / "context_occupied_description.json"
path.write_text(json.dumps(report, indent=2, default=str))
print(f"\ndescription -> {path.resolve()}")
