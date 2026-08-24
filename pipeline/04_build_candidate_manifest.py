"""Map every candidate site to its structure, before any matching happens.

The earlier manifest was built from `matched_pairs.csv`, so a site had to be
matched before anyone asked whether it could be scored. That ordering lets
unscoreable sites into matched sets and removes them afterwards, which quietly
unbalances the sets that matching had just balanced. This runner takes the
whole candidate pool instead and knows nothing about matching.

Usage:  build_candidate_manifest.py [dataset|controls|secretory]
"""
import csv, gzip, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.input_paths import resolve_input
from experimental_glycosylation_sites.manifest import build_manifest

WHICH = sys.argv[1] if len(sys.argv) > 1 else "dataset"
OUT = Path(f"results/manifests/candidate_manifest_{WHICH}.csv")

seqs = {}
with gzip.open(resolve_input("raw/uniprot/uniprot_reviewed_glycoproteins_2026-04-27.tsv.gz"),
               "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        seqs[row["Entry"]] = row.get("Sequence", "")
for cache in ("data/cache/negative_control_proteins.csv.gz",
              "data/cache/secretory_unannotated_proteins.csv.gz"):
    if Path(cache).exists():
        ctrl = pd.read_csv(cache, low_memory=False)
        seqs.update(dict(zip(ctrl.Entry, ctrl.Sequence.fillna(""))))

paths = {}
for directory in ("data/cache/pdb",
                  "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    d = Path(directory)
    for p in list(d.glob("*.pdb")) + list(d.glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

if WHICH == "dataset":
    rows = pd.read_csv("results/datasets/site_structural_features.csv", low_memory=False)
    rows = rows[rows.features_available &
                rows.occupancy_status.isin(["occupied_supported", "observed_unmodified"])].copy()
    rows["control_set"] = pd.NA
    # ortholog cluster travels with the site: it is the resampling unit later
    assoc = pd.read_csv("results/datasets/site_pair_associations.csv", low_memory=False)
    cluster = (assoc.groupby(["accession", "position"]).cluster_id
                    .agg(lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x)})[:3]))
                    .rename("ortholog_clusters").reset_index())
    rows = rows.merge(cluster, on=["accession", "position"], how="left")
elif WHICH == "secretory":
    rows = pd.read_csv("results/datasets/secretory_unannotated_features.csv", low_memory=False)
    rows = rows[rows.features_available].copy()
    rows["ortholog_clusters"] = pd.NA
else:
    rows = pd.read_csv("results/datasets/negative_control_features.csv", low_memory=False)
    rows = rows[rows.features_available].copy()
    rows["ortholog_clusters"] = pd.NA

keep = ["accession", "position", "occupancy_status", "control_set", "ortholog_clusters",
        "structure_pdb_id", "structure_chain_id", "structure_resseq",
        "rsa", "n_neighbours_8a", "hydrophobic_fraction_8a"]
rows = rows[[c for c in keep if c in rows.columns]]

manifest, exclusions = build_manifest(rows, seqs, paths)
manifest.to_csv(OUT, index=False)
exclusions.to_csv(f"results/manifests/candidate_manifest_{WHICH}_exclusions.csv", index=False)

print(f"{WHICH}: {len(rows)} candidate sites -> {len(manifest)} mapped, {len(exclusions)} excluded")
assert len(manifest) + len(exclusions) == len(rows), "partition incomplete"
if len(exclusions):
    print(exclusions.exclusion_reason.value_counts().to_string())
if len(manifest) and "occupancy_status" in manifest:
    print("\nby class:"); print(manifest.occupancy_status.value_counts(dropna=False).to_string())
print("\nsubtype:"); print(manifest.subtype.value_counts().to_string())
