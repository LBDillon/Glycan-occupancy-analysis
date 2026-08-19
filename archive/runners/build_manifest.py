import csv, gzip, sys, json, pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.manifest import build_manifest, within_structure_pairs

pairs = pd.read_csv("results/matched_pairs.csv", low_memory=False)
dfeat = pd.read_csv("results/site_structural_features.csv", low_memory=False)
cfeat = pd.read_csv("results/negative_control_features.csv", low_memory=False)
assoc = pd.read_csv("results/site_pair_associations.csv", low_memory=False)

# sequences: dataset from the UniProt snapshot, controls from the control cache
seqs = {}
with gzip.open("../../data/raw/uniprot/uniprot_reviewed_glycoproteins_2026-04-27.tsv.gz",
               "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        seqs[row["Entry"]] = row.get("Sequence", "")
ctrl_prot = pd.read_csv("data/cache/negative_control_proteins.csv.gz", low_memory=False)
seqs.update(dict(zip(ctrl_prot.Entry, ctrl_prot.Sequence.fillna(""))))

# Dataset sites were mapped against the ortholog database's structure cache;
# control sites against this module's. Both must be indexed or the occupied side
# vanishes into structure_not_cached.
paths = {}
for directory in ("data/cache/pdb",
                  "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    d = Path(directory)
    for p in list(d.glob("*.pdb")) + list(d.glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)
print(f"structure files indexed: {len(paths)}")

cluster = (assoc.groupby(["accession", "position"]).cluster_id
                .agg(lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x)})[:3]))
                .rename("ortholog_clusters").reset_index())

# one row per (comparison, matched set, role)
rows = []
for r in pairs.itertuples(index=False):
    rows.append({"comparison": r.comparison,
                 "matched_set_id": f"{r.comparison}:{r.case_accession}:{r.case_position}",
                 "role": "occupied", "accession": r.case_accession, "position": r.case_position,
                 "match_distance": None, "match_rank": 0})
    rows.append({"comparison": r.comparison,
                 "matched_set_id": f"{r.comparison}:{r.case_accession}:{r.case_position}",
                 "role": "control", "accession": r.control_accession, "position": r.control_position,
                 "match_distance": r.distance, "match_rank": r.match_rank})
rows = pd.DataFrame(rows).drop_duplicates(["comparison", "matched_set_id", "role", "accession", "position"])

feat = pd.concat([
    dfeat.assign(control_set=pd.NA),
    cfeat,
], ignore_index=True)[["accession", "position", "occupancy_status", "control_set",
                       "structure_pdb_id", "structure_chain_id", "structure_resseq"]]
rows = rows.merge(feat, on=["accession", "position"], how="left")
rows = rows.merge(cluster, on=["accession", "position"], how="left")

manifest, exclusions = build_manifest(rows, seqs, paths)
manifest.to_csv("results/scoring_manifest.csv", index=False)
exclusions.to_csv("results/scoring_manifest_exclusions.csv", index=False)

print(f"input rows: {len(rows)}   manifest: {len(manifest)}   excluded: {len(exclusions)}")
assert len(manifest) + len(exclusions) == len(rows), "partition incomplete"
print("\nexclusion reasons:")
if len(exclusions):
    print(exclusions.exclusion_reason.value_counts().to_string())
print("\nmanifest by comparison and role:")
print(manifest.groupby(["comparison", "role"]).size().to_string())
print("\nsequon subtype:")
print(manifest.subtype.value_counts().to_string())
print("\nmiddle X residue (top 8):")
print(manifest.plus1_aa.value_counts().head(8).to_string())

sub = within_structure_pairs(manifest)
sub.to_csv("results/within_structure_subset.csv", index=False)
print(f"\nwithin-protein/structure occupied-vs-unmodified pairs: {len(sub)}")
if len(sub):
    print(sub.groupby("shared_by").size().to_string())
