"""Score the dataset sites that never entered a matched set.

The frozen config defines the reference-SD population as all structurally
scoreable dataset sites, not just the matched ones, so the remainder must be
scored for the equivalence margin to follow the rule as written.
"""
import csv, gzip, json, sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.manifest import build_manifest
from experimental_glycosylation_sites.mpnn_scoring import (
    DEFAULT_MODEL, load_model, conditional_probabilities, sequon_score)

rest = pd.read_csv("/tmp/unmatched_dataset_sites.csv", low_memory=False)
seqs = {}
with gzip.open("../../data/raw/uniprot/uniprot_reviewed_glycoproteins_2026-04-27.tsv.gz",
               "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        seqs[row["Entry"]] = row.get("Sequence", "")

paths = {}
for d in ("data/cache/pdb", "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

man, exc = build_manifest(rest, seqs, paths)
print(f"{len(rest)} unmatched sites -> {len(man)} mappable, {len(exc)} excluded")
if len(exc):
    print(exc.exclusion_reason.value_counts().to_string())

model = load_model(Path("../../ProteinMPNN"))
rows = []
for r in man.itertuples(index=False):
    path = paths.get(str(r.structure_pdb_id).upper())
    idx = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
    try:
        probs = conditional_probabilities(path, r.structure_chain_id, model,
                                          n_decoding_orders=8, seed=0, positions=list(idx))
    except Exception as exc_:
        print(f"  {r.accession} {r.structure_pdb_id}: {type(exc_).__name__}", flush=True)
        continue
    sc = sequon_score(probs, *idx)
    vec = {k: json.dumps([round(x, 6) for x in sc.pop(k)])
           for k in ("probs_n", "probs_plus1", "probs_plus2")}
    rows.append({"accession": r.accession, "position": r.position,
                 "structure_pdb_id": r.structure_pdb_id,
                 "structure_chain_id": r.structure_chain_id,
                 "triplet": r.triplet, "subtype": r.subtype, "model": DEFAULT_MODEL,
                 "conditioning": "conditional", "n_orders": 8, "seed": 0, **sc, **vec})

out = pd.DataFrame(rows)
out.to_csv("results/mpnn_conditional_scores_unmatched.csv", index=False)
print(f"scored {len(out)} additional dataset sites")
