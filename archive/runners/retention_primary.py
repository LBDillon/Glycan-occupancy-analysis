"""Retention for the primary comparison only: the 22 matched occupied/unmodified pairs.

The main sweep walks chains alphabetically, so these scattered sites would not
surface for hours. They are the only ones the primary retention contrast needs.
Written to a separate file and merged later, so neither run can clobber the other.
"""
import sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.mpnn_scoring import load_model
from experimental_glycosylation_sites.retention import (
    PREPRINT_CONDITION, STANDARD_CONDITION, classify_retention, design_sequences)

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
OUT = Path("results/mpnn_retention_primary.csv")
N_STD, N_PRE = STANDARD_CONDITION["n_designs"], PREPRINT_CONDITION["n_designs"]

manifest = pd.read_csv("results/scoring_manifest.csv", low_memory=False)
pairs = pd.read_csv("results/matched_pairs.csv", low_memory=False)
prim = pairs[pairs.comparison == "vs_observed_unmodified"]
wanted = set(zip(prim.case_accession, prim.case_position)) | \
         set(zip(prim.control_accession, prim.control_position))
sites = manifest[[(a, p) in wanted for a, p in zip(manifest.accession, manifest.position)]]
sites = sites.drop_duplicates(KEY).reset_index(drop=True)

already = set()
for f in ("results/mpnn_retention.csv", str(OUT)):
    if Path(f).exists():
        already |= set(map(tuple, pd.read_csv(f, low_memory=False)[KEY].astype(str).values))
sites = sites[[tuple(map(str, r)) not in already for r in sites[KEY].values]]

paths = {}
for d in ("data/cache/pdb", "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

print(f"{len(sites)} primary-comparison sites still needing retention", flush=True)
model = load_model(Path("../../ProteinMPNN"))
rows, t0 = [], time.time()
for (pdb_id, chain_id), group in sites.groupby(["structure_pdb_id", "structure_chain_id"]):
    path = paths.get(str(pdb_id).upper())
    if path is None:
        continue
    try:
        designs = design_sequences(path, chain_id, model, n_designs=N_STD,
                                   temperature=0.1, seed=0)
    except Exception as exc:
        print(f"  {pdb_id}/{chain_id}: {type(exc).__name__}", flush=True)
        continue
    for r in group.itertuples(index=False):
        idx = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
        std, pre = classify_retention(designs, *idx), classify_retention(designs[:N_PRE], *idx)
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, "temperature": 0.1, "seed": 0,
                     **{f"std_{k}": v for k, v in std.items()},
                     **{f"pre_{k}": v for k, v in pre.items()}})
    print(f"  {len(rows)}/{len(sites)} ({time.time()-t0:.0f}s)", flush=True)

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"\nwrote {len(rows)} sites in {(time.time()-t0)/60:.0f} min")
