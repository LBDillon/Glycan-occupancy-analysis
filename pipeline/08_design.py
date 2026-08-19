import sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.mpnn_scoring import load_model
from experimental_glycosylation_sites.retention import (
    PREPRINT_CONDITION, STANDARD_CONDITION, classify_retention, design_sequences)

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

# Manifest and destination are arguments so the same sweep serves each control
# population without a second copy of this file.
MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/manifests/scoring_manifest.csv")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results/designs/mpnn_retention.csv")
N_STD, N_PRE, TEMP, SEED = STANDARD_CONDITION["n_designs"], PREPRINT_CONDITION["n_designs"], 0.1, 0

manifest = pd.read_csv(MANIFEST, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
sites = manifest.drop_duplicates(KEY).reset_index(drop=True)
paths = {}
for d in ("data/cache/pdb", "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

done = set()
if OUT.exists():
    done = set(map(tuple, pd.read_csv(OUT, low_memory=False)[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already done", flush=True)

model = load_model(Path("../../ProteinMPNN"))
groups = list(sites.groupby(["structure_pdb_id", "structure_chain_id"]))
print(f"{len(sites)} sites in {len(groups)} chains; {N_STD} designs each at T={TEMP}", flush=True)

rows, t0 = [], time.time()
for gi, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    todo = [r for r in group.itertuples(index=False)
            if tuple(str(getattr(r, k)) for k in KEY) not in done]
    if not todo:
        continue
    path = paths.get(str(pdb_id).upper())
    if path is None:
        continue
    try:
        designs = design_sequences(path, chain_id, model, n_designs=N_STD,
                                   temperature=TEMP, seed=SEED)
    except Exception as exc:
        print(f"  {pdb_id}/{chain_id}: {type(exc).__name__}", flush=True)
        continue
    for r in todo:
        idx = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
        std = classify_retention(designs, *idx)
        pre = classify_retention(designs[:N_PRE], *idx)
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, "temperature": TEMP, "seed": SEED,
                     **{f"std_{k}": v for k, v in std.items()},
                     **{f"pre_{k}": v for k, v in pre.items()}})
    if gi % 25 == 0:
        el = time.time() - t0
        print(f"  {gi}/{len(groups)} chains, {len(rows)} sites "
              f"({el/gi:.1f}s/chain, eta {(len(groups)-gi)*el/gi/60:.0f} min)", flush=True)
        if rows:
            pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
            rows = []

if rows:
    pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
frame = pd.read_csv(OUT, low_memory=False).drop_duplicates(KEY)
frame.to_csv(OUT, index=False)
print(f"\nretention recorded for {len(frame)} sites; elapsed {(time.time()-t0)/60:.0f} min")
