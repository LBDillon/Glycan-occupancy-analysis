import json, sys, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.mpnn_scoring import (
    DEFAULT_MODEL, IncompleteBackboneError, InvalidProbabilityVector,
    load_model, conditional_probabilities, sequon_score)

N_ORDERS, SEED = 8, 0
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

# Manifest and destination are arguments so the same scorer serves the dataset
# sites and the diagnostic control pools without a second copy of this file.
MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/scoring_manifest.csv")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results/mpnn_conditional_scores.csv")
FAILURES = OUT.with_name(OUT.stem + "_failures.csv")

manifest = pd.read_csv(MANIFEST, low_memory=False)
if "scoreable" in manifest.columns:
    # Unscoreable sites are excluded before matching, so reaching the scorer at
    # all would mean the pipeline ran out of order.
    manifest = manifest[manifest.scoreable.astype(bool)]
sites = manifest.drop_duplicates(KEY).reset_index(drop=True)

paths = {}
for d in ("data/cache/pdb", "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

done = set()
if OUT.exists():                       # resumable
    prev = pd.read_csv(OUT, low_memory=False)
    done = set(map(tuple, prev[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already scored", flush=True)

model = load_model(Path("../../ProteinMPNN"), model_name=DEFAULT_MODEL)
groups = list(sites.groupby(["structure_pdb_id", "structure_chain_id"]))
print(f"{len(sites)} unique sites in {len(groups)} structure-chain groups", flush=True)

rows, failures, t0 = [], [], time.time()
for gi, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    todo = [r for r in group.itertuples(index=False)
            if tuple(str(getattr(r, k)) for k in KEY) not in done]
    if not todo:
        continue
    path = paths.get(str(pdb_id).upper())
    if path is None:
        failures += [{"accession": r.accession, "position": r.position,
                      "structure_pdb_id": pdb_id, "reason": "structure_not_cached"} for r in todo]
        continue

    # every sequon residue needed from this chain, computed in one pass
    wanted = sorted({int(v) for r in todo for v in
                     (r.n_model_index, r.plus1_model_index, r.plus2_model_index)})
    try:
        probs, computed = conditional_probabilities(path, chain_id, model,
                                                    n_decoding_orders=N_ORDERS, seed=SEED,
                                                    positions=wanted)
    except Exception as exc:
        failures += [{"accession": r.accession, "position": r.position,
                      "structure_pdb_id": pdb_id,
                      "reason": f"{type(exc).__name__}: {str(exc)[:80]}"} for r in todo]
        continue

    for r in todo:
        # Per site, not per chain: one residue with an incomplete backbone must
        # not discard the other sites decoded from the same chain.
        try:
            sc = sequon_score(probs, int(r.n_model_index), int(r.plus1_model_index),
                              int(r.plus2_model_index), computed=computed)
        except (IncompleteBackboneError, InvalidProbabilityVector) as exc:
            failures.append({"accession": r.accession, "position": r.position,
                             "structure_pdb_id": pdb_id,
                             "structure_chain_id": chain_id,
                             "reason": f"{type(exc).__name__}: {str(exc)[:120]}"})
            continue
        vectors = {k: json.dumps([round(x, 6) for x in sc.pop(k)])
                   for k in ("probs_n", "probs_plus1", "probs_plus2")}
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, "model": DEFAULT_MODEL,
                     "conditioning": "conditional", "n_orders": N_ORDERS, "seed": SEED,
                     **sc, **vectors})

    if gi % 25 == 0:
        el = time.time() - t0
        print(f"  {gi}/{len(groups)} groups, {len(rows)} sites "
              f"({el/gi:.2f}s/group, eta {(len(groups)-gi)*el/gi/60:.0f} min)", flush=True)
        if rows:
            # append-only checkpoint: the run is genuinely resumable only if
            # completed work reaches disk before the end.
            pd.DataFrame(rows).to_csv(
                OUT, mode="a", header=not OUT.exists(), index=False)
            rows = []

if rows:
    pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
frame = pd.read_csv(OUT, low_memory=False) if OUT.exists() else pd.DataFrame()
frame = frame.drop_duplicates(KEY)
frame.to_csv(OUT, index=False)
pd.DataFrame(failures).to_csv(FAILURES, index=False)

print(f"\nscored {len(frame)} of {len(sites)} unique sites; {len(failures)} failures")
if failures:
    print(pd.DataFrame(failures).reason.value_counts().to_string())
print(f"elapsed {(time.time()-t0)/60:.0f} min")
