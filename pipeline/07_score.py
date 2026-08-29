"""Conditional sequon scores for one manifest, under one model.

The manifest and destination are arguments so the same scorer serves the dataset
sites and every control pool without a second copy of this file; `--model`
serves a second model the same way. Scoring is grouped by chain because both
supported models do their expensive work per chain: ProteinMPNN runs one decoder
pass per requested position, ESM-IF one teacher-forced pass for the whole chain.

Usage:  07_score.py [manifest] [out] [--model NAME] [--device DEV]
"""
import json, sys, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.runner_support import (
    apply_shard, build_adapter, needs_header, parse_args, read_resumable_csv,
    resolve_device, structure_paths)

args = parse_args(sys.argv[1:],
                  "results/manifests/scoring_manifest.csv",
                  "results/scores/mpnn_conditional_scores.csv",
                  description=__doc__)
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

MANIFEST, OUT = Path(args.manifest), Path(args.out)
FAILURES = OUT.with_name(OUT.stem + "_failures.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(MANIFEST, low_memory=False)
if "scoreable" in manifest.columns:
    # Unscoreable sites are excluded before matching, so reaching the scorer at
    # all would mean the pipeline ran out of order.
    manifest = manifest[manifest.scoreable.astype(bool)]
sites = manifest.drop_duplicates(KEY).reset_index(drop=True)

paths = structure_paths(tuple(args.structure_dir))

done = set()
if OUT.exists():                       # resumable
    prev = read_resumable_csv(OUT, empty_columns=KEY)
    done = set(map(tuple, prev[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already scored", flush=True)

device = resolve_device(args.device)
adapter = build_adapter(args.model, device, mask_mode=args.mask_mode,
                        structure_mode=args.structure_mode)
provenance = adapter.describe()
print(f"model {args.model} ({provenance['model']}, {provenance['conditioning']}) "
      f"on {device}", flush=True)

groups = list(sites.groupby(["structure_pdb_id", "structure_chain_id"]))
groups = apply_shard(groups, args.shard)
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
                      "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                      "reason": "structure_not_cached"} for r in todo]
        continue

    # every sequon residue needed from this chain, prepared in one go
    wanted = sorted({int(v) for r in todo for v in
                     (r.n_model_index, r.plus1_model_index, r.plus2_model_index)})
    try:
        context = adapter.prepare_chain(path, chain_id, wanted)
    except Exception as exc:
        failures += [{"accession": r.accession, "position": r.position,
                      "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                      "reason": f"{type(exc).__name__}: {str(exc)[:80]}"} for r in todo]
        continue

    for r in todo:
        # Per site, not per chain: one residue the model declined to evaluate
        # must not discard the other sites read from the same chain.
        try:
            sc = adapter.score_from(context,
                                    (int(r.n_model_index), int(r.plus1_model_index),
                                     int(r.plus2_model_index)),
                                    expected_triplet=r.triplet)
        except Exception as exc:
            failures.append({"accession": r.accession, "position": r.position,
                             "structure_pdb_id": pdb_id,
                             "structure_chain_id": chain_id,
                             "reason": f"{type(exc).__name__}: {str(exc)[:120]}"})
            continue
        vectors = {k: json.dumps([round(x, 6) for x in sc.pop(k)])
                   for k in ("probs_n", "probs_plus1", "probs_plus2")}
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, **provenance, **sc, **vectors})

    if gi % 25 == 0:
        el = time.time() - t0
        print(f"  {gi}/{len(groups)} groups, {len(rows)} sites "
              f"({el/gi:.2f}s/group, eta {(len(groups)-gi)*el/gi/60:.0f} min)", flush=True)
        if rows:
            # append-only checkpoint: the run is genuinely resumable only if
            # completed work reaches disk before the end.
            pd.DataFrame(rows).to_csv(
                OUT, mode="a", header=needs_header(OUT), index=False)
            rows = []

if rows:
    pd.DataFrame(rows).to_csv(OUT, mode="a", header=needs_header(OUT), index=False)
frame = read_resumable_csv(OUT, empty_columns=KEY)
frame = frame.drop_duplicates(KEY)
frame.to_csv(OUT, index=False)
# Reindexed so the columns are the same whether there were failures, one
# kind of failure, or several. A table whose schema depends on which errors
# occurred cannot be joined against reliably.
pd.DataFrame(failures, columns=KEY + ["reason"]).to_csv(FAILURES, index=False)

print(f"\nscored {len(frame)} of {len(sites)} unique sites; {len(failures)} failures")
if failures:
    print(pd.DataFrame(failures).reason.value_counts().to_string())
print(f"elapsed {(time.time()-t0)/60:.0f} min")
