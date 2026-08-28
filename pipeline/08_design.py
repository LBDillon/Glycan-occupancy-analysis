"""Sequon retention in unconstrained designs, under one model.

The conditional score asks what probability a model holds at a site; this asks
what it writes when it actually generates a sequence. Nothing is fixed or
biased, so the sequon is free to disappear.

Designs are produced per chain and read at every sequon on that chain. Both
adapters return sequences indexed as the manifest indexes the chain, so
`classify_retention` is model-agnostic.

Usage:  08_design.py [manifest] [out] [--model NAME] [--device DEV]
"""
import sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.retention import (
    PREPRINT_CONDITION, STANDARD_CONDITION, classify_retention)
from experimental_glycosylation_sites.runner_support import (
    apply_shard, build_adapter, parse_args, resolve_device, structure_paths)

args = parse_args(sys.argv[1:],
                  "results/manifests/scoring_manifest.csv",
                  "results/designs/mpnn_retention.csv",
                  description=__doc__)
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
N_STD, N_PRE = STANDARD_CONDITION["n_designs"], PREPRINT_CONDITION["n_designs"]
TEMP, SEED = 0.1, 0

# How the designs were produced. Autoregressive sampling for the models that
# decode a residue at a time; an adapter whose procedure differs says so itself,
# so a retention table always records what generated it.
MANIFEST, OUT = Path(args.manifest), Path(args.out)
OUT.parent.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(MANIFEST, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
sites = manifest.drop_duplicates(KEY).reset_index(drop=True)
paths = structure_paths(tuple(args.structure_dir))

done = set()
if OUT.exists():
    done = set(map(tuple, pd.read_csv(OUT, low_memory=False)[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already done", flush=True)

device = resolve_device(args.device)
adapter = build_adapter(args.model, device, max_batch=args.max_batch)
provenance = adapter.describe()
GENERATION = (adapter.describe_generation()
              if hasattr(adapter, "describe_generation")
              else {"generation": "autoregressive_sampling"})
print(f"model {args.model} ({provenance['model']}) on {device}", flush=True)

groups = list(sites.groupby(["structure_pdb_id", "structure_chain_id"]))
groups = apply_shard(groups, args.shard)
print(f"{len(sites)} sites in {len(groups)} chains; {N_STD} designs each at T={TEMP}", flush=True)

rows, t0, failures = [], time.time(), []
for gi, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    todo = [r for r in group.itertuples(index=False)
            if tuple(str(getattr(r, k)) for k in KEY) not in done]
    if not todo:
        continue
    path = paths.get(str(pdb_id).upper())
    if path is None:
        continue
    try:
        designs = adapter.design(path, chain_id, n_designs=N_STD,
                                 temperature=TEMP, seed=SEED)
    except Exception as exc:
        print(f"  {pdb_id}/{chain_id}: {type(exc).__name__}: {str(exc)[:80]}", flush=True)
        failures += [{"accession": r.accession, "position": r.position,
                      "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                      "reason": f"{type(exc).__name__}: {str(exc)[:120]}"} for r in todo]
        continue
    for r in todo:
        idx = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
        std = classify_retention(designs, *idx)
        pre = classify_retention(designs[:N_PRE], *idx)
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, "temperature": TEMP, "seed": SEED,
                     "model": provenance["model"], **GENERATION,
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
if failures:
    pd.DataFrame(failures).to_csv(OUT.with_name(OUT.stem + "_failures.csv"), index=False)
print(f"\nretention recorded for {len(frame)} sites; {len(failures)} chain failures; "
      f"elapsed {(time.time()-t0)/60:.0f} min")
