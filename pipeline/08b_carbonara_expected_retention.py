"""CARBonAra sequon retention: the exact figure, and 32 samples beside it.

CARBonAra is one-shot, so its designs are independent across positions and the
retention rate has a closed form. That closed form is not an approximation of
the sampler — it is exactly what sampling 32 sequences estimates, with none of
the Monte Carlo error. At a rate near 0.1, 32 draws carry roughly +/-0.05.

So the primary outcome here is `exact_full_sequon_retained`. The 32-design
columns are produced from the same forward pass and kept as a compatibility
check: they should agree with the exact value to within sampling error, and a
disagreement means the sampling path and the closed form have diverged.

**The rate is not comparable to ProteinMPNN's or ESM-IF's.** Those decode a
residue at a time, so their designs carry correlations this model cannot
express. What this supports is a within-CARBonAra occupancy-associated
difference in independent-marginal retention.

Usage:  08b_carbonara_expected_retention.py [manifest] [out] [--device DEV]
"""
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.carbonara_scoring import (
    design_distribution, design_sequences, expected_retention)
from experimental_glycosylation_sites.retention import (STANDARD_CONDITION,
                                                        classify_retention)
from experimental_glycosylation_sites.runner_support import (
    apply_shard, build_adapter, parse_args, resolve_device, structure_paths)

args = parse_args(sys.argv[1:],
                  "results/manifests/manifest_matched_secretory.csv",
                  "results/designs/retention_secretory_carbonara.csv",
                  description=__doc__)
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
N_DESIGNS, TEMP, SEED = STANDARD_CONDITION["n_designs"], 0.1, 0

MANIFEST, OUT = Path(args.manifest), Path(args.out)
FAILURES = OUT.with_name(OUT.stem + "_failures.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(MANIFEST, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
sites = manifest.drop_duplicates(KEY).reset_index(drop=True)
paths = structure_paths(tuple(args.structure_dir))

done = set()
if OUT.exists():                       # resumable
    prev = pd.read_csv(OUT, low_memory=False)
    done = set(map(tuple, prev[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already done", flush=True)

device = resolve_device(args.device)
adapter = build_adapter("carbonara", device)
provenance = adapter.describe()
generation = adapter.describe_generation()
print(f"carbonara ({provenance['model']}, {generation['generation']}) on {device}",
      flush=True)

groups = apply_shard(list(sites.groupby(["structure_pdb_id", "structure_chain_id"])),
                     args.shard)
print(f"{len(sites)} sites in {len(groups)} chains", flush=True)

rows, failures, t0 = [], [], time.time()
for gi, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    todo = [r for r in group.itertuples(index=False)
            if tuple(str(getattr(r, k)) for k in KEY) not in done]
    if not todo:
        continue
    path = paths.get(str(pdb_id).upper())
    if path is None:
        failures += [{**{k: getattr(r, k) for k in KEY},
                      "reason": "structure_not_cached"} for r in todo]
        continue

    # One unconditioned pass serves the exact figure and all 32 designs.
    try:
        mapping = adapter._mapping(path, chain_id)
        sharpened, usable = design_distribution(mapping, adapter.model, TEMP,
                                                carbonara_dir=adapter.dir)
        designs = design_sequences(mapping, adapter.model, n_designs=N_DESIGNS,
                                   temperature=TEMP, seed=SEED,
                                   distribution=(sharpened, usable))
    except Exception as exc:
        failures += [{**{k: getattr(r, k) for k in KEY},
                      "reason": f"{type(exc).__name__}: {str(exc)[:120]}"}
                     for r in todo]
        continue

    for r in todo:
        idx = (int(r.n_model_index), int(r.plus1_model_index),
               int(r.plus2_model_index))
        try:
            exact = expected_retention(sharpened, usable, mapping, idx)
        except Exception as exc:
            failures.append({**{k: getattr(r, k) for k in KEY},
                             "reason": f"{type(exc).__name__}: {str(exc)[:120]}"})
            continue
        sampled = classify_retention(designs, *idx)
        rows.append({**{k: getattr(r, k) for k in KEY}, "triplet": r.triplet,
                     "subtype": r.subtype, "temperature": TEMP, "seed": SEED,
                     "model": provenance["model"], **generation,
                     **exact,
                     **{f"sampled_{k}": v for k, v in sampled.items()}})

    if gi % 25 == 0:
        el = time.time() - t0
        print(f"  {gi}/{len(groups)} chains, {len(rows)} sites "
              f"({el/gi:.1f}s/chain, eta {(len(groups)-gi)*el/gi/60:.0f} min)",
              flush=True)
        if rows:
            pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(),
                                      index=False)
            rows = []

if rows:
    pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
frame = pd.read_csv(OUT, low_memory=False) if OUT.exists() else pd.DataFrame()
frame = frame.drop_duplicates(KEY)
frame.to_csv(OUT, index=False)
pd.DataFrame(failures).to_csv(FAILURES, index=False)

print(f"\n{len(frame)} of {len(sites)} sites; {len(failures)} failures")
if len(frame):
    exact = frame.exact_full_sequon_retained
    sampled = frame.sampled_frac_full_sequon_retained
    print(f"exact retention   : {exact.mean():.4f}")
    print(f"sampled ({N_DESIGNS} designs): {sampled.mean():.4f}   "
          f"mean |difference| {abs(exact - sampled).mean():.4f}")
    print("  the two should agree to within sampling error; a large gap means "
          "the sampler and the closed form have diverged")
print(f"elapsed {(time.time()-t0)/60:.1f} min")
