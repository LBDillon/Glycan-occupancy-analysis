"""Which manifest sites the chosen model can actually decode.

Run before matching, not after. A site whose backbone is incomplete at any of
its three sequon residues cannot be scored, and ProteinMPNN signals this only
through its own mask: it returns a row of zeros that exponentiates to
twenty-one ones and scores near +13.8. Establishing this up front keeps
unscoreable sites out of the matching pool, so matched sets do not lose members
after they are balanced.

No model forward pass is needed — for either model, scoreability is a property
of the coordinates alone. ESM-IF adds a second reason a site can be unscoreable:
its parser (biotite) need not agree residue-for-residue with the parser the
manifest's indices came from (Biopython), so a manifest index with no ESM-IF
counterpart is reported here rather than discovered during scoring.

Usage:  05_scoreability.py [manifest] [out] [--model NAME] [--device DEV]
"""
import sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.runner_support import (
    build_adapter, parse_args, resolve_device, structure_paths)

args = parse_args(sys.argv[1:],
                  "results/manifests/scoring_manifest.csv",
                  "results/manifests/scoreability.csv",
                  description=__doc__)
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

MANIFEST, OUT = Path(args.manifest), Path(args.out)
sites = pd.read_csv(MANIFEST, low_memory=False).drop_duplicates(KEY).reset_index(drop=True)
paths = structure_paths(tuple(args.structure_dir))

device = resolve_device(args.device)
adapter = build_adapter(args.model, device, mask_mode=args.mask_mode)
print(f"model {args.model} on {device}", flush=True)

rows, cache, t0 = [], {}, time.time()
groups = list(sites.groupby(["structure_pdb_id", "structure_chain_id"]))
print(f"{len(sites)} sites in {len(groups)} chains", flush=True)

for gi, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    path = paths.get(str(pdb_id).upper())
    key = (str(pdb_id).upper(), chain_id)
    if key not in cache:
        if path is None:
            cache[key] = ("structure_not_cached", None)
        else:
            try:
                cache[key] = (None, adapter.decodable_positions(path, chain_id))
            except Exception as exc:
                cache[key] = (f"{type(exc).__name__}: {str(exc)[:60]}", None)
    reason, decodable = cache[key]

    for r in group.itertuples(index=False):
        idx = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
        if decodable is None:
            rows.append({**{k: getattr(r, k) for k in KEY},
                         "scoreable": False, "reason": reason})
            continue
        bad = [i for i in idx if i >= len(decodable) or not bool(decodable[i])]
        rows.append({**{k: getattr(r, k) for k in KEY},
                     "scoreable": not bad,
                     "reason": "" if not bad else
                               f"incomplete_backbone_at_model_index_{bad}"})
    if gi % 200 == 0:
        el = time.time() - t0
        print(f"  {gi}/{len(groups)} chains ({el/gi:.2f}s/chain, "
              f"eta {(len(groups)-gi)*el/gi/60:.0f} min)", flush=True)

out = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)
print(f"\n{int(out.scoreable.sum())} of {len(out)} sites scoreable")
print(out[~out.scoreable].reason.str.replace(r"_\[.*", "", regex=True).value_counts().to_string())
print(f"elapsed {(time.time()-t0)/60:.1f} min")
