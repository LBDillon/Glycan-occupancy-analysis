"""Which manifest sites ProteinMPNN can actually decode.

Run before matching, not after. A site whose backbone is incomplete at any of
its three sequon residues cannot be scored, and ProteinMPNN signals this only
through its own mask: it returns a row of zeros that exponentiates to
twenty-one ones and scores near +13.8. Establishing this up front keeps
unscoreable sites out of the matching pool, so matched sets do not lose members
after they are balanced.

No model forward pass is needed — the mask is a property of the coordinates.
"""
import sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.mpnn_scoring import decodable_positions

MPNN_DIR = Path("../../ProteinMPNN")

MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/scoring_manifest.csv")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results/scoreability.csv")
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

sites = pd.read_csv(MANIFEST, low_memory=False).drop_duplicates(KEY).reset_index(drop=True)

paths = {}
for d in ("data/cache/pdb",
          "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

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
                cache[key] = (None, decodable_positions(path, chain_id, MPNN_DIR))
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
out.to_csv(OUT, index=False)
print(f"\n{int(out.scoreable.sum())} of {len(out)} sites scoreable")
print(out[~out.scoreable].reason.str.replace(r"_\[.*", "", regex=True).value_counts().to_string())
print(f"elapsed {(time.time()-t0)/60:.1f} min")
