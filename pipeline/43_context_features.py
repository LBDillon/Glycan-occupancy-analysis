"""Extract context features v2 over the manifest's populations.

Resumable and shardable, like the model stages: `--shard K/N` splits by chain so
a chain's DSSP runs once in one task, and rerunning skips what is already on
disk. Grouped by chain because DSSP and the SASA calculation are per structure.

`dssp_ok` and `dssp_reason` are written for every site, never dropped. DSSP
coverage is not equal between arms — about 95% for occupied sites against 89%
for the broader control sets, because failures concentrate in large structures
and the control sets draw more of them — so any analysis using secondary
structure has to report coverage per population rather than filtering silently.

Usage:
    43_context_features.py [--populations occupied,internal_control,...]
                           [--out results/datasets/context_features.csv]
                           [--shard K/N]
"""
import argparse, sys, time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_features import (FEATURE_COLUMNS,
                                                               QC_COLUMNS,
                                                               sequon_context)
from experimental_glycosylation_sites.runner_support import (apply_shard,
                                                             structure_paths)

KEY = ["accession", "position", "population"]
DEFAULT_POPULATIONS = "occupied,internal_control,secretory_unannotated"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", default="results/datasets/context_manifest.csv")
parser.add_argument("--out", default="results/datasets/context_features.csv")
parser.add_argument("--populations", default=DEFAULT_POPULATIONS,
                    help="comma-separated; 'all' for every population")
parser.add_argument("--shard", default=None, metavar="K/N")
parser.add_argument("--structure-dir", action="append", default=[])
args = parser.parse_args()

manifest = pd.read_csv(args.manifest, low_memory=False)
manifest = manifest[manifest.features_available.astype(bool)]
if args.populations != "all":
    wanted = [p.strip() for p in args.populations.split(",")]
    manifest = manifest[manifest.population.isin(wanted)]
manifest = manifest.dropna(subset=["structure_pdb_id", "structure_chain_id",
                                   "structure_resseq"])

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
done = set()
if out.exists():
    previous = pd.read_csv(out, low_memory=False)
    done = set(map(tuple, previous[KEY].astype(str).values))
    print(f"resuming: {len(done)} sites already extracted", flush=True)

paths = structure_paths(tuple(args.structure_dir))
groups = list(manifest.groupby(["structure_pdb_id", "structure_chain_id"]))
groups = apply_shard(groups, args.shard)
print(f"{len(manifest)} sites in {len(groups)} chains "
      f"({manifest.population.value_counts().to_dict()})", flush=True)

rows, failures, t0 = [], [], time.time()
for index, ((pdb_id, chain_id), group) in enumerate(groups, 1):
    todo = [r for r in group.itertuples(index=False)
            if tuple(str(getattr(r, k)) for k in KEY) not in done]
    if not todo:
        continue
    path = paths.get(str(pdb_id).upper())
    if path is None:
        failures += [{**{k: getattr(r, k) for k in KEY},
                      "reason": "structure_not_cached"} for r in todo]
        continue
    for r in todo:
        try:
            features = sequon_context(path, str(chain_id), int(r.structure_resseq),
                                      icode=getattr(r, "structure_icode", ""),
                                      pdb_id=str(pdb_id))
        except Exception as exc:
            features = None
            failures.append({**{k: getattr(r, k) for k in KEY},
                             "reason": f"{type(exc).__name__}: {str(exc)[:80]}"})
        if features is None:
            failures.append({**{k: getattr(r, k) for k in KEY},
                             "reason": "site_not_locatable"})
            continue
        # The manifest's triplet is what the site is SUPPOSED to be. A mismatch
        # means the structure mapping is wrong, so it is recorded rather than
        # scored -- the same guard the model adapters use.
        features["triplet_expected"] = getattr(r, "sequon_triplet", "")
        features["triplet_matches"] = (
            features.get("triplet_observed") == features["triplet_expected"])
        rows.append({**{k: getattr(r, k) for k in KEY},
                     "subtype": getattr(r, "subtype", ""), **features})

    if index % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {index}/{len(groups)} chains, {len(rows)} sites "
              f"({elapsed/index:.2f}s/chain, eta {(len(groups)-index)*elapsed/index/60:.0f} min)",
              flush=True)
        if rows:
            pd.DataFrame(rows).to_csv(out, mode="a", header=not out.exists(), index=False)
            rows = []

if rows:
    pd.DataFrame(rows).to_csv(out, mode="a", header=not out.exists(), index=False)
frame = pd.read_csv(out, low_memory=False).drop_duplicates(KEY) if out.exists() else pd.DataFrame()
if len(frame):
    # Biological panel first, QC and provenance after it, so a reader cannot
    # mistake a technical field for a predictor by column position.
    lead = KEY + ["subtype"]
    ordered = ([c for c in lead if c in frame.columns]
               + [c for c in FEATURE_COLUMNS if c in frame.columns]
               + [c for c in QC_COLUMNS if c in frame.columns])
    frame = frame[ordered + [c for c in frame.columns if c not in ordered]]
frame.to_csv(out, index=False)
if failures:
    pd.DataFrame(failures).to_csv(out.with_name(out.stem + "_failures.csv"), index=False)

print(f"\nextracted {len(frame)} sites; {len(failures)} failures; "
      f"elapsed {(time.time()-t0)/60:.1f} min")
if len(frame):
    # Chain-level dssp_ok says DSSP ran; it does not say all three sequon
    # positions carry a call. The second number is the one an analysis using
    # secondary structure can actually use.
    print("\nDSSP coverage by population (report this, do not filter on it):")
    print(f"  {'population':24}{'chain ran':>12}{'all 3 positions':>18}")
    complete = frame[[f"{p}_dssp_ok" for p in ("n", "plus1", "plus2")]].all(axis=1)
    for name, group in frame.groupby("population"):
        ok = int(group.dssp_ok.sum())
        full = int(complete[group.index].sum())
        print(f"  {name:24}{ok:6d}/{len(group):<5d}{full:11d}/{len(group):<6d}"
              f" ({100*full/len(group):5.1f}%)")
    broken = int((~frame.mapping_continuous.astype(bool)).sum())
    print(f"\nsequon mapping discontinuous (a +1 or +2 was never resolved): {broken}")
    bad = frame[~frame.triplet_matches.astype(bool)]
    print(f"triplet mismatches (structure mapping wrong): {len(bad)}")
    if len(bad):
        print("  ", bad.head(3)[["accession", "position", "triplet_expected",
                                 "triplet_observed"]].to_dict("records"))
