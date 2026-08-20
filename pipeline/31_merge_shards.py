"""Combine SLURM array shards into the single table the analysis expects.

A sharded run is only as trustworthy as its merge. This does three things the
analysis cannot do for itself:

  * checks every shard of the array is present, because a missing one is a
    silently short table rather than an error -- exactly the failure mode that
    lost 21% of the first ESM-IF retention run;
  * checks no site appears twice, which would mean the shards overlapped;
  * reports how many sites were expected from the manifest and how many arrived.

Usage:
    31_merge_shards.py 'results/retention_scoring_manifest_esm_if.shard*.csv' \
                       results/designs/retention_scoring_manifest_esm_if.csv \
                       [--expect-manifest results/manifests/scoring_manifest.csv] \
                       [--shards 16]
"""
import argparse, glob, re, sys
from pathlib import Path

import pandas as pd

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("pattern", help="glob matching the shard files (quote it)")
parser.add_argument("out")
parser.add_argument("--expect-manifest", default=None,
                    help="manifest the run came from, to report coverage")
parser.add_argument("--shards", type=int, default=None,
                    help="array width; warns if a shard index is missing")
args = parser.parse_args()

# The failures side-files sit next to the shards and match the same glob, so
# exclude them explicitly rather than relying on the caller's pattern.
files = sorted(f for f in glob.glob(args.pattern) if "_failures" not in Path(f).name)
if not files:
    raise SystemExit(f"no files match {args.pattern!r}")

frames, empty = [], []
for path in files:
    if Path(path).stat().st_size < 50:
        empty.append(path)
        continue
    frames.append(pd.read_csv(path, low_memory=False))
if not frames:
    raise SystemExit(f"every file matching {args.pattern!r} is empty")

merged = pd.concat(frames, ignore_index=True)
before = len(merged)
merged = merged.drop_duplicates(KEY)

print(f"shard files found : {len(files)}" + (f" ({len(empty)} empty)" if empty else ""))
print(f"rows              : {before} -> {len(merged)} after de-duplication")

if before != len(merged):
    print(f"  WARNING: {before - len(merged)} duplicate sites across shards. Shards "
          "should partition the chains, so this means two tasks ran the same one.")

if args.shards:
    seen = {int(m.group(1)) for p in files
            for m in [re.search(r"shard(\d+)", Path(p).name)] if m}
    missing = sorted(set(range(args.shards)) - seen)
    if missing:
        print(f"  WARNING: shards {missing} are absent. The merged table is "
              "INCOMPLETE -- resubmit those array indices before analysing it.")
    else:
        print(f"all {args.shards} shards present")

if args.expect_manifest and Path(args.expect_manifest).exists():
    man = pd.read_csv(args.expect_manifest, low_memory=False)
    if "scoreable" in man.columns:
        man = man[man.scoreable.astype(bool)]
    expected = len(man.drop_duplicates(KEY))
    pct = 100 * len(merged) / expected if expected else 0
    print(f"coverage          : {len(merged)} of {expected} manifest sites ({pct:.1f}%)")
    if pct < 95:
        print("  WARNING: under 95% coverage. Check the *_failures.csv files "
              "before treating this as a complete run.")

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(args.out, index=False)
print(f"wrote {args.out}")

failures = sorted(glob.glob(args.pattern.replace(".csv", "_failures.csv")))
rows = [pd.read_csv(f, low_memory=False) for f in failures
        if Path(f).stat().st_size > 50]
if rows:
    fail = pd.concat(rows, ignore_index=True)
    out = Path(args.out).with_name(Path(args.out).stem + "_failures.csv")
    fail.to_csv(out, index=False)
    print(f"\n{len(fail)} failures -> {out}")
    print(fail.reason.str.replace(r"\[.*", "", regex=True).str.slice(0, 70)
              .value_counts().head(5).to_string())
