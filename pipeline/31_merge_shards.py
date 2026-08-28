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
parser.add_argument("--allow-incomplete", action="store_true",
                    help="write the merged table even when shards are missing. For inspecting a partial run, never for analysis.")
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

# A shard writes its failures side-file only after its final flush, so a shard
# whose data file exists without that sibling is one that died partway through.
# Its output looks like any other shard's: shorter, and in no way marked.
truncated = [Path(f).name for f in files
             if not Path(f.replace(".csv", "_failures.csv")).exists()]
if truncated:
    message = (f"{truncated} have no _failures.csv beside them, which means "
               "those shards were killed before finishing. Their rows are a "
               "partial run, not a short one.")
    if not args.allow_incomplete:
        raise SystemExit(f"MERGE REFUSED: {message}\n"
                         "Pass --allow-incomplete to write it anyway.")
    print(f"  WARNING: {message}")

failure_files = sorted(glob.glob(args.pattern.replace(".csv", "_failures.csv")))
failure_rows = [pd.read_csv(f, low_memory=False) for f in failure_files
                if Path(f).stat().st_size > 50]
fail = pd.concat(failure_rows, ignore_index=True) if failure_rows else None

if args.shards:
    seen = {int(m.group(1)) for p in files
            for m in [re.search(r"shard(\d+)", Path(p).name)] if m}
    missing = sorted(set(range(args.shards)) - seen)
    if missing:
        # Exit non-zero rather than warn. A warning scrolls past, and what is
        # left behind is a file that looks exactly like a complete one -- which
        # is the failure this stage exists to prevent, not to narrate.
        # --allow-incomplete is for deliberately inspecting a partial run.
        message = (f"shards {missing} are absent of {args.shards}. The merged "
                   "table is INCOMPLETE. Resubmit those array indices, and "
                   "delete the partial shards first: a resubmitted array "
                   "resumes from existing output, so stale shards are mistaken "
                   "for finished work.")
        if not args.allow_incomplete:
            raise SystemExit(f"MERGE REFUSED: {message}\n"
                             "Pass --allow-incomplete to write it anyway.")
        print(f"  WARNING: {message}")
    else:
        print(f"all {args.shards} shards present")

if args.expect_manifest and Path(args.expect_manifest).exists():
    man = pd.read_csv(args.expect_manifest, low_memory=False)
    if "scoreable" in man.columns:
        man = man[man.scoreable.astype(bool)]
    expected = len(man.drop_duplicates(KEY))
    pct = 100 * len(merged) / expected if expected else 0
    n_failed = 0 if fail is None else len(fail.drop_duplicates(
        [c for c in KEY if c in fail.columns]))
    print(f"coverage          : {len(merged)} of {expected} manifest sites "
          f"({pct:.1f}%), {n_failed} recorded as failures")

    # Every manifest site must be either scored or recorded as having failed.
    # A site that is neither was never attempted, which means a shard stopped
    # early -- and a percentage threshold cannot see that: the run this check
    # was written for sat at 97.1%, comfortably above any threshold worth
    # setting, while sixteen sites had simply never been tried.
    unaccounted = expected - len(merged) - n_failed
    if unaccounted > 0:
        message = (f"{unaccounted} manifest sites are neither scored nor "
                   "recorded as failures, so they were never attempted. The "
                   "merged table is INCOMPLETE.")
        if not args.allow_incomplete:
            raise SystemExit(f"MERGE REFUSED: {message}\n"
                             "Pass --allow-incomplete to write it anyway.")
        print(f"  WARNING: {message}")
    elif unaccounted < 0:
        print(f"  WARNING: {-unaccounted} more sites than the manifest lists; "
              "the manifest and the run may not correspond.")

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(args.out, index=False)
print(f"wrote {args.out}")

if fail is not None:
    out = Path(args.out).with_name(Path(args.out).stem + "_failures.csv")
    fail.to_csv(out, index=False)
    print(f"\n{len(fail)} failures -> {out}")
    print(fail.reason.str.replace(r"\[.*", "", regex=True).str.slice(0, 70)
              .value_counts().head(5).to_string())
