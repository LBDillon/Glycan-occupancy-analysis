"""Merge the context extraction shards, refusing anything short of complete.

Separate from 31_merge_shards because the context table's key includes
`population`: the same (accession, position) is legitimately both an occupied
site and, in another arm, a control, so the scoring merge's key would delete
real rows and report a clean deduplication.

Every condition that could reduce coverage is fatal here -- missing shards,
duplicate keys, recorded extraction failures, manifest rows that never arrived.

Usage:
    43b_merge_context_shards.py [--pattern 'glyco_context/results/datasets/context_features.shard*.csv']
                                [--manifest results/datasets/context_manifest.csv]
                                [--populations occupied,internal_control,secretory_unannotated]
                                [--shards 4] [--out results/datasets/context_features.csv]
"""
import argparse, glob, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.context_merge import (KEY, ShardError,
                                                            merge_context_shards)

DEFAULT_POPULATIONS = "occupied,internal_control,secretory_unannotated"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pattern", default="glyco_context/results/datasets/context_features.shard*.csv")
parser.add_argument("--manifest", default="glyco_context/results/datasets/context_manifest.csv")
parser.add_argument("--populations", default=DEFAULT_POPULATIONS)
parser.add_argument("--shards", type=int, default=4)
parser.add_argument("--out", default="glyco_context/results/datasets/context_features.csv")
args = parser.parse_args()

paths = sorted(p for p in glob.glob(args.pattern) if "_failures" not in Path(p).name)
print(f"shard files: {len(paths)}")
for path in paths:
    print(f"  {path}")

manifest = pd.read_csv(args.manifest, low_memory=False)
manifest = manifest[manifest.features_available.astype(bool)]
manifest = manifest[manifest.population.isin(
    [p.strip() for p in args.populations.split(",")])]
manifest = manifest.dropna(subset=["structure_pdb_id", "structure_chain_id",
                                   "structure_resseq"])
print(f"eligible manifest rows: {len(manifest)}")

try:
    merged = merge_context_shards(paths, manifest, expected_shards=args.shards)
except ShardError as error:
    raise SystemExit(f"MERGE FAILED: {error}")

merged.to_csv(args.out, index=False)
print(f"\nmerged {len(merged)} rows -> {args.out}")
print(merged.population.value_counts().to_string())
