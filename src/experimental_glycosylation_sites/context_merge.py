"""Merging extraction shards, where the dangerous outcome is silence.

A missing shard yields a short table rather than an error, and a short context
atlas looks exactly like a complete one. Every condition here that could reduce
coverage raises instead: missing shards, duplicate keys, recorded extraction
failures, and manifest rows that never arrived.

The key includes `population`. A site can legitimately appear under more than
one population -- the same (accession, position) is both an occupied site and,
in another arm, a control -- so keying on (accession, position) alone would
delete real rows while reporting a clean deduplication.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

KEY = ["accession", "position", "population"]


class ShardError(RuntimeError):
    """A shard set that cannot be merged without losing or inventing rows."""


def _shard_index(path: Path) -> "int | None":
    match = re.search(r"shard(\d+)", Path(path).name)
    return int(match.group(1)) if match else None


def merge_context_shards(paths, manifest: pd.DataFrame,
                         expected_shards: "int | None" = None) -> pd.DataFrame:
    """The merged table, or ShardError explaining why it would be incomplete."""
    paths = [Path(p) for p in paths]
    if not paths:
        raise ShardError("no shard files given")

    if expected_shards is not None:
        found = {i for i in (_shard_index(p) for p in paths) if i is not None}
        missing = sorted(set(range(expected_shards)) - found)
        if missing:
            raise ShardError(
                f"shard index {missing} missing; expected {expected_shards} shards, "
                f"found {sorted(found)}. A missing shard is a short table, not an error.")

    # Recorded failures mean sites the extractor could not process at all.
    failures = []
    for path in paths:
        side = path.with_name(path.stem + "_failures.csv")
        if side.exists() and side.stat().st_size > 5:
            failures.append(pd.read_csv(side, low_memory=False))
    if failures:
        combined = pd.concat(failures, ignore_index=True)
        if len(combined):
            reasons = combined.get("reason", pd.Series(dtype=str)).value_counts().to_dict()
            raise ShardError(
                f"{len(combined)} extraction failure(s) recorded: {reasons}. "
                "Resolve them rather than merging an incomplete table.")

    merged = pd.concat([pd.read_csv(p, low_memory=False) for p in paths],
                       ignore_index=True)

    duplicated = merged.duplicated(KEY, keep=False)
    if duplicated.any():
        example = merged[duplicated].iloc[0][KEY].to_dict()
        raise ShardError(
            f"{int(duplicated.sum())} rows share a context key, e.g. {example}. "
            "Shards must partition the sites; overlapping shards are a sharding "
            "defect and are not deduplicated here.")

    expected_rows = manifest[KEY].drop_duplicates()
    arrived = expected_rows.merge(merged[KEY].assign(_seen=1), on=KEY, how="left")
    absent = arrived[arrived._seen.isna()]
    if len(absent):
        example = absent.iloc[0][KEY].to_dict()
        raise ShardError(
            f"{len(absent)} manifest row(s) missing from the merged table, "
            f"e.g. {example}.")
    return merged
