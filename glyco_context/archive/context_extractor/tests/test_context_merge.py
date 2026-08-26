"""Merging the extraction shards, where silence is the dangerous outcome.

A missing shard produces a short table, not an error; the first ESM-IF
retention run lost 21% that way. These checks turn every such case into a
failure, because a context atlas built on 90% of its sites looks exactly like
one built on all of them.
"""
from __future__ import annotations

import pandas as pd
import pytest

from glyco_context.context_merge import (ShardError,
                                                            merge_context_shards)

MANIFEST = pd.DataFrame([
    {"accession": "P1", "position": 1, "population": "occupied"},
    {"accession": "P1", "position": 1, "population": "internal_control"},
    {"accession": "P2", "position": 5, "population": "occupied"},
])


def _shard(tmp_path, name, rows):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_merges_every_shard(tmp_path):
    a = _shard(tmp_path, "f.shard0.csv", [MANIFEST.iloc[0].to_dict(),
                                          MANIFEST.iloc[1].to_dict()])
    b = _shard(tmp_path, "f.shard1.csv", [MANIFEST.iloc[2].to_dict()])
    merged = merge_context_shards([a, b], MANIFEST, expected_shards=2)
    assert len(merged) == 3


def test_population_is_part_of_the_key(tmp_path):
    """P1/1 appears twice under different populations and both must survive."""
    a = _shard(tmp_path, "f.shard0.csv", [MANIFEST.iloc[0].to_dict(),
                                          MANIFEST.iloc[1].to_dict()])
    b = _shard(tmp_path, "f.shard1.csv", [MANIFEST.iloc[2].to_dict()])
    merged = merge_context_shards([a, b], MANIFEST, expected_shards=2)
    assert set(merged.population) == {"occupied", "internal_control"}


def test_missing_shard_is_fatal(tmp_path):
    a = _shard(tmp_path, "f.shard0.csv", [MANIFEST.iloc[0].to_dict()])
    with pytest.raises(ShardError, match="shard"):
        merge_context_shards([a], MANIFEST, expected_shards=2)


def test_duplicate_key_is_fatal_not_deduplicated(tmp_path):
    row = MANIFEST.iloc[0].to_dict()
    a = _shard(tmp_path, "f.shard0.csv", [row])
    b = _shard(tmp_path, "f.shard1.csv", [row])
    with pytest.raises(ShardError, match="duplicate"):
        merge_context_shards([a, b], MANIFEST, expected_shards=2)


def test_incomplete_manifest_coverage_is_fatal(tmp_path):
    a = _shard(tmp_path, "f.shard0.csv", [MANIFEST.iloc[0].to_dict()])
    b = _shard(tmp_path, "f.shard1.csv", [MANIFEST.iloc[1].to_dict()])
    with pytest.raises(ShardError, match="missing from the merged table"):
        merge_context_shards([a, b], MANIFEST, expected_shards=2)


def test_extraction_failures_are_fatal(tmp_path):
    a = _shard(tmp_path, "f.shard0.csv", [MANIFEST.iloc[0].to_dict(),
                                          MANIFEST.iloc[1].to_dict()])
    b = _shard(tmp_path, "f.shard1.csv", [MANIFEST.iloc[2].to_dict()])
    pd.DataFrame([{"accession": "P9", "position": 1, "population": "occupied",
                   "reason": "structure_not_cached"}]).to_csv(
        tmp_path / "f.shard1_failures.csv", index=False)
    with pytest.raises(ShardError, match="failure"):
        merge_context_shards([a, b], MANIFEST, expected_shards=2)
