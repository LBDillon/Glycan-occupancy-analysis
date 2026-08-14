from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from experimental_glycosylation_sites.config import load_config
from experimental_glycosylation_sites.pipeline import run_uniprot_baseline

MODULE_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = MODULE_ROOT / "tests" / "snapshots" / "uniprot_baseline_2026-04-27.json"
ENRICHED = MODULE_ROOT / "tests" / "snapshots" / "enriched_2026-08-06.json"


def fingerprint(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()[:16]


# The UniProt baseline reads only these three. Gating on every configured path
# would let an unrelated missing file silently skip the primary regression guard.
BASELINE_INPUTS = ("pairs_master", "homology_qc", "uniprot_tsv")


@pytest.fixture(scope="module")
def config():
    cfg = load_config(MODULE_ROOT / "config" / "default.toml")
    missing = [
        f"{key}: {cfg.paths[key]}"
        for key in BASELINE_INPUTS
        if not cfg.paths[key].exists()
    ]
    if missing:
        pytest.skip("canonical inputs unavailable: " + "; ".join(missing))
    return cfg


def test_uniprot_baseline_matches_snapshot(config):
    expected = json.loads(SNAPSHOT.read_text())["counts"]
    actual = run_uniprot_baseline(config)["counts"]

    deltas = {
        key: {"expected": value, "actual": actual[key]}
        for key, value in expected.items()
        if actual[key] != value
    }
    if deltas:
        prints = {
            key: fingerprint(path)
            for key, path in sorted(config.paths.items())
            if path.is_file()
        }
        pytest.fail(
            "Baseline counts drifted from the 2026-04-27 snapshot.\n"
            f"Deltas: {json.dumps(deltas, indent=2)}\n"
            f"Current input fingerprints: {json.dumps(prints, indent=2)}\n"
            "Investigate whether UniProt or the canonical tables changed "
            "before editing the snapshot file."
        )


def test_enriched_counts_match_snapshot(config):
    if not ENRICHED.exists():
        pytest.skip("enriched fixture not yet frozen")
    cache = Path(config.paths["cache_dir"]) / "glygen_protein_detail.jsonl"
    if not cache.exists():
        pytest.skip("GlyGen cache absent; run fetch-glygen first")

    from experimental_glycosylation_sites.pipeline import run_full

    expected = json.loads(ENRICHED.read_text())["counts"]
    actual = run_full(config, fetch=False)
    deltas = {
        key: {"expected": value, "actual": actual[key]}
        for key, value in expected.items()
        if actual[key] != value
    }
    if deltas:
        pytest.fail(
            "Enriched counts drifted from the 2026-08-06 snapshot.\n"
            f"Deltas: {json.dumps(deltas, indent=2)}\n"
            "Check whether the GlyGen cache or cached structures changed."
        )
