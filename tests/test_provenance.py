from __future__ import annotations

import hashlib
from pathlib import Path

from experimental_glycosylation_sites.config import load_config
from experimental_glycosylation_sites.provenance import build_manifest, hash_file


def make_config(tmp_path: Path):
    (tmp_path / "config").mkdir()
    (tmp_path / "input.csv").write_text("a,b\n1,2\n")
    cfg = tmp_path / "config" / "test.toml"
    cfg.write_text(
        "[paths]\n"
        'pairs_master = "../input.csv"\n'
        'results_dir = "../results"\n'
        "[layers]\n"
        "uniprot = true\n"
        "glygen = false\n"
        "[policy]\n"
        'qualifying_uniprot_tiers = ["manual_experimental"]\n'
        "[api]\n"
        "delay_seconds = 0.4\n"
    )
    return load_config(cfg)


def test_hash_file_matches_hashlib(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(b"hello")
    assert hash_file(path) == hashlib.sha256(b"hello").hexdigest()


def test_manifest_records_input_hashes(tmp_path):
    config = make_config(tmp_path)
    manifest = build_manifest(config, {"all_sites": 5}, {})
    entry = manifest["inputs"]["pairs_master"]
    assert entry["sha256"] == hash_file(tmp_path / "input.csv")
    assert entry["size_bytes"] == 8
    assert "modified" in entry


def test_manifest_records_config_hash_and_layers(tmp_path):
    config = make_config(tmp_path)
    manifest = build_manifest(config, {}, {})
    assert manifest["config"]["hash"] == config.config_hash
    assert manifest["config"]["layers"] == {"uniprot": True, "glygen": False}


def test_manifest_records_counts_and_extra(tmp_path):
    manifest = build_manifest(make_config(tmp_path), {"all_sites": 5}, {"glygen_fetched": 3})
    assert manifest["counts"]["all_sites"] == 5
    assert manifest["extra"]["glygen_fetched"] == 3


def test_manifest_skips_output_directories(tmp_path):
    manifest = build_manifest(make_config(tmp_path), {}, {})
    assert "results_dir" not in manifest["inputs"]
