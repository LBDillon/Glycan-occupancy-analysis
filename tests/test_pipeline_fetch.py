from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from experimental_glycosylation_sites import glyconnect as glyconnect_layer
from experimental_glycosylation_sites import glygen as glygen_layer
from experimental_glycosylation_sites.config import load_config
from experimental_glycosylation_sites.pipeline import run_full


def _write_config(tmp_path: Path) -> Path:
    """A hermetic mini-database exercising the whole run_full wiring.

    Four analysis-ready candidates with deliberately different cross-reference
    coverage, so a fetch list that ignores the cross-references is visibly wrong.
    """
    (tmp_path / "config").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sibling").mkdir(parents=True, exist_ok=True)

    pairs = pd.DataFrame([
        {"pair_id": "pair1", "pos_accession": "P00001",
         "best_loss_positive_position": 10, "analysis_ready": True},
        {"pair_id": "pair2", "pos_accession": "P00002",
         "best_loss_positive_position": 20, "analysis_ready": True},
        {"pair_id": "pair3", "pos_accession": "P00003",
         "best_loss_positive_position": 30, "analysis_ready": True},
        # No cross-reference at all: still a candidate site, never fetched.
        {"pair_id": "pair4", "pos_accession": "P00004",
         "best_loss_positive_position": 40, "analysis_ready": True},
    ])
    pairs.to_csv(tmp_path / "sibling" / "pairs.csv", index=False)

    pd.DataFrame([
        {"pair_id": f"pair{n}", "homology_qc_bucket": "strict_ortholog_like"}
        for n in range(1, 5)
    ]).to_csv(tmp_path / "sibling" / "homology.csv", index=False)

    tsv_path = tmp_path / "sibling" / "uniprot.tsv.gz"
    with gzip.open(tsv_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("Entry\tGlycosylation\tSequence\tGlyGen\tGlyConnect\n")
        handle.write("P00001\t\tMNTS\tG12345\tC12345\n")
        handle.write("P00002\t\tMNTS\t\tC99999\n")
        handle.write("P00003\t\tMNTS\tG99999\t\n")
        handle.write("P00004\t\tMNTS\t\t\n")

    cfg = tmp_path / "config" / "test.toml"
    cfg.write_text(
        "[paths]\n"
        'pairs_master = "../sibling/pairs.csv"\n'
        'homology_qc = "../sibling/homology.csv"\n'
        'uniprot_tsv = "../sibling/uniprot.tsv.gz"\n'
        'cache_dir = "../cache"\n'
        'results_dir = "../results"\n'
        "[layers]\n"
        "uniprot = true\n"
        "glygen = true\n"
        "glyconnect = true\n"
        "structure = false\n"
        "[policy]\n"
        "require_analysis_ready = true\n"
        'qualifying_uniprot_tiers = ["manual_experimental"]\n'
        'strict_buckets = ["strict_ortholog_like"]\n'
        "[api]\n"
        "delay_seconds = 0.0\n"
    )
    return cfg


def _record_fetches(monkeypatch) -> dict[str, list[str]]:
    """Capture what each layer would request, without touching the network."""
    calls: dict[str, list[str]] = {}

    def fake_glygen(accessions, config):
        calls["glygen"] = list(accessions)
        return Path(config.paths["cache_dir"]) / "glygen_protein_detail.jsonl"

    def fake_glyconnect(accessions, config):
        calls["glyconnect"] = list(accessions)
        return Path(config.paths["cache_dir"]) / "glyconnect_protein_detail.jsonl"

    monkeypatch.setattr(glygen_layer, "fetch_details", fake_glygen)
    monkeypatch.setattr(glyconnect_layer, "fetch_details", fake_glyconnect)
    return calls


def test_run_fetch_requests_only_cross_referenced_accessions(tmp_path, monkeypatch):
    # Regression: run --fetch previously built its own unfiltered accession list,
    # so it requested every candidate from GlyGen. Accessions with no GlyGen
    # cross-reference answer HTTP 500, never 404, so those requests can only ever
    # burn the run's time budget.
    config = load_config(_write_config(tmp_path))
    calls = _record_fetches(monkeypatch)

    run_full(config, fetch=True)

    assert calls["glygen"] == ["P00001", "P00003"]
    assert calls["glyconnect"] == ["P00001", "P00002"]


def test_run_without_fetch_requests_nothing(tmp_path, monkeypatch):
    config = load_config(_write_config(tmp_path))
    calls = _record_fetches(monkeypatch)

    run_full(config, fetch=False)

    assert calls == {}


def test_fetch_filtering_does_not_narrow_the_candidate_universe(tmp_path, monkeypatch):
    """Only fetching is filtered; evidence still covers every candidate site.

    P00004 has no cross-reference in either source and so is never requested,
    but it must still appear as a candidate and be partitioned into a result
    table - otherwise the filter would be silently dropping sites.
    """
    config = load_config(_write_config(tmp_path))
    _record_fetches(monkeypatch)

    counts = run_full(config, fetch=True)

    assert counts["candidate_sites"] == 4
    assert counts["enriched_all_sites"] + counts["excluded_sites"] == 4

    candidates = pd.read_csv(Path(config.paths["results_dir"]) / "candidate_sites.csv")
    assert set(candidates["accession"]) == {"P00001", "P00002", "P00003", "P00004"}
