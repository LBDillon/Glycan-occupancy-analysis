from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from experimental_glycosylation_sites.config import Config
from experimental_glycosylation_sites.glygen import (
    build_site_evidence,
    classify_glygen_entry,
    extract_sites,
    fetch_details,
)

FIXTURE = Path(__file__).parent / "fixtures" / "glygen_O43570.json"


def _config(tmp_path: Path, retries: int = 3) -> Config:
    return Config(
        paths={"cache_dir": tmp_path / "cache"},
        layers={},
        policy={},
        api={
            "glygen_url": "http://example.invalid/{accession}/",
            "delay_seconds": 0.0,
            "timeout_seconds": 1,
            "max_retries": retries,
        },
        source_path=tmp_path / "config.toml",
        config_hash="test",
    )


def entry(category, databases, n_type="N-linked", start=10, end=10):
    return {
        "type": n_type,
        "site_category": category,
        "site_category_dict": {category: True},
        "start_pos": start,
        "end_pos": end,
        "evidence": [{"database": db, "id": "x"} for db in databases],
    }


def test_predicted_is_the_uniprot_rule_and_never_qualifies():
    assert classify_glygen_entry(entry("predicted", ["UniProtKB"])) == "glygen_predicted"


def test_reported_with_glycan_is_independent():
    result = classify_glygen_entry(entry("reported_with_glycan", ["PubMed", "GlyConnect"]))
    assert result == "glygen_reported_with_glycan"


def test_reported_citing_only_uniprot_is_marked_circular():
    result = classify_glygen_entry(entry("reported", ["UniProtKB"]))
    assert result == "glygen_reported_uniprot_derived"


def test_reported_citing_pubmed_is_independent():
    result = classify_glygen_entry(entry("reported", ["UniProtKB", "PubMed"]))
    assert result == "glygen_reported_independent"


def test_o_linked_entries_are_ignored():
    sites = extract_sites({"glycosylation": [entry("reported_with_glycan", ["PubMed"], n_type="O-linked")]})
    assert sites == {}


def test_range_entries_are_ignored():
    sites = extract_sites({"glycosylation": [entry("reported_with_glycan", ["PubMed"], start=10, end=12)]})
    assert sites == {}


def test_entry_without_position_is_ignored():
    broken = entry("reported_with_glycan", ["PubMed"])
    del broken["start_pos"]
    assert extract_sites({"glycosylation": [broken]}) == {}


def test_strongest_tier_wins_when_a_position_has_several_entries():
    detail = {"glycosylation": [
        entry("predicted", ["UniProtKB"], start=10, end=10),
        entry("reported_with_glycan", ["PubMed"], start=10, end=10),
    ]}
    assert extract_sites(detail)[10]["tier"] == "glygen_reported_with_glycan"


def test_real_fixture_yields_expected_tiers():
    detail = json.loads(FIXTURE.read_text())
    sites = extract_sites(detail)
    assert sites[162]["tier"] == "glygen_predicted"
    assert sites[80]["tier"] == "glygen_reported_with_glycan"
    assert "29741879" in sites[80]["pubmed_ids"]


def test_build_site_evidence_covers_only_candidate_positions():
    detail = json.loads(FIXTURE.read_text())
    candidates = pd.DataFrame([
        {"accession": "O43570", "position": 80},
        {"accession": "O43570", "position": 999},
    ])
    frame = build_site_evidence(candidates, {"O43570": detail})
    assert len(frame) == 2
    by_position = frame.set_index("position")
    assert by_position.loc[80, "glygen_tier"] == "glygen_reported_with_glycan"
    assert by_position.loc[999, "glygen_tier"] == ""


def test_fetch_details_does_not_retry_a_permanent_http_error(tmp_path, monkeypatch):
    # GlyGen answers HTTP 500 (not 404) for an accession it has no entry for.
    # That is permanent, so a single request should be made, not max_retries.
    calls = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        raise urllib.error.HTTPError(url, 500, "Internal Server Error", None, None)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    config = _config(tmp_path, retries=3)
    cache_path = fetch_details(["P00001"], config)

    assert len(calls) == 1
    record = json.loads(cache_path.read_text().strip())
    assert record["detail"] is None
    assert "HTTPError" in record["error"]


def test_fetch_details_retries_a_transient_url_error(tmp_path, monkeypatch):
    calls = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        raise urllib.error.URLError("connection refused")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    config = _config(tmp_path, retries=3)
    cache_path = fetch_details(["P00001"], config)

    assert len(calls) == 3
    record = json.loads(cache_path.read_text().strip())
    assert record["detail"] is None
    assert "URLError" in record["error"]
