from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from experimental_glycosylation_sites.glygen import (
    build_site_evidence,
    classify_glygen_entry,
    extract_sites,
)

FIXTURE = Path(__file__).parent / "fixtures" / "glygen_O43570.json"


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
