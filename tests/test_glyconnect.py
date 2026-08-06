from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.glyconnect import (
    build_site_evidence,
    extract_positions,
    parse_site_label,
)


def detail() -> dict:
    return {
        "id": 1064,
        "structures": [
            {"type": "N-Linked", "composition_string": "Hex:6 HexNAc:5 NeuAc:3",
             "sites": [{"site": "Asn-80", "reviewed": False}]},
            {"type": "N-Linked", "composition_string": "Hex:5 HexNAc:2",
             "sites": [{"site": "Asn-80"}, {"site": "Asn-120"}]},
            {"type": "O-Linked", "composition_string": "HexNAc:1",
             "sites": [{"site": "Thr-40"}]},
        ],
    }


def test_parses_standard_site_label():
    assert parse_site_label("Asn-80") == 80


def test_rejects_non_asparagine_label():
    assert parse_site_label("Thr-40") is None


def test_rejects_unparseable_label():
    assert parse_site_label("") is None
    assert parse_site_label("Asn-?") is None
    assert parse_site_label(None) is None


def test_extracts_only_n_linked_asparagine_positions():
    positions = extract_positions(detail())
    assert set(positions) == {80, 120}


def test_counts_distinct_glycan_structures_per_site():
    positions = extract_positions(detail())
    assert positions[80]["glycan_count"] == 2
    assert positions[120]["glycan_count"] == 1


def test_compositions_are_recorded():
    positions = extract_positions(detail())
    assert "Hex:6 HexNAc:5 NeuAc:3" in positions[80]["compositions"]


def test_build_site_evidence_marks_only_matching_candidates():
    candidates = pd.DataFrame([
        {"accession": "O43570", "position": 80},
        {"accession": "O43570", "position": 999},
    ])
    frame = build_site_evidence(candidates, {"O43570": detail()})
    by_position = frame.set_index("position")
    assert bool(by_position.loc[80, "glyconnect_supported"]) is True
    assert bool(by_position.loc[999, "glyconnect_supported"]) is False
