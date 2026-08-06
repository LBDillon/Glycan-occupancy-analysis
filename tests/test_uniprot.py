from __future__ import annotations

from experimental_glycosylation_sites.uniprot import parse_glycosylation_column


def parse_one(text: str):
    features = parse_glycosylation_column(text, "P00001")
    assert len(features) == 1
    return features[0]


def test_parses_exact_n_linked_asparagine_site():
    text = 'CARBOHYD 64; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000269|PubMed:123"'
    feature = parse_one(text)
    assert feature.position == 64
    assert feature.glyco_type == "N-linked"
    assert feature.evidence_codes == frozenset({"ECO:0000269"})
    assert feature.parse_status == "ok"


def test_glycation_is_not_n_linked():
    text = 'CARBOHYD 10; /note="N-linked (Glc) (glycation) lysine"; /evidence="ECO:0000269"'
    assert parse_one(text).glyco_type == "glycation"


def test_o_linked_is_separate_type():
    text = 'CARBOHYD 5; /note="O-linked (GalNAc...) threonine"'
    assert parse_one(text).glyco_type == "O-linked"


def test_range_feature_is_rejected_not_truncated():
    text = 'CARBOHYD 10..12; /note="N-linked (GlcNAc...) asparagine"'
    feature = parse_one(text)
    assert feature.position is None
    assert feature.parse_status == "uncertain_or_range_position"


def test_uncertain_position_is_rejected():
    text = 'CARBOHYD ?; /note="N-linked (GlcNAc...) asparagine"'
    feature = parse_one(text)
    assert feature.position is None
    assert feature.parse_status == "uncertain_or_range_position"


def test_multiple_features_parsed_independently():
    text = (
        'CARBOHYD 10; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000269"; '
        'CARBOHYD 20; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000255"'
    )
    features = parse_glycosylation_column(text, "P00001")
    assert [f.position for f in features] == [10, 20]
    assert features[0].evidence_codes == frozenset({"ECO:0000269"})
    assert features[1].evidence_codes == frozenset({"ECO:0000255"})


def test_multiple_evidence_codes_on_one_feature():
    text = (
        'CARBOHYD 30; /note="N-linked (GlcNAc...) asparagine"; '
        '/evidence="ECO:0000269|PubMed:1, ECO:0007744|PDB:1ABC"'
    )
    assert parse_one(text).evidence_codes == frozenset({"ECO:0000269", "ECO:0007744"})


def test_empty_and_null_input():
    assert parse_glycosylation_column("", "P00001") == []
    assert parse_glycosylation_column(None, "P00001") == []


def test_n_linked_without_asparagine_is_other():
    text = 'CARBOHYD 15; /note="N-linked (GlcNAc...) tryptophan"'
    assert parse_one(text).glyco_type == "N-linked-other"
