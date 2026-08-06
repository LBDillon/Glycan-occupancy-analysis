from __future__ import annotations

import gzip

from experimental_glycosylation_sites.uniprot import (
    load_uniprot_features,
    parse_glycosylation_column,
)


def parse_one(text: str):
    features = parse_glycosylation_column(text, "P00001")
    assert len(features) == 1
    return features[0]


def _write_tsv(path, rows, opener=open, **open_kwargs):
    """Write a minimal UniProt-style TSV: an Entry/Glycosylation header plus rows."""
    with opener(path, "wt", encoding="utf-8", newline="", **open_kwargs) as handle:
        handle.write("Entry\tGlycosylation\n")
        for entry, glycosylation in rows:
            handle.write(f"{entry}\t{glycosylation}\n")


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


def test_load_uniprot_features_parses_plain_tsv(tmp_path):
    tsv_path = tmp_path / "snapshot.tsv"
    _write_tsv(
        tsv_path,
        [("P00001", 'CARBOHYD 64; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000269"')],
    )

    features, missing = load_uniprot_features(tsv_path, {"P00001"})

    assert len(features) == 1
    assert features[0].accession == "P00001"
    assert features[0].position == 64
    assert missing == set()


def test_load_uniprot_features_parses_gzipped_tsv(tmp_path):
    tsv_path = tmp_path / "snapshot.tsv.gz"
    _write_tsv(
        tsv_path,
        [("P00002", 'CARBOHYD 10; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000269"')],
        opener=gzip.open,
    )

    features, missing = load_uniprot_features(tsv_path, {"P00002"})

    assert len(features) == 1
    assert features[0].accession == "P00002"
    assert features[0].position == 10
    assert missing == set()


def test_load_uniprot_features_skips_unrequested_accessions(tmp_path):
    tsv_path = tmp_path / "snapshot.tsv"
    _write_tsv(
        tsv_path,
        [
            ("P00001", 'CARBOHYD 64; /note="N-linked (GlcNAc...) asparagine"'),
            ("P99999", 'CARBOHYD 5; /note="N-linked (GlcNAc...) asparagine"'),
        ],
    )

    features, missing = load_uniprot_features(tsv_path, {"P00001"})

    assert [f.accession for f in features] == ["P00001"]
    assert missing == set()


def test_load_uniprot_features_reports_missing_accessions(tmp_path):
    tsv_path = tmp_path / "snapshot.tsv"
    _write_tsv(
        tsv_path,
        [("P00001", 'CARBOHYD 64; /note="N-linked (GlcNAc...) asparagine"')],
    )

    features, missing = load_uniprot_features(tsv_path, {"P00001", "P00002"})

    assert [f.accession for f in features] == ["P00001"]
    assert missing == {"P00002"}


def test_load_uniprot_features_empty_glycosylation_marks_seen_not_missing(tmp_path):
    tsv_path = tmp_path / "snapshot.tsv"
    _write_tsv(tsv_path, [("P00001", "")])

    features, missing = load_uniprot_features(tsv_path, {"P00001"})

    assert features == []
    assert missing == set()
