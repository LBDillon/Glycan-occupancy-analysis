from __future__ import annotations

from experimental_glycosylation_sites.evidence import (
    UNIPROT_TIER_ORDER,
    classify_uniprot_tier,
)


def test_literature_code_is_manual_experimental():
    assert classify_uniprot_tier(frozenset({"ECO:0000269"})) == "manual_experimental"


def test_combinatorial_code_is_not_called_pdb_backed():
    assert classify_uniprot_tier(frozenset({"ECO:0007744"})) == "manual_combinatorial"


def test_strongest_tier_wins_when_multiple_codes_present():
    codes = frozenset({"ECO:0000255", "ECO:0000305", "ECO:0000269"})
    assert classify_uniprot_tier(codes) == "manual_experimental"


def test_curator_inference_outranks_sequence_model():
    codes = frozenset({"ECO:0000255", "ECO:0000305"})
    assert classify_uniprot_tier(codes) == "manual_curator_inference"


def test_sequence_similarity_code_classified():
    assert classify_uniprot_tier(frozenset({"ECO:0000250"})) == "sequence_similarity"


def test_automatic_codes_classified():
    assert classify_uniprot_tier(frozenset({"ECO:0000256"})) == "automatic_sequence_model"
    assert classify_uniprot_tier(frozenset({"ECO:0000259"})) == "automatic_sequence_model"


def test_no_codes_is_annotation_without_qualifying_evidence():
    assert classify_uniprot_tier(frozenset()) == "annotation_without_qualifying_evidence"


def test_unknown_code_is_annotation_without_qualifying_evidence():
    assert classify_uniprot_tier(frozenset({"ECO:9999999"})) == "annotation_without_qualifying_evidence"


def test_tier_order_is_strongest_first():
    assert UNIPROT_TIER_ORDER[0] == "manual_experimental"
    assert UNIPROT_TIER_ORDER[1] == "manual_combinatorial"
