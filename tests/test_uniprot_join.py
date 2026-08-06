from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.evidence import join_uniprot_evidence
from experimental_glycosylation_sites.models import UniProtFeature

QUALIFYING = ["manual_experimental", "manual_combinatorial"]


def candidates() -> pd.DataFrame:
    return pd.DataFrame([
        {"accession": "P1", "position": 10},   # experimental
        {"accession": "P1", "position": 20},   # sequence model only
        {"accession": "P1", "position": 30},   # no feature at this position
        {"accession": "P2", "position": 40},   # curator inference
        {"accession": "P3", "position": 50},   # accession absent from snapshot
        {"accession": "P4", "position": 60},   # O-linked feature at the position
    ])


def features() -> list[UniProtFeature]:
    def feature(acc, pos, kind, codes):
        return UniProtFeature(acc, pos, kind, frozenset(codes), "", "ok")

    return [
        feature("P1", 10, "N-linked", {"ECO:0000269"}),
        feature("P1", 20, "N-linked", {"ECO:0000255"}),
        feature("P1", 99, "N-linked", {"ECO:0000269"}),   # different position
        feature("P2", 40, "N-linked", {"ECO:0000305"}),
        feature("P4", 60, "O-linked", {"ECO:0000269"}),
    ]


def joined() -> pd.DataFrame:
    return join_uniprot_evidence(candidates(), features(), {"P3"}, QUALIFYING)


def reason_for(accession: str, position: int) -> str:
    row = joined()
    row = row[(row.accession == accession) & (row.position == position)].iloc[0]
    return row.exclusion_reason


def test_experimental_site_qualifies():
    row = joined()
    row = row[(row.accession == "P1") & (row.position == 10)].iloc[0]
    assert bool(row.uniprot_qualifies) is True
    assert row.uniprot_tier == "manual_experimental"
    assert row.exclusion_reason == ""


def test_nearby_feature_does_not_satisfy_a_different_position():
    assert reason_for("P1", 30) == "exact_feature_absent"


def test_sequence_model_only_is_excluded_with_reason():
    assert reason_for("P1", 20) == "sequence_model_only"


def test_curator_inference_is_excluded_with_reason():
    assert reason_for("P2", 40) == "curator_inference_only"


def test_missing_accession_gets_its_own_reason():
    assert reason_for("P3", 50) == "accession_absent_from_snapshot"


def test_non_n_linked_feature_does_not_count():
    assert reason_for("P4", 60) == "exact_feature_absent"


def test_every_candidate_appears_exactly_once():
    result = joined()
    assert len(result) == len(candidates())
    assert not result.duplicated(["accession", "position"]).any()


def test_partition_is_complete():
    result = joined()
    qualifying = result[result.uniprot_qualifies]
    excluded = result[~result.uniprot_qualifies]
    assert len(qualifying) + len(excluded) == len(result)
    assert (qualifying.exclusion_reason == "").all()
    assert (excluded.exclusion_reason != "").all()


def test_evidence_codes_are_recorded_sorted_and_joined():
    result = join_uniprot_evidence(
        pd.DataFrame([{"accession": "P1", "position": 10}]),
        [UniProtFeature("P1", 10, "N-linked",
                        frozenset({"ECO:0007744", "ECO:0000269"}), "", "ok")],
        set(), QUALIFYING,
    )
    assert result.iloc[0].uniprot_evidence_codes == "ECO:0000269|ECO:0007744"
