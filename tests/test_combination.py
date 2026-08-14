from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.evidence import combine_layers

POLICY = {
    "qualifying_uniprot_tiers": ["manual_experimental", "manual_combinatorial"],
    "qualifying_glygen_tiers": ["glygen_reported_with_glycan", "glygen_reported_independent"],
    "qualifying_structure_tiers": ["structure_linked_glycan"],
}

SITES = [("P1", 10), ("P2", 20), ("P3", 30), ("P4", 40), ("P5", 50)]


def uniprot_frame() -> pd.DataFrame:
    tiers = {
        ("P1", 10): ("manual_experimental", True, ""),
        ("P2", 20): ("manual_sequence_model", False, "sequence_model_only"),
        ("P3", 30): ("manual_sequence_model", False, "sequence_model_only"),
        ("P4", 40): ("exact_feature_absent", False, "exact_feature_absent"),
        ("P5", 50): ("manual_curator_inference", False, "curator_inference_only"),
    }
    return pd.DataFrame([
        {"accession": a, "position": p, "uniprot_tier": tiers[(a, p)][0],
         "uniprot_evidence_codes": "", "uniprot_qualifies": tiers[(a, p)][1],
         "exclusion_reason": tiers[(a, p)][2]}
        for a, p in SITES
    ])


def glygen_frame() -> pd.DataFrame:
    tiers = {("P2", 20): "glygen_reported_with_glycan", ("P3", 30): "glygen_predicted"}
    return pd.DataFrame([
        {"accession": a, "position": p, "glygen_tier": tiers.get((a, p), ""),
         "glygen_categories": "", "glygen_evidence_databases": "",
         "glygen_pubmed_ids": "", "glygen_glytoucan_ids": ""}
        for a, p in SITES
    ])


def structure_frame() -> pd.DataFrame:
    tiers = {("P4", 40): "structure_linked_glycan", ("P3", 30): "structure_residue_resolved"}
    return pd.DataFrame([
        {"accession": a, "position": p,
         "structure_tier": tiers.get((a, p), "structure_not_assessed"),
         "structure_pdb_id": "", "structure_chain_id": "", "structure_resseq": None,
         "structure_icode": "", "structure_detail": ""}
        for a, p in SITES
    ])


def combined() -> pd.DataFrame:
    return combine_layers(
        uniprot_frame(), glygen_frame(), None, structure_frame(), POLICY
    ).set_index(["accession", "position"])


def test_uniprot_only_support():
    row = combined().loc[("P1", 10)]
    assert bool(row.experimental_positive) is True
    assert row.support_sources == "uniprot"
    assert row.occupancy_status == "occupied_supported"


def test_glygen_can_promote_a_site_uniprot_rejects():
    row = combined().loc[("P2", 20)]
    assert bool(row.experimental_positive) is True
    assert row.support_sources == "glygen"


def test_structure_link_can_promote_a_site_with_no_annotation():
    row = combined().loc[("P4", 40)]
    assert bool(row.experimental_positive) is True
    assert row.support_sources == "structure"


def test_unsupported_site_is_unknown_not_negative():
    row = combined().loc[("P3", 30)]
    assert bool(row.experimental_positive) is False
    assert row.occupancy_status == "unknown"
    assert row.support_sources == ""


def test_resolved_residue_without_glycan_does_not_confer_support():
    row = combined().loc[("P3", 30)]
    assert row.structure_tier == "structure_residue_resolved"
    assert bool(row.experimental_positive) is False


def test_glygen_predicted_never_confers_support():
    row = combined().loc[("P3", 30)]
    assert row.glygen_tier == "glygen_predicted"
    assert bool(row.experimental_positive) is False


def test_observed_unmodified_not_emitted_without_an_internal_control():
    """A bare residue is a silence unless the structure proves glycans were modelled."""
    assert "observed_unmodified" not in set(combined().occupancy_status)


def _structure_frame_with_control(controlled: set) -> pd.DataFrame:
    frame = structure_frame()
    frame["structure_internal_control"] = [
        (a, p) in controlled for a, p in zip(frame.accession, frame.position)
    ]
    return frame


def test_observed_unmodified_emitted_only_with_an_internal_control():
    result = combine_layers(
        uniprot_frame(), glygen_frame(), None,
        _structure_frame_with_control({("P3", 30)}), POLICY,
    ).set_index(["accession", "position"])
    assert result.loc[("P3", 30)].occupancy_status == "observed_unmodified"
    # P5 has no control and stays unknown rather than being swept along
    assert result.loc[("P5", 50)].occupancy_status == "unknown"


def test_supported_site_is_occupied_even_with_an_internal_control():
    """Occupied outranks observed-unmodified; a glycan seen anywhere wins."""
    result = combine_layers(
        uniprot_frame(), glygen_frame(), None,
        _structure_frame_with_control({("P1", 10), ("P3", 30)}), POLICY,
    ).set_index(["accession", "position"])
    assert result.loc[("P1", 10)].occupancy_status == "occupied_supported"


def test_observed_unmodified_can_be_disabled_by_policy():
    policy = dict(POLICY, observed_unmodified_from_internal_control=False)
    result = combine_layers(
        uniprot_frame(), glygen_frame(), None,
        _structure_frame_with_control({("P3", 30)}), policy,
    )
    assert "observed_unmodified" not in set(result.occupancy_status)


def test_support_count_reflects_independent_layers():
    result = combine_layers(
        uniprot_frame(), glygen_frame(), None, structure_frame(), POLICY
    ).set_index(["accession", "position"])
    assert int(result.loc[("P1", 10)].support_count) == 1
    assert int(result.loc[("P3", 30)].support_count) == 0


def test_multiple_layers_are_listed_in_stable_order():
    uniprot = uniprot_frame()
    glygen = glygen_frame()
    glygen.loc[glygen.accession == "P1", "glygen_tier"] = "glygen_reported_with_glycan"
    result = combine_layers(uniprot, glygen, None, structure_frame(), POLICY)
    row = result.set_index(["accession", "position"]).loc[("P1", 10)]
    assert row.support_sources == "uniprot|glygen"
    assert int(row.support_count) == 2


def test_exclusion_reason_cleared_when_another_layer_supports():
    row = combined().loc[("P2", 20)]
    assert row.exclusion_reason == ""


def test_exclusion_reason_retained_when_nothing_supports():
    row = combined().loc[("P3", 30)]
    assert row.exclusion_reason == "sequence_model_only"


def test_disabled_layers_are_tolerated():
    result = combine_layers(uniprot_frame(), None, None, None, POLICY)
    assert int(result.experimental_positive.sum()) == 1
    assert set(result.support_sources) <= {"uniprot", ""}
