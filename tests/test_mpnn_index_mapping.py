"""Mapping manifest indices onto ProteinMPNN's own parse.

`model_index` is an ordinal into the chain as `structures._parse_chains` reads
it: observed residues only, in order. ProteinMPNN's `parse_PDB` instead walks
the contiguous author residue-number range, inserting a placeholder for every
absent number and expanding insertion codes:

    for resn in range(min_resn, max_resn+1):
        if resn in seq:
            for k in sorted(seq[resn]): seq_.append(...)
        else: seq_.append(20)

So the two indices coincide only when the chain is numbered without gaps. The
mapping is arithmetic on residue numbers, not a sequence alignment -- and it is
verified against ProteinMPNN's own sequence before being trusted.
"""
from __future__ import annotations

from experimental_glycosylation_sites.mpnn_scoring import build_index_map


def test_contiguous_numbering_maps_one_to_one():
    ids = [(1, ""), (2, ""), (3, "")]
    assert build_index_map(ids, "ACD", "ACD") == {0: 0, 1: 1, 2: 2}


def test_a_numbering_gap_shifts_everything_after_it():
    """Residue 3 is absent, so ProteinMPNN puts an X where it would have been."""
    ids = [(1, ""), (2, ""), (4, "")]
    assert build_index_map(ids, "ACD", "ACXD") == {0: 0, 1: 1, 2: 3}


def test_a_wide_gap_shifts_by_its_full_width():
    ids = [(1, ""), (10, "")]
    model = "A" + "X" * 8 + "C"
    assert build_index_map(ids, "AC", model) == {0: 0, 1: 9}


def test_insertion_codes_occupy_their_own_slots():
    """36, 36A, 36B are three positions in both parses."""
    ids = [(36, ""), (36, "A"), (36, "B"), (37, "")]
    assert build_index_map(ids, "LKND", "LKND") == {0: 0, 1: 1, 2: 2, 3: 3}


def test_mapping_is_rejected_when_it_disagrees_with_the_model_sequence():
    """The arithmetic is a hypothesis; the model's own sequence settles it."""
    ids = [(1, ""), (2, ""), (3, "")]
    assert build_index_map(ids, "ACD", "WWW") == {}


def test_empty_inputs_give_an_empty_map():
    assert build_index_map([], "", "") == {}


def test_designs_are_projected_into_manifest_index_space():
    """`design()` must return sequences the manifest's indices address.

    ESM-IF's adapter documents this contract -- "returned in the manifest's
    index space" -- and honours it. ProteinMPNN returned raw model-space
    sequences, so retention was read at unmapped positions.
    """
    from experimental_glycosylation_sites.mpnn_scoring import to_manifest_space
    # manifest 0,1,2 -> model 0,1,3 (a placeholder sits at model index 2)
    mapping = {0: 0, 1: 1, 2: 3}
    assert to_manifest_space("ACXD", mapping) == "ACD"


def test_projection_is_identity_for_a_gapless_chain():
    from experimental_glycosylation_sites.mpnn_scoring import to_manifest_space
    assert to_manifest_space("ACDE", {i: i for i in range(4)}) == "ACDE"


def test_unmapped_manifest_positions_become_X():
    """An unresolvable position must not silently borrow a neighbour."""
    from experimental_glycosylation_sites.mpnn_scoring import to_manifest_space
    assert to_manifest_space("AC", {0: 0, 2: 1}) == "AXC"


def test_projection_of_an_empty_mapping_is_empty():
    from experimental_glycosylation_sites.mpnn_scoring import to_manifest_space
    assert to_manifest_space("ACDE", {}) == ""
