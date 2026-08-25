"""Holding positions fixed during ProteinMPNN design.

ProteinMPNN decides what to design from `chain_M_pos`: 1 means design this
position, 0 means keep the native residue. Fixing a sequon is therefore a mask,
not a post-hoc repair of the output -- repairing afterwards would let the model
condition on residues it was going to change, which is a different experiment.

The failure mode these guard against is the one that cost a day already: an
off-by-one or a wrong index silently protects the wrong residue, and the design
still looks perfectly reasonable.
"""
from __future__ import annotations

import numpy as np
import pytest

from glyco_context.fixed_design import design_mask, sequon_positions


def test_mask_is_all_design_when_nothing_is_fixed():
    assert list(design_mask(5, [])) == [1.0, 1.0, 1.0, 1.0, 1.0]


def test_mask_zeroes_exactly_the_fixed_positions():
    assert list(design_mask(5, [1, 3])) == [1.0, 0.0, 1.0, 0.0, 1.0]


def test_mask_rejects_a_position_outside_the_chain():
    """Silently ignoring it would protect nothing and report success."""
    with pytest.raises(ValueError, match="outside"):
        design_mask(5, [7])


def test_mask_rejects_a_negative_position():
    with pytest.raises(ValueError, match="outside"):
        design_mask(5, [-1])


def test_mask_tolerates_duplicate_positions():
    assert list(design_mask(4, [2, 2])) == [1.0, 1.0, 0.0, 1.0]


def test_full_sequon_policy_fixes_all_three_residues():
    """SugarFix's `full_sequon`: the Asn, the X and the hydroxyl."""
    assert sequon_positions(10, policy="full_sequon") == [10, 11, 12]


def test_functional_preserve_policy_leaves_x_free():
    """SugarFix's `functional_preserve`: the motif stays valid, X may change."""
    assert sequon_positions(10, policy="functional_preserve") == [10, 12]


def test_unknown_policy_is_rejected():
    with pytest.raises(ValueError, match="policy"):
        sequon_positions(10, policy="whatever")


def test_sequon_index_verifies_against_the_sequence():
    from glyco_context.fixed_design import verify_sequon_index
    seq = "AAANGTAAA"
    assert verify_sequon_index(seq, 3, "NGT")
    assert not verify_sequon_index(seq, 4, "NGT"), "off by one must be caught"
    assert not verify_sequon_index(seq, 3, "NAT"), "wrong expected triplet"
    assert not verify_sequon_index("AAANPTAAA", 3, "NPT"), "N-P-T is not a sequon"
    assert not verify_sequon_index(seq, 100, "NGT"), "past the end"
