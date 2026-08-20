from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from experimental_glycosylation_sites.mpnn_scoring import (
    AA_INDEX,
    ALPHABET,
    EPSILON,
    logit,
    sequon_score,
)


def probs_array(orders: int, length: int, assignments: dict) -> np.ndarray:
    """Uniform probabilities with specific (position, residue) values forced in.

    The remainder of each row is spread evenly so rows still sum to one, which is
    what the scorer receives from the model.
    """
    a = np.zeros((orders, length, 21))
    for o in range(orders):
        for pos in range(length):
            forced = assignments.get(pos, {})
            spare = 1.0 - sum(forced.values())
            free = 21 - len(forced)
            a[o, pos, :] = spare / free
            for aa, p in forced.items():
                a[o, pos, AA_INDEX[aa]] = p
    return a


def test_alphabet_matches_proteinmpnn_order():
    """Indexing with a hand-written ordering would silently read the wrong residue."""
    # Assert against ProteinMPNN's own source, not a copy of it. The previous
    # version of this test hard-coded the string the module happened to hold,
    # so it locked the wrong alphabet in place instead of catching it.
    assert ALPHABET == "ACDEFGHIKLMNPQRSTVWYX"
    assert AA_INDEX["N"] == 11 and AA_INDEX["S"] == 15 and AA_INDEX["T"] == 16
    assert AA_INDEX["P"] == 12
    assert [ALPHABET[i] for i in (0, 20)] == ["A", "X"]


def test_logit_is_symmetric_and_clamped():
    assert logit(0.5) == 0.0
    assert logit(0.75) == pytest.approx(-logit(0.25))
    assert math.isfinite(logit(0.0)) and math.isfinite(logit(1.0))
    assert logit(0.0) == pytest.approx(math.log(EPSILON / (1 - EPSILON)))


def test_score_averages_log_odds_of_n_and_s_or_t():
    a = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    result = sequon_score(a, n_index=1, plus1_index=2, plus2_index=3)
    # both terms are logit(0.5) = 0, so the mean is 0
    assert result["conditional_sequon_score"] == pytest.approx(0.0, abs=1e-9)
    assert result["p_asn_at_n"] == pytest.approx(0.5)
    assert result["p_ser_or_thr_at_plus2"] == pytest.approx(0.5)


def test_serine_and_threonine_are_pooled_at_plus2():
    """Either residue completes a sequon, so the score uses their sum."""
    only_s = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.4, "T": 0.0}})
    split = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.2, "T": 0.2}})
    a = sequon_score(only_s, 1, 2, 3)["conditional_sequon_score"]
    b = sequon_score(split, 1, 2, 3)["conditional_sequon_score"]
    assert a == pytest.approx(b)


def test_middle_residue_does_not_enter_the_primary_score():
    """Any residue but proline permits a sequon, so +1 is diagnostic only."""
    low = probs_array(1, 5, {1: {"N": 0.5}, 2: {"P": 0.01}, 3: {"S": 0.5}})
    high = probs_array(1, 5, {1: {"N": 0.5}, 2: {"P": 0.90}, 3: {"S": 0.5}})
    assert (sequon_score(low, 1, 2, 3)["conditional_sequon_score"]
            == pytest.approx(sequon_score(high, 1, 2, 3)["conditional_sequon_score"]))
    # but it is still reported
    assert sequon_score(high, 1, 2, 3)["p_pro_at_plus1"] == pytest.approx(0.90)


def test_variation_across_decoding_orders_is_retained():
    """Conditional probabilities depend on decoding order; the spread is evidence."""
    # Each order is built as a whole distribution rather than by overwriting one
    # entry of a finished array: that shortcut leaves the row summing to 1.3, and
    # the scorer now rejects rows that are not distributions.
    a = np.concatenate([
        probs_array(1, 5, {1: {"N": p_n}, 3: {"S": 0.5}})
        for p_n in (0.5, 0.8, 0.2)
    ])
    result = sequon_score(a, 1, 2, 3)
    assert result["n_decoding_orders"] == 3
    assert result["conditional_sequon_score_sd"] > 0


def test_full_probability_vectors_are_returned():
    a = probs_array(2, 5, {1: {"N": 0.5}, 3: {"S": 0.5}})
    result = sequon_score(a, 1, 2, 3)
    for key in ("probs_n", "probs_plus1", "probs_plus2"):
        assert len(result[key]) == 21
        assert sum(result[key]) == pytest.approx(1.0, abs=1e-6)


def test_higher_motif_probability_gives_a_higher_score():
    weak = probs_array(1, 5, {1: {"N": 0.05}, 3: {"S": 0.05}})
    strong = probs_array(1, 5, {1: {"N": 0.80}, 3: {"S": 0.80}})
    assert (sequon_score(strong, 1, 2, 3)["conditional_sequon_score"]
            > sequon_score(weak, 1, 2, 3)["conditional_sequon_score"])


# --------------------------------------------------------------------------
# Incomplete backbones
#
# ProteinMPNN's conditional_probs fills only the positions surviving its own
# chain_M * mask and returns zeros elsewhere. Exponentiating a zero row gives
# twenty-one ones, P(S)+P(T) = 2, and a score near +13.8. This corrupted 105 of
# 2,564 scored sites before it was caught, and inflated the reference SD from
# 1.33 to 2.62, so the guard is tested from both directions: a row that is not a
# distribution, and a position the model declined to decode.
# --------------------------------------------------------------------------

def test_all_ones_row_is_rejected():
    """The exact corruption: an unfilled log-prob row exponentiated to ones."""
    from experimental_glycosylation_sites.mpnn_scoring import InvalidProbabilityVector

    corrupted = np.ones((1, 5, 21))
    with pytest.raises(InvalidProbabilityVector, match="sums to 21"):
        sequon_score(corrupted, 1, 2, 3)


def test_one_corrupted_row_among_valid_ones_is_rejected():
    """A site fails if any of its three residues is bad, not only the first."""
    from experimental_glycosylation_sites.mpnn_scoring import InvalidProbabilityVector

    a = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    a[0, 3, :] = 1.0                       # the +2 residue alone is unfilled
    with pytest.raises(InvalidProbabilityVector, match="model index 3"):
        sequon_score(a, 1, 2, 3)


def test_corruption_in_a_later_decoding_order_is_rejected():
    """Averaging hides a bad order; the check runs per order, before averaging."""
    from experimental_glycosylation_sites.mpnn_scoring import InvalidProbabilityVector

    a = probs_array(8, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    a[6, 1, :] = 1.0
    with pytest.raises(InvalidProbabilityVector, match="decoding order 6"):
        sequon_score(a, 1, 2, 3)


def test_undecoded_position_is_rejected_even_when_rows_look_valid():
    """The two guards are independent: `computed` catches what the sums cannot."""
    from experimental_glycosylation_sites.mpnn_scoring import IncompleteBackboneError

    a = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    computed = np.array([True, True, True, False, True])
    with pytest.raises(IncompleteBackboneError, match=r"\[3\]"):
        sequon_score(a, 1, 2, 3, computed=computed)


def test_valid_site_still_scores_with_both_guards_active():
    """The guard must not reject good data."""
    a = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    result = sequon_score(a, 1, 2, 3, computed=np.ones(5, dtype=bool))
    assert result["conditional_sequon_score"] == pytest.approx(0.0, abs=1e-9)


def test_probabilities_slightly_off_one_are_tolerated():
    """float32 softmax noise must not be mistaken for corruption."""
    a = probs_array(1, 5, {1: {"N": 0.5}, 3: {"S": 0.25, "T": 0.25}})
    a *= 1.0 + 5e-5
    sequon_score(a, 1, 2, 3)               # does not raise


# --------------------------------------------------------------------------
# Integration: a real structure with a real missing backbone atom.
#
# The unit tests above construct the corrupted array by hand, which only proves
# the guard rejects what we think the model returns. This runs ProteinMPNN on a
# 21-residue window of real backbone geometry whose sequon +2 residue is missing
# its O, and asserts the model does in fact produce the all-ones row there.
# --------------------------------------------------------------------------

FIXTURE = Path(__file__).parent / "fixtures" / "incomplete_backbone.pdb"
MPNN_DIR = Path(__file__).resolve().parents[3] / "ProteinMPNN"
SEQUON_INDICES = (9, 10, 11)               # +2 residue at 11 has no O

requires_mpnn = pytest.mark.skipif(
    not (MPNN_DIR / "vanilla_model_weights" / "v_48_020.pt").exists(),
    reason="ProteinMPNN checkpoint not available",
)


@requires_mpnn
def test_proteinmpnn_leaves_incomplete_backbone_positions_undecoded():
    """Reproduces the defect end to end, then confirms the scorer refuses it."""
    from experimental_glycosylation_sites.mpnn_scoring import (
        IncompleteBackboneError, conditional_probabilities, load_model)

    model = load_model(MPNN_DIR)
    probs, computed = conditional_probabilities(
        FIXTURE, "A", model, n_decoding_orders=2, seed=0,
        positions=list(SEQUON_INDICES))

    n_index, plus1_index, plus2_index = SEQUON_INDICES

    # the model's own mask excludes only the residue missing its O
    assert computed[n_index] and computed[plus1_index]
    assert not computed[plus2_index]

    # and the row it returns there is the all-ones vector, not a distribution
    assert probs[0, n_index].sum() == pytest.approx(1.0, abs=1e-3)
    assert probs[0, plus2_index].sum() == pytest.approx(21.0, abs=1e-3)
    p_ser_or_thr = (probs[0, plus2_index, AA_INDEX["S"]]
                    + probs[0, plus2_index, AA_INDEX["T"]])
    assert p_ser_or_thr == pytest.approx(2.0, abs=1e-3)

    with pytest.raises(IncompleteBackboneError):
        sequon_score(probs, *SEQUON_INDICES, computed=computed)


@requires_mpnn
def test_complete_backbone_positions_score_normally():
    """The same structure scores where the backbone is intact."""
    from experimental_glycosylation_sites.mpnn_scoring import (
        conditional_probabilities, load_model)

    model = load_model(MPNN_DIR)
    intact = (2, 3, 4)
    probs, computed = conditional_probabilities(
        FIXTURE, "A", model, n_decoding_orders=2, seed=0, positions=list(intact))

    assert all(computed[i] for i in intact)
    result = sequon_score(probs, *intact, computed=computed)
    assert math.isfinite(result["conditional_sequon_score"])
    assert 0.0 <= result["p_ser_or_thr_at_plus2"] <= 1.0
