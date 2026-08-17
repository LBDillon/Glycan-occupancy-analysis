from __future__ import annotations

import math

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
    assert ALPHABET == "ARNDCQEGHILKMFPSTWYVX"
    assert AA_INDEX["N"] == 2 and AA_INDEX["S"] == 15 and AA_INDEX["T"] == 16
    assert AA_INDEX["P"] == 14


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
    a = probs_array(3, 5, {1: {"N": 0.5}, 3: {"S": 0.5}})
    a[1, 1, AA_INDEX["N"]] = 0.8
    a[2, 1, AA_INDEX["N"]] = 0.2
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
