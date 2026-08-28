"""ESM3's two knobs: the structure track, and what is masked at the sequon.

Hermetic — a fake model and tokeniser, so none of this needs EvolutionaryScale's
`esm`, a gated checkpoint, a GPU or the network.

The tests that matter are the ones asserting the manipulation really happens:
that `seq_only` withholds the structure tokens rather than merely being labelled
that way, and that `single` and `joint` mask different things. A mode that is
named but not applied would produce a full set of plausible numbers.
"""
from __future__ import annotations

import numpy as np
import pytest

from experimental_glycosylation_sites import esm3_scoring as e3

VOCAB = ["<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E",
         "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y", "M", "H", "W",
         "C", "X", "B", "U", "Z", "O", ".", "-", "|", "<mask>"]
MASK_ID = VOCAB.index("<mask>")


class FakeTokenizer:
    mask_token_id = MASK_ID

    def encode(self, text):
        return ([VOCAB.index("<cls>")] + [VOCAB.index(c) for c in text]
                + [VOCAB.index("<eos>")])

    def convert_ids_to_tokens(self, ids):
        return [VOCAB[i] for i in ids]

    def convert_tokens_to_ids(self, token):
        return VOCAB.index(token)


class FakeModel:
    """Records every forward call so the manipulation can be inspected."""

    def __init__(self):
        self.calls = []

    def forward(self, sequence_tokens=None, structure_tokens=None):
        import torch

        self.calls.append({"sequence": sequence_tokens.clone(),
                           "structure": None if structure_tokens is None
                           else structure_tokens.clone()})
        b, t = sequence_tokens.shape
        logits = torch.zeros((b, t, len(VOCAB)))
        return type("Out", (), {"sequence_logits": logits})()


def context_for(sequence="MANKSTV", with_structure=True):
    import torch

    tokenizer = FakeTokenizer()
    tokens = torch.tensor(tokenizer.encode(sequence))
    structure = torch.arange(len(tokens)) if with_structure else None
    return (sequence, tokens, structure), tokenizer


# --------------------------------------------------------------------------
# Provenance: the four modes must be distinguishable in the output.
# --------------------------------------------------------------------------

def test_every_mode_combination_has_its_own_conditioning_string():
    seen = {e3.conditioning(s, m)
            for s in e3.STRUCTURE_MODES for m in e3.MASK_MODES}
    assert len(seen) == 4
    assert e3.conditioning("struct_cond", "single") == \
        "masked_structure_conditioned_single"
    assert e3.conditioning("seq_only", "joint") == "masked_sequence_only_joint"


def test_the_structure_track_is_named_in_the_conditioning():
    """A structure-conditioned and a sequence-only run of the SAME model would
    otherwise be indistinguishable once written to disk."""
    assert "structure_conditioned" in e3.conditioning("struct_cond", "single")
    assert "sequence_only" in e3.conditioning("seq_only", "single")


@pytest.mark.parametrize("structure_mode,mask_mode",
                         [("nonsense", "single"), ("struct_cond", "nonsense")])
def test_an_unknown_mode_is_refused(structure_mode, mask_mode):
    with pytest.raises(ValueError):
        e3._check(structure_mode, mask_mode)


# --------------------------------------------------------------------------
# The manipulation actually happens.
# --------------------------------------------------------------------------

def test_seq_only_withholds_the_structure_tokens():
    """The whole value of this model is that the track can be switched off."""
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 structure_mode="seq_only", mask_mode="single")
    assert all(c["structure"] is None for c in model.calls)


def test_struct_cond_passes_the_structure_tokens():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 structure_mode="struct_cond", mask_mode="single")
    assert all(c["structure"] is not None for c in model.calls)
    # one structure row per masked sequence row
    assert model.calls[0]["structure"].shape[0] == model.calls[0]["sequence"].shape[0]


def test_single_masks_one_position_per_row():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 mask_mode="single")
    batch = model.calls[0]["sequence"]
    assert batch.shape[0] == 3, "one row per scored position"
    for row, index in enumerate((2, 3, 4)):
        masked = (batch[row] == MASK_ID).nonzero().flatten().tolist()
        assert masked == [index + e3.TOKEN_OFFSET]


def test_joint_masks_all_three_at_once_in_a_single_row():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 mask_mode="joint")
    batch = model.calls[0]["sequence"]
    assert batch.shape[0] == 1, "joint needs one pass, not three"
    masked = (batch[0] == MASK_ID).nonzero().flatten().tolist()
    assert masked == [2 + e3.TOKEN_OFFSET, 3 + e3.TOKEN_OFFSET,
                      4 + e3.TOKEN_OFFSET]


def test_the_token_offset_is_where_the_masks_land():
    """Position i is token i+1 because <cls> is prepended. An error here masks
    a neighbouring residue and reads a plausible number for the wrong site."""
    context, tokenizer = context_for("MANKSTV")
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (0,), mask_mode="single")
    batch = model.calls[0]["sequence"]
    assert int(batch[0][0]) == VOCAB.index("<cls>")
    assert int(batch[0][1]) == MASK_ID, "position 0 must mask token 1"


def test_a_wrong_token_offset_is_caught_at_load():
    class Shifted(FakeTokenizer):
        def encode(self, text):
            return [VOCAB.index(c) for c in text]      # no <cls>
    with pytest.raises(RuntimeError, match="token offset"):
        e3._assert_token_offset(Shifted())


def test_the_correct_token_offset_is_accepted():
    e3._assert_token_offset(FakeTokenizer())


# --------------------------------------------------------------------------
# Score and guards.
# --------------------------------------------------------------------------

def rows(**overrides):
    base = {i: np.full(len(VOCAB), 1.0 / len(VOCAB)) for i in (0, 1, 2)}
    base.update(overrides)
    return base


def test_the_score_matches_a_hand_calculation():
    tokenizer = FakeTokenizer()
    n_row = np.zeros(len(VOCAB)); n_row[VOCAB.index("N")] = 0.5
    n_row[VOCAB.index("A")] = 0.5
    p1 = np.zeros(len(VOCAB)); p1[VOCAB.index("P")] = 0.1; p1[VOCAB.index("G")] = 0.9
    p2 = np.zeros(len(VOCAB)); p2[VOCAB.index("S")] = 0.5
    p2[VOCAB.index("T")] = 0.25; p2[VOCAB.index("A")] = 0.25

    score = e3.sequon_score({0: n_row, 1: p1, 2: p2}, tokenizer, 0, 1, 2)
    assert score["p_asn_at_n"] == pytest.approx(0.5)
    assert score["p_ser_or_thr_at_plus2"] == pytest.approx(0.75)
    assert score["p_pro_at_plus1"] == pytest.approx(0.1)
    assert score["conditional_sequon_score"] == pytest.approx(np.log(3.0) / 2)
    assert score["conditional_sequon_score_sd"] == 0.0
    assert score["n_decoding_orders"] == 1


def test_an_unnormalised_row_is_refused():
    bad = rows(); bad[1] = bad[1] * 3
    with pytest.raises(e3.InvalidProbabilityVector):
        e3.check_scoreable(bad, (0, 1, 2))


def test_a_missing_position_is_refused():
    with pytest.raises(e3.InvalidProbabilityVector):
        e3.check_scoreable({0: np.full(len(VOCAB), 1.0 / len(VOCAB))}, (0, 1, 2))


def test_a_triplet_mismatch_is_refused():
    e3.check_triplet("MANKSTV", (2, 3, 4), "NKS")
    with pytest.raises(e3.SequonMismatchError):
        e3.check_triplet("MANKSTV", (2, 3, 4), "NQT")
