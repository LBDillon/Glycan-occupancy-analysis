"""ProGen2's tokenisation, alignment and score.

Hermetic: a fake tokeniser and a fake model, so none of this needs
`transformers`, a checkpoint download, a GPU or the network.

The alignment test is the one that matters. Everything downstream assumes
`logits[i]` predicts the residue at manifest index `i`, which holds only because
BOS is prepended and the final row dropped. That is the same class of assumption
as ESMC's token offset and ProteinMPNN's residue numbering, both of which were
wrong until checked.
"""
from __future__ import annotations

import numpy as np
import pytest

from experimental_glycosylation_sites import progen2_scoring as pg

VOCAB = ["<|pad|>", "<|bos|>", "<|eos|>", "1", "2"] + list(pg.STANDARD_AA)


class FakeTokenizer:
    """One token per residue, as ProGen2 does."""

    bos_token_id = 1

    def __init__(self, vocab=None):
        self.vocab = list(vocab if vocab is not None else VOCAB)
        self.index = {t: i for i, t in enumerate(self.vocab)}

    def __call__(self, text, add_special_tokens=True, **kwargs):
        return {"input_ids": [self.index[c] for c in text]}

    def convert_ids_to_tokens(self, ids):
        return [self.vocab[i] for i in ids]


class MergingTokenizer(FakeTokenizer):
    """A BPE-style tokeniser that merges residues, as ProtGPT2's does."""

    def __call__(self, text, add_special_tokens=True, **kwargs):
        return {"input_ids": [0] * ((len(text) + 1) // 2)}


class FakeModel:
    """Emits a chosen residue at each position, so alignment is checkable."""

    def __init__(self, favours, vocab_size=len(VOCAB)):
        self.favours = favours          # what row i should predict
        self.vocab_size = vocab_size
        self.seen = None

    def __call__(self, input_ids=None, **kwargs):
        import torch

        self.seen = input_ids
        length = input_ids.shape[1]
        logits = torch.full((1, length, self.vocab_size), -10.0)
        # row i predicts token i+1, so favour the residue at manifest index i
        for i, aa in enumerate(self.favours):
            if i < length:
                logits[0, i, VOCAB.index(aa)] = 10.0
        return type("Out", (), {"logits": logits})()


def probabilities_for(sequence, favours=None):
    tokenizer = FakeTokenizer()
    model = FakeModel(favours if favours is not None else sequence)
    return pg.conditional_probabilities(sequence, model, tokenizer,
                                        tokenizer.bos_token_id), tokenizer, model


# --------------------------------------------------------------------------
# Tokenisation — the ProtGPT2 failure mode, caught rather than assumed away.
# --------------------------------------------------------------------------

def test_a_per_residue_tokeniser_is_accepted():
    index = pg.verify_tokenisation(FakeTokenizer())
    assert set(index) == set(pg.STANDARD_AA)
    assert index["N"] == VOCAB.index("N")


def test_a_merging_tokeniser_is_refused():
    """ProtGPT2's BPE vocabulary merges residues, so positions cannot align.

    Refused loudly here rather than producing an off-by-some that still returns
    a plausible number.
    """
    with pytest.raises(pg.TokenisationError, match="one token per residue"):
        pg.verify_tokenisation(MergingTokenizer())


def test_a_tokeniser_whose_ids_do_not_decode_back_is_refused():
    scrambled = FakeTokenizer()
    scrambled.vocab[scrambled.index["N"]] = "Q"
    with pytest.raises(pg.TokenisationError, match="decodes to"):
        pg.verify_tokenisation(scrambled)


# --------------------------------------------------------------------------
# Alignment.
# --------------------------------------------------------------------------

def test_row_i_predicts_the_residue_at_manifest_index_i():
    """With BOS prepended, no further shift is needed. This is the measurement."""
    sequence = "MANKSTV"
    probabilities, tokenizer, _ = probabilities_for(sequence)
    assert probabilities.shape[0] == len(sequence)
    recovered = "".join(VOCAB[int(row.argmax())] for row in probabilities)
    assert recovered == sequence


def test_bos_is_prepended_and_the_last_row_dropped():
    sequence = "MANKS"
    probabilities, _, model = probabilities_for(sequence)
    assert model.seen.shape[1] == len(sequence) + 1, "BOS was not prepended"
    assert int(model.seen[0, 0]) == FakeTokenizer.bos_token_id
    # the final logit row predicts a residue past the chain end
    assert probabilities.shape[0] == len(sequence)


def test_rows_are_probability_distributions():
    probabilities, _, _ = probabilities_for("MANKSTV")
    assert np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-5)


def test_an_off_by_one_would_be_caught():
    """A model whose rows are shifted must not reproduce the sequence."""
    sequence = "MANKSTV"
    shifted = sequence[1:] + "A"
    probabilities, _, _ = probabilities_for(sequence, favours=shifted)
    recovered = "".join(VOCAB[int(row.argmax())] for row in probabilities)
    assert recovered != sequence


# --------------------------------------------------------------------------
# Score and guards.
# --------------------------------------------------------------------------

def uniform_rows(n=5):
    row = np.full(len(VOCAB), 1.0 / len(VOCAB))
    return np.tile(row, (n, 1))


def test_the_score_matches_a_hand_calculation():
    index = {aa: VOCAB.index(aa) for aa in pg.STANDARD_AA}
    rows = np.zeros((3, len(VOCAB)))
    rows[0, index["N"]] = 0.5
    rows[0, index["A"]] = 0.5
    rows[1, index["P"]] = 0.1
    rows[1, index["G"]] = 0.9
    rows[2, index["S"]] = 0.5
    rows[2, index["T"]] = 0.25
    rows[2, index["A"]] = 0.25

    score = pg.sequon_score(rows, index, 0, 1, 2)
    assert score["p_asn_at_n"] == pytest.approx(0.5)
    assert score["p_ser_or_thr_at_plus2"] == pytest.approx(0.75)
    assert score["p_pro_at_plus1"] == pytest.approx(0.1)
    assert score["conditional_sequon_score"] == pytest.approx(np.log(3.0) / 2)
    assert score["conditional_sequon_score_sd"] == 0.0
    assert score["n_decoding_orders"] == 1


def test_an_unnormalised_row_is_refused():
    rows = uniform_rows()
    rows[1] *= 2
    with pytest.raises(pg.InvalidProbabilityVector):
        pg.check_scoreable(rows, (0, 1, 2))


def test_an_out_of_range_index_is_refused():
    with pytest.raises(IndexError):
        pg.check_scoreable(uniform_rows(3), (0, 1, 99))


def test_a_triplet_mismatch_is_refused():
    pg.check_triplet("MANKSTV", (2, 3, 4), "NKS")
    with pytest.raises(pg.SequonMismatchError):
        pg.check_triplet("MANKSTV", (2, 3, 4), "NQT")


def test_the_conditioning_is_esm_ifs_string_not_esmcs():
    """Causal and prefix-only, so it must never be pooled with a masked LM."""
    from experimental_glycosylation_sites.esmif_scoring import CONDITIONING as ESMIF
    assert pg.CONDITIONING == ESMIF == "autoregressive_prefix"
