"""ESMC's guards, tested without the 300M checkpoint.

As with ESM-IF, what is tested here is the code that decides *which residue gets
scored* and *what is read out of the distribution* — the parts with a history in
this project. The model itself is exercised on Colab.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from experimental_glycosylation_sites.esmc_scoring import (
    DEFAULT_MASK_MODE, MASK_MODES, N_ORDERS, SCORE_SD, TOKEN_OFFSET,
    InvalidProbabilityVector, check_scoreable, conditioning, sequon_score)
from experimental_glycosylation_sites.mpnn_scoring import \
    sequon_score as mpnn_sequon_score

def _have_esm_sdk() -> bool:
    """Is EvolutionaryScale's `esm` importable, rather than fair-esm's?

    `find_spec` on a submodule imports its parent, so with fair-esm installed
    this RAISES rather than returning None -- `esm` exists but `esm.models` does
    not. Both packages claim the import name `esm`, so the question is genuinely
    "which one is on the path", and it has to be asked defensively.
    """
    try:
        return importlib.util.find_spec("esm.models.esmc") is not None
    except (ImportError, AttributeError, ValueError):
        return False


HAVE_ESM_SDK = _have_esm_sdk()


class FakeTokenizer:
    """A vocabulary in a deliberately non-alphabetical order.

    Different from both ProteinMPNN's ordering and ESM-IF's, so a test passing
    here cannot be passing because two orderings coincide.
    """

    TOKENS = ["<cls>", "<pad>", "<eos>"] + list("LAGVSERTIDPKQNFYMHWC") + ["<mask>"]

    def convert_tokens_to_ids(self, token):
        return self.TOKENS.index(token)

    def convert_ids_to_tokens(self, ids):
        return [self.TOKENS[i] for i in ids]

    @property
    def mask_token_id(self):
        return self.TOKENS.index("<mask>")

    def __len__(self):
        return len(self.TOKENS)


def vectors(tokenizer, assignments):
    """Rows summing to one, with specific (position, residue) values forced in."""
    width = len(tokenizer)
    out = {}
    for position, forced in assignments.items():
        row = np.zeros(width)
        spare = 1.0 - sum(forced.values())
        free = width - len(forced)
        row[:] = spare / free
        for residue, probability in forced.items():
            row[tokenizer.convert_tokens_to_ids(residue)] = probability
        out[position] = row
    return out


# --- what gets read out of the distribution ------------------------------

def test_score_reads_the_residues_it_claims_to():
    t = FakeTokenizer()
    v = vectors(t, {3: {"N": 0.7}, 4: {"P": 0.25}, 5: {"S": 0.3, "T": 0.4}})
    r = sequon_score(v, t, 3, 4, 5)
    assert r["p_asn_at_n"] == pytest.approx(0.7)
    assert r["p_ser_at_plus2"] == pytest.approx(0.3)
    assert r["p_thr_at_plus2"] == pytest.approx(0.4)
    assert r["p_ser_or_thr_at_plus2"] == pytest.approx(0.7)
    assert r["p_pro_at_plus1"] == pytest.approx(0.25)


def test_score_is_the_mean_of_the_two_log_odds():
    t = FakeTokenizer()
    v = vectors(t, {3: {"N": 0.5}, 4: {}, 5: {"S": 0.25, "T": 0.25}})
    assert sequon_score(v, t, 3, 4, 5)["conditional_sequon_score"] == pytest.approx(0.0, abs=1e-9)


def test_column_names_match_proteinmpnn_exactly():
    """Downstream stays model-agnostic only while every scorer emits one schema."""
    t = FakeTokenizer()
    esmc = sequon_score(vectors(t, {3: {}, 4: {}, 5: {}}), t, 3, 4, 5)
    mpnn = mpnn_sequon_score(np.full((2, 6, 21), 1.0 / 21), 3, 4, 5)
    assert set(esmc) == set(mpnn)


def test_deterministic_pass_reports_one_order_and_no_spread():
    t = FakeTokenizer()
    r = sequon_score(vectors(t, {3: {}, 4: {}, 5: {}}), t, 3, 4, 5)
    assert r["n_decoding_orders"] == N_ORDERS == 1
    assert r["conditional_sequon_score_sd"] == SCORE_SD == 0.0


# --- the distribution guard ----------------------------------------------

def test_rejects_unnormalised_rows():
    t = FakeTokenizer()
    v = vectors(t, {3: {}, 4: {}, 5: {}})
    v[4] = np.ones(len(t))
    with pytest.raises(InvalidProbabilityVector, match="not 1"):
        check_scoreable(v)


def test_rejects_negative_values_that_still_sum_to_one():
    t = FakeTokenizer()
    v = vectors(t, {3: {}, 4: {}, 5: {}})
    row = np.zeros(len(t)); row[0] = -0.5; row[1] = 1.5
    v[5] = row
    with pytest.raises(InvalidProbabilityVector, match=r"outside \[0, 1\]"):
        check_scoreable(v)


# --- masking modes --------------------------------------------------------

def test_conditioning_names_the_mask_mode():
    assert conditioning("single") == "masked_sequence_single"
    assert conditioning("joint") == "masked_sequence_joint"
    assert DEFAULT_MASK_MODE == "single" and set(MASK_MODES) == {"single", "joint"}


def test_adapter_rejects_an_unknown_mask_mode():
    from experimental_glycosylation_sites.adapters.esmc import ESMCAdapter

    with pytest.raises(ValueError, match="mask_mode"):
        ESMCAdapter(mask_mode="everything")


def test_adapter_is_a_scorer_but_not_a_designer():
    """ESMC is a masked LM: sampling from it would not be backbone-conditioned."""
    from experimental_glycosylation_sites.adapters.base import (SequenceDesigner,
                                                                SequonScorer)
    from experimental_glycosylation_sites.adapters.esmc import ESMCAdapter

    a = ESMCAdapter()
    assert isinstance(a, SequonScorer)
    assert not isinstance(a, SequenceDesigner)


def test_describe_reports_the_mask_mode():
    from experimental_glycosylation_sites.adapters.esmc import ESMCAdapter

    assert ESMCAdapter(mask_mode="joint").describe()["conditioning"] == "masked_sequence_joint"


# --- model-dependent ------------------------------------------------------

@pytest.mark.skipif(not HAVE_ESM_SDK, reason="EvolutionaryScale esm not installed")
def test_token_offset_round_trips_against_the_real_tokenizer():
    """The check the ProteinMPNN alphabet defect went undetected for."""
    from esm.models.esmc import ESMC

    from experimental_glycosylation_sites.esmc_scoring import _assert_token_offset

    tokenizer = ESMC.from_pretrained("esmc_300m").tokenizer
    _assert_token_offset(tokenizer)                      # raises if wrong
    probe = "MKTAYIAKQRNLTSHFSRQ"
    ids = tokenizer(probe, return_tensors="pt")["input_ids"][0].tolist()
    assert tokenizer.convert_ids_to_tokens([ids[TOKEN_OFFSET]])[0] == probe[0]
