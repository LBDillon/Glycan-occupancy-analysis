"""The ESM-IF adapter's guards, tested without the 1.7 GB checkpoint.

Everything here exercises the code that decides *which residue gets scored* —
the index mapping, the triplet check, the distribution check. That is the part
with a history: the manifest's indices come from Biopython and ESM-IF's come
from biotite, and a silent disagreement scores the wrong residue while still
returning a plausible number.

Tests needing the real model are marked and skip when it is absent.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from experimental_glycosylation_sites.esmif_scoring import (
    CONDITIONING, N_ORDERS, SCORE_SD, ChainMapping, IncompleteBackboneError,
    InvalidProbabilityVector, SequonMismatchError, check_scoreable, sequon_score)
from experimental_glycosylation_sites.mpnn_scoring import \
    sequon_score as mpnn_sequon_score

HAVE_ESM = importlib.util.find_spec("esm") is not None
requires_esm = pytest.mark.skipif(not HAVE_ESM, reason="fair-esm not installed")


class FakeAlphabet:
    """Stands in for ESM-IF's dictionary: a distinct, non-alphabetical order.

    Deliberately not the same ordering as ProteinMPNN's, so a test that passes
    here cannot be passing because two orderings happen to coincide.
    """

    ORDER = "<pad><eos>LAGVSERTIDPKQNFYMHWC"

    def __init__(self):
        self.tokens = ["<pad>", "<eos>"] + list("LAGVSERTIDPKQNFYMHWC")

    def get_idx(self, token):
        return self.tokens.index(token)

    def __len__(self):
        return len(self.tokens)


def uniform_probs(length, alphabet, assignments):
    """Rows that sum to one, with specific (position, residue) values forced in."""
    width = len(alphabet)
    array = np.zeros((length, width))
    for position in range(length):
        forced = assignments.get(position, {})
        spare = 1.0 - sum(forced.values())
        free = width - len(forced)
        array[position, :] = spare / free
        for residue, probability in forced.items():
            array[position, alphabet.get_idx(residue)] = probability
    return array


def mapping(to_esm, esm_seq, finite=None, manifest_length=None):
    length = len(esm_seq)
    if finite is None:
        finite = np.ones(length, dtype=bool)
    return ChainMapping(to_esm, esm_seq, finite, np.zeros((length, 3, 3)),
                        manifest_length if manifest_length is not None else length)


# --- the index mapping ----------------------------------------------------

def test_identity_mapping_returns_the_same_indices():
    m = mapping({i: i for i in range(10)}, "AAANLTAAAA")
    assert m.map_indices((3, 4, 5)) == [3, 4, 5]


def test_offset_mapping_is_applied():
    """A chain the two parsers read differently must be translated, not trusted."""
    # ESM-IF sees two extra residues before the region the manifest indexes.
    m = mapping({i: i + 2 for i in range(10)}, "XXAAANLTAA")
    assert m.map_indices((3, 4, 5)) == [5, 6, 7]
    m.check_triplet(m.map_indices((3, 4, 5)), "NLT")


def test_unmapped_index_raises_rather_than_guessing():
    m = mapping({0: 0, 1: 1, 2: 2}, "NLTAA")
    with pytest.raises(IncompleteBackboneError, match="no ESM-IF counterpart"):
        m.map_indices((0, 1, 7))


def test_incomplete_backbone_raises():
    finite = np.ones(10, dtype=bool)
    finite[4] = False                      # +1 residue is missing an atom
    m = mapping({i: i for i in range(10)}, "AAANLTAAAA", finite=finite)
    with pytest.raises(IncompleteBackboneError, match="backbone incomplete"):
        m.map_indices((3, 4, 5))


def test_triplet_mismatch_raises():
    """The guard that stops a parser disagreement being scored as biology."""
    m = mapping({i: i for i in range(10)}, "AAAKLTAAAA")     # K where N was expected
    with pytest.raises(SequonMismatchError, match="manifest records"):
        m.check_triplet(m.map_indices((3, 4, 5)), "NLT")


def test_triplet_check_passes_when_identities_agree():
    m = mapping({i: i for i in range(10)}, "AAANLTAAAA")
    m.check_triplet(m.map_indices((3, 4, 5)), "NLT")


# --- the distribution guard ----------------------------------------------

def test_check_scoreable_rejects_unnormalised_rows():
    alphabet = FakeAlphabet()
    probs = uniform_probs(6, alphabet, {})
    probs[3] = np.ones(len(alphabet))                # the failure mode by shape
    with pytest.raises(InvalidProbabilityVector, match="not 1"):
        check_scoreable(probs, (3, 4, 5))


def test_check_scoreable_rejects_negative_values():
    alphabet = FakeAlphabet()
    probs = uniform_probs(6, alphabet, {})
    # Sums to exactly 1, so it slips past the sum check and only the range
    # check can catch it — which is why both checks exist.
    probs[4, :] = 0.0
    probs[4, 0] = -0.5
    probs[4, 1] = 1.5
    with pytest.raises(InvalidProbabilityVector, match=r"outside \[0, 1\]"):
        check_scoreable(probs, (3, 4, 5))


def test_check_scoreable_accepts_real_distributions():
    alphabet = FakeAlphabet()
    check_scoreable(uniform_probs(6, alphabet, {}), (3, 4, 5))


# --- the score itself -----------------------------------------------------

def test_score_reads_the_residues_it_claims_to():
    """A wrong alphabet ordering is exactly how this went wrong for ProteinMPNN."""
    alphabet = FakeAlphabet()
    probs = uniform_probs(6, alphabet,
                          {3: {"N": 0.7}, 4: {"P": 0.25}, 5: {"S": 0.3, "T": 0.4}})
    result = sequon_score(probs, alphabet, 3, 4, 5)
    assert result["p_asn_at_n"] == pytest.approx(0.7)
    assert result["p_ser_at_plus2"] == pytest.approx(0.3)
    assert result["p_thr_at_plus2"] == pytest.approx(0.4)
    assert result["p_ser_or_thr_at_plus2"] == pytest.approx(0.7)
    assert result["p_pro_at_plus1"] == pytest.approx(0.25)


def test_score_is_the_mean_of_the_two_log_odds():
    alphabet = FakeAlphabet()
    probs = uniform_probs(6, alphabet, {3: {"N": 0.5}, 5: {"S": 0.25, "T": 0.25}})
    result = sequon_score(probs, alphabet, 3, 4, 5)
    # both terms are logit(0.5) = 0, so the mean is 0
    assert result["conditional_sequon_score"] == pytest.approx(0.0, abs=1e-9)


def test_deterministic_pass_reports_one_order_and_no_spread():
    alphabet = FakeAlphabet()
    result = sequon_score(uniform_probs(6, alphabet, {}), alphabet, 3, 4, 5)
    assert result["n_decoding_orders"] == N_ORDERS == 1
    assert result["conditional_sequon_score_sd"] == SCORE_SD == 0.0


def test_column_names_match_proteinmpnn_exactly():
    """Downstream is model-agnostic only while both scorers emit one schema."""
    alphabet = FakeAlphabet()
    esmif = sequon_score(uniform_probs(6, alphabet, {}), alphabet, 3, 4, 5)

    mpnn_probs = np.full((2, 6, 21), 1.0 / 21)
    mpnn = mpnn_sequon_score(mpnn_probs, 3, 4, 5)
    assert set(esmif) == set(mpnn)


def test_conditioning_is_named_for_what_it_is():
    assert CONDITIONING == "autoregressive_prefix"


# --- model-dependent -------------------------------------------------------

@requires_esm
def test_biotite_patch_is_idempotent_and_tolerant():
    from biotite.sequence import ProteinSequence

    from experimental_glycosylation_sites.esmif_scoring import patch_biotite

    patch_biotite()
    patch_biotite()                                   # must not double-wrap
    assert ProteinSequence.convert_letter_3to1("ALA") == "A"
    # the failure that aborted whole chains: an unknown residue name
    assert ProteinSequence.convert_letter_3to1("PCA") == "X"


@requires_esm
def test_inverse_folding_imports_after_the_patch():
    from experimental_glycosylation_sites.esmif_scoring import patch_biotite

    patch_biotite()
    import esm.inverse_folding  # noqa: F401
