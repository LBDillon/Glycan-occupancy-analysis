"""ProGen2 adapter — the causal sequence-only model.

Scorer only. ProGen2 generates sequences, but generation here is unconditioned
by any backbone, so "redesigning this chain" is not a thing it can be asked to
do — a design would be a fresh protein rather than a redesign of this one, and
retention would measure nothing about the site. `isinstance(adapter,
SequenceDesigner)` is False, as it is for ESMC and for the same reason.

Read `progen2_scoring`'s module docstring before comparing its numbers. In
short: the conditional is prefix-only, so ESM-IF is the like-for-like comparison
and ESMC is not — the two differ in conditioning as well as in input, and
`conditioning` is recorded as `autoregressive_prefix` so they can never be
pooled.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..progen2_scoring import (CONDITIONING, DEFAULT_MODEL, N_ORDERS,
                               chain_sequence, check_triplet,
                               conditional_probabilities, decodable_positions,
                               load_model, sequon_score, verify_tokenisation)


class ProGen2Adapter:
    """Implements SequonScorer. Not a SequenceDesigner."""

    name = "progen2"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 seed: int = 0):
        self.device = device
        self.model_name = model_name
        self.seed = seed
        self._model = None
        self._tokenizer = None
        self._bos = None
        self._aa_index = None
        self._cached_key = None
        self._cached = None

    def _load(self):
        if self._model is None:
            self._model, self._tokenizer, self._bos = load_model(
                self.model_name, self.device)
            self._aa_index = verify_tokenisation(self._tokenizer)
        return self._model, self._tokenizer, self._bos

    @property
    def model(self):
        return self._load()[0]

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write.

        `conditioning` is ESM-IF's string, not ESMC's: this model is causal and
        prefix-only, so pooling it with a masked language model would average
        two different estimands.
        """
        return {"model": self.model_name, "conditioning": CONDITIONING,
                "n_orders": N_ORDERS, "seed": self.seed}

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """All True: no backbone requirement. Needs no model pass."""
        return decodable_positions(structure_path, chain_id)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """One teacher-forced pass for the whole chain.

        `positions` is ignored: the pass yields every position at once, so
        restricting it would save nothing and make the second sequon on a chain
        cost as much as the first.
        """
        key = (str(structure_path), str(chain_id))
        if key != self._cached_key:
            model, tokenizer, bos = self._load()
            sequence = chain_sequence(structure_path, chain_id)
            probabilities = conditional_probabilities(
                sequence, model, tokenizer, bos, device=self.device)
            self._cached = (sequence, probabilities)
            self._cached_key = key
        return self._cached

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon, checking the residues really are the ones meant.

        No cross-parser mapping is needed here — the sequence scored is the one
        `model_index` indexes — but the triplet is still checked, because "no
        mapping needed" is itself an assumption about which sequence was loaded.
        """
        sequence, probabilities = context
        if expected_triplet:
            check_triplet(sequence, indices, expected_triplet)
        return sequon_score(probabilities, self._aa_index, *indices)

    def score_site(self, structure_path: Path, chain_id: str, indices,
                   expected_triplet: "str | None" = None):
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices, expected_triplet)
