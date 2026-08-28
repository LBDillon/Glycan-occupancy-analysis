"""ESM3 adapter — the model that can withhold its own structure track.

Scorer only, for ESMC's reason: generation from a masked language model would
condition on sequence rather than on a backbone, so "retention" would not mean
what it means for the inverse-folding models.

Two knobs, and the first is what makes this model worth having:

    structure_mode  struct_cond | seq_only    the structure track, on or off
    mask_mode       single | joint            what is hidden at the sequon

Crossed, they give a 2x2 inside one model with one tokeniser and one masking
scheme. Every other structure-versus-sequence comparison in this benchmark is
between models and so confounds the question with architecture and training.

Read `esm3_scoring`'s module docstring before quoting a number.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..esm3_scoring import (DEFAULT_MASK_MODE, DEFAULT_MODEL,
                            DEFAULT_STRUCTURE_MODE, N_ORDERS, chain_context,
                            check_triplet, conditional_probabilities,
                            conditioning, decodable_positions, load_model,
                            sequon_score)


class ESM3Adapter:
    """Implements SequonScorer. Not a SequenceDesigner."""

    name = "esm3"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 structure_mode: str = DEFAULT_STRUCTURE_MODE,
                 mask_mode: str = DEFAULT_MASK_MODE, seed: int = 0):
        from ..esm3_scoring import _check

        _check(structure_mode, mask_mode)
        self.device = device
        self.model_name = model_name
        self.structure_mode = structure_mode
        self.mask_mode = mask_mode
        self.seed = seed
        self._model = None
        self._tokenizer = None
        self._cached_key = None
        self._cached = None

    def _load(self):
        if self._model is None:
            self._model, self._tokenizer = load_model(self.device, self.model_name)
        return self._model, self._tokenizer

    @property
    def model(self):
        return self._load()[0]

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write.

        The conditioning string names the structure track, so a
        structure-conditioned run and a sequence-only run of the SAME model can
        never be pooled -- which is the one confusion this model exists to avoid.
        """
        return {"model": f"{self.model_name}/{self.structure_mode}",
                "conditioning": conditioning(self.structure_mode, self.mask_mode),
                "n_orders": N_ORDERS, "seed": self.seed}

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Needs no model pass: it is whether ESM3's parse matches the manifest's."""
        return decodable_positions(structure_path, chain_id)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """Encode the chain once; the masked passes happen per sequon.

        Each scored position needs its own mask, so unlike ESM-IF there is no
        single pass that serves the whole chain -- but only the three sequon
        positions are ever masked, not all L of them as the source notebook does.
        """
        model, _ = self._load()
        key = (str(structure_path), str(chain_id))
        if key != self._cached_key:
            self._cached = chain_context(structure_path, chain_id, model,
                                         device=self.device)
            self._cached_key = key
        return self._cached

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon, checking the residues really are the ones meant."""
        model, tokenizer = self._load()
        sequence, _, _ = context
        if expected_triplet:
            check_triplet(sequence, indices, expected_triplet)
        probabilities = conditional_probabilities(
            context, model, tokenizer, indices,
            structure_mode=self.structure_mode, mask_mode=self.mask_mode,
            device=self.device)
        return sequon_score(probabilities, tokenizer, *indices)

    def score_site(self, structure_path: Path, chain_id: str, indices,
                   expected_triplet: "str | None" = None):
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices, expected_triplet)
