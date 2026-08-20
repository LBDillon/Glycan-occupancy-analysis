"""ESMC adapter — the sequence-only arm of the benchmark.

Implements `SequonScorer` and deliberately not `SequenceDesigner`: ESMC is a
masked language model, so sampling from it would condition on sequence rather
than on a backbone, and "retention" would not mean what it means for the
inverse-folding models.

Read `esmc_scoring`'s module docstring before comparing its numbers with the
others. In short: it sees no structure, it is scored on the same chain sequence
the manifest indexes, and `mask_mode` chooses between conditioning on every other
native residue (`single`, primary) and hiding the whole sequon (`joint`,
sensitivity).

Requires EvolutionaryScale's `esm`, which collides on the import name `esm` with
`fair-esm` (ESM-IF). Only one can be installed at a time.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..esmc_scoring import (DEFAULT_MASK_MODE, DEFAULT_MODEL, MASK_MODES,
                            N_ORDERS, SequonMismatchError, chain_sequence,
                            conditioning, decodable_positions, load_model,
                            masked_distributions, sequon_score)


class ESMCAdapter:
    """Implements SequonScorer only."""

    name = "esmc"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 mask_mode: str = DEFAULT_MASK_MODE, seed: int = 0):
        if mask_mode not in MASK_MODES:
            raise ValueError(f"mask_mode must be one of {MASK_MODES}, got {mask_mode!r}")
        self.device = device
        self.model_name = model_name
        self.mask_mode = mask_mode
        self.seed = seed
        self._model = None
        self._tokenizer = None
        self._cached_key = None
        self._cached_sequence = None

    @property
    def model(self):
        if self._model is None:
            self._model, self._tokenizer = load_model(self.device, self.model_name)
        return self._model

    @property
    def tokenizer(self):
        if self._tokenizer is None:
            self.model
        return self._tokenizer

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write."""
        return {"model": self.model_name, "conditioning": conditioning(self.mask_mode),
                "n_orders": N_ORDERS, "seed": self.seed}

    def _sequence(self, structure_path: Path, chain_id: str) -> str:
        key = (str(structure_path), str(chain_id))
        if key != self._cached_key:
            self._cached_sequence = chain_sequence(structure_path, chain_id)
            self._cached_key = key
        return self._cached_sequence

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Every residue: a sequence model has no backbone requirement."""
        return decodable_positions(structure_path, chain_id)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """Mask and forward once per group, batched across the chain.

        Under `single` the groups are the individual positions, which are
        independent, so every masked variant the chain needs goes through in one
        batched sweep. Under `joint` a sequon must be masked as a unit and the
        caller's flat `positions` no longer identifies which triplet is which, so
        the work is deferred to `score_from`.
        """
        sequence = self._sequence(structure_path, chain_id)
        if self.mask_mode == "joint" or positions is None:
            return sequence, None
        groups = [(int(p),) for p in sorted({int(p) for p in positions})]
        return sequence, masked_distributions(sequence, groups, self.model,
                                              self.tokenizer, self.device)

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon, checking the residues are the ones the manifest meant."""
        sequence, cache = context
        indices = tuple(int(i) for i in indices)

        if expected_triplet:
            observed = "".join(sequence[i] for i in indices
                               if 0 <= i < len(sequence))
            if observed != expected_triplet:
                raise SequonMismatchError(
                    f"ESMC reads {observed!r} where the manifest records "
                    f"{expected_triplet!r}")

        if cache is None:                      # joint: mask the sequon as a unit
            vectors = masked_distributions(sequence, [indices], self.model,
                                           self.tokenizer, self.device)[indices]
        else:
            vectors = {i: cache[(i,)][i] for i in indices}

        return sequon_score(vectors, self.tokenizer, *indices)

    def score_site(self, structure_path: Path, chain_id: str, indices,
                   expected_triplet: "str | None" = None):
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices, expected_triplet)
