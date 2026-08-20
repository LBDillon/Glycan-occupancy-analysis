"""ESM-IF1 adapter — the benchmark's second model.

Implements both SequonScorer and SequenceDesigner over
`esm_if1_gvp4_t16_142M_UR50`. All of the substance lives in `esmif_scoring`;
this is the registration surface and the per-chain cache.

Read `esmif_scoring`'s module docstring before comparing its numbers with
ProteinMPNN's. In short: ESM-IF is autoregressive, so its conditional is
prefix-only rather than ProteinMPNN's full bidirectional one, there is no
decoding-order distribution to average over, and the two models' raw score
magnitudes are therefore not comparable. The matched-pair contrast within each
model is.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..esmif_scoring import (CONDITIONING, DEFAULT_MODEL, N_ORDERS,
                             ChainUnreadableError, chain_mapping,
                             conditional_probabilities, decodable_positions,
                             design_sequences, load_model, sequon_score)


class ESMIFAdapter:
    """Implements both SequonScorer and SequenceDesigner."""

    name = "esm_if"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 seed: int = 0, max_batch: "int | None" = None):
        self.device = device
        self.max_batch = max_batch
        self.model_name = model_name
        self.seed = seed
        self._model = None
        self._alphabet = None
        # One chain at a time. The runners group their manifests by chain, so a
        # single slot removes every repeat forward pass; holding more would just
        # pin activations for chains nobody is going to ask about again.
        self._cached_key = None
        self._cached_mapping = None

    @property
    def model(self):
        if self._model is None:
            self._model, self._alphabet = load_model(self.device, self.model_name)
        return self._model

    @property
    def alphabet(self):
        if self._alphabet is None:
            self.model
        return self._alphabet

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write.

        `n_orders` is 1 and the conditioning is named for what it is: a single
        deterministic left-to-right pass, not an average over sampled decoding
        orders. The column exists so both models share a schema, and its value
        here says plainly that no averaging happened.
        """
        return {"model": self.model_name, "conditioning": CONDITIONING,
                "n_orders": N_ORDERS, "seed": self.seed}

    def _mapping(self, structure_path: Path, chain_id: str):
        key = (str(structure_path), str(chain_id))
        if key != self._cached_key:
            self._cached_mapping = chain_mapping(structure_path, chain_id)
            self._cached_key = key
        return self._cached_mapping

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Needs no model pass: it is the coordinate mask plus the index mapping."""
        return decodable_positions(structure_path, chain_id)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """One teacher-forced pass for the whole chain.

        `positions` is ignored: ESM-IF decodes the chain in a single pass, so
        restricting it would save nothing and would only make the second sequon
        on a chain cost as much as the first.
        """
        mapping = self._mapping(structure_path, chain_id)
        probabilities = conditional_probabilities(
            mapping, self.model, self.alphabet, device=self.device)
        return mapping, probabilities

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon, checking the residues really are the ones meant.

        The manifest's indices come from Biopython's parse; ESM-IF's come from
        biotite's. `map_indices` translates between them and `check_triplet`
        refuses anything whose residue identities do not reproduce the
        manifest's — the guard that stops a parser disagreement being scored as
        a biological observation.
        """
        mapping, probabilities = context
        mapped = mapping.map_indices(indices)
        if expected_triplet:
            mapping.check_triplet(mapped, expected_triplet)
        return sequon_score(probabilities, self.alphabet, *mapped)

    def score_site(self, structure_path: Path, chain_id: str, indices,
                   expected_triplet: "str | None" = None):
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices, expected_triplet)

    # --- SequenceDesigner ---------------------------------------------------
    def design(self, structure_path: Path, chain_id: str, n_designs: int,
               temperature: float, seed: int = 0) -> "list[str]":
        """Unconstrained designs, returned in the manifest's index space."""
        mapping = self._mapping(structure_path, chain_id)
        return design_sequences(mapping, self.model, self.alphabet,
                                n_designs=n_designs, temperature=temperature,
                                device=self.device, seed=seed,
                                max_batch=self.max_batch)
