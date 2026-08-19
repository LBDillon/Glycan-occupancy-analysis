"""ProteinMPNN adapter — the reference implementation of the interface.

Thin wrapper over `mpnn_scoring` and `retention`, which hold the actual logic.
It exists so that a second model has something concrete to copy, and so the
pipeline never imports a specific model directly.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..mpnn_scoring import (DEFAULT_DECODING_ORDERS, DEFAULT_MODEL,
                            conditional_probabilities, decodable_positions,
                            load_model, sequon_score)
from ..retention import design_sequences


class ProteinMPNNAdapter:
    """Implements both SequonScorer and SequenceDesigner."""

    name = "proteinmpnn"

    def __init__(self, proteinmpnn_dir: Path = Path("../../ProteinMPNN"),
                 checkpoint: str = DEFAULT_MODEL, device: str = "cpu",
                 n_decoding_orders: int = DEFAULT_DECODING_ORDERS, seed: int = 0):
        self.dir = Path(proteinmpnn_dir)
        self.checkpoint = checkpoint
        self.device = device
        self.n_decoding_orders = n_decoding_orders
        self.seed = seed
        self._model = None

    @property
    def model(self):
        if self._model is None:
            self._model = load_model(self.dir, self.checkpoint, self.device)
        return self._model

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Backbone completeness; needs no model pass, so it can precede matching."""
        return decodable_positions(structure_path, chain_id, self.dir)

    def score_site(self, structure_path: Path, chain_id: str, indices):
        """Score one sequon. Raises if the model did not evaluate all three."""
        probabilities, computed = conditional_probabilities(
            structure_path, chain_id, self.model, device=self.device,
            n_decoding_orders=self.n_decoding_orders, seed=self.seed,
            positions=list(indices))
        return sequon_score(probabilities, *indices, computed=computed)

    # --- SequenceDesigner ---------------------------------------------------
    def design(self, structure_path: Path, chain_id: str, n_designs: int,
               temperature: float, seed: int = 0) -> list[str]:
        return design_sequences(structure_path, chain_id, self.model,
                                n_designs=n_designs, temperature=temperature,
                                device=self.device, seed=seed)
