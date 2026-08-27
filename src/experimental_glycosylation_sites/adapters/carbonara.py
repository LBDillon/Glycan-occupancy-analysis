"""CARBonAra adapter — a scorer only.

Implements `SequonScorer` and deliberately not `SequenceDesigner`. Sequence
generation upstream goes through `imprint_sampling`, which samples stochastically
from raw confidences; neither belongs in this benchmark, and a `design` method
here would let stage 08 produce retention numbers that look like ProteinMPNN's
without meaning the same thing.

All of the substance lives in `carbonara_scoring`. Read its module docstring
before putting these numbers beside another model's: the conditional is
`P(residue | backbone, all other native residues)` from one deterministic pass,
the probability vectors are twenty entries in CARBonAra's abundance-sorted order
rather than an alphabetical one, and the input is one protein chain with every
heteroatom removed — not the molecular context CARBonAra is known for.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..carbonara_scoring import (CONDITIONING, DEFAULT_MODEL, N_ORDERS,
                                 chain_mapping, conditional_probabilities,
                                 decodable_positions, load_model, sequon_score)


class CARBonAraAdapter:
    """Implements SequonScorer. Not a SequenceDesigner."""

    name = "carbonara"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 carbonara_dir: "Path | str | None" = None, seed: int = 0):
        self.device = device
        self.model_name = model_name
        # None means "discover it when first needed", so constructing the
        # adapter — and importing the package — never requires the checkout.
        self.dir = Path(carbonara_dir) if carbonara_dir is not None else None
        self.seed = seed
        self._model = None
        # One chain at a time. The runners group their manifests by chain, so a
        # single slot removes every repeat parse; holding more would pin
        # structures for chains nobody will ask about again.
        self._cached_key = None
        self._cached_mapping = None

    @property
    def model(self):
        if self._model is None:
            self._model = load_model(self.dir, self.model_name, self.device)
        return self._model

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write.

        `n_orders` is 1 and the conditioning is named for what it is: a single
        deterministic pass per position, not an average over sampled decoding
        orders. `seed` is carried so the models share a schema, but nothing here
        is sampled and no value of it changes a score.
        """
        return {"model": self.model_name, "conditioning": CONDITIONING,
                "n_orders": N_ORDERS, "seed": self.seed}

    def _mapping(self, structure_path: Path, chain_id: str):
        key = (str(structure_path), str(chain_id))
        if key != self._cached_key:
            self._cached_mapping = chain_mapping(structure_path, chain_id,
                                                 carbonara_dir=self.dir)
            self._cached_key = key
        return self._cached_mapping

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Needs no model pass: it is the index mapping plus backbone completeness."""
        return decodable_positions(structure_path, chain_id,
                                   carbonara_dir=self.dir)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """One forward pass per requested position.

        CARBonAra pays per position — each needs its own imprint, with that
        residue's identity withheld — so `positions` genuinely restricts the
        work and passing it matters. Passing None means every residue on the
        chain, which is one pass per residue and should be reserved for a chain
        you actually want in full.

        Positions with no verified counterpart, or with an incomplete backbone,
        are dropped here rather than computed. `score_from` raises for them
        individually, so one bad residue costs its own site and not the chain.
        """
        mapping = self._mapping(structure_path, chain_id)
        if positions is None:
            wanted = [index for index in sorted(mapping.to_model.values())
                      if bool(mapping.backbone_ok[index])]
        else:
            wanted = sorted({
                mapping.to_model[int(i)] for i in positions
                if int(i) in mapping.to_model
                and bool(mapping.backbone_ok[mapping.to_model[int(i)]])})
        probabilities = conditional_probabilities(
            mapping, self.model, wanted, device=self.device,
            carbonara_dir=self.dir)
        return mapping, probabilities

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon, checking the residues really are the ones meant.

        The manifest's indices come from Biopython's parse; CARBonAra's come
        from gemmi's, after a renumbering that discards insertion codes.
        `map_indices` translates and `check_triplet` refuses anything whose
        residue identities do not reproduce the manifest's — the guard that
        stops a parser disagreement being scored as a biological observation.
        """
        mapping, probabilities = context
        mapped = mapping.map_indices(indices)
        if expected_triplet:
            mapping.check_triplet(mapped, expected_triplet)
        return sequon_score(probabilities, *mapped)

    def score_site(self, structure_path: Path, chain_id: str, indices,
                   expected_triplet: "str | None" = None):
        """Score one sequon. Raises if the model did not evaluate all three."""
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices, expected_triplet)
