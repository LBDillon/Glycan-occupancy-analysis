"""ProteinMPNN adapter — the reference implementation of the interface.

Thin wrapper over `mpnn_scoring` and `retention`, which hold the actual logic.
It exists so that a second model has something concrete to copy, and so the
pipeline never imports a specific model directly.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..mpnn_scoring import (ALPHABET, DEFAULT_DECODING_ORDERS, DEFAULT_MODEL,
                            chain_mapping, to_manifest_space,
                            conditional_probabilities, decodable_positions,
                            load_model, sequon_score)
from ..retention import design_sequences


class ProteinMPNNAdapter:
    """Implements both SequonScorer and SequenceDesigner."""

    name = "proteinmpnn"

    def __init__(self, proteinmpnn_dir: Path = Path("../../ProteinMPNN"),
                 checkpoint: str = DEFAULT_MODEL, device: str = "cpu",
                 n_decoding_orders: int = DEFAULT_DECODING_ORDERS, seed: int = 0,
                 max_batch: "int | None" = None, mask_mode: str = "single"):
        self.dir = Path(proteinmpnn_dir)
        self.checkpoint = checkpoint
        self.device = device
        self.n_decoding_orders = n_decoding_orders
        self.seed = seed
        self.max_batch = max_batch
        if mask_mode not in ("single", "joint"):
            raise ValueError(f"mask_mode must be 'single' or 'joint', got {mask_mode!r}")
        self.mask_mode = mask_mode
        self._model = None

    @property
    def model(self):
        if self._model is None:
            self._model = load_model(self.dir, self.checkpoint, self.device)
        return self._model

    def describe(self) -> dict:
        """Provenance the runners stamp onto every row they write."""
        conditioning = ("conditional" if self.mask_mode == "single"
                        else "conditional_sequon_masked")
        return {"model": self.checkpoint, "conditioning": conditioning,
                "n_orders": self.n_decoding_orders, "seed": self.seed}

    # --- SequonScorer -------------------------------------------------------
    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Backbone completeness; needs no model pass, so it can precede matching."""
        return decodable_positions(structure_path, chain_id, self.dir)

    def prepare_chain(self, structure_path: Path, chain_id: str, positions=None):
        """One decoder sweep for every position wanted on this chain.

        ProteinMPNN's `conditional_probs` runs a full decoder pass per position,
        so a 300-residue chain costs 300 passes when three residues are wanted.
        Passing `positions` restricts it to those, and the per-position result is
        unchanged by which others were requested.
        """
        # The manifest's indices count observed residues; ProteinMPNN's parse
        # walks the author numbering and inserts a placeholder per absent
        # number, so the two coincide only for gapless chains. Translate here,
        # once, rather than trusting them to agree.
        mapping, model_sequence = chain_mapping(structure_path, chain_id)
        if positions is not None:
            positions = [mapping[int(i)] for i in positions if int(i) in mapping]
        if self.mask_mode == "joint":
            # Nothing can be shared across sites: each residue is scored with a
            # different set hidden. Carry the identifiers so score_from can work.
            return (None, None, structure_path, chain_id, mapping, model_sequence)
        probabilities, computed = conditional_probabilities(
            structure_path, chain_id, self.model, device=self.device,
            n_decoding_orders=self.n_decoding_orders, seed=self.seed,
            positions=None if positions is None else list(positions))
        return (probabilities, computed, structure_path, chain_id, mapping, model_sequence)

    @staticmethod
    def _mapped(mapping, model_sequence, indices, expected_triplet):
        """Manifest indices translated into ProteinMPNN's, and checked.

        The same guard ESM-IF applies: refuse anything whose residue identities
        do not reproduce the manifest's, so a parser disagreement cannot be
        scored as a biological observation.
        """
        mapped = []
        for index in indices:
            target = mapping.get(int(index))
            if target is None:
                raise KeyError(f"model_index {int(index)} has no ProteinMPNN counterpart")
            mapped.append(int(target))
        if expected_triplet:
            observed = "".join(model_sequence[i] if i < len(model_sequence) else "?"
                               for i in mapped)
            if observed != expected_triplet:
                raise ValueError(
                    f"residues at the mapped indices are {observed!r}, "
                    f"not the manifest's {expected_triplet!r}")
        return tuple(mapped)

    def score_from(self, context, indices, expected_triplet=None):
        """Score one sequon from a prepared chain. Raises if any row is unfilled.

        Under `joint` the prepared chain is not used: each sequon residue needs
        its own pass, with the other two pushed later in the decoding order so it
        cannot condition on them. That is three passes per site rather than one
        shared sweep, which is the price of the sensitivity.
        """
        mapping, model_sequence = context[4], context[5]
        indices = self._mapped(mapping, model_sequence, indices, expected_triplet)

        if self.mask_mode == "single":
            probabilities, computed = context[0], context[1]
            return sequon_score(probabilities, *indices, computed=computed)

        structure_path, chain_id = context[2], context[3]
        rows, computed = None, None
        for position in indices:
            hide = [i for i in indices if i != position]
            probs, comp = conditional_probabilities(
                structure_path, chain_id, self.model, device=self.device,
                n_decoding_orders=self.n_decoding_orders, seed=self.seed,
                positions=[position], hide_positions=hide)
            if rows is None:
                rows, computed = np.zeros_like(probs), np.zeros_like(comp)
            rows[:, position, :] = probs[:, position, :]
            computed[position] = comp[position]
        return sequon_score(rows, *indices, computed=computed)

    def score_site(self, structure_path: Path, chain_id: str, indices):
        """Score one sequon. Raises if the model did not evaluate all three."""
        context = self.prepare_chain(structure_path, chain_id, list(indices))
        return self.score_from(context, indices)

    # --- SequenceDesigner ---------------------------------------------------
    def design(self, structure_path: Path, chain_id: str, n_designs: int,
               temperature: float, seed: int = 0) -> list[str]:
        """Unconstrained designs, returned in the manifest's index space.

        ProteinMPNN decodes in its own indexing, which inserts a placeholder for
        every absent author residue number. Returning that unprojected would
        have `classify_retention` read the manifest's indices against a sequence
        they do not address -- so the projection happens here, once, rather than
        being left to each caller to remember.
        """
        mapping, _ = chain_mapping(structure_path, chain_id)
        designs = design_sequences(structure_path, chain_id, self.model,
                                   n_designs=n_designs, temperature=temperature,
                                   device=self.device, seed=seed,
                                   max_batch=self.max_batch)
        if not mapping:
            raise KeyError(
                f"chain {chain_id} of {Path(structure_path).name} cannot be "
                "mapped onto ProteinMPNN's parse; its designs are not readable "
                "at manifest indices")
        return [to_manifest_space(d, mapping) for d in designs]
