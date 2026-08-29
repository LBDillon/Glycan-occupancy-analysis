"""ESM3 adapter — the model that can withhold its own structure track.

Scorer and designer. Design is available only under `struct_cond`: retention
asks what a model writes when it redesigns a BACKBONE, and with the structure
track withheld there is no backbone to redesign against -- generation would be
unconditional, which is why ProGen2 has no retention row either. So the 2x2
belongs to scoring; design is one arm of it.

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
                            DESIGN_MAX_STEPS, DESIGN_MIN_STEPS,
                            DESIGN_STEP_DIVISOR, conditioning,
                            decodable_positions, design_sequences, design_steps,
                            load_model, sequon_score)


class ESM3Adapter:
    """Implements SequonScorer, and SequenceDesigner under struct_cond."""

    name = "esm3"

    def __init__(self, device: str = "cpu", model_name: str = DEFAULT_MODEL,
                 structure_mode: str = DEFAULT_STRUCTURE_MODE,
                 mask_mode: str = DEFAULT_MASK_MODE, seed: int = 0,
                 num_steps: "int | None" = None,
                 max_batch: "int | None" = None, use_batch: bool = True):
        from ..esm3_scoring import _check

        _check(structure_mode, mask_mode)
        self.device = device
        self.model_name = model_name
        self.structure_mode = structure_mode
        self.mask_mode = mask_mode
        self.seed = seed
        self.num_steps = num_steps
        # Designs are generated in batches, so the batch-times-length blow-up
        # that the slot budget exists to prevent applies here exactly as it does
        # to ProteinMPNN's decoding, and max_batch caps it as it does there.
        self.max_batch = max_batch
        # An escape hatch to the one-at-a-time path, which is slower by roughly
        # the batch size but needs only `generate`. Kept because the batched
        # path depends on `batch_generate` existing on the loaded model.
        self.use_batch = use_batch
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
        described = {"model": f"{self.model_name}/{self.structure_mode}",
                     "conditioning": conditioning(self.structure_mode,
                                                  self.mask_mode),
                     "n_orders": N_ORDERS, "seed": self.seed}
        if self.num_steps:
            described["design_num_steps"] = self.num_steps
        return described

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

    def describe_generation(self) -> dict:
        """Provenance for the design rows, which is not the scoring provenance.

        Stage 08 labels designs `autoregressive_sampling` unless an adapter
        says otherwise, and that label would be wrong here: ESM3 is a masked
        diffusion model that fills a fully masked track over several passes,
        not a left-to-right decoder. Retention measured under this procedure
        and under an autoregressive one are not the same quantity, and the
        column is what keeps them from being pooled.
        """
        return {"generation": "masked_diffusion_unmasking",
                "num_steps": self.num_steps or "scaled_by_length",
                "step_divisor": DESIGN_STEP_DIVISOR,
                "step_bounds": [DESIGN_MIN_STEPS, DESIGN_MAX_STEPS],
                "batched": bool(self.use_batch),
                "native_procedure": True}

    # --- SequenceDesigner ---------------------------------------------------
    def design(self, structure_path: Path, chain_id: str, n_designs: int,
               temperature: float, seed: "int | None" = None) -> "list[str]":
        """Unconstrained designs, in the manifest's index space.

        Refuses under `seq_only` rather than quietly generating without a
        backbone. A sequence produced with the structure track withheld is not
        a redesign of this chain, and scoring retention on it would compare an
        unconditional sample against the inverse-folding models.
        """
        if self.structure_mode != "struct_cond":
            raise ValueError(
                f"design needs the structure track, but structure_mode is "
                f"{self.structure_mode!r}. Retention measures what the model "
                "writes for a given backbone; with no backbone there is "
                "nothing to redesign.")
        model, _ = self._load()
        return design_sequences(
            structure_path, chain_id, model, n_designs=n_designs,
            temperature=temperature,
            seed=self.seed if seed is None else seed,
            device=self.device, num_steps=self.num_steps,
            max_batch=self.max_batch, use_batch=self.use_batch)

    def design_steps_for(self, length: int) -> int:
        """Unmasking passes this adapter would use for a chain of this length."""
        return design_steps(length, self.num_steps)
