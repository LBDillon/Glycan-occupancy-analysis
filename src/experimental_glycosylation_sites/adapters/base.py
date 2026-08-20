"""The interface a model must implement to enter this benchmark.

Everything downstream of scoring — matching, contrasts, the cluster bootstrap,
significance testing, the figures — works on model-agnostic tables keyed by
(accession, position). So adding a second model means writing one adapter and
nothing else.

An adapter supplies two capabilities, and may supply only one:

    SequonScorer  what probability the model holds at the three sequon residues,
                  given the native backbone and the rest of the native sequence
    SequenceDesigner  what the model actually writes when asked to redesign the
                  chain, from which sequon retention is read

The two answer different questions and are analysed separately. A model that
cannot generate sequences can still be scored; a model that offers no per-residue
probabilities can still be measured on retention.

## Adding a model

1. Create `adapters/<model>.py` with a class implementing one or both protocols.
2. Register it in `adapters/__init__.py`.
3. Run the pipeline with `--model <name>`. Scores land in
   `results/scores/scores_<model>.csv`, designs in `results/retention_<model>.csv`.

Two invariants any adapter must honour, both learned the hard way:

**Never return a score for a residue the model did not actually evaluate.**
ProteinMPNN silently returns a zero row for residues with incomplete backbones,
which exponentiates to twenty-one ones and scores +13.8. Report such sites as
unscoreable rather than scoring them. `decodable_positions`-style prescreening,
settled before matching, is the pattern to copy.

**Probabilities must be probabilities.** Rows are checked to sum to one; an
adapter that returns logits or unnormalised values will be rejected downstream
rather than silently mis-scored.
"""
from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class SequonScorer(Protocol):
    """Per-site conditional probabilities at the three sequon residues."""

    name: str

    def decodable_positions(self, structure_path: Path, chain_id: str) -> np.ndarray:
        """Boolean array over the chain: which residues the model can evaluate.

        Must be answerable WITHOUT running the model, so scoreability can be
        settled before matching rather than after. If a model has no such
        restriction, return all True.
        """

    def score_site(self, structure_path: Path, chain_id: str,
                   indices: "tuple[int, int, int]") -> dict:
        """Score one sequon, given the 0-based model indices of its three residues.

        Returns a dict carrying at least `conditional_sequon_score`. Must RAISE
        rather than return a value if the model did not genuinely evaluate all
        three positions — the caller records that as an exclusion. Scoring is
        per site rather than per manifest so that a chain's probabilities can be
        computed once and read for every sequon on it.
        """

    # --- the batched path -------------------------------------------------
    # `score_site` is the honest unit of the interface, but scoring a manifest
    # one site at a time is wasteful in opposite ways for different models, and
    # neither waste is the caller's business:
    #
    #   ProteinMPNN runs one decoder pass PER POSITION, so it wants to be told
    #   every position on a chain at once and compute only those.
    #   ESM-IF runs ONE teacher-forced pass for the whole chain, so the second
    #   and subsequent sequons on a chain should cost nothing at all.
    #
    # Both are served by splitting the work: `prepare_chain` does whatever the
    # model does once per chain and returns an opaque context, `score_from`
    # reads a sequon out of it. `score_site` is then just the two composed, and
    # a runner that groups its manifest by chain uses the pair directly.

    def prepare_chain(self, structure_path: Path, chain_id: str,
                      positions: "list[int] | None" = None):
        """Do the once-per-chain work; return an opaque context for `score_from`.

        `positions` is every model index the caller will ask about on this
        chain, so a model that pays per position can restrict itself to them. A
        model that computes the whole chain regardless may ignore it.
        """

    def score_from(self, context, indices: "tuple[int, int, int]",
                   expected_triplet: "str | None" = None) -> dict:
        """Score one sequon from a prepared chain context.

        Same contract as `score_site`: raise rather than return a value the
        model did not genuinely produce. `expected_triplet` is the residue
        identity the manifest recorded; an adapter whose parser may disagree
        with the manifest's MUST check it and raise on a mismatch, because a
        silent off-by-some scores the wrong residue and still looks plausible.
        """


@runtime_checkable
class SequenceDesigner(Protocol):
    """Unconstrained sequence generation, for retention measurement."""

    name: str

    def design(self, structure_path: Path, chain_id: str, n_designs: int,
               temperature: float, seed: int) -> list[str]:
        """`n_designs` full-length sequences for the chain, as one-letter strings.

        Nothing may be fixed or biased: the sequon must be free to disappear, or
        retention measures the constraint rather than the model.
        """
