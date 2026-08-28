"""ESM3 scoring — one model with the structure track switched on and off.

Adapted from `score_proteins_esm3_colab.ipynb`. The notebook's masking scheme is
kept; its statistic is not. It reports a whole-chain mean pseudo-log-likelihood,
which is the protein-level score this module already argues cannot resolve a
three-residue question. Here the same masked pass is read at the three sequon
positions instead.

## Why this model is worth having

Every other "does structure matter?" comparison in this benchmark is between
models — ESMC against ProteinMPNN, ESM-IF, CARBonAra — and so confounds the
structure question with architecture, training data and tokenisation. ESM3
carries a structure track that can simply be withheld:

    struct_cond   VQ-VAE structure tokens from the backbone, intact
    seq_only      the same model, same tokeniser, same masking, no structure

The difference between those two is the structure contribution with nothing else
varying. Combined with `mask_mode`, it gives a 2x2 inside one model — structure
on/off crossed with motif visible/hidden — which no cross-model comparison here
can match.

## What the conditional is

Masked, bidirectional, exactly as for ESMC: one position is replaced by the mask
token and the distribution read there. `conditional_sequon_score_sd` is
structurally zero and `n_decoding_orders` is one.

`mask_mode` selects what is hidden:

    single   the scored position only, so the other two sequon residues are
             visible and P(S/T) can be read off the upstream asparagine
    joint    all three at once, which removes that shortcut at the cost of
             conditioning on strictly less context

The notebook masks every position in turn to build a whole-chain PLL. Only three
positions per site are needed here, so a site costs one or three forward passes
rather than L of them.

## Which sequence

The chain as `_parse_chains` reads it, so `model_index` indexes it directly --
the same choice made for ESMC. ESM3's own parser (`ProteinChain.from_pdb`) is
checked against it per chain rather than assumed to agree: it did agree exactly
on the three chains used for the smoke test, but ESM-IF's parser disagreed on
about 5% of sites and that was also invisible until checked.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .mpnn_scoring import EPSILON, PROBABILITY_SUM_TOLERANCE, logit  # noqa: F401

DEFAULT_MODEL = "esm3-open"

STRUCTURE_MODES = ("struct_cond", "seq_only")
MASK_MODES = ("single", "joint")
DEFAULT_STRUCTURE_MODE = "struct_cond"
DEFAULT_MASK_MODE = "single"

N_ORDERS = 1
SCORE_SD = 0.0

# The tokeniser prepends <cls>, so sequence position i is token index i + 1.
# Asserted against the tokeniser at load rather than trusted.
TOKEN_OFFSET = 1


class ChainUnreadableError(ValueError):
    """ESM3 could not be given a chain whose indices we can justify."""


class SequonMismatchError(ValueError):
    """ESM3 read different residues than the manifest recorded."""


class InvalidProbabilityVector(ValueError):
    """A row that is not a probability distribution."""


def conditioning(structure_mode: str, mask_mode: str) -> str:
    """The provenance string. Names the structure track, because that is the
    variable this model exists to manipulate."""
    track = ("structure_conditioned" if structure_mode == "struct_cond"
             else "sequence_only")
    return f"masked_{track}_{mask_mode}"


def _check(structure_mode: str, mask_mode: str) -> None:
    if structure_mode not in STRUCTURE_MODES:
        raise ValueError(f"structure_mode must be one of {STRUCTURE_MODES}, "
                         f"got {structure_mode!r}")
    if mask_mode not in MASK_MODES:
        raise ValueError(f"mask_mode must be one of {MASK_MODES}, got {mask_mode!r}")


def load_model(device: str = "cpu", model_name: str = DEFAULT_MODEL):
    """Load ESM3 in eval mode and verify the token offset. Returns `(model, tokenizer)`.

    Needs EvolutionaryScale's `esm`, which collides on the import name with
    `fair-esm` (ESM-IF); the two cannot share an environment. The checkpoint is
    gated on HuggingFace and must be accepted there before first use.
    """
    import torch
    from esm.models.esm3 import ESM3

    model = ESM3.from_pretrained(model_name).to(device).eval()
    tokenizer = model.tokenizers.sequence
    _assert_token_offset(tokenizer)
    return model, tokenizer


def _assert_token_offset(tokenizer) -> None:
    """Round-trip a probe through the tokeniser and back.

    The same class of check the ProteinMPNN alphabet defect went undetected for:
    an offset here reads a neighbouring residue and returns a number that looks
    entirely reasonable.
    """
    probe = "MNKTA"
    ids = tokenizer.encode(probe)
    for offset, residue in enumerate(probe):
        token = tokenizer.convert_ids_to_tokens([ids[offset + TOKEN_OFFSET]])[0]
        if token != residue:
            raise RuntimeError(
                f"token offset {TOKEN_OFFSET} is wrong: position {offset} of "
                f"{probe!r} decodes to {token!r}, not {residue!r}")


def chain_context(structure_path, chain_id: str, model,
                  pdb_id: "str | None" = None, device: str = "cpu"):
    """Encode one chain, checking ESM3's parse against the manifest's.

    Returns `(sequence, sequence_tokens, structure_tokens)`. A chain whose
    sequence ESM3 reads differently is refused rather than scored at indices
    that address a different residue.
    """
    from esm.sdk.api import ESMProtein
    from esm.utils.structure.protein_chain import ProteinChain

    from .structures import _parse_chains

    path = Path(structure_path)
    identifier = str(pdb_id or path.stem)

    chains = _parse_chains(path, identifier)
    native = next((c for c in chains if c.chain_id == str(chain_id)), None)
    if native is None or not native.sequence:
        raise ChainUnreadableError(f"chain {chain_id!r} absent from {path.name}")

    try:
        chain = ProteinChain.from_pdb(str(path), chain_id=str(chain_id))
    except Exception as exc:
        raise ChainUnreadableError(f"{type(exc).__name__}: {str(exc)[:120]}") from exc

    if chain.sequence != native.sequence:
        raise ChainUnreadableError(
            f"ESM3 reads {len(chain.sequence)} residues where the manifest's "
            f"parse lists {len(native.sequence)}; the manifest's indices would "
            "not address the same residues")

    encoded = model.encode(ESMProtein.from_protein_chain(chain))
    sequence_tokens = encoded.sequence.to(device)
    structure_tokens = (encoded.structure.to(device)
                        if encoded.structure is not None else None)
    return native.sequence, sequence_tokens, structure_tokens


def conditional_probabilities(context, model, tokenizer, indices,
                              structure_mode: str = DEFAULT_STRUCTURE_MODE,
                              mask_mode: str = DEFAULT_MASK_MODE,
                              device: str = "cpu") -> "dict[int, np.ndarray]":
    """P(residue at each requested index), with the chosen tracks and masking.

    Under `single` each position is masked alone, so the passes are independent
    and batched together. Under `joint` all three are masked at once and one pass
    serves all three, which is both cheaper and a different estimand.
    """
    import torch

    _check(structure_mode, mask_mode)
    _, sequence_tokens, structure_tokens = context
    mask_id = tokenizer.mask_token_id
    wanted = sorted({int(i) for i in indices})

    if structure_mode == "seq_only":
        structure_tokens = None

    if mask_mode == "joint":
        batch = sequence_tokens.unsqueeze(0).clone()
        for index in wanted:
            batch[0, index + TOKEN_OFFSET] = mask_id
        rows = _forward(model, batch, structure_tokens, device)
        return {index: rows[0, index + TOKEN_OFFSET] for index in wanted}

    batch = sequence_tokens.unsqueeze(0).repeat(len(wanted), 1).clone()
    for row, index in enumerate(wanted):
        batch[row, index + TOKEN_OFFSET] = mask_id
    rows = _forward(model, batch, structure_tokens, device)
    return {index: rows[row, index + TOKEN_OFFSET]
            for row, index in enumerate(wanted)}


def _forward(model, sequence_batch, structure_tokens, device) -> np.ndarray:
    """One forward pass; returns softmaxed logits as `[batch, tokens, vocab]`."""
    import torch

    structure_batch = None
    if structure_tokens is not None:
        structure_batch = structure_tokens.unsqueeze(0).repeat(
            sequence_batch.shape[0], 1)
    with torch.no_grad():
        out = model.forward(sequence_tokens=sequence_batch,
                            structure_tokens=structure_batch)
        probabilities = torch.softmax(out.sequence_logits.float(), dim=-1)
    return np.asarray(probabilities.cpu().numpy(), dtype=float)


def decodable_positions(structure_path, chain_id: str,
                        pdb_id: "str | None" = None) -> np.ndarray:
    """Which manifest indices ESM3 can evaluate, in MANIFEST index space.

    A sequence position is decodable when ESM3's parse reproduces the manifest's
    chain; the structure track adds no further restriction, because a residue
    with no backbone still receives a structure token. Needs no model pass.
    """
    from esm.utils.structure.protein_chain import ProteinChain

    from .structures import _parse_chains

    path = Path(structure_path)
    try:
        chains = _parse_chains(path, str(pdb_id or path.stem))
        native = next((c for c in chains if c.chain_id == str(chain_id)), None)
        if native is None:
            return np.zeros(0, dtype=bool)
        chain = ProteinChain.from_pdb(str(path), chain_id=str(chain_id))
    except Exception:
        return np.zeros(0, dtype=bool)

    if chain.sequence != native.sequence:
        return np.zeros(len(native.sequence), dtype=bool)
    return np.ones(len(native.sequence), dtype=bool)


def check_triplet(sequence: str, indices, expected: str) -> None:
    observed = "".join(sequence[i] if i < len(sequence) else "?" for i in indices)
    if observed != expected:
        raise SequonMismatchError(
            f"the chain reads {observed!r} where the manifest records {expected!r}")


def check_scoreable(probabilities: "dict[int, np.ndarray]", indices) -> None:
    for index in indices:
        if index not in probabilities:
            raise InvalidProbabilityVector(f"index {index} was not evaluated")
        row = np.asarray(probabilities[index], dtype=float)
        if not np.all(np.isfinite(row)):
            raise InvalidProbabilityVector(
                f"row at index {index} has non-finite entries")
        total = float(row.sum())
        if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at index {index} sums to {total:.6g}, not 1")


def sequon_score(probabilities: "dict[int, np.ndarray]", tokenizer,
                 n_index: int, plus1_index: int, plus2_index: int) -> dict:
    """Score one sequon, with the column names every other adapter uses."""
    check_scoreable(probabilities, (n_index, plus1_index, plus2_index))
    index_of = tokenizer.convert_tokens_to_ids

    n_row = np.asarray(probabilities[n_index], dtype=float)
    plus1_row = np.asarray(probabilities[plus1_index], dtype=float)
    plus2_row = np.asarray(probabilities[plus2_index], dtype=float)

    p_n = float(n_row[index_of("N")])
    p_s = float(plus2_row[index_of("S")])
    p_t = float(plus2_row[index_of("T")])
    p_pro = float(plus1_row[index_of("P")])

    return {
        "conditional_sequon_score": 0.5 * (logit(p_n) + logit(p_s + p_t)),
        "conditional_sequon_score_sd": SCORE_SD,
        "n_decoding_orders": N_ORDERS,
        "p_asn_at_n": p_n,
        "p_ser_at_plus2": p_s,
        "p_thr_at_plus2": p_t,
        "p_ser_or_thr_at_plus2": p_s + p_t,
        "p_pro_at_plus1": p_pro,
        "logit_p_asn": logit(p_n),
        "logit_p_ser_or_thr": logit(p_s + p_t),
        "probs_n": n_row.tolist(),
        "probs_plus1": plus1_row.tolist(),
        "probs_plus2": plus2_row.tolist(),
    }
