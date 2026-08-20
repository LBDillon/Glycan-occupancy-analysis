"""Masked-language-model probabilities at sequon positions, from ESMC.

The third model in the benchmark, and the first that sees **no structure at
all**. ProteinMPNN and ESM-IF both condition on a backbone; ESMC conditions only
on surrounding sequence. So it answers a question neither of them can:

    does sequence context alone distinguish occupied sequons from structurally
    matched sequons carrying no glycan?

If it does, the effect the structure-conditioned models report need not be
structural. If it does not, the structure-conditioned result is doing work that
sequence alone cannot.

## What is scored, and on which sequence

The chain sequence as `structures._parse_chains` reads it — the same string the
manifest's `model_index` is an ordinal into. That is a deliberate choice over the
full UniProt sequence: it keeps model indices, scoreability and the matched pairs
byte-identical to the structure-conditioned models, so "sequence alone versus
sequence plus structure" is a like-for-like comparison rather than a comparison
that also changes the context window. The cost is that ESMC sees only the
resolved chain — unresolved loops and truncated termini are absent, exactly as
they are absent for the other two models.

## Masking

`P(N)` at the first sequon position and `P(S)+P(T)` at the third, read from the
model's distribution at a masked position. Two schemes:

`single` (default) masks one position at a time, so each score is
`P(residue at i | every other native residue)`. This is the closest sequence-only
analogue of ProteinMPNN's conditional.

`joint` masks all three sequon positions at once. It exists because `single`
leaves a confound: masking only the +2 residue still shows the model a native
asparagine two positions upstream, and N-X-S/T is a heavily learned motif, so the
model can infer S/T from the N rather than from anything about this site. That
inflates both arms of a matched pair and therefore largely cancels in the paired
contrast, but it compresses the dynamic range and could hide a real difference.
`joint` removes it, at the cost of conditioning on strictly less context.

Report `single` as primary and `joint` as a sensitivity; they are different
estimands, not competing estimates of one.

## Why there is no SequenceDesigner

ESMC is a masked language model, not an inverse-folding model. Sampling from it
would produce sequences conditioned on sequence context rather than on a
backbone, so "retention" would measure something other than what it measures for
the other two models. It implements `SequonScorer` only.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .mpnn_scoring import (EPSILON, PROBABILITY_SUM_TOLERANCE,  # noqa: F401
                           InvalidProbabilityVector, _prepare_environment, logit)

DEFAULT_MODEL = "esmc_300m"

# A masked forward pass is deterministic: one "order", zero spread. Named so the
# score dict cannot drift from what this docstring claims.
N_ORDERS = 1
SCORE_SD = 0.0

MASK_MODES = ("single", "joint")
DEFAULT_MASK_MODE = "single"

# The tokenizer prepends <cls>, so sequence position i is token index i + 1.
# Asserted against the tokenizer at load rather than trusted -- a silent offset
# here scores the neighbouring residue and still returns a plausible number.
TOKEN_OFFSET = 1

# Batch of masked variants per forward. Each variant is a full copy of the
# chain, so this bounds memory on long chains rather than sequence count.
DEFAULT_BATCH = 16


class SequonMismatchError(ValueError):
    """The scored residues do not reproduce the manifest's triplet."""


def conditioning(mask_mode: str) -> str:
    return f"masked_sequence_{mask_mode}"


def load_model(device: str = "cpu", model_name: str = DEFAULT_MODEL):
    """Load ESMC in eval mode and verify the token offset. Returns (model, tokenizer).

    Requires EvolutionaryScale's `esm` package, which collides on the import name
    `esm` with `fair-esm` (ESM-IF). The two cannot be installed together; see
    docs/third_model_esmc.md.
    """
    _prepare_environment()
    import torch

    try:
        from esm.models.esmc import ESMC
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise ImportError(
            "ESMC needs EvolutionaryScale's `esm` package (pip install esm). "
            "It installs a top-level module named `esm` and will shadow "
            "`fair-esm`, breaking the ESM-IF adapter -- use a separate "
            "environment for each.") from exc

    model = ESMC.from_pretrained(model_name, device=torch.device(device)).eval()
    _assert_token_offset(model.tokenizer)
    return model, model.tokenizer


def _assert_token_offset(tokenizer) -> None:
    """Round-trip a probe sequence through the tokenizer and back.

    The same class of check the ProteinMPNN alphabet defect went undetected for:
    an assumption about how a model indexes its own input is not a fact until
    something reproduces the input from it.
    """
    probe = "MKTAYIAKQRNLTSHFSRQ"
    ids = tokenizer(probe, return_tensors="pt")["input_ids"][0].tolist()
    recovered = "".join(
        tokenizer.convert_ids_to_tokens([ids[i + TOKEN_OFFSET]])[0]
        for i in range(len(probe)))
    if recovered != probe:
        raise RuntimeError(
            f"ESMC token offset {TOKEN_OFFSET} does not round-trip: "
            f"{recovered!r} != {probe!r}. Scoring would read the wrong residue.")


def chain_sequence(structure_path: Path, chain_id: str,
                   pdb_id: "str | None" = None) -> str:
    """The chain as the manifest indexes it, so no index mapping is needed."""
    from .structures import _parse_chains

    path = Path(structure_path)
    chains = _parse_chains(path, str(pdb_id or path.stem))
    chain = next((c for c in chains if c.chain_id == str(chain_id)), None)
    if chain is None:
        raise ValueError(f"chain {chain_id!r} absent from {path.name}")
    return chain.sequence


def decodable_positions(structure_path: Path, chain_id: str,
                        pdb_id: "str | None" = None) -> np.ndarray:
    """Every residue of the chain: a sequence model has no backbone requirement.

    ESMC's scoreable set is therefore a superset of the structure-conditioned
    models'. Sites they cannot score are still excluded from the comparison,
    because the matched pairs were frozen on their scoreability -- so this being
    permissive never widens a matched set.
    """
    try:
        return np.ones(len(chain_sequence(structure_path, chain_id, pdb_id)), dtype=bool)
    except Exception:
        return np.zeros(0, dtype=bool)


def masked_distributions(sequence: str, groups: "list[tuple[int, ...]]", model,
                         tokenizer, device: str = "cpu",
                         batch_size: int = DEFAULT_BATCH) -> "dict[tuple, dict]":
    """Distributions at masked positions, one forward per group.

    Each group is a tuple of sequence indices masked together: `(i,)` for the
    `single` scheme, `(n, plus1, plus2)` for `joint`. Returns
    `{group: {position: probability vector}}`.
    """
    _prepare_environment()
    import torch

    encoded = tokenizer(sequence, return_tensors="pt")["input_ids"]
    length = len(sequence)

    results: "dict[tuple, dict]" = {}
    for start in range(0, len(groups), batch_size):
        chunk = groups[start:start + batch_size]
        batch = encoded.repeat(len(chunk), 1).clone()
        for row, group in enumerate(chunk):
            for index in group:
                if not 0 <= index < length:
                    raise IndexError(
                        f"position {index} outside chain of length {length}")
                batch[row, index + TOKEN_OFFSET] = tokenizer.mask_token_id

        with torch.no_grad():
            logits = model(batch.to(device)).sequence_logits

        probabilities = torch.softmax(logits.float(), dim=-1).cpu().numpy()
        for row, group in enumerate(chunk):
            results[group] = {
                index: probabilities[row, index + TOKEN_OFFSET] for index in group
            }
    return results


def check_scoreable(vectors) -> None:
    """Raise unless every row is a genuine probability distribution."""
    for index, row in vectors.items():
        total = float(row.sum())
        if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at sequence index {index} sums to {total:.6g}, not 1")
        if float(row.min()) < 0.0 or float(row.max()) > 1.0 + PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at sequence index {index} has values outside [0, 1]")


def sequon_score(vectors: dict, tokenizer, n_index: int, plus1_index: int,
                 plus2_index: int) -> dict:
    """Score one sequon, with the column names the other scorers use."""
    check_scoreable(vectors)
    index_of = tokenizer.convert_tokens_to_ids

    p_n = float(vectors[n_index][index_of("N")])
    p_s = float(vectors[plus2_index][index_of("S")])
    p_t = float(vectors[plus2_index][index_of("T")])
    p_pro = float(vectors[plus1_index][index_of("P")])

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
        "probs_n": vectors[n_index].tolist(),
        "probs_plus1": vectors[plus1_index].tolist(),
        "probs_plus2": vectors[plus2_index].tolist(),
    }
