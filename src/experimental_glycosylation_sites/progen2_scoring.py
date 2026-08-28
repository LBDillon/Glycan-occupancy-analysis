"""ProGen2 scoring — the causal sequence-only model.

Adapted from `score_proteins_progen2_colab.ipynb`, which computes a whole-chain
mean log-likelihood. **That statistic is not what this benchmark measures.** A
per-residue mean over several hundred positions cannot resolve a three-residue
question, which is the objection this module already raises against protein-level
scores. The forward pass is the notebook's; the summation is replaced by reading
the three sequon positions out of it.

## Where it sits

ProGen2 is autoregressive and sequence-only, which fills the empty cell in the
conditioning grid:

                        masked / bidirectional        causal
    structure+sequence  ProteinMPNN, CARBonAra        ESM-IF
    sequence only       ESMC                          ProGen2

That makes ESM-IF the sharp comparison rather than ESMC: both are causal and
prefix-only, so the difference between them isolates what the backbone adds under
identical conditioning. Against ESMC the conditioning also changes, so the two
must never be pooled — `conditioning` is recorded as `autoregressive_prefix`,
the same string ESM-IF uses and deliberately not ESMC's.

## What the conditional is

    P(residue at i | residues 1..i-1)

from one teacher-forced pass, so `conditional_sequon_score_sd` is structurally
zero and `n_decoding_orders` is one, exactly as for ESM-IF.

A consequence worth stating: the asparagine term is **already** motif-blind. The
prefix for position i ends at i-1, so the model has not seen N-X-S/T when it
predicts the N. Only the +1 and +2 terms see any of the motif, through the
residue at the asparagine position.

`mask_mode="joint"` therefore integrates that residue out of those two terms and
leaves the asparagine row untouched — there is nothing to hide from it. **This is
a smaller manipulation than ESMC's or ESM-IF's**, which hide all three positions.
Two of three terms change here, and the first cannot. A near-zero change for this
model means correspondingly less.

## Which sequence

The chain as `structures._parse_chains` reads it, so `model_index` is an ordinal
into it directly and no cross-parser mapping is needed — the same choice made for
ESMC, and for the same reason: sequence-alone versus sequence-plus-structure
stays like-for-like rather than also changing the context window. The cost is
that unresolved loops and truncated termini are absent, as they are for every
other model here.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .mpnn_scoring import EPSILON, PROBABILITY_SUM_TOLERANCE, logit  # noqa: F401

# `hugohrban/*` are the HuggingFace ports the notebook uses; the originals are
# not packaged for `transformers`. base is 764M and runs on CPU; xlarge is 6.4B
# and needs a GPU.
DEFAULT_MODEL = "hugohrban/progen2-base"

CONDITIONING = "autoregressive_prefix"
N_ORDERS = 1
SCORE_SD = 0.0

# The twenty this benchmark scores. ProGen2's vocabulary is larger -- it carries
# direction-control and special tokens -- so these are located in it at load
# time rather than assumed to sit at any particular index.
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


class ChainUnreadableError(ValueError):
    """No usable chain sequence for this structure."""


class SequonMismatchError(ValueError):
    """The sequence does not carry the residues the manifest recorded."""


class InvalidProbabilityVector(ValueError):
    """A row that is not a probability distribution."""


class ContextTooLongError(ValueError):
    """The chain exceeds ProGen2's context window, so it cannot be scored."""


class TokenisationError(RuntimeError):
    """ProGen2 did not tokenise one residue per token, so positions do not align."""


def verify_tokenisation(tokenizer) -> "dict[str, int]":
    """Confirm one token per residue, and locate the twenty. Raise otherwise.

    Everything downstream assumes `logits[i]` predicts the residue at manifest
    index `i`, which holds only if each amino acid is exactly one token. ProGen2
    is documented as per-residue — unlike ProtGPT2, whose BPE vocabulary merges
    residues and would break this — but "documented" is not "checked", and an
    off-by-one here reads a different residue and returns a plausible number.
    """
    ids = tokenizer(STANDARD_AA, add_special_tokens=False)["input_ids"]
    if len(ids) != len(STANDARD_AA):
        raise TokenisationError(
            f"{len(STANDARD_AA)} residues tokenised to {len(ids)} tokens; this "
            "model does not use one token per residue and its positions cannot "
            "be aligned with the manifest's indices")

    index_of = {aa: int(i) for aa, i in zip(STANDARD_AA, ids)}
    if len(set(index_of.values())) != len(STANDARD_AA):
        raise TokenisationError("two residues share a token id")

    # And the map really is the identity it claims to be, read back through the
    # tokeniser rather than trusted.
    for aa, token_id in index_of.items():
        back = tokenizer.convert_ids_to_tokens([token_id])[0]
        if back.strip().upper() != aa:
            raise TokenisationError(
                f"token {token_id} decodes to {back!r}, not {aa!r}")
    return index_of


def load_model(model_name: str = DEFAULT_MODEL, device: str = "cpu"):
    """Load ProGen2 and verify its tokenisation. Returns `(model, tokenizer, bos)`.

    `trust_remote_code` is required: ProGen2's architecture ships with the
    checkpoint rather than living in `transformers`.
    """
    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    verify_tokenisation(tokenizer)

    dtype = torch.bfloat16 if str(device).startswith("cuda") else torch.float32
    model = AutoModelForCausalLM.from_pretrained(
        model_name, trust_remote_code=True, torch_dtype=dtype)
    model = model.to(device).eval()

    # ProGen2 uses '1' to mark N-to-C direction, which the published code and
    # the notebook both prepend. Taken from the tokeniser when it exposes one,
    # because the HuggingFace ports do not all agree with the original.
    bos = tokenizer.bos_token_id
    if bos is None:
        bos = tokenizer(("1"), add_special_tokens=False)["input_ids"][0]
    return model, tokenizer, int(bos)


def chain_sequence(structure_path, chain_id: str,
                   pdb_id: "str | None" = None) -> str:
    """The chain as `_parse_chains` reads it -- what `model_index` indexes."""
    from .structures import _parse_chains

    path = Path(structure_path)
    chains = _parse_chains(path, str(pdb_id or path.stem))
    native = next((c for c in chains if c.chain_id == str(chain_id)), None)
    if native is None or not native.sequence:
        raise ChainUnreadableError(
            f"chain {chain_id!r} absent from {path.name}")
    return native.sequence


def decodable_positions(structure_path, chain_id: str,
                        pdb_id: "str | None" = None) -> np.ndarray:
    """All True: a sequence model has no backbone requirement.

    Returned in manifest index space, and a superset of the structure models'
    scoreable sets — which never widens a matched set, because the pairs were
    frozen on ProteinMPNN's scoreability.
    """
    try:
        return np.ones(len(chain_sequence(structure_path, chain_id, pdb_id)),
                       dtype=bool)
    except ChainUnreadableError:
        return np.zeros(0, dtype=bool)


def context_limit(model) -> "int | None":
    """ProGen2's context window, from the config. None when it cannot be read."""
    config = getattr(model, "config", None)
    for attribute in ("n_positions", "max_position_embeddings", "n_ctx"):
        value = getattr(config, attribute, None)
        if isinstance(value, int) and value > 0:
            return value
    return None


def check_context_length(sequence: str, model) -> None:
    """Refuse a chain longer than the context window, with a legible reason.

    Without this the failure is a raw tensor-size mismatch from inside the
    model, which says nothing about which chain or why. ProGen2-base carries
    2048 positions and the manifest contains chains beyond that -- 3JAV:A is
    2328 residues -- so this is a real exclusion, not a defensive one.

    Truncating the prefix instead was rejected: for a causal model the prefix
    *is* the conditioning, so a truncated one answers a different question, and
    silently answering a different question is worse than declining.
    """
    limit = context_limit(model)
    if limit is not None and len(sequence) + 1 > limit:
        raise ContextTooLongError(
            f"chain is {len(sequence)} residues and needs {len(sequence) + 1} "
            f"tokens with BOS, over ProGen2's {limit}-token context")


def conditional_probabilities(sequence: str, model, tokenizer, bos: int,
                              device: str = "cpu") -> np.ndarray:
    """`[L, vocab]` of P(residue at i | residues 1..i-1), one forward pass.

    The alignment, which is the whole measurement:

        tokens   [BOS, t_1, ..., t_L]        length L+1
        logits[j] predicts token j+1
        logits[i] therefore predicts t_{i+1}, the residue at manifest index i

    so row `i` of the result is the distribution for manifest index `i` with no
    further shift. Prepending BOS is what makes that hold; without it row 0 would
    be undefined and everything after it off by one.
    """
    import torch

    check_context_length(sequence, model)
    ids = tokenizer(sequence, add_special_tokens=False)["input_ids"]
    tokens = torch.tensor([bos] + list(ids), device=device).unsqueeze(0)

    with torch.no_grad():
        logits = model(input_ids=tokens).logits[0]
    # Drop the final row: it predicts a residue past the end of the chain.
    log_probs = torch.log_softmax(logits[:-1].float(), dim=-1)
    return np.asarray(log_probs.exp().cpu().numpy(), dtype=float)


def check_scoreable(probabilities: np.ndarray, indices) -> None:
    """Raise unless all three rows are genuine probability distributions."""
    length = probabilities.shape[0]
    for index in indices:
        if not 0 <= index < length:
            raise IndexError(
                f"index {index} outside the chain's {length} residues")
        row = probabilities[index]
        if not np.all(np.isfinite(row)):
            raise InvalidProbabilityVector(
                f"row at index {index} has non-finite entries")
        total = float(row.sum())
        if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at index {index} sums to {total:.6g}, not 1")


def check_triplet(sequence: str, indices, expected: str) -> None:
    observed = "".join(sequence[i] if i < len(sequence) else "?" for i in indices)
    if observed != expected:
        raise SequonMismatchError(
            f"the chain reads {observed!r} where the manifest records {expected!r}")


def sequon_score(probabilities: np.ndarray, aa_index: "dict[str, int]",
                 n_index: int, plus1_index: int, plus2_index: int) -> dict:
    """Score one sequon, with the column names every other adapter uses.

    The same statistic: the mean of the log odds of asparagine at the first
    position and of serine-or-threonine at the third, the middle residue excluded
    because any residue but proline permits a sequon.

    Vectors are returned over ProGen2's full vocabulary rather than restricted to
    the twenty, matching what ESM-IF does with its 35 tokens: the row is the
    model's actual output, and renormalising it onto a subset would report a
    distribution the model never produced.
    """
    check_scoreable(probabilities, (n_index, plus1_index, plus2_index))

    p_n = float(probabilities[n_index, aa_index["N"]])
    p_s = float(probabilities[plus2_index, aa_index["S"]])
    p_t = float(probabilities[plus2_index, aa_index["T"]])
    p_pro = float(probabilities[plus1_index, aa_index["P"]])

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
        "probs_n": probabilities[n_index].tolist(),
        "probs_plus1": probabilities[plus1_index].tolist(),
        "probs_plus2": probabilities[plus2_index].tolist(),
    }


# --------------------------------------------------------------------------
# Joint masking, for the motif-hidden arm.
# --------------------------------------------------------------------------

MASK_MODES = ("single", "joint")
DEFAULT_MASK_MODE = "single"

# The same string ESM-IF uses for the same operation, so the two are pooled or
# separated on the merits rather than by accident of naming.
CONDITIONING_JOINT = "autoregressive_prefix_marginalised"

# 20 variants of one chain is a small batch, but chains here run to 1287
# residues and activation memory is batch x length. Capped, and halved on OOM,
# for the reason recorded in esmif_scoring: a host OOM is delivered by the
# kernel and cannot be caught, so the cap has to be conservative up front.
DEFAULT_MARGINAL_BATCH = 8


def conditioning(mask_mode: str) -> str:
    if mask_mode not in MASK_MODES:
        raise ValueError(f"mask_mode must be one of {MASK_MODES}, got {mask_mode!r}")
    return CONDITIONING if mask_mode == "single" else CONDITIONING_JOINT


def marginalised_probabilities(sequence: str, model, tokenizer, bos: int,
                               aa_index: "dict[str, int]", indices,
                               device: str = "cpu",
                               max_batch: int = DEFAULT_MARGINAL_BATCH
                               ) -> np.ndarray:
    """`conditional_probabilities`, with the sequon asparagine integrated out.

    What this hides, and what it cannot: a causal model's prefix for position i
    ends at i-1, so **the asparagine term is already motif-blind** — there is
    nothing to hide from it, and its row is returned unchanged. Only the +1 and
    +2 terms see any of the motif, through the residue at the asparagine
    position, and those are the rows this replaces with

        P(x at j) = sum_a P(a at i | 1..i-1) * P(x at j | 1..i-1, a at i, ...)

    weighting each substitution by the model's own predictive distribution at i
    rather than uniformly, so the result is a genuine marginal rather than an
    average over an arbitrary set.

    That makes ProGen2's masking manipulation **smaller than ESMC's or
    ESM-IF's**, which hide all three positions. Here two of the three terms
    change and the first cannot. A near-zero change for this model therefore
    means less than it would for them, and should not be read as "masking does
    not matter".

    The X position is left native in each variant. Marginalising it too would be
    400 sequences rather than 20, and any residue but proline permits a sequon,
    so it carries little of the motif to hide.
    """
    import torch

    n_index, plus1_index, plus2_index = (int(i) for i in indices)
    base = conditional_probabilities(sequence, model, tokenizer, bos, device=device)
    if max(n_index, plus1_index, plus2_index) >= len(sequence):
        return base

    letters = list(STANDARD_AA)
    weights = np.array([base[n_index, aa_index[aa]] for aa in letters], dtype=float)
    total = weights.sum()
    if total <= 0:
        return base
    weights /= total                      # renormalised over the twenty

    variants = []
    for aa in letters:
        altered = sequence[:n_index] + aa + sequence[n_index + 1:]
        variants.append([bos] + list(
            tokenizer(altered, add_special_tokens=False)["input_ids"]))
    tokens = torch.tensor(variants, device=device)

    size = max(1, min(len(letters), int(max_batch)))
    while True:
        try:
            pieces = []
            with torch.no_grad():
                for start in range(0, tokens.shape[0], size):
                    logits = model(input_ids=tokens[start:start + size]).logits
                    pieces.append(torch.softmax(logits[:, :-1].float(), dim=-1).cpu())
            stacked = torch.cat(pieces, dim=0).numpy()
            break
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            if size > 1 and "out of memory" in str(exc).lower():
                if str(device).startswith("cuda"):
                    torch.cuda.empty_cache()
                size = max(1, size // 2)
                print(f"    OOM in marginalisation; retrying with batch {size}",
                      flush=True)
                continue
            raise

    marginal = base.copy()
    for index in (plus1_index, plus2_index):
        marginal[index] = np.tensordot(weights, stacked[:, index, :], axes=(0, 0))
    return marginal
