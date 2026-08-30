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


class DesignFailedError(RuntimeError):
    """ESM3 produced no usable sequence for a chain it agreed to read.

    Separate from ChainUnreadableError: that one means the chain was refused
    before any model ran, this one means generation itself came back wrong.
    Conflating them would hide a generation problem among the parse mismatches,
    which are expected and numerous.
    """


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


def checked_chain(structure_path, chain_id: str,
                  pdb_id: "str | None" = None):
    """ESM3's own parse of one chain, checked against the manifest's.

    Returns `(native_sequence, ProteinChain)`. Shared by scoring and design so
    that a chain refused for one is refused for the other on identical grounds:
    the whole point of the guard is that ESM3's indices and the manifest's
    address the same residues, and design reads the sequon back out by index
    exactly as scoring reads probabilities at it.
    """
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
    return native.sequence, chain


def chain_context(structure_path, chain_id: str, model,
                  pdb_id: "str | None" = None, device: str = "cpu"):
    """Encode one chain, checking ESM3's parse against the manifest's.

    Returns `(sequence, sequence_tokens, structure_tokens)`. A chain whose
    sequence ESM3 reads differently is refused rather than scored at indices
    that address a different residue.
    """
    from esm.sdk.api import ESMProtein

    native_sequence, chain = checked_chain(structure_path, chain_id, pdb_id)
    encoded = model.encode(ESMProtein.from_protein_chain(chain))
    sequence_tokens = encoded.sequence.to(device)
    structure_tokens = (encoded.structure.to(device)
                        if encoded.structure is not None else None)
    return native_sequence, sequence_tokens, structure_tokens


# How many unmasking steps ESM3 takes to write one sequence. ESM3 is a masked
# diffusion model: it fills a fully masked track over several passes, each
# conditioned on what the previous ones committed to. Sampling every position
# from ONE pass would draw from independent marginals rather than the joint --
# the same distinction `independent_calibrated_sampling` names for CARBonAra --
# so the number of steps is a real parameter, not a speed knob.
#
# One step per residue is the highest-fidelity schedule and costs L forward
# passes per design, which at 32 designs a chain is not affordable here. This
# scales with length and is bounded at both ends, and is recorded in the
# provenance of every row so a run made under a different schedule cannot be
# pooled with one made under this.
DESIGN_STEP_DIVISOR = 8
DESIGN_MIN_STEPS, DESIGN_MAX_STEPS = 8, 64

# ESM3 needs its own residue-slot budget, not ProteinMPNN's. That one is 6000,
# calibrated against a small decoder; ESM3-open carries 1.4B parameters and far
# larger activations per token, and reusing 6000 here put a median chain at 25
# designs at once and lost 320 of 427 chains to CUDA OOM. 2000 puts a
# 234-residue chain at 8 and a 1172-residue one at 1.
#
# The budget is a starting point, not a guarantee: it is a linear rule for a
# cost that is not linear, so the batch is halved and retried on OOM as
# ProteinMPNN's has always been.
ESM3_DESIGN_SLOT_BUDGET = 2000


def design_batch(length: int, n_designs: int,
                 max_batch: "int | None" = None) -> int:
    """How many designs to generate at once for a chain of this length."""
    if max_batch:
        return max(1, min(int(max_batch), int(n_designs)))
    return max(1, min(int(n_designs),
                      ESM3_DESIGN_SLOT_BUDGET // max(int(length), 1)))


def free_device_memory(device: str = "cpu") -> None:
    """Release cached blocks so the next chain starts from a clean allocator.

    Without this the failures were dominated by attempts to allocate two
    MEGAbytes on a 40 GB card: once one long chain had filled and fragmented
    the cache, every later chain died regardless of its own size, so a single
    bad chain took the rest of the run with it.
    """
    import torch

    if str(device).startswith("cuda") and torch.cuda.is_available():
        torch.cuda.empty_cache()


def design_steps(length: int, num_steps: "int | None" = None) -> int:
    """Unmasking passes for a chain of this length."""
    if num_steps:
        return max(1, min(int(num_steps), int(length)))
    scaled = int(length) // DESIGN_STEP_DIVISOR
    return max(DESIGN_MIN_STEPS, min(DESIGN_MAX_STEPS, max(scaled, 1)))


def design_sequences(structure_path, chain_id: str, model, n_designs: int,
                     temperature: float, seed: int = 0, device: str = "cpu",
                     pdb_id: "str | None" = None,
                     num_steps: "int | None" = None,
                     max_batch: "int | None" = None,
                     use_batch: bool = True) -> "list[str]":
    """Unconstrained inverse-folding designs for one chain.

    The whole sequence track is masked and rewritten from the structure tokens.
    Nothing is fixed and nothing is biased, so the sequon is free to disappear
    -- otherwise retention would measure the constraint rather than the model.

    Generation goes through ESM3's own `generate` rather than a decoding loop
    written here. That is the vendor's schedule for its own masked-diffusion
    model, and getting it subtly wrong would change every design while still
    producing plausible sequences.

    Designs come back in the manifest's index space, because `checked_chain`
    has already refused any chain whose parse disagrees with it -- so position i
    of a returned string is position i of the manifest.
    """
    import torch
    from esm.sdk.api import ESMProtein, GenerationConfig

    native_sequence, chain = checked_chain(structure_path, chain_id, pdb_id)
    free_device_memory(device)
    length = len(native_sequence)
    steps = design_steps(length, num_steps)
    n_designs = int(n_designs)

    def masked():
        """A fresh protein with the sequence track cleared."""
        protein = ESMProtein.from_protein_chain(chain)
        protein.sequence = None          # the whole track is rewritten
        protein.function_annotations = None
        return protein

    def config():
        return GenerationConfig(track="sequence", num_steps=steps,
                                temperature=float(temperature))

    def checked(sequence, which):
        if not sequence:
            raise DesignFailedError(
                f"ESM3 returned no sequence for chain {chain_id!r} of "
                f"{Path(structure_path).name} on design {which}")
        if len(sequence) != length:
            raise DesignFailedError(
                f"ESM3 returned {len(sequence)} residues for chain {chain_id!r} "
                f"of {Path(structure_path).name} where the chain has {length}; "
                "the manifest's indices would not address the same residues")
        return str(sequence)

    # One design at a time leaves an A100 almost idle: a pilot measured ~13
    # unmasking steps per second whatever the chain length, so the cost was
    # steps x designs and barely touched by length. Batching the designs fills
    # the device, and the batch is bounded by the same residue-slot budget that
    # stops long chains exhausting memory during ProteinMPNN's decoding.
    #
    # The batched and sequential paths draw from the same distribution under the
    # same schedule, but consume randomness differently, so a run is
    # reproducible against itself and not across the two. `describe_generation`
    # records which was used.
    batched = use_batch and hasattr(model, "batch_generate")
    designs: "list[str]" = []

    if batched:
        size = design_batch(length, n_designs, max_batch)
        while len(designs) < n_designs:
            take = min(size, n_designs - len(designs))
            torch.manual_seed(int(seed) + len(designs))
            try:
                with torch.no_grad():
                    generated = model.batch_generate(
                        [masked() for _ in range(take)],
                        [config() for _ in range(take)])
            except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
                # Halve and retry rather than lose the chain, exactly as
                # ProteinMPNN's decoding has always done. A CUDA OOM is a
                # catchable Python exception, unlike the host OOM that the slot
                # budget exists to prevent.
                if size > 1 and "out of memory" in str(exc).lower():
                    free_device_memory(device)
                    size = max(1, size // 2)
                    print(f"    OOM; retrying with batch {size}", flush=True)
                    continue
                raise
            if len(generated) != take:
                raise DesignFailedError(
                    f"ESM3 returned {len(generated)} designs for a batch of "
                    f"{take} on chain {chain_id!r} of "
                    f"{Path(structure_path).name}")
            designs += [checked(getattr(g, "sequence", None), len(designs) + i)
                        for i, g in enumerate(generated)]
    else:
        for k in range(n_designs):
            # Seeded per design, so a design is reproducible individually and a
            # rerun that resumes midway cannot silently draw a different
            # sequence for a site it already has.
            torch.manual_seed(int(seed) + k)
            with torch.no_grad():
                generated = model.generate(masked(), config())
            designs.append(checked(getattr(generated, "sequence", None), k))

    free_device_memory(device)
    if len(designs) != n_designs:
        raise DesignFailedError(
            f"{len(designs)} designs for a request of {n_designs} on chain "
            f"{chain_id!r} of {Path(structure_path).name}")
    return designs


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
    """One forward pass; returns softmaxed logits as `[batch, tokens, vocab]`.

    Autocast on CUDA, as the reference notebook does. ESM3 loads bfloat16
    weights on a GPU, and running it without autocast raises "mat1 and mat2
    must have the same dtype, but got Float and BFloat16" the moment an fp32
    activation meets a bf16 Linear. Omitting it to keep GPU and CPU numerically
    identical does not work: it does not run at all.

    So a GPU run and a CPU run of this model will differ numerically. That is
    the same situation CARBonAra is in for a different reason, and the same rule
    applies -- every arm of a comparison must come from one environment. The
    softmax is taken in fp32 either way, which keeps the difference well below
    the effects being measured, but it is not zero.
    """
    import torch

    structure_batch = None
    if structure_tokens is not None:
        structure_batch = structure_tokens.unsqueeze(0).repeat(
            sequence_batch.shape[0], 1)

    on_cuda = str(device).startswith("cuda")
    with torch.no_grad():
        if on_cuda:
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                out = model.forward(sequence_tokens=sequence_batch,
                                    structure_tokens=structure_batch)
        else:
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
