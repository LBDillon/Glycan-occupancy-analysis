"""Sequon retention in unconstrained ProteinMPNN designs.

The conditional score asks what probability the model holds at a site. This asks
what it actually does when generating a sequence, which is a different question:
sampling combines many residue decisions, and a motif can survive or vanish for
reasons that have little to do with the model's preference at that one position.

Retention is therefore a secondary, operational outcome. It is descriptive until
it has a justified margin of its own, and it is not on its own evidence that the
model does or does not recognise N-X-S/T as a motif.

Nothing is fixed or biased during generation: the designs are unconstrained, so
the sequon is free to disappear.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .mpnn_scoring import ALPHABET, _prepare_environment

# The generation settings used by the earlier scoping analysis - a quick baseline
# run over a handful of proteins, not a published preprint - recovered from its
# repository history so retention here is comparable with what it saw. The key
# name is kept as "preprint" only because it is already written into result files.
SCOPING_CONDITION = {"name": "preprint", "temperature": 0.1, "n_designs": 8}
PREPRINT_CONDITION = SCOPING_CONDITION  # backwards-compatible alias

# Eight designs give a per-site retention estimate with a standard error near
# 0.18 at p=0.5 — far too coarse for a site-level analysis. The standardised
# condition keeps that temperature and raises only the sample count, to
# 32, which brings that standard error to about 0.09 at a cost the full corpus
# can absorb. Designs are independent draws, so the first eight of a 32-design
# run are themselves a valid scoping-condition sample and are reported as such
# rather than generated separately.
STANDARD_CONDITION = {"name": "standard", "temperature": 0.1, "n_designs": 32}

# Activation memory during generation scales with batch x chain length, and these
# chains vary six-fold: 299 residues at the median, 1287 at the longest. A fixed
# batch of 32 is fine for a short chain and 41,000 residue-slots for a long one,
# which is what killed 63 of 64 ARC array tasks with OUT_OF_MEMORY.
#
# That failure cannot be caught and retried, either: a host OOM is delivered by
# the kernel's OOM killer or the SLURM cgroup, so the process dies without
# Python seeing an exception. The batch has to be bounded BEFORE allocating,
# not reduced after failing.
#
# So the batch is chosen from a budget on the product. 6000 slots keeps a
# 1287-residue chain at batch 4 and still lets a 200-residue chain run 30 at
# once, which is where the batching pays off anyway.
DESIGN_SLOT_BUDGET = 6000


def batch_for_length(length: int, n_designs: int,
                     max_batch: "int | None" = None) -> int:
    """How many designs to decode at once for a chain of this length."""
    if max_batch:
        return max(1, min(int(max_batch), n_designs))
    return max(1, min(n_designs, DESIGN_SLOT_BUDGET // max(int(length), 1)))


RETENTION_CATEGORIES = (
    "full_sequon_retained",
    "asn_retained_motif_lost",
    "ser_thr_retained_motif_lost",
    "proline_introduced_at_x",
    "complete_motif_loss",
)


def design_mask(length: int, fixed_positions) -> np.ndarray:
    """1.0 where the model may design, 0.0 where the native residue is kept.

    An out-of-range position raises rather than being ignored: silently dropping
    it would protect nothing while the run still reports success, and the design
    would look entirely reasonable.
    """
    mask = np.ones(int(length), dtype=np.float32)
    for position in fixed_positions:
        index = int(position)
        if not 0 <= index < length:
            raise ValueError(
                f"fixed position {index} is outside a chain of length {length}")
        mask[index] = 0.0
    return mask


def design_sequences(
    pdb_path: Path,
    chain_id: str,
    model,
    n_designs: int,
    temperature: float,
    device: str = "cpu",
    seed: int = 0,
    max_batch: "int | None" = None,
    fixed_positions: "list[int] | None" = None,
) -> list[str]:
    """Designs for one chain, as one-letter sequences.

    `fixed_positions` holds chain indices at their native residue, via
    ProteinMPNN's own `chain_M_pos` mask. Passing None keeps the historical
    behaviour exactly: nothing is fixed and the sequon is free to disappear,
    which is what the retention measurement requires.

    All `n_designs` are decoded as one batch: `sample()` already carries a batch
    dimension, and the featurised tensors describe a single backbone, so they are
    tiled rather than looped over.

    **How much this buys depends entirely on the device**, and an earlier version
    of this note claimed a speedup that does not hold on CPU. Measured on one
    124-residue chain, CPU, 2026-08-26: 1 design 0.20 s, 8 designs 1.12 s, 32
    designs 4.28 s. That is about 1.5x over sampling one at a time, not the order
    of magnitude batching gives on a GPU, because a single sequence already
    saturates the cores. The batch is still worth having, and it is not what
    makes a corpus run tractable.

    What does is loading the model once and reusing it -- 0.85 s paid once rather
    than per chain -- and keeping everything in memory instead of round-tripping
    through files. Per chain the cost tracks length: 4.2 s for 124 residues,
    7.0 s for 200, against a corpus median of 332.

    Each batch row draws its own `randn`, so the decoding orders stay independent
    and the designs are independent samples rather than 32 copies.

    `max_batch` caps the batch when memory is tight; on OOM the batch is halved
    and retried rather than losing the chain.
    """
    _prepare_environment()
    import torch

    from protein_mpnn_utils import StructureDatasetPDB, parse_PDB, tied_featurize

    parsed = parse_PDB(str(pdb_path), input_chain_list=[chain_id])
    protein = StructureDatasetPDB(parsed, truncate=None, max_length=20000)[0]

    out = tied_featurize(
        [protein], device, {protein["name"]: ([chain_id], [])},
        None, None, None, None, None, ca_only=False,
    )
    X, S, mask = out[0], out[1], out[2]
    chain_M, chain_encoding_all, chain_M_pos, residue_idx = out[4], out[5], out[10], out[12]

    if fixed_positions:
        # chain_M_pos is ProteinMPNN's own design mask: 1 designs the position,
        # 0 keeps the native residue. Masking before sampling is not the same as
        # repairing the output afterwards -- repairing would let the model
        # condition on residues it was about to overwrite.
        keep = design_mask(chain_M_pos.shape[1], fixed_positions)
        chain_M_pos = chain_M_pos * torch.from_numpy(keep).to(chain_M_pos.device)
    # sample() dereferences these regardless of their defaults, so they must be
    # the real featurised tensors: omit_AA_mask at 11, bias_by_res at 18.
    omit_AA_mask, bias_by_res = out[11], out[18]

    # 'X' is the unknown-residue token and is never a design output
    omit = np.array([aa == "X" for aa in ALPHABET], dtype=np.float32)
    bias = np.zeros(len(ALPHABET), dtype=np.float32)
    generator = torch.Generator(device="cpu").manual_seed(seed)

    def tile(tensor, size):
        return tensor.repeat(size, *([1] * (tensor.dim() - 1)))

    def decode(size: int) -> list[str]:
        randn = torch.randn((size,) + tuple(chain_M.shape[1:]),
                            generator=generator).to(device)
        with torch.no_grad():
            sample = model.sample(
                tile(X, size), randn, tile(S, size), tile(chain_M, size),
                tile(chain_encoding_all, size), tile(residue_idx, size),
                mask=tile(mask, size), temperature=temperature, omit_AAs_np=omit,
                bias_AAs_np=bias, chain_M_pos=tile(chain_M_pos, size),
                omit_AA_mask=tile(omit_AA_mask, size),
                bias_by_res=tile(bias_by_res, size),
            )
        return ["".join(ALPHABET[i] for i in row)
                for row in sample["S"].cpu().numpy()]

    sequences: list[str] = []
    remaining = n_designs
    size = batch_for_length(X.shape[1], n_designs, max_batch)
    while remaining > 0:
        try:
            chunk = decode(min(size, remaining))
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            if size > 1 and "out of memory" in str(exc).lower():
                if device.startswith("cuda"):
                    torch.cuda.empty_cache()
                size = max(1, size // 2)
                print(f"    OOM; retrying with batch {size}", flush=True)
                continue
            raise
        sequences.extend(chunk)
        remaining -= len(chunk)
    return sequences


def classify_retention(designs: list[str], n_index: int, plus1_index: int,
                       plus2_index: int) -> dict:
    """Per-design outcomes at one sequon, summarised as fractions.

    The categories are deliberately not mutually exclusive at the residue level —
    a design can keep the asparagine and lose the hydroxyl — so each is reported
    as its own fraction rather than forced into a single partition. Only
    `full_sequon_retained` and `complete_motif_loss` are mutually exclusive.
    """
    counts = dict.fromkeys(RETENTION_CATEGORIES, 0)
    n_total = 0

    for design in designs:
        if max(n_index, plus1_index, plus2_index) >= len(design):
            continue
        n_total += 1
        asn = design[n_index] == "N"
        hydroxyl = design[plus2_index] in ("S", "T")
        proline = design[plus1_index] == "P"
        full = asn and hydroxyl and not proline

        if full:
            counts["full_sequon_retained"] += 1
        if asn and not full:
            counts["asn_retained_motif_lost"] += 1
        if hydroxyl and not full:
            counts["ser_thr_retained_motif_lost"] += 1
        if proline:
            counts["proline_introduced_at_x"] += 1
        if not asn and not hydroxyl:
            counts["complete_motif_loss"] += 1

    if n_total == 0:
        return {"n_designs_scored": 0, **{f"frac_{k}": None for k in RETENTION_CATEGORIES}}

    result = {"n_designs_scored": n_total}
    result.update({f"frac_{k}": round(counts[k] / n_total, 4) for k in RETENTION_CATEGORIES})
    result.update({f"n_{k}": counts[k] for k in RETENTION_CATEGORIES})
    return result


# --------------------------------------------------------------------------
# Incomplete backbones affect generation too, by a different route.
#
# sample() computes chain_mask = chain_mask * chain_M_pos * mask and then commits
#
#     S_t = S_t * chain_mask + S_true * (1 - chain_mask)
#
# so a residue whose backbone is incomplete is never redesigned: it is written
# back at its native identity in every design. On a real fixture whose sequon
# +2 oxygen is missing, all three designs return native T at that position
# while the asparagine two residues earlier is redesigned away.
#
# The distortion therefore depends on which residue is affected rather than
# pushing retention one way. A frozen +2 inflates ser_thr_retained_motif_lost;
# a frozen N inflates asn_retained_motif_lost; all three frozen would report a
# fully retained sequon that the model never had the chance to alter. In every
# case the number describes the parser, not the model, so these sites are
# excluded rather than corrected.
# --------------------------------------------------------------------------

RETENTION_REQUIRES_SCOREABLE = True
