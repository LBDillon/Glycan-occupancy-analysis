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

RETENTION_CATEGORIES = (
    "full_sequon_retained",
    "asn_retained_motif_lost",
    "ser_thr_retained_motif_lost",
    "proline_introduced_at_x",
    "complete_motif_loss",
)


def design_sequences(
    pdb_path: Path,
    chain_id: str,
    model,
    n_designs: int,
    temperature: float,
    device: str = "cpu",
    seed: int = 0,
    max_batch: "int | None" = None,
) -> list[str]:
    """Unconstrained designs for one chain, as one-letter sequences.

    All `n_designs` are decoded as one batch: `sample()` already carries a batch
    dimension, and the featurised tensors describe a single backbone, so tiling
    them turns 32 sequential decodes into one. Sampling one at a time projected
    to roughly 32 hours over this corpus, which no Colab session survives.

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
    size = min(max_batch or n_designs, n_designs)
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
