"""Conditional ProteinMPNN probabilities at sequon positions.

The question is how strongly ProteinMPNN prefers to keep the motif-forming
residues where they already are, so the measurement is the conditional
probability

    P(amino acid at position | backbone and the rest of the native sequence)

read off the native structure without generating anything and without altering
the sequon. This is not the whole-protein score: three residues contribute
almost nothing to an average over a few hundred, so a protein-level number
cannot answer a site-level question.

Two properties of the model matter for the benchmark's validity.

ProteinMPNN reads only ATOM records — glycans are HETATM and are invisible to
it — so the occupancy label cannot reach the model through its input. That is
what makes an occupied-versus-unoccupied comparison meaningful rather than
circular.

Its conditional probabilities depend on the decoding order, which is sampled.
A single order is therefore one draw from a distribution, so scores are averaged
over a fixed panel of seeded orders and the spread across them is retained.

One failure mode has to be guarded explicitly. `conditional_probs` allocates a
zero tensor of log-probabilities and fills only the positions that survive
`chain_M * mask`, where `mask` is zero wherever the backbone is incomplete. A
residue with a missing N, CA, C or O is therefore returned as a row of zeros,
and exponentiating that row yields twenty-one ones — not a distribution, but a
value that propagates silently into a score of about +13.8. A manifest can
confirm that a residue has coordinates without confirming that ProteinMPNN
accepted it, so both the computed set and the returned vectors are checked here
rather than trusted.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

# ProteinMPNN's own alphabet, in its own order. Index into probability vectors
# with this, never with a hand-written ordering.
ALPHABET = "ARNDCQEGHILKMFPSTWYVX"
AA_INDEX = {aa: i for i, aa in enumerate(ALPHABET)}

# The checkpoint used by the preprint, so the conditional scores and the later
# design-retention analysis describe the same model.
DEFAULT_MODEL = "v_48_020"
DEFAULT_DECODING_ORDERS = 8

# Probabilities are clamped before the logit so a saturated position yields a
# large finite score rather than an infinity that would poison every average
# downstream.
EPSILON = 1e-6

# A softmax in float32, averaged over eight orders, stays within a few times
# 1e-6 of one. A tolerance of 1e-3 is far looser than that and still four
# orders of magnitude tighter than the 21.0 an unfilled row produces.
PROBABILITY_SUM_TOLERANCE = 1e-3


class IncompleteBackboneError(ValueError):
    """A requested residue was rejected by ProteinMPNN's own backbone mask."""


class InvalidProbabilityVector(ValueError):
    """A returned row is not a probability distribution."""


def _prepare_environment() -> None:
    """Single-threaded, and tolerant of the duplicate libomp on macOS.

    torch and miniforge each ship a libomp. Pinning one thread removes the
    threading hazard the duplicate warning is actually about, and makes the
    numbers reproducible.
    """
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def _ensure_importable(proteinmpnn_dir: "Path | None") -> None:
    """Put ProteinMPNN's own modules on the import path."""
    import sys

    if proteinmpnn_dir is None:
        return
    resolved = str(Path(proteinmpnn_dir).resolve())
    if resolved not in sys.path:
        sys.path.insert(0, resolved)


def load_model(proteinmpnn_dir: Path, model_name: str = DEFAULT_MODEL, device: str = "cpu"):
    """Load a ProteinMPNN checkpoint with backbone noise disabled."""
    _prepare_environment()
    import torch

    proteinmpnn_dir = Path(proteinmpnn_dir).resolve()
    _ensure_importable(proteinmpnn_dir)
    from protein_mpnn_utils import ProteinMPNN

    checkpoint_path = proteinmpnn_dir / "vanilla_model_weights" / f"{model_name}.pt"
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)

    model = ProteinMPNN(
        num_letters=21,
        node_features=128,
        edge_features=128,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        # Inference is on the native backbone: the training-time noise level is a
        # property of the checkpoint, not something to reapply here.
        augment_eps=0.0,
        k_neighbors=checkpoint["num_edges"],
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()
    return model


def conditional_probabilities(
    pdb_path: Path,
    chain_id: str,
    model,
    device: str = "cpu",
    n_decoding_orders: int = DEFAULT_DECODING_ORDERS,
    seed: int = 0,
    positions: "list[int] | None" = None,
) -> np.ndarray:
    """Per-order conditional probabilities for one chain, shape [orders, L, 21].

    `positions` restricts which residues are computed. ProteinMPNN's
    conditional_probs runs a full decoder pass for every position it is asked
    about, so a 300-residue chain costs 300 passes when three are wanted. The
    per-position result depends only on that position and the decoding order, not
    on which other positions were requested, so restricting the set returns
    identical values for them — verified against the unrestricted computation.
    Rows outside `positions` are left at zero and must not be read.

    Returns `(probabilities, computed)`. `computed` is a boolean array of length
    L that is true only where ProteinMPNN actually ran a decoder pass, which is
    where the request and the model's own backbone mask agree. Callers must
    consult it: the rows it marks false are zeros that exponentiate to ones, and
    are indistinguishable from a real result by shape or dtype alone.
    """
    _prepare_environment()
    import torch

    from protein_mpnn_utils import StructureDatasetPDB, parse_PDB, tied_featurize

    parsed = parse_PDB(str(pdb_path), input_chain_list=[chain_id])
    dataset = StructureDatasetPDB(parsed, truncate=None, max_length=20000)
    protein = dataset[0]

    # tied_featurize returns a 20-tuple. Index it explicitly: positional unpacking
    # with a star silently shifts every field if the signature ever changes, and
    # residue_idx sits at 12, well past where a short unpack would reach.
    out = tied_featurize(
        [protein], device, {protein["name"]: ([chain_id], [])},
        None, None, None, None, None, ca_only=False,
    )
    X, S, mask = out[0], out[1], out[2]
    chain_M, chain_encoding_all, chain_M_pos, residue_idx = out[4], out[5], out[10], out[12]

    design_mask = chain_M * chain_M_pos
    if positions is not None:
        length = design_mask.shape[1]
        wanted = torch.zeros_like(design_mask)
        for index in positions:
            index = int(index)
            if not 0 <= index < length:
                raise IndexError(
                    f"position {index} outside chain {chain_id} of length {length}")
            wanted[0, index] = 1.0
        design_mask = design_mask * wanted

    # conditional_probs recomputes this internally as chain_M * mask; mirroring it
    # here is what lets a caller distinguish an unfilled row from a real one.
    computed = (design_mask * mask)[0].cpu().numpy() > 0

    per_order = []
    generator = torch.Generator(device="cpu").manual_seed(seed)
    with torch.no_grad():
        for _ in range(n_decoding_orders):
            randn = torch.randn(chain_M.shape, generator=generator).to(device)
            log_probs = model.conditional_probs(
                X, S, mask, design_mask, residue_idx, chain_encoding_all, randn, False
            )
            per_order.append(torch.exp(log_probs)[0].cpu().numpy())
    return np.stack(per_order), computed


def logit(probability: float) -> float:
    p = float(np.clip(probability, EPSILON, 1.0 - EPSILON))
    return float(np.log(p / (1.0 - p)))


def check_scoreable(
    probabilities: np.ndarray,
    indices: "tuple[int, int, int]",
    computed: "np.ndarray | None" = None,
) -> None:
    """Raise unless all three sequon rows are genuine probability distributions.

    Two independent checks, because either alone can be passed by a bad result.
    The `computed` set catches a position ProteinMPNN declined to decode; the
    sum check catches anything else that leaves a row malformed, including a
    row of zeros that has been exponentiated into ones. A site that fails is
    excluded and recorded, never scored — an unfilled row scores near +13.8 and
    would otherwise sit in the data looking like an unusually confident site.
    """
    if computed is not None:
        rejected = [i for i in indices if not bool(computed[i])]
        if rejected:
            raise IncompleteBackboneError(
                "ProteinMPNN did not decode model "
                f"{'index' if len(rejected) == 1 else 'indices'} {rejected}: "
                "incomplete backbone (missing N, CA, C or O)")

    for order in range(probabilities.shape[0]):
        for index in indices:
            row = probabilities[order, index]
            total = float(row.sum())
            if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
                raise InvalidProbabilityVector(
                    f"row at model index {index}, decoding order {order} sums to "
                    f"{total:.6g}, not 1")
            if float(row.min()) < 0.0 or float(row.max()) > 1.0 + PROBABILITY_SUM_TOLERANCE:
                raise InvalidProbabilityVector(
                    f"row at model index {index}, decoding order {order} has values "
                    f"outside [0, 1]: min {float(row.min()):.6g}, max {float(row.max()):.6g}")


def sequon_score(
    probabilities: np.ndarray,
    n_index: int,
    plus1_index: int,
    plus2_index: int,
    computed: "np.ndarray | None" = None,
) -> dict:
    """Score one sequon from a [orders, L, 21] probability array.

    The primary score averages the log odds of asparagine at the first position
    and of serine-or-threonine at the third. The middle residue is excluded
    because any residue except proline permits a sequon, so a preference there
    is not a preference for the motif; proline is recorded separately as the
    residue whose presence would abolish it.

    Validation is unconditional. `computed` is optional only because a caller
    may hold probabilities from elsewhere, but the distribution check always
    runs, so an invalid row cannot be scored by omitting an argument.
    """
    check_scoreable(probabilities, (n_index, plus1_index, plus2_index), computed)

    per_order = []
    for order in range(probabilities.shape[0]):
        p_n = probabilities[order, n_index, AA_INDEX["N"]]
        p_st = (probabilities[order, plus2_index, AA_INDEX["S"]]
                + probabilities[order, plus2_index, AA_INDEX["T"]])
        per_order.append(0.5 * (logit(p_n) + logit(p_st)))

    mean_probs = probabilities.mean(axis=0)
    p_n = float(mean_probs[n_index, AA_INDEX["N"]])
    p_s = float(mean_probs[plus2_index, AA_INDEX["S"]])
    p_t = float(mean_probs[plus2_index, AA_INDEX["T"]])
    p_pro = float(mean_probs[plus1_index, AA_INDEX["P"]])

    return {
        "conditional_sequon_score": float(np.mean(per_order)),
        "conditional_sequon_score_sd": float(np.std(per_order, ddof=1)) if len(per_order) > 1 else 0.0,
        "n_decoding_orders": int(probabilities.shape[0]),
        "p_asn_at_n": p_n,
        "p_ser_at_plus2": p_s,
        "p_thr_at_plus2": p_t,
        "p_ser_or_thr_at_plus2": p_s + p_t,
        "p_pro_at_plus1": p_pro,
        "logit_p_asn": logit(p_n),
        "logit_p_ser_or_thr": logit(p_s + p_t),
        "probs_n": mean_probs[n_index].tolist(),
        "probs_plus1": mean_probs[plus1_index].tolist(),
        "probs_plus2": mean_probs[plus2_index].tolist(),
    }


def decodable_positions(pdb_path: Path, chain_id: str,
                        proteinmpnn_dir: "Path | None" = None) -> np.ndarray:
    """Which residues of a chain ProteinMPNN will decode, without running it.

    The model's backbone mask is computed by its own featuriser from the parsed
    coordinates, so scoreability is a property of the structure alone and can be
    settled before any scoring. That ordering matters for the study design: if
    scoreability were established afterwards, sites would be matched and then
    silently lost, leaving matched sets that no longer balance.

    Returns a boolean array of length L, true where the backbone is complete.
    """
    _prepare_environment()
    _ensure_importable(proteinmpnn_dir)

    from protein_mpnn_utils import StructureDatasetPDB, parse_PDB, tied_featurize

    parsed = parse_PDB(str(pdb_path), input_chain_list=[chain_id])
    protein = StructureDatasetPDB(parsed, truncate=None, max_length=20000)[0]
    out = tied_featurize(
        [protein], "cpu", {protein["name"]: ([chain_id], [])},
        None, None, None, None, None, ca_only=False,
    )
    mask, chain_M_pos = out[2], out[10]
    return ((mask * chain_M_pos)[0].cpu().numpy() > 0)
