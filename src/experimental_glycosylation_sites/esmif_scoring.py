"""Conditional ESM-IF probabilities at sequon positions.

The counterpart to `mpnn_scoring`, answering the same site-level question with a
second inverse-folding model so that a null is not a property of one model's
inductive bias. The output columns are deliberately identical, so every stage
downstream of scoring is untouched.

## What ESM-IF's conditional actually is, and how it differs

ProteinMPNN's `conditional_probs` gives

    P(aa at i | backbone, ALL other native residues)

averaged over a panel of sampled decoding orders. ESM-IF is an autoregressive
encoder-decoder trained strictly left to right, so the only conditional it can
honestly produce is the prefix one

    P(aa at i | backbone, native residues 1..i-1)

read off a single teacher-forced pass. Two consequences, both of which belong in
any write-up that puts the two models side by side.

**The conditioning is asymmetric where ProteinMPNN's is not.** Scoring the +2
residue, both models see the native asparagine upstream; only ProteinMPNN also
sees the sequence C-terminal to the site. So the two numbers answer neighbouring
questions, not the same one, and their raw magnitudes are not comparable. What
*is* comparable is the within-model matched-pair contrast, which is what the
analysis actually rests on.

**There is no decoding-order distribution.** A teacher-forced pass is
deterministic, so `conditional_sequon_score_sd` is structurally 0 and
`n_decoding_orders` is 1. Those columns are retained rather than dropped so the
score tables from the two models share a schema; they are not evidence that
ESM-IF is more precise than ProteinMPNN.

## Why indices cannot be taken on trust

The manifest's `model_index` is an ordinal into the chain as `structures._parse_chains`
reads it — Biopython, keeping residues that satisfy `is_aa(standard=False)` and
carry a CA. ESM-IF reads structures through biotite, whose residue set is not
guaranteed to be the same one. On this corpus the two agree for the large
majority of chains and disagree for a handful, which is precisely the dangerous
case: a silent off-by-some scores the wrong residue and still returns a
plausible-looking number.

So the mapping is established explicitly, and two defects are handled:

- biotite raises `KeyError` converting a non-standard residue name (PCA and
  friends) to one letter. The residue count is unaffected — `coords` is built by
  `get_atom_coords_residuewise`, independently of the sequence string — so
  mapping unknown names to `X` recovers the chain *without shifting any index*.
- where the two sequences still differ, the module's own aligner decides the
  correspondence, and a site is scored only if its three residues map to
  positions whose identities reproduce the manifest's triplet.

Anything that fails is reported unscoreable, never scored — the same invariant
`mpnn_scoring` enforces against ProteinMPNN's all-zero rows.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .mpnn_scoring import (EPSILON, PROBABILITY_SUM_TOLERANCE,  # noqa: F401
                           IncompleteBackboneError, InvalidProbabilityVector,
                           _prepare_environment, logit)

DEFAULT_MODEL = "esm_if1_gvp4_t16_142M_UR50"

# A teacher-forced pass is deterministic: one "order", zero spread. Named
# constants so the score dict cannot drift from what the docstring claims.
N_ORDERS = 1
SCORE_SD = 0.0

CONDITIONING = "autoregressive_prefix"


class ChainUnreadableError(ValueError):
    """ESM-IF could not load or align this chain."""


class SequonMismatchError(ValueError):
    """The mapped residues do not reproduce the manifest's triplet."""


_PATCHED = False


def patch_biotite() -> None:
    """Make fair-esm 2.0.0 work with biotite >= 1.0, without shifting indices.

    Two independent breakages. `filter_backbone` was renamed
    `filter_peptide_backbone`, which stops `esm.inverse_folding` importing at
    all. And `ProteinSequence.convert_letter_3to1` raises on residue names it
    does not know, which aborts an entire chain over one modified residue; since
    the coordinate array is built separately from the sequence string, returning
    `X` keeps every residue in place and costs only that residue's identity.
    """
    global _PATCHED
    if _PATCHED:
        return

    import biotite.structure as bs

    if not hasattr(bs, "filter_backbone") and hasattr(bs, "filter_peptide_backbone"):
        bs.filter_backbone = bs.filter_peptide_backbone  # type: ignore[attr-defined]

    from biotite.sequence import ProteinSequence

    original = ProteinSequence.convert_letter_3to1

    def tolerant(symbol):
        try:
            return original(symbol)
        except Exception:
            return "X"

    ProteinSequence.convert_letter_3to1 = staticmethod(tolerant)
    _PATCHED = True


def load_model(device: str = "cpu", model_name: str = DEFAULT_MODEL):
    """Load ESM-IF1 in eval mode. Returns `(model, alphabet)`."""
    _prepare_environment()
    patch_biotite()
    if model_name != DEFAULT_MODEL:
        raise ValueError(f"unsupported ESM-IF weights: {model_name!r}")

    import esm.pretrained

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    return model.eval().to(device), alphabet


class ChainMapping:
    """The correspondence between manifest indices and ESM-IF's residue array.

    `to_esm` maps a manifest `model_index` to an ESM-IF index. `finite` marks
    ESM-IF residues whose N, CA and C are all present — the model's own notion of
    a usable backbone, and the analogue of ProteinMPNN's mask.
    """

    __slots__ = ("to_esm", "esm_seq", "finite", "coords", "manifest_length")

    def __init__(self, to_esm, esm_seq, finite, coords, manifest_length):
        self.to_esm = to_esm
        self.esm_seq = esm_seq
        self.finite = finite
        self.coords = coords
        self.manifest_length = manifest_length

    def map_indices(self, indices) -> "list[int]":
        """Manifest indices -> ESM-IF indices, or raise."""
        mapped = [self.to_esm.get(int(i)) for i in indices]
        missing = [int(i) for i, m in zip(indices, mapped) if m is None]
        if missing:
            raise IncompleteBackboneError(
                f"manifest {'index' if len(missing) == 1 else 'indices'} {missing} "
                "have no ESM-IF counterpart")
        rejected = [m for m in mapped if not bool(self.finite[m])]
        if rejected:
            raise IncompleteBackboneError(
                f"ESM-IF backbone incomplete (missing N, CA or C) at "
                f"{'index' if len(rejected) == 1 else 'indices'} {rejected}")
        return mapped

    def check_triplet(self, mapped, expected: str) -> None:
        observed = "".join(self.esm_seq[m] for m in mapped)
        if observed != expected:
            raise SequonMismatchError(
                f"ESM-IF reads {observed!r} where the manifest records {expected!r}")


def chain_mapping(structure_path: Path, chain_id: str, pdb_id: "str | None" = None) -> ChainMapping:
    """Build the manifest-index -> ESM-IF-index correspondence for one chain.

    Needs no model pass, so scoreability stays answerable before matching.
    """
    patch_biotite()

    from esm.inverse_folding.util import (extract_coords_from_structure,
                                          load_structure)

    from .structures import _alignment_pairs, _parse_chains

    path = Path(structure_path)
    identifier = pdb_id or path.stem

    try:
        structure = load_structure(str(path), str(chain_id))
        coords, esm_seq = extract_coords_from_structure(structure)
    except Exception as exc:
        raise ChainUnreadableError(f"{type(exc).__name__}: {str(exc)[:120]}") from exc

    chains = _parse_chains(path, str(identifier))
    native = next((c for c in chains if c.chain_id == str(chain_id)), None)
    if native is None:
        raise ChainUnreadableError(f"chain {chain_id!r} absent from the Biopython parse")

    if native.sequence == esm_seq:
        # The common case: the two parsers agree residue for residue.
        to_esm = {i: i for i in range(len(esm_seq))}
    else:
        # They do not. Let the module's own aligner decide, and convert its
        # 1-based pairs to the 0-based indices everything else here speaks.
        to_esm = {a - 1: b - 1 for a, b in _alignment_pairs(native.sequence, esm_seq)}

    finite = np.all(np.isfinite(coords), axis=(-1, -2))
    return ChainMapping(to_esm, esm_seq, finite, coords, len(native.sequence))


def decodable_positions(structure_path: Path, chain_id: str,
                        pdb_id: "str | None" = None) -> np.ndarray:
    """Which manifest indices ESM-IF can evaluate, in the MANIFEST index space.

    Returned in manifest space rather than ESM-IF's because that is what
    `05_scoreability.py` indexes with, and because scoreability has to be
    comparable across models. A chain ESM-IF cannot read is entirely
    unscoreable rather than an exception, so one bad structure does not stop a
    sweep.
    """
    try:
        mapping = chain_mapping(structure_path, chain_id, pdb_id)
    except ChainUnreadableError:
        from .structures import _parse_chains

        try:
            chains = _parse_chains(Path(structure_path), pdb_id or Path(structure_path).stem)
            native = next((c for c in chains if c.chain_id == str(chain_id)), None)
            length = len(native.sequence) if native else 0
        except Exception:
            length = 0
        return np.zeros(length, dtype=bool)

    decodable = np.zeros(mapping.manifest_length, dtype=bool)
    for manifest_index, esm_index in mapping.to_esm.items():
        if 0 <= manifest_index < mapping.manifest_length and bool(mapping.finite[esm_index]):
            decodable[manifest_index] = True
    return decodable


def conditional_probabilities(mapping: ChainMapping, model, alphabet,
                              device: str = "cpu") -> np.ndarray:
    """Teacher-forced per-residue distributions for one chain, shape [L, vocab].

    One forward pass over the whole chain. Column i of the decoder output
    predicts the token at ESM-IF residue i, conditioned on the backbone and on
    native residues 0..i-1 — so the residue being scored is never among its own
    inputs, and no label can leak through the sequence.
    """
    _prepare_environment()
    import torch
    import torch.nn.functional as F

    from esm.inverse_folding.util import CoordBatchConverter

    converter = CoordBatchConverter(alphabet)
    batch_coords, confidence, _, tokens, padding_mask = converter(
        [(mapping.coords, None, mapping.esm_seq)], device=device)

    with torch.no_grad():
        logits, _ = model.forward(batch_coords, padding_mask, confidence, tokens[:, :-1])

    # logits are [B, vocab, L]; transpose to [L, vocab] and normalise in float32
    # so a half-precision run cannot produce rows that fail the sum check.
    return F.softmax(logits[0].transpose(0, 1).float(), dim=-1).cpu().numpy()


def check_scoreable(probabilities: np.ndarray, indices) -> None:
    """Raise unless all three rows are genuine probability distributions.

    ESM-IF has no all-zero-row failure mode of ProteinMPNN's kind, but the check
    is kept because it is cheap and because it is the invariant the interface
    promises: a row that is not a distribution is never scored.
    """
    for index in indices:
        row = probabilities[index]
        total = float(row.sum())
        if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at ESM-IF index {index} sums to {total:.6g}, not 1")
        if float(row.min()) < 0.0 or float(row.max()) > 1.0 + PROBABILITY_SUM_TOLERANCE:
            raise InvalidProbabilityVector(
                f"row at ESM-IF index {index} has values outside [0, 1]: "
                f"min {float(row.min()):.6g}, max {float(row.max()):.6g}")


def sequon_score(probabilities: np.ndarray, alphabet, n_index: int,
                 plus1_index: int, plus2_index: int) -> dict:
    """Score one sequon, with the column names `mpnn_scoring.sequon_score` uses.

    The same statistic: the mean of the log odds of asparagine at the first
    position and of serine-or-threonine at the third, with the middle residue
    excluded because any residue but proline permits a sequon.
    """
    check_scoreable(probabilities, (n_index, plus1_index, plus2_index))

    index_of = alphabet.get_idx
    p_n = float(probabilities[n_index, index_of("N")])
    p_s = float(probabilities[plus2_index, index_of("S")])
    p_t = float(probabilities[plus2_index, index_of("T")])
    p_pro = float(probabilities[plus1_index, index_of("P")])

    score = 0.5 * (logit(p_n) + logit(p_s + p_t))

    return {
        "conditional_sequon_score": score,
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
# Generation, for the retention outcome.
# --------------------------------------------------------------------------

# ProteinMPNN's design pass omits the unknown-residue token, so ESM-IF's is
# restricted to the same twenty letters. Without this the decoder can emit a
# special token at a sequon position, which `classify_retention` would score as
# a lost motif — a statement about the vocabulary, not about the model.
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def design_sequences(mapping: ChainMapping, model, alphabet, n_designs: int,
                     temperature: float, device: str = "cpu",
                     seed: int = 0) -> "list[str]":
    """`n_designs` unconstrained sequences for one chain, in MANIFEST index space.

    All `n_designs` are decoded as one batch. ESM-IF is autoregressive, so a
    chain of length L costs L sequential decoder steps however many sequences
    are wanted; sampling them one at a time pays that latency `n_designs` times
    over for no benefit. Batching is what makes the retention sweep affordable.

    Sequences are returned indexed as the manifest indexes the chain, with `X`
    wherever ESM-IF has no counterpart residue, so a caller can read
    `design[n_model_index]` exactly as it does for ProteinMPNN. Those `X`
    positions are unscoreable by construction — `decodable_positions` excludes
    them — so no scored site ever reads one.

    Nothing is fixed or biased: the sequon is free to disappear.
    """
    _prepare_environment()
    import torch
    import torch.nn.functional as F

    from esm.inverse_folding.util import CoordBatchConverter

    dictionary = model.decoder.dictionary
    converter = CoordBatchConverter(dictionary)
    # One entry per design, all the same backbone: the batch dimension IS the
    # sample dimension, so the encoder runs once and every decoder step advances
    # all designs together.
    batch_coords, confidence, _, _, padding_mask = converter(
        [(mapping.coords, None, None)] * n_designs, device=device)

    length = len(mapping.coords)
    mask_idx = dictionary.get_idx("<mask>")
    tokens = torch.full((n_designs, 1 + length), mask_idx,
                        dtype=torch.long, device=device)
    tokens[:, 0] = dictionary.get_idx("<cath>")

    allowed = torch.tensor([dictionary.get_idx(aa) for aa in STANDARD_AA],
                           dtype=torch.long, device=device)

    generator = torch.Generator(device="cpu").manual_seed(seed)
    incremental_state = {}
    with torch.no_grad():
        encoder_out = model.encoder(batch_coords, padding_mask, confidence)
        for step in range(1, length + 1):
            logits, _ = model.decoder(tokens[:, :step], encoder_out,
                                      incremental_state=incremental_state)
            # [B, vocab, 1] -> [B, vocab] for the single new position
            step_logits = logits[:, :, -1].float()
            probabilities = F.softmax(step_logits[:, allowed] / temperature, dim=-1)
            # Sampling is driven from a CPU generator so a run is reproducible
            # from its seed whether it executes on CPU or GPU.
            choice = torch.multinomial(probabilities.cpu(), 1, generator=generator)
            tokens[:, step] = allowed[choice.squeeze(-1).to(device)]

    esm_designs = [
        "".join(dictionary.get_tok(t) for t in row)
        for row in tokens[:, 1:].cpu().tolist()
    ]

    to_manifest = {esm: manifest for manifest, esm in mapping.to_esm.items()}
    rebuilt = []
    for design in esm_designs:
        chars = ["X"] * mapping.manifest_length
        for esm_index, residue in enumerate(design):
            manifest_index = to_manifest.get(esm_index)
            if manifest_index is not None and 0 <= manifest_index < mapping.manifest_length:
                chars[manifest_index] = residue
        rebuilt.append("".join(chars))
    return rebuilt
