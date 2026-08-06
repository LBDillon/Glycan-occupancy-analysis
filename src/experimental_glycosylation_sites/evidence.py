from __future__ import annotations

# Evidence and Conclusion Ontology tier names.
#
# NOTE: ECO:0007744 is "combinatorial computational and experimental evidence
# used in manual assertion". It does NOT mean the site was observed in a
# structure. Do not describe it as structural or PDB evidence.
ECO_TIERS: dict[str, str] = {
    "ECO:0000269": "manual_experimental",
    "ECO:0007744": "manual_combinatorial",
    "ECO:0000305": "manual_curator_inference",
    "ECO:0000250": "sequence_similarity",
    "ECO:0000255": "manual_sequence_model",
    "ECO:0000256": "automatic_sequence_model",
    "ECO:0000259": "automatic_sequence_model",
}

UNIPROT_TIER_ORDER: tuple[str, ...] = (
    "manual_experimental",
    "manual_combinatorial",
    "manual_curator_inference",
    "manual_sequence_model",
    "sequence_similarity",
    "automatic_sequence_model",
    "annotation_without_qualifying_evidence",
)

NO_QUALIFYING_EVIDENCE = "annotation_without_qualifying_evidence"


def classify_uniprot_tier(codes: frozenset[str]) -> str:
    """Return the strongest tier implied by a feature's ECO codes."""
    tiers = {ECO_TIERS[code] for code in codes if code in ECO_TIERS}
    for tier in UNIPROT_TIER_ORDER:
        if tier in tiers:
            return tier
    return NO_QUALIFYING_EVIDENCE
