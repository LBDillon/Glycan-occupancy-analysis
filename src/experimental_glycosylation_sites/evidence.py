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


EXCLUSION_BY_TIER = {
    # A site excluded only because the policy did not select its tier has good
    # evidence that is simply not being counted. Saying "no qualifying evidence"
    # would be false, and the sensitivity workflow in docs/evidence_sources.md
    # depends on this distinction being reported truthfully.
    "manual_experimental": "experimental_tier_not_selected",
    "manual_combinatorial": "combinatorial_tier_not_selected",
    "manual_curator_inference": "curator_inference_only",
    "manual_sequence_model": "sequence_model_only",
    "automatic_sequence_model": "sequence_model_only",
    "sequence_similarity": "sequence_similarity_only",
    "annotation_without_qualifying_evidence": "annotation_without_qualifying_evidence",
}


def join_uniprot_evidence(
    candidates: "pd.DataFrame",
    features: "list[UniProtFeature]",
    missing_accessions: set[str],
    qualifying_tiers: list[str],
) -> "pd.DataFrame":
    """Attach exact-position UniProt evidence to each candidate site.

    Only an N-linked asparagine feature at exactly the candidate position
    counts. Nearby features are never accepted. Every candidate appears in the
    output exactly once, either qualifying or with an exclusion reason.
    """
    import pandas as pd

    exact: dict[tuple[str, int], frozenset[str]] = {}
    for feature in features:
        if feature.parse_status != "ok" or feature.position is None:
            continue
        if feature.glyco_type != "N-linked":
            continue
        key = (feature.accession, feature.position)
        exact[key] = exact.get(key, frozenset()) | feature.evidence_codes

    qualifying = set(qualifying_tiers)
    rows = []
    for accession, position in zip(candidates["accession"], candidates["position"]):
        accession, position = str(accession), int(position)
        codes = exact.get((accession, position))

        if accession in missing_accessions:
            tier, reason, qualifies = "", "accession_absent_from_snapshot", False
        elif codes is None:
            tier, reason, qualifies = "exact_feature_absent", "exact_feature_absent", False
        else:
            tier = classify_uniprot_tier(codes)
            qualifies = tier in qualifying
            reason = "" if qualifies else EXCLUSION_BY_TIER.get(tier, "annotation_without_qualifying_evidence")

        rows.append({
            "accession": accession,
            "position": position,
            "uniprot_tier": tier,
            "uniprot_evidence_codes": "|".join(sorted(codes)) if codes else "",
            "uniprot_qualifies": qualifies,
            "exclusion_reason": reason,
        })

    return pd.DataFrame(rows).sort_values(["accession", "position"]).reset_index(drop=True)


LAYER_ORDER: tuple[str, ...] = ("uniprot", "glygen", "glyconnect", "structure")

OCCUPIED = "occupied_supported"
UNKNOWN = "unknown"
# Reserved. Populating it requires evidence that a site was examined and found
# bare; no current source provides that. Phase 1 must never emit it.
OBSERVED_UNMODIFIED = "observed_unmodified"


def combine_layers(
    uniprot: "pd.DataFrame",
    glygen: "pd.DataFrame | None",
    glyconnect: "pd.DataFrame | None",
    structure: "pd.DataFrame | None",
    policy: dict,
) -> "pd.DataFrame":
    """Merge evidence layers into per-site positivity and occupancy status.

    Layers are independent: any one of them may support a site on its own, and a
    site rejected by UniProt can still be supported by GlyGen or by a structural
    glycan linkage. A site with no supporting layer is `unknown` - never a
    biological negative, because absence of annotation overwhelmingly means
    nobody looked.
    """
    import pandas as pd

    merged = uniprot.copy()
    for frame in (glygen, glyconnect, structure):
        if frame is not None:
            merged = merged.merge(frame, on=["accession", "position"], how="left")

    qualifying_glygen = set(policy.get("qualifying_glygen_tiers", []))
    qualifying_structure = set(policy.get("qualifying_structure_tiers", []))

    support = {
        "uniprot": merged["uniprot_qualifies"].fillna(False).astype(bool),
        "glygen": (
            merged["glygen_tier"].fillna("").isin(qualifying_glygen)
            if "glygen_tier" in merged.columns
            else pd.Series(False, index=merged.index)
        ),
        "glyconnect": (
            merged["glyconnect_supported"].fillna(False).astype(bool)
            if "glyconnect_supported" in merged.columns and policy.get("glyconnect_qualifies", False)
            else pd.Series(False, index=merged.index)
        ),
        "structure": (
            merged["structure_tier"].fillna("").isin(qualifying_structure)
            if "structure_tier" in merged.columns
            else pd.Series(False, index=merged.index)
        ),
    }

    merged["support_sources"] = [
        "|".join(layer for layer in LAYER_ORDER if support[layer].iloc[i])
        for i in range(len(merged))
    ]
    merged["support_count"] = sum(series.astype(int) for series in support.values())
    merged["experimental_positive"] = merged["support_count"] > 0
    # Precedence: supported sites are occupied; otherwise a site whose structure
    # carries an internal control (a glycan modelled at another residue, in a host
    # that can glycosylate) is observed-unmodified; everything else is unknown.
    #
    # This is the ONLY route by which observed_unmodified may be emitted. Without
    # that control a bare asparagine is a silence, not an observation, and calling
    # it a negative would reintroduce the error the whole package guards against.
    merged["occupancy_status"] = UNKNOWN
    if "structure_internal_control" in merged.columns and policy.get(
        "observed_unmodified_from_internal_control", True
    ):
        control = merged["structure_internal_control"].fillna(False).astype(bool)
        merged.loc[control, "occupancy_status"] = OBSERVED_UNMODIFIED
    merged.loc[merged["experimental_positive"], "occupancy_status"] = OCCUPIED
    merged.loc[merged["experimental_positive"], "exclusion_reason"] = ""

    return merged.sort_values(["accession", "position"]).reset_index(drop=True)
