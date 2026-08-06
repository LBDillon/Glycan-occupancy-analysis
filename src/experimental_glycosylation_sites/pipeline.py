from __future__ import annotations

import pandas as pd

from .config import Config
from .evidence import join_uniprot_evidence
from .orthologs import (
    assign_subsets,
    build_associations,
    build_candidate_sites,
    count_null_positions,
)
from .uniprot import load_uniprot_features


def run_uniprot_baseline(config: Config) -> dict:
    """UniProt-only evidence pass. This is the frozen regression baseline."""
    policy = config.policy
    require_ready = bool(policy.get("require_analysis_ready", True))

    pairs = pd.read_csv(config.paths["pairs_master"], low_memory=False)
    homology = pd.read_csv(config.paths["homology_qc"], low_memory=False)

    candidates = build_candidate_sites(pairs, require_ready)
    associations = build_associations(pairs, homology, require_ready)
    subsets = assign_subsets(
        associations,
        list(policy.get("strict_buckets", [])),
        list(policy.get("plausible_buckets", [])),
    )

    features, missing = load_uniprot_features(
        config.paths["uniprot_tsv"], set(candidates["accession"])
    )
    evidence = join_uniprot_evidence(
        candidates, features, missing, list(policy.get("qualifying_uniprot_tiers", []))
    )

    merged = evidence.merge(subsets, on=["accession", "position"], how="left")
    qualifying = merged[merged["uniprot_qualifies"]]
    strict = qualifying[qualifying["in_strict"].fillna(False)]
    plausible = qualifying[qualifying["in_strict_plus_plausible"].fillna(False)]

    def counts(frame: pd.DataFrame) -> tuple[int, int]:
        return len(frame), frame["accession"].nunique()

    all_sites, all_proteins = counts(qualifying)
    sp_sites, sp_proteins = counts(plausible)
    s_sites, s_proteins = counts(strict)

    return {
        "candidates": candidates,
        "associations": associations,
        "evidence": merged,
        "all": qualifying,
        "strict_plus_plausible": plausible,
        "strict": strict,
        "excluded": merged[~merged["uniprot_qualifies"]],
        "counts": {
            "candidate_sites": len(candidates),
            "excluded_sites": int((~merged["uniprot_qualifies"]).sum()),
            "null_position_rows": count_null_positions(pairs, require_ready),
            "all_sites": all_sites,
            "all_proteins": all_proteins,
            "strict_plus_plausible_sites": sp_sites,
            "strict_plus_plausible_proteins": sp_proteins,
            "strict_sites": s_sites,
            "strict_proteins": s_proteins,
        },
    }
