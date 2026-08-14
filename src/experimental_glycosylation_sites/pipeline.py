from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .config import Config
from .evidence import join_uniprot_evidence
from .orthologs import (
    assign_subsets,
    build_associations,
    build_candidate_sites,
    count_null_positions,
)
from .uniprot import accessions_with_xref, load_uniprot_features

from . import glyconnect as glyconnect_layer
from . import glygen as glygen_layer
from . import structures as structure_layer
from .evidence import combine_layers
from .provenance import build_manifest
from .table_io import write_table


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


def _uniprot_sequences(tsv_path, accessions: set[str]) -> dict[str, str]:
    import csv
    import gzip

    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    sequences = {}
    with opener(tsv_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("Entry") or "").strip()
            if accession in accessions:
                sequences[accession] = (row.get("Sequence") or "").strip()
    return sequences


def _feature_inventory(config: Config, accessions: set[str]) -> pd.DataFrame:
    """Every parsed CARBOHYD feature, including ones no candidate site uses.

    Provenance and intermediate evidence, not a headline result: it records what
    the snapshot actually contained, so exclusions can be audited against the
    source rather than taken on trust.
    """
    from .evidence import classify_uniprot_tier
    from .uniprot import load_uniprot_features

    features, _ = load_uniprot_features(config.paths["uniprot_tsv"], accessions)
    return pd.DataFrame([
        {
            "accession": feature.accession,
            "position": feature.position,
            "glyco_type": feature.glyco_type,
            "parse_status": feature.parse_status,
            "evidence_codes": "|".join(sorted(feature.evidence_codes)),
            "tier": classify_uniprot_tier(feature.evidence_codes),
            "note": feature.raw_note,
        }
        for feature in features
    ])


def run_full(config: Config, fetch: bool = False) -> dict:
    """Run every enabled evidence layer and write all result tables."""
    baseline = run_uniprot_baseline(config)
    candidates = baseline["candidates"]
    accessions = sorted(set(candidates["accession"]))
    results = Path(config.paths["results_dir"])
    cache_dir = Path(config.paths["cache_dir"])
    extra: dict = {}

    def to_fetch(xref_column: str) -> list[str]:
        """Accessions worth requesting from a source that indexes `xref_column`.

        Narrows the fetch list only. Evidence below is still built from the full
        candidate set, so filtering here cannot change any published count.
        """
        return accessions_with_xref(
            config.paths["uniprot_tsv"], set(accessions), xref_column
        )

    glygen_frame = None
    if config.layers.get("glygen"):
        if fetch:
            glygen_layer.fetch_details(to_fetch("GlyGen"), config)
        cache = glygen_layer.load_cache(cache_dir / "glygen_protein_detail.jsonl")
        extra["glygen_cached_accessions"] = len(cache)
        glygen_frame = glygen_layer.build_site_evidence(candidates, cache)
        write_table(glygen_frame, results / "glygen_site_evidence.csv")

    glyconnect_frame = None
    if config.layers.get("glyconnect"):
        if fetch:
            glyconnect_layer.fetch_details(to_fetch("GlyConnect"), config)
        cache = glyconnect_layer.load_cache(cache_dir / "glyconnect_protein_detail.jsonl")
        extra["glyconnect_cached_accessions"] = len(cache)
        glyconnect_frame = glyconnect_layer.build_site_evidence(candidates, cache)
        write_table(glyconnect_frame, results / "glyconnect_site_evidence.csv")

    structure_frame = None
    if config.layers.get("structure"):
        manifest = structure_layer.load_manifest(
            config.paths["structure_manifest"], config.paths.get("structure_dir")
        )
        sequences = _uniprot_sequences(config.paths["uniprot_tsv"], set(accessions))
        extra["accessions_with_cached_structure"] = sum(
            1 for accession in accessions if accession in manifest
        )
        structure_frame = structure_layer.build_site_evidence(candidates, sequences, manifest)
        write_table(structure_frame, results / "structure_site_evidence.csv")

    combined = combine_layers(
        baseline["evidence"].drop(columns=["in_strict", "in_strict_plus_plausible", "n_associations"],
                                  errors="ignore"),
        glygen_frame, glyconnect_frame, structure_frame, config.policy,
    )
    subsets = baseline["evidence"][
        ["accession", "position", "in_strict", "in_strict_plus_plausible", "n_associations"]
    ]
    combined = combined.merge(subsets, on=["accession", "position"], how="left")

    positive = combined[combined["experimental_positive"]]
    strict = positive[positive["in_strict"].fillna(False)]
    plausible = positive[positive["in_strict_plus_plausible"].fillna(False)]
    excluded = combined[~combined["experimental_positive"]]
    curator = combined[combined["uniprot_tier"] == "manual_curator_inference"]

    write_table(_feature_inventory(config, set(accessions)), results / "uniprot_exact_n_linked_sites.csv")
    write_table(candidates, results / "candidate_sites.csv")
    write_table(baseline["associations"], results / "site_pair_associations.csv")
    write_table(baseline["all"], results / "experimental_sites_uniprot_baseline.csv")
    write_table(positive, results / "experimental_sites_all.csv")
    write_table(plausible, results / "experimental_sites_strict_plus_plausible.csv")
    write_table(strict, results / "experimental_sites_strict.csv")
    write_table(excluded, results / "excluded_sites.csv")
    if config.policy.get("curator_inferred_sensitivity", True):
        write_table(curator, results / "curator_inferred_sensitivity_sites.csv")

    counts = dict(baseline["counts"])
    counts.update({
        "enriched_all_sites": len(positive),
        "enriched_all_proteins": positive["accession"].nunique(),
        "enriched_strict_plus_plausible_sites": len(plausible),
        "enriched_strict_plus_plausible_proteins": plausible["accession"].nunique(),
        "enriched_strict_sites": len(strict),
        "enriched_strict_proteins": strict["accession"].nunique(),
        "excluded_sites": len(excluded),
        "support_uniprot": int(combined["support_sources"].str.contains("uniprot").sum()),
        "support_glygen": int(combined["support_sources"].str.contains("glygen").sum()),
        "support_structure": int(combined["support_sources"].str.contains("structure").sum()),
    })

    results.mkdir(parents=True, exist_ok=True)
    (results / "summary.json").write_text(json.dumps(counts, indent=2))
    (results / "provenance.json").write_text(
        json.dumps(build_manifest(config, counts, extra), indent=2)
    )
    return counts
