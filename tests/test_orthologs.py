from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.orthologs import (
    assign_subsets,
    build_associations,
    build_candidate_sites,
    count_null_positions,
)


def pairs_frame() -> pd.DataFrame:
    return pd.DataFrame([
        # same site repeated across three pair rows
        {"pair_id": "P1:X", "pos_accession": "P1", "best_loss_positive_position": 10.0,
         "analysis_ready": True, "cluster_id": "c1", "neg_accession": "X", "source": "direct"},
        {"pair_id": "P1:Y", "pos_accession": "P1", "best_loss_positive_position": 10.0,
         "analysis_ready": True, "cluster_id": "c1", "neg_accession": "Y", "source": "direct"},
        {"pair_id": "P1:Z", "pos_accession": "P1", "best_loss_positive_position": 10.0,
         "analysis_ready": True, "cluster_id": "c1", "neg_accession": "Z", "source": "two_hop"},
        # a second distinct site on the same protein
        {"pair_id": "P1b:X", "pos_accession": "P1", "best_loss_positive_position": 55.0,
         "analysis_ready": True, "cluster_id": "c1", "neg_accession": "X", "source": "direct"},
        # excluded: not analysis ready
        {"pair_id": "P2:X", "pos_accession": "P2", "best_loss_positive_position": 20.0,
         "analysis_ready": False, "cluster_id": "c2", "neg_accession": "X", "source": "direct"},
        # excluded: null position
        {"pair_id": "P3:X", "pos_accession": "P3", "best_loss_positive_position": None,
         "analysis_ready": True, "cluster_id": "c3", "neg_accession": "X", "source": "direct"},
    ])


def homology_frame() -> pd.DataFrame:
    return pd.DataFrame([
        {"pair_id": "P1:X", "homology_qc_bucket": "weak_or_confounded"},
        {"pair_id": "P1:Y", "homology_qc_bucket": "plausible_ortholog_like"},
        {"pair_id": "P1:Z", "homology_qc_bucket": "strict_ortholog_like"},
        {"pair_id": "P1b:X", "homology_qc_bucket": "family_level_homolog"},
    ])


def test_repeated_pair_rows_collapse_to_one_site():
    sites = build_candidate_sites(pairs_frame(), require_analysis_ready=True)
    assert len(sites) == 2
    assert list(zip(sites.accession, sites.position)) == [("P1", 10), ("P1", 55)]


def test_position_is_integer_not_float():
    sites = build_candidate_sites(pairs_frame(), require_analysis_ready=True)
    assert sites.position.dtype.kind == "i"


def test_analysis_ready_filter_excludes_rows():
    sites = build_candidate_sites(pairs_frame(), require_analysis_ready=True)
    assert "P2" not in set(sites.accession)
    relaxed = build_candidate_sites(pairs_frame(), require_analysis_ready=False)
    assert "P2" in set(relaxed.accession)


def test_null_positions_are_dropped_and_counted():
    assert count_null_positions(pairs_frame(), require_analysis_ready=True) == 1
    sites = build_candidate_sites(pairs_frame(), require_analysis_ready=True)
    assert "P3" not in set(sites.accession)


def test_associations_preserve_every_pair_row():
    assoc = build_associations(pairs_frame(), homology_frame(), require_analysis_ready=True)
    site_rows = assoc[(assoc.accession == "P1") & (assoc.position == 10)]
    assert len(site_rows) == 3
    assert set(site_rows.pair_id) == {"P1:X", "P1:Y", "P1:Z"}


def test_subset_membership_uses_best_bucket_across_associations():
    assoc = build_associations(pairs_frame(), homology_frame(), require_analysis_ready=True)
    subsets = assign_subsets(assoc, ["strict_ortholog_like"], ["plausible_ortholog_like"])
    row = subsets[(subsets.accession == "P1") & (subsets.position == 10)].iloc[0]
    assert bool(row.in_strict) is True
    assert bool(row.in_strict_plus_plausible) is True
    assert int(row.n_associations) == 3


def test_site_with_only_weak_buckets_is_in_neither_subset():
    assoc = build_associations(pairs_frame(), homology_frame(), require_analysis_ready=True)
    subsets = assign_subsets(assoc, ["strict_ortholog_like"], ["plausible_ortholog_like"])
    row = subsets[(subsets.accession == "P1") & (subsets.position == 55)].iloc[0]
    assert bool(row.in_strict) is False
    assert bool(row.in_strict_plus_plausible) is False


def test_subsets_have_one_row_per_site_not_per_pair():
    assoc = build_associations(pairs_frame(), homology_frame(), require_analysis_ready=True)
    subsets = assign_subsets(assoc, ["strict_ortholog_like"], ["plausible_ortholog_like"])
    assert len(subsets) == 2
