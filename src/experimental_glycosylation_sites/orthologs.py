from __future__ import annotations

import pandas as pd

POSITION_COLUMN = "best_loss_positive_position"
ACCESSION_COLUMN = "pos_accession"

ASSOCIATION_COLUMNS = [
    "accession", "position", "pair_id", "cluster_id",
    "neg_accession", "source", "homology_qc_bucket",
]


def _ready(pairs: pd.DataFrame, require_analysis_ready: bool) -> pd.DataFrame:
    if not require_analysis_ready:
        return pairs
    if "analysis_ready" not in pairs.columns:
        # Falling through to the unfiltered frame here inflates the candidate
        # universe from 4,307 to 5,836 with no indication anything went wrong.
        # A missing column is a broken input, not a reason to skip the filter.
        raise KeyError(
            "analysis_ready filtering was requested but the column is absent from "
            "pairs_master; refusing to silently return the unfiltered frame"
        )
    return pairs[pairs["analysis_ready"].astype("boolean").fillna(False)]


def count_null_positions(pairs: pd.DataFrame, require_analysis_ready: bool) -> int:
    """Rows that pass the readiness filter but carry no positive-side position."""
    frame = _ready(pairs, require_analysis_ready)
    return int(frame[POSITION_COLUMN].isna().sum())


def build_candidate_sites(pairs: pd.DataFrame, require_analysis_ready: bool) -> pd.DataFrame:
    """Deduplicated (accession, position) candidate universe.

    A site occupied on one protein appears in as many pair rows as it has
    ortholog comparisons; dedup here is what stops those rows inflating counts.
    """
    frame = _ready(pairs, require_analysis_ready).dropna(subset=[POSITION_COLUMN])
    sites = pd.DataFrame({
        "accession": frame[ACCESSION_COLUMN].astype(str),
        "position": frame[POSITION_COLUMN].astype(int),
    })
    return (
        sites.drop_duplicates()
        .sort_values(["accession", "position"])
        .reset_index(drop=True)
    )


def build_associations(
    pairs: pd.DataFrame, homology: pd.DataFrame, require_analysis_ready: bool
) -> pd.DataFrame:
    """One row per (site, pair). Kept separate from site tables by design."""
    frame = _ready(pairs, require_analysis_ready).dropna(subset=[POSITION_COLUMN]).copy()
    frame["accession"] = frame[ACCESSION_COLUMN].astype(str)
    frame["position"] = frame[POSITION_COLUMN].astype(int)

    buckets = homology[["pair_id", "homology_qc_bucket"]].drop_duplicates("pair_id")
    merged = frame.merge(buckets, on="pair_id", how="left")

    for column in ASSOCIATION_COLUMNS:
        if column not in merged.columns:
            merged[column] = pd.NA

    return (
        merged[ASSOCIATION_COLUMNS]
        .sort_values(["accession", "position", "pair_id"])
        .reset_index(drop=True)
    )


def assign_subsets(
    associations: pd.DataFrame,
    strict_buckets: list[str],
    plausible_buckets: list[str],
) -> pd.DataFrame:
    """Collapse associations to one row per site with subset membership.

    A site joins a subset if ANY of its associations qualifies. Homology
    decides subset membership only; it is never occupancy evidence.
    """
    frame = associations.copy()
    strict = set(strict_buckets)
    plausible = strict | set(plausible_buckets)
    frame["_strict"] = frame["homology_qc_bucket"].isin(strict)
    frame["_plausible"] = frame["homology_qc_bucket"].isin(plausible)

    grouped = frame.groupby(["accession", "position"], as_index=False).agg(
        in_strict=("_strict", "any"),
        in_strict_plus_plausible=("_plausible", "any"),
        n_associations=("pair_id", "count"),
    )
    return grouped.sort_values(["accession", "position"]).reset_index(drop=True)
