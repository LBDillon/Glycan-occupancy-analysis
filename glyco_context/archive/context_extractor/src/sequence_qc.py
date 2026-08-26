"""Row-level validation of the sequence context behind each site.

Checking a whole population at once -- "does this population have any sequences
at all?" -- catches a cache that failed to load and misses one that loaded
half. These checks are per row, so partial loss is visible, and each failure
names a reason rather than reporting a count.

The sequon rule is N-X-S/T with X != P. A row failing it has not been filtered
badly; it means the position was mapped to the wrong residue somewhere upstream,
because a site that is not a sequon cannot be an N-linked glycosylation site.
"""
from __future__ import annotations

import pandas as pd

KEY = ["accession", "position", "population"]


def _length(value) -> "int | None":
    try:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _reason(row) -> "str | None":
    length = _length(row.get("uniprot_length"))
    if length is None or length <= 0:
        return "no_uniprot_sequence"

    triplet = "" if pd.isna(row.get("sequon_triplet")) else str(row.get("sequon_triplet", ""))
    if len(triplet) != 3:
        return "incomplete_expected_triplet"

    if not (triplet[0] == "N" and triplet[1] != "P" and triplet[2] in ("S", "T")):
        return "expected_triplet_not_a_sequon"

    position = _length(row.get("position"))
    if position is None or position < 1 or position + 2 > length:
        return "sequon_outside_sequence"
    return None


def sequence_context_failures(frame: pd.DataFrame) -> pd.DataFrame:
    """One row per site failing a sequence-context check, with its reason."""
    if not len(frame):
        return pd.DataFrame(columns=KEY + ["reason"])
    reasons = frame.apply(lambda row: _reason(row.to_dict()), axis=1)
    failed = frame[reasons.notna()].copy()
    failed["reason"] = reasons[reasons.notna()]
    columns = [c for c in KEY if c in failed.columns] + ["reason"]
    return failed[columns].reset_index(drop=True)


def add_sequence_distances(frame: pd.DataFrame) -> pd.DataFrame:
    """Where the site sits in the full-length protein, from the manifest.

    Derived here rather than in the structural extractor because they are
    properties of the UniProt sequence: a structure that stops short of the
    C-terminus must not shorten them.

        uniprot_residues_after_asn     residues following the asparagine
        uniprot_residues_after_sequon  residues following the whole N-X-S/T

    A negative value means the sequon runs past the end of the sequence, which
    the row-level checks exclude -- so it is a defect rather than data, and
    raises rather than propagating into the views.
    """
    out = frame.copy()
    length = pd.to_numeric(out.get("uniprot_length"), errors="coerce")
    position = pd.to_numeric(out.get("position"), errors="coerce")
    out["uniprot_residues_after_asn"] = length - position
    out["uniprot_residues_after_sequon"] = length - (position + 2)
    for column in ("uniprot_residues_after_asn", "uniprot_residues_after_sequon"):
        negative = out[column].notna() & (out[column] < 0)
        if negative.any():
            example = out[negative].iloc[0]
            raise ValueError(
                f"{column} is negative for {negative.sum()} row(s), e.g. "
                f"{example.get('accession')}/{example.get('position')} with "
                f"uniprot_length={example.get('uniprot_length')}")
    return out
