"""Invariants the context table must satisfy.

Every check here corresponds to a defect that reached a finished table once. A
report that only counted rows and coverage passed all of them, which is why
they are asserted rather than described.
"""
from __future__ import annotations

import pandas as pd

KEY = ("accession", "position")


def _violation(name: str, description: str, offending: pd.DataFrame) -> dict:
    example = {}
    if len(offending):
        row = offending.iloc[0]
        example = {k: row[k] for k in KEY if k in offending.columns}
    return {"invariant": name, "description": description,
            "rows": int(len(offending)), "example": example}


def check_invariants(frame: pd.DataFrame) -> "list[dict]":
    """Every invariant the table breaks, with a count and one example each."""
    violations, checks = [], []

    def has(*columns):
        return all(c in frame.columns for c in columns)

    if has("distance_to_n_terminus_resolved", "distance_to_c_terminus_resolved",
           "chain_length_resolved"):
        total = (frame.distance_to_n_terminus_resolved
                 + frame.distance_to_c_terminus_resolved + 1)
        checks.append((
            "terminal_distances_sum_to_chain_length",
            "d_N + d_C + 1 must equal the resolved chain length; author-number "
            "arithmetic against a residue count does not satisfy it",
            total.ne(frame.chain_length_resolved) & total.notna()))

    if has("has_nd2", "nd2_atoms_8a_same_chain"):
        nd2 = frame.has_nd2.fillna(False).astype(bool)
        present = frame.nd2_atoms_8a_same_chain.notna()
        checks.append((
            "nd2_counts_present_iff_nd2",
            "ND2-centred counts must be present exactly when ND2 is resolved; "
            "a missing side chain must not read as an uncrowded site",
            nd2 != present))

    for prefix in ("n", "plus1", "plus2"):
        if has(f"{prefix}_dssp_ok", f"{prefix}_ss"):
            ok = frame[f"{prefix}_dssp_ok"].fillna(False).astype(bool)
            checks.append((
                "ss_requires_dssp_ok",
                "a secondary-structure call requires a DSSP entry at that position",
                frame[f"{prefix}_ss"].notna() & ~ok))

    if has("triplet_observed"):
        observed = frame.triplet_observed.fillna("").astype(str)
        checks.append((
            "triplet_has_three_positions",
            "the observed triplet always has three positions, '?' where unresolved",
            observed.str.len().ne(3)))

    if has("loop_run_length"):
        length = frame.loop_run_length
        checks.append((
            "loop_run_length_is_positive",
            "a loop run containing the Asn has at least one residue",
            length.notna() & (length < 1)))

    for column in [c for c in frame.columns if c.endswith("_fraction_8a")]:
        values = frame[column]
        checks.append((
            "fractions_within_unit_interval",
            "neighbour fractions are proportions and must lie in [0, 1]",
            values.notna() & ((values < 0) | (values > 1))))

    for name, description, mask in checks:
        mask = mask.fillna(False)
        if mask.any():
            violations.append(_violation(name, description, frame[mask]))
    return violations
