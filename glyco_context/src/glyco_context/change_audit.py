"""Attribute every difference between the original and corrected tables.

The purpose is not to summarise the diff but to account for it. Each changed
row is assigned to the named correction that explains it; a row that changed
for a reason we cannot name is reported as UNEXPLAINED rather than absorbed
into the new table, because that is the signal that something else moved.

Categories, in the order they are tested:

    insertion_code_propagation      the site sits in an insertion block and the
                                    Asn was previously read from the wrong one
    invisible_gap_jump              matched the expected triplet before, but the
                                    +1/+2 came from across a gap
    gap_aware_mapping               +1/+2 changed because the walk no longer
                                    crosses unobserved residues
    dssp_multichain_repair          DSSP now runs on a multi-character chain id
    dihedral_no_longer_crosses_gap  phi/psi withdrawn across a gap
    terminal_distance_units         residue counts instead of author numbers
"""
from __future__ import annotations

import pandas as pd

KEY = ["accession", "position", "population"]
DIHEDRALS = [f"{p}_{a}" for p in ("n", "plus1", "plus2") for a in ("phi", "psi")]
TERMINALS = ["distance_to_n_terminus_resolved", "distance_to_c_terminus_resolved",
             "chain_length_resolved"]


def _differs(left: pd.Series, right: pd.Series) -> pd.Series:
    both_null = left.isna() & right.isna()
    return (left.astype(object).ne(right.astype(object)) & ~both_null).fillna(False)


def attribute_changes(old: pd.DataFrame, new: pd.DataFrame,
                      compare: "list[str] | None" = None) -> pd.DataFrame:
    """One row per changed site: its category and the columns that moved."""
    keys = [k for k in KEY if k in old.columns and k in new.columns]
    merged = old.merge(new, on=keys, how="inner", suffixes=("_old", "_new"))
    if not len(merged):
        return pd.DataFrame(columns=keys + ["category", "changed_columns"])

    shared = [c for c in old.columns if c in new.columns and c not in keys]
    if compare is not None:
        shared = [c for c in shared if c in compare]

    changed = {c: _differs(merged[f"{c}_old"], merged[f"{c}_new"]) for c in shared}
    any_change = pd.DataFrame(changed).any(axis=1) if changed else pd.Series(
        False, index=merged.index)

    def column(name, side):
        # pandas only suffixes overlapping names, so a column present in just
        # one table keeps its bare name after the merge.
        for candidate in (f"{name}_{side}", name):
            if candidate in merged:
                return merged[candidate]
        return pd.Series(None, index=merged.index)

    def flag(name, side, default=False):
        return column(name, side).fillna(default).astype(bool)

    icode = column("n_icode", "new").fillna("").astype(str).str.strip().ne("")
    continuous_new = flag("mapping_continuous", "new", True)
    matched_old = flag("triplet_matches", "old")
    triplet_moved = _differs(column("triplet_observed", "old"),
                             column("triplet_observed", "new"))
    dssp_gained = ~flag("dssp_ok", "old") & flag("dssp_ok", "new")
    dihedral_moved = pd.DataFrame(
        {c: changed[c] for c in DIHEDRALS if c in changed}).any(axis=1) \
        if any(c in changed for c in DIHEDRALS) else pd.Series(False, index=merged.index)
    terminal_moved = pd.DataFrame(
        {c: changed[c] for c in TERMINALS if c in changed}).any(axis=1) \
        if any(c in changed for c in TERMINALS) else pd.Series(False, index=merged.index)

    category = pd.Series("UNEXPLAINED", index=merged.index)
    # Later assignments do not override earlier ones: precedence is top down.
    assigned = pd.Series(False, index=merged.index)
    for name, mask in (
        ("insertion_code_propagation", icode & triplet_moved),
        ("invisible_gap_jump", matched_old & ~continuous_new),
        ("gap_aware_mapping", triplet_moved & ~continuous_new),
        ("dssp_multichain_repair", dssp_gained),
        ("dihedral_no_longer_crosses_gap", dihedral_moved),
        ("terminal_distance_units", terminal_moved),
    ):
        take = mask & ~assigned & any_change
        category[take] = name
        assigned |= take

    result = merged.loc[any_change, keys].copy()
    result["category"] = category[any_change]
    result["changed_columns"] = [
        ",".join(c for c in shared if changed[c].iloc[i])
        for i in range(len(merged)) if any_change.iloc[i]]
    return result.reset_index(drop=True)
