"""Three views of the context table, each answering a different question.

A single "clean dataset" cannot serve the atlas, because the reasons a site is
imperfect are not interchangeable. A crystallographer's N->Q knockout, a +1 that
differs between isoform and construct, and a +2 that was never resolved are
three different facts, and an analysis that lumps them either throws away good
Asn measurements or quietly measures the wrong residue.

    triplet_core      every feature in the row describes the sequon it names
    asn_centred       the Asn was measured; +1 and +2 may not have been
    construct_review  everything excluded from triplet_core, kept for inspection

`triplet_core` requires the triplet to match *and* the mapping to be
continuous. The second condition is not redundant: the triplet check compares
residue identities, so a +1 taken from the far side of a numbering gap passes
whenever the landing residue happens to have the same identity. Nine sites in
the first extraction did exactly that.

`asn_centred` is deliberately wider and deliberately restricted in use. It is
valid for features centred on the asparagine -- its RSA, its dihedrals, its ND2
neighbourhood -- and not for +1 or +2 RSA, secondary structure or dihedrals.
"""
from __future__ import annotations

import pandas as pd


def _observed(frame: pd.DataFrame) -> pd.Series:
    return frame.get("triplet_observed", pd.Series("", index=frame.index)).fillna("").astype(str)


def _expected(frame: pd.DataFrame) -> pd.Series:
    return frame.get("triplet_expected", pd.Series("", index=frame.index)).fillna("").astype(str)


def asn_matches(frame: pd.DataFrame) -> pd.Series:
    """Whether the residue actually measured at the site is the expected Asn."""
    observed, expected = _observed(frame), _expected(frame)
    return (observed.str[:1] == expected.str[:1]) & observed.str[:1].ne("")


def is_core(frame: pd.DataFrame) -> pd.Series:
    matches = frame.get("triplet_matches", pd.Series(False, index=frame.index))
    matches = matches.fillna(False).astype(bool)
    continuous = frame.get("mapping_continuous", pd.Series(False, index=frame.index))
    continuous = continuous.fillna(False).astype(bool)
    return matches & continuous


def split_views(frame: pd.DataFrame) -> "dict[str, pd.DataFrame]":
    """The three views, as copies so downstream edits cannot alias each other."""
    core = is_core(frame)
    return {
        "triplet_core": frame[core].copy(),
        "asn_centred": frame[asn_matches(frame)].copy(),
        "construct_review": frame[~core].copy(),
    }
