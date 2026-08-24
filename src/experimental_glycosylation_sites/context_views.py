"""Three views of the context table, each answering a different question.

A single "clean dataset" cannot serve the atlas, because the reasons a site is
imperfect are not interchangeable. A crystallographer's N->Q knockout, a +1 that
differs between isoform and construct, and a +2 that was never resolved are
three different facts, and an analysis that lumps them either throws away good
Asn measurements or quietly measures the wrong residue.

    triplet_core      every feature in the row describes the sequon it names
    asn_core          the Asn was measured; +1 and +2 may not have been
    construct_review  everything excluded from triplet_core, kept for inspection

`triplet_core` requires the triplet to match *and* the mapping to be
continuous. The second condition is not redundant: the triplet check compares
residue identities, so a +1 taken from the far side of a numbering gap passes
whenever the landing residue happens to have the same identity. Nine sites in
the first extraction did exactly that.

`asn_core` is deliberately wider and deliberately restricted in use. It is
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


def _flag(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.get(column, pd.Series(False, index=frame.index)).fillna(False).astype(bool)


def is_core(frame: pd.DataFrame) -> pd.Series:
    """Triplet agrees, all three residues located, and the mapping continuous."""
    return (_flag(frame, "triplet_matches")
            & _flag(frame, "triplet_resolved")
            & _flag(frame, "mapping_continuous"))


def exclusion_reason(frame: pd.DataFrame) -> pd.Series:
    """Why each site is not in `triplet_core`, empty string where it is.

    The reasons are ordered by what a reader would act on first: a substituted
    asparagine is a construct question, an unresolved position is a density
    question, and a discontinuous mapping with a matching triplet is the silent
    case that motivated the continuity requirement.
    """
    observed, expected = _observed(frame), _expected(frame)
    reason = pd.Series("", index=frame.index)
    unset = ~is_core(frame)

    no_expectation = unset & expected.str.len().ne(3)
    asn_substituted = unset & ~no_expectation & ~asn_matches(frame)
    unresolved = unset & ~no_expectation & ~asn_substituted & observed.str.contains(r"\?")
    discontinuous = (unset & ~no_expectation & ~asn_substituted & ~unresolved
                     & ~_flag(frame, "mapping_continuous"))
    substitution = unset & ~(no_expectation | asn_substituted | unresolved | discontinuous)

    reason[no_expectation] = "no_expected_triplet"
    reason[asn_substituted] = "asn_substituted"
    reason[unresolved] = "unresolved_position"
    reason[discontinuous] = "discontinuous_mapping"
    reason[substitution] = "sequence_substitution"
    return reason


def split_views(frame: pd.DataFrame) -> "dict[str, pd.DataFrame]":
    """The three views, as copies so downstream edits cannot alias each other."""
    core = is_core(frame)
    return {
        "triplet_core": frame[core].copy(),
        "asn_core": frame[asn_matches(frame)].copy(),
        "construct_review": frame[~core].copy(),
    }
