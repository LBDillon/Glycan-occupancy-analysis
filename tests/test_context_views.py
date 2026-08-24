"""The three analysis views, and why a site belongs to each.

The strict view exists so that every feature in a row describes the sequon it
claims to. That needs two conditions, not one: the triplet must match *and* the
mapping must be continuous. Nine rows in the first run matched the expected
triplet while the +1 or +2 came from the far side of a numbering gap, so a
triplet check alone readmits exactly the defect the view exists to exclude.
"""
from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.context_views import split_views


def _frame(rows):
    return pd.DataFrame(rows)


BASE = {"accession": "P1", "position": 1, "triplet_expected": "NAS",
        "triplet_observed": "NAS", "triplet_matches": True,
        "mapping_continuous": True}


def test_clean_site_is_in_the_strict_view():
    views = split_views(_frame([BASE]))
    assert len(views["triplet_core"]) == 1
    assert len(views["construct_review"]) == 0


def test_matching_triplet_with_a_discontinuous_mapping_is_excluded():
    """The nine invisible rows: right triplet, wrong residues."""
    row = {**BASE, "mapping_continuous": False}
    views = split_views(_frame([row]))
    assert len(views["triplet_core"]) == 0, "a gap-jumped +1 is not a measured +1"
    assert len(views["construct_review"]) == 1


def test_asn_only_view_keeps_sites_where_just_the_asn_matches():
    """N->Q knockouts are excluded; a differing +1 or +2 is not."""
    row = {**BASE, "triplet_observed": "NKS", "triplet_matches": False}
    views = split_views(_frame([row]))
    assert len(views["triplet_core"]) == 0
    assert len(views["asn_centred"]) == 1
    assert len(views["construct_review"]) == 1


def test_substituted_asn_is_not_in_the_asn_centred_view():
    row = {**BASE, "triplet_observed": "QAS", "triplet_matches": False}
    views = split_views(_frame([row]))
    assert len(views["asn_centred"]) == 0
    assert len(views["construct_review"]) == 1


def test_unresolved_positions_do_not_count_as_a_matching_asn():
    row = {**BASE, "triplet_observed": "N??", "triplet_matches": False,
           "mapping_continuous": False}
    views = split_views(_frame([row]))
    assert len(views["triplet_core"]) == 0
    assert len(views["asn_centred"]) == 1, "the Asn itself was still measured"


def test_every_site_lands_in_exactly_one_of_core_or_review():
    rows = [BASE,
            {**BASE, "mapping_continuous": False},
            {**BASE, "triplet_observed": "QAS", "triplet_matches": False}]
    views = split_views(_frame(rows))
    assert len(views["triplet_core"]) + len(views["construct_review"]) == len(rows)
