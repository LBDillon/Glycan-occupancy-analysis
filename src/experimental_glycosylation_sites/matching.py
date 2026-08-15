"""Match occupied sites to unoccupied ones on structural context.

Occupied sites are more solvent-exposed than any of the negative sets — median
relative accessibility 0.43 against 0.31-0.33 — because an
oligosaccharyltransferase has to reach the residue. Structure-based models can
see exposure. So an unmatched comparison of model scores would mostly restate
that occupied sites are exposed, which is already known and is not a claim about
whether a model understands glycosylation.

Matching removes that. Each occupied site is paired with unoccupied sites of
comparable local environment, so a residual difference in score cannot be
attributed to accessibility or packing.

The balance report is the part to read first. Matching that leaves a large
standardised mean difference has not worked, and a comparison built on it is not
interpretable however clean the downstream statistics look.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Local environment only. Anything describing the protein's identity, taxonomy or
# compartment is deliberately excluded: those differ between the groups by
# construction, and matching them away would remove the comparison itself.
MATCH_FEATURES = ("rsa", "n_neighbours_8a", "hydrophobic_fraction_8a")

# Maximum distance, in pooled standard deviations across all matching features,
# between a case and an acceptable control. Pairs further apart than this are not
# matched at all rather than matched badly — an unmatched case is visible in the
# report, whereas a bad match silently weakens the comparison.
DEFAULT_CALIPER = 0.25


def standardised_mean_difference(case: pd.Series, control: pd.Series) -> float:
    """Group difference in pooled standard deviations.

    The conventional diagnostic for matching. Below about 0.1 is usually taken as
    adequate balance; the sign says which group sits higher.
    """
    case, control = case.dropna(), control.dropna()
    if len(case) < 2 or len(control) < 2:
        return float("nan")
    pooled = np.sqrt((case.var(ddof=1) + control.var(ddof=1)) / 2)
    if pooled == 0:
        return 0.0
    return float((case.mean() - control.mean()) / pooled)


def _scaled(frame: pd.DataFrame, features: tuple[str, ...], scales: dict) -> np.ndarray:
    return np.column_stack([
        frame[f].astype(float).to_numpy() / (scales[f] or 1.0) for f in features
    ])


def match_controls(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    features: tuple[str, ...] = MATCH_FEATURES,
    k: int = 5,
    caliper: float = DEFAULT_CALIPER,
    seed: int = 0,
) -> pd.DataFrame:
    """Greedy nearest-neighbour matching without replacement.

    Without replacement, because reusing one convenient control across many cases
    would make the comparison look better powered than it is. Cases are processed
    in a seeded random order so no systematic advantage goes to whichever site
    happens to sort first.

    Returns one row per matched pair, long format. Cases that find no control
    inside the caliper simply produce no rows, and the count is reported.
    """
    columns = [
        "case_accession", "case_position", "control_accession",
        "control_position", "distance", "match_rank",
    ]
    usable_case = cases.dropna(subset=list(features))
    usable_control = controls.dropna(subset=list(features)).reset_index(drop=True)
    if usable_case.empty or usable_control.empty:
        return pd.DataFrame(columns=columns)

    pooled = pd.concat([usable_case[list(features)], usable_control[list(features)]])
    scales = {f: float(pooled[f].astype(float).std(ddof=1)) for f in features}

    case_points = _scaled(usable_case, features, scales)
    control_points = _scaled(usable_control, features, scales)

    available = np.ones(len(usable_control), dtype=bool)
    order = np.random.default_rng(seed).permutation(len(usable_case))
    rows = []

    for index in order:
        distances = np.linalg.norm(control_points - case_points[index], axis=1)
        distances[~available] = np.inf
        nearest = np.argsort(distances)[:k]
        case_row = usable_case.iloc[index]
        for rank, position in enumerate(nearest, start=1):
            if not np.isfinite(distances[position]) or distances[position] > caliper:
                break
            control_row = usable_control.iloc[position]
            available[position] = False
            rows.append({
                "case_accession": case_row.accession,
                "case_position": int(case_row.position),
                "control_accession": control_row.accession,
                "control_position": int(control_row.position),
                "distance": round(float(distances[position]), 4),
                "match_rank": rank,
            })

    # An empty result is a legitimate outcome — every case fell outside the
    # caliper — so it must still carry the schema rather than an empty frame with
    # no columns, which would fail downstream instead of reporting zero matches.
    return pd.DataFrame(rows, columns=columns)


def balance_report(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    pairs: pd.DataFrame,
    features: tuple[str, ...] = MATCH_FEATURES,
) -> dict:
    """Standardised mean differences before and after matching, per feature."""
    matched_cases = cases.merge(
        pairs[["case_accession", "case_position"]].drop_duplicates(),
        left_on=["accession", "position"],
        right_on=["case_accession", "case_position"],
    )
    matched_controls = controls.merge(
        pairs[["control_accession", "control_position"]].drop_duplicates(),
        left_on=["accession", "position"],
        right_on=["control_accession", "control_position"],
    )

    # Sites, not proteins. Counting accessions here while counting sites below
    # would make the two fail to reconcile whenever a protein carries several
    # sites — which most of them do.
    matched_keys = (
        set(zip(pairs.case_accession, pairs.case_position)) if len(pairs) else set()
    )
    report = {
        "cases_total": int(len(cases)),
        "cases_matched": len(matched_keys),
        "cases_unmatched": int(len(cases) - len(matched_keys)),
        "case_proteins_matched": int(pairs.case_accession.nunique() if len(pairs) else 0),
        "controls_used": int(len(matched_controls)),
        "controls_available": int(len(controls)),
        "mean_pairs_per_case": (
            round(len(pairs) / max(len(matched_cases), 1), 2) if len(pairs) else 0.0
        ),
        "features": {},
    }
    for feature in features:
        report["features"][feature] = {
            "smd_before": round(standardised_mean_difference(cases[feature], controls[feature]), 4),
            "smd_after": round(
                standardised_mean_difference(matched_cases[feature], matched_controls[feature]), 4
            ),
            "case_median": round(float(cases[feature].median()), 4),
            "control_median_before": round(float(controls[feature].median()), 4),
            "control_median_after": (
                round(float(matched_controls[feature].median()), 4)
                if len(matched_controls) else None
            ),
        }
    return report
