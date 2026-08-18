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

# Sequon subtype is matched exactly rather than as a distance. NXT is occupied
# more often than NXS and the two are chemically distinct, so a comparison that
# pairs an occupied NXS against an unoccupied NXT confounds subtype with
# occupancy. Exact matching removes that at the cost of pairs, which is the
# right trade when the alternative is an uninterpretable difference.
DEFAULT_EXACT = ("subtype",)


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
    exact: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Greedy nearest-neighbour matching without replacement.

    `exact` names columns a control must equal the case on before distance is
    considered at all. Requiring exact agreement is not the same as adding the
    column to `features`: a distance-based match will accept a near neighbour of
    the wrong subtype whenever the structural fit is good enough, which is
    precisely the confound the constraint exists to remove.

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
    usable_case = cases.dropna(subset=list(features) + list(exact))
    usable_control = controls.dropna(subset=list(features) + list(exact)).reset_index(drop=True)
    if usable_case.empty or usable_control.empty:
        return pd.DataFrame(columns=columns)

    pooled = pd.concat([usable_case[list(features)], usable_control[list(features)]])
    scales = {f: float(pooled[f].astype(float).std(ddof=1)) for f in features}

    case_points = _scaled(usable_case, features, scales)
    control_points = _scaled(usable_control, features, scales)

    available = np.ones(len(usable_control), dtype=bool)
    order = np.random.default_rng(seed).permutation(len(usable_case))
    rows = []

    exact_case = {c: usable_case[c].to_numpy() for c in exact}
    exact_control = {c: usable_control[c].to_numpy() for c in exact}

    for index in order:
        distances = np.linalg.norm(control_points - case_points[index], axis=1)
        for column in exact:
            distances = np.where(
                exact_control[column] == exact_case[column][index], distances, np.inf)
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


def match_controls_optimal(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    features: tuple[str, ...] = MATCH_FEATURES,
    caliper: float = DEFAULT_CALIPER,
    exact: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Deterministic 1:1 matching: most pairs first, then least total distance.

    Greedy matching processes cases in a random order and takes the nearest
    available control each time. That is a heuristic, and with a scarce control
    pool it is a consequential one: an early case can take a control that was
    the only admissible partner for a later case, so both the number of pairs
    and which pairs form depend on the seed. The result then carries a
    dependence on an arbitrary choice that no one would defend on its merits.

    This solves the assignment directly. Ineligible pairings — different
    subtype, or further apart than the caliper — are given a cost far larger
    than any achievable total of real distances, so a minimum-cost assignment
    first maximises the number of admissible pairs and only then minimises the
    distance across them. There is no seed and no ordering.

    Ties are possible in principle: two controls can sit at identical distance
    from a case. Rows are therefore sorted by accession and position before the
    assignment, so the tie is broken by a stable property of the data rather
    than by whatever order the caller happened to supply. This makes the output
    reproducible across callers; it is not a claim that the chosen optimum is
    the only one.
    """
    from scipy.optimize import linear_sum_assignment

    columns = [
        "case_accession", "case_position", "control_accession",
        "control_position", "distance", "match_rank",
    ]
    usable_case = (cases.dropna(subset=list(features) + list(exact))
                        .sort_values(["accession", "position"], kind="mergesort")
                        .reset_index(drop=True))
    usable_control = (controls.dropna(subset=list(features) + list(exact))
                              .sort_values(["accession", "position"], kind="mergesort")
                              .reset_index(drop=True))
    if usable_case.empty or usable_control.empty:
        return pd.DataFrame(columns=columns)

    pooled = pd.concat([usable_case[list(features)], usable_control[list(features)]])
    scales = {f: float(pooled[f].astype(float).std(ddof=1)) for f in features}
    case_points = _scaled(usable_case, features, scales)
    control_points = _scaled(usable_control, features, scales)

    # [n_cases, n_controls]
    cost = np.linalg.norm(case_points[:, None, :] - control_points[None, :, :], axis=2)
    eligible = cost <= caliper
    for column in exact:
        eligible &= (usable_case[column].to_numpy()[:, None]
                     == usable_control[column].to_numpy()[None, :])

    # Any total of admissible distances is bounded by caliper * n_pairs, so a
    # penalty of this size can never be worth paying to shorten real distances.
    penalty = 1.0 + caliper * (min(len(usable_case), len(usable_control)) + 1)
    work = np.where(eligible, cost, penalty)

    rows_i, cols_i = linear_sum_assignment(work)
    rows = []
    for i, j in zip(rows_i, cols_i):
        if not eligible[i, j]:
            continue                    # the solver had to place it somewhere
        rows.append({
            "case_accession": usable_case.iloc[i].accession,
            "case_position": int(usable_case.iloc[i].position),
            "control_accession": usable_control.iloc[j].accession,
            "control_position": int(usable_control.iloc[j].position),
            "distance": round(float(cost[i, j]), 4),
            "match_rank": 1,
        })
    frame = pd.DataFrame(rows, columns=columns)
    return frame.sort_values(["case_accession", "case_position"]).reset_index(drop=True)


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


def _weighted_moments(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    total = weights.sum()
    if total <= 0 or len(values) == 0:
        return float("nan"), float("nan")
    mean = float((values * weights).sum() / total)
    var = float((weights * (values - mean) ** 2).sum() / total)
    return mean, var


def weighted_balance_report(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    pairs: pd.DataFrame,
    features: tuple[str, ...] = MATCH_FEATURES,
) -> dict:
    """Balance using matched-set weights.

    Each occupied site carries total weight one, and the controls matched to it
    share weight one between them. Without this a case matched to five controls
    would contribute five times as much to the control mean as a case matched to
    one, so the balance statistic would describe the matching's bookkeeping rather
    than the comparison being made.
    """
    if pairs.empty:
        return {"cases_matched": 0, "controls_used": 0, "features": {}}

    per_case = pairs.groupby(["case_accession", "case_position"]).size()
    weights = pairs.join(
        per_case.rename("k"), on=["case_accession", "case_position"]
    ).assign(weight=lambda f: 1.0 / f.k)

    case_keys = pairs[["case_accession", "case_position"]].drop_duplicates()
    matched_cases = cases.merge(
        case_keys, left_on=["accession", "position"],
        right_on=["case_accession", "case_position"],
    )
    control_side = weights.merge(
        controls, left_on=["control_accession", "control_position"],
        right_on=["accession", "position"],
    )

    report = {
        "cases_matched": int(len(case_keys)),
        "controls_used": int(len(pairs)),
        "control_weight_total": round(float(weights.weight.sum()), 6),
        "features": {},
    }
    for feature in features:
        case_values = matched_cases[feature].astype(float).to_numpy()
        control_values = control_side[feature].astype(float).to_numpy()
        control_weights = control_side.weight.to_numpy()
        case_mean, case_var = _weighted_moments(case_values, np.ones(len(case_values)))
        ctl_mean, ctl_var = _weighted_moments(control_values, control_weights)
        pooled = np.sqrt((case_var + ctl_var) / 2)
        report["features"][feature] = {
            "smd_weighted": round(float((case_mean - ctl_mean) / pooled) if pooled else 0.0, 4),
            "case_mean": round(case_mean, 4),
            "control_weighted_mean": round(ctl_mean, 4),
        }
    return report
