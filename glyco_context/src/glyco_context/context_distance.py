"""How far a site's local chemistry sits from natural occupied context.

One number per site, so a wild type and its designs can be compared directly:

    D = median over the panel of | (x - mu_ref) / sigma_ref |

The median rather than the mean, because a single feature with a small reference
spread would otherwise dominate the score. The absolute value because either
direction is a departure.

`mu_ref` and `sigma_ref` are computed with the scored site's **protein** held
out. Holding out only the row would let a protein with several sequons leak into
its own reference, and wild-type sites would look more natural than they are.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def reference_moments(reference: pd.DataFrame, features, exclude_accession=None):
    """Per-feature mean and standard deviation over the held-out reference."""
    frame = reference
    if exclude_accession is not None and "accession" in frame.columns:
        frame = frame[frame.accession != exclude_accession]
    mu, sigma = {}, {}
    for feature in features:
        values = pd.to_numeric(frame.get(feature), errors="coerce").dropna() \
            if feature in frame.columns else pd.Series(dtype=float)
        mu[feature] = float(values.mean()) if len(values) else float("nan")
        sigma[feature] = float(values.std(ddof=1)) if len(values) > 1 else float("nan")
    return mu, sigma


def context_distance(row, mu: dict, sigma: dict, features) -> float:
    """Median absolute standardised departure, or NaN if nothing is measurable.

    A feature whose reference has no spread is skipped rather than treated as an
    infinite departure: it carries no information about typicality.
    """
    departures = []
    for feature in features:
        value = row.get(feature)
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        spread, centre = sigma.get(feature), mu.get(feature)
        if spread is None or centre is None:
            continue
        if not np.isfinite(spread) or spread <= 0 or not np.isfinite(centre):
            continue
        departures.append(abs((float(value) - centre) / spread))
    return float(np.median(departures)) if departures else float("nan")
