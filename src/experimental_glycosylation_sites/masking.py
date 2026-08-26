"""How much of a model's occupancy preference survives hiding the motif.

Both arms of every comparison contain the sequon, so motif recognition alone
cannot separate them. Hiding the whole motif and rescoring asks what the
surroundings say on their own. The quantity is a difference of differences:

    change = contrast(motif visible) - contrast(motif hidden)

per matched pair, where each contrast is already occupied minus control. A large
positive change means the preference depended on seeing the motif in context; a
change near zero means it did not.

Only pairs both schemes scored can contribute, and the resample unit is the same
connected group the rest of the analysis uses -- ortholog clusters linked through
shared control proteins.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEY = ["case_accession", "case_position"]


def masking_change(visible: pd.DataFrame, hidden: pd.DataFrame,
                   n_boot: int = 20000, seed: int = 0, alpha: float = 0.05) -> dict:
    """Paired change in contrast between the two masking schemes."""
    joined = visible[KEY + ["contrast", "resample_unit"]].merge(
        hidden[KEY + ["contrast"]], on=KEY, suffixes=("_visible", "_hidden"))
    out = {"n_pairs": int(len(joined)),
           "n_units": int(joined.resample_unit.nunique()) if len(joined) else 0,
           "mean_visible": float(joined.contrast_visible.mean()) if len(joined) else np.nan,
           "mean_hidden": float(joined.contrast_hidden.mean()) if len(joined) else np.nan,
           "mean": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    if not len(joined):
        return out

    change = joined.contrast_visible - joined.contrast_hidden
    out["mean"] = float(change.mean())
    groups = [g.to_numpy(float) for _, g in change.groupby(joined.resample_unit)]
    rng = np.random.default_rng(seed)
    draws = np.array([np.concatenate([groups[i] for i in
                                      rng.integers(0, len(groups), len(groups))]).mean()
                      for _ in range(n_boot)])
    out["ci_low"] = float(np.percentile(draws, 100 * alpha / 2))
    out["ci_high"] = float(np.percentile(draws, 100 * (1 - alpha / 2)))
    below, above = float((draws <= 0).mean()), float((draws >= 0).mean())
    out["p"] = float(min(1.0, max(2 * min(below, above), 1.0 / n_boot)))
    return out
