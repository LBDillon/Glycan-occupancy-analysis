from __future__ import annotations

from pathlib import Path

import pandas as pd

SORT_KEYS = ["accession", "position"]


def write_table(frame: pd.DataFrame, path: Path) -> None:
    """Write a CSV with deterministic row order so runs diff cleanly."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = [key for key in SORT_KEYS if key in frame.columns]
    ordered = frame.sort_values(keys).reset_index(drop=True) if keys else frame
    ordered.to_csv(path, index=False)
