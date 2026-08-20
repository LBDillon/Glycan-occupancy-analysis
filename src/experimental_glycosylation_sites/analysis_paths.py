"""One naming convention for the analysis stages, so a second model is a flag.

Stages 09, 10, 10b and 11 used to hard-code every input and output filename.
That was fine while one model existed and actively dangerous once a second did:
the corrected and per-model score files land under new names, so the stages went
on reading the old ones and reported stale numbers **without failing**. A silent
wrong answer is the worst failure mode this repository has, and this module
exists to remove that particular one.

## The convention

A *variant* is a short tag naming which run produced a set of numbers —
`alphabet_corrected`, `esm_if`, `esmc_single`. Every path is then the base name
with the variant appended:

    results/scores/scores_dataset.csv                    variant ""  (legacy)
    results/scores/scores_dataset_esm_if.csv             variant "esm_if"
    results/analysis/analysis_optimal_esm_if.json        variant "esm_if"

The empty variant reproduces every pre-existing filename exactly, so running any
stage without `--variant` behaves as it always did.

Retention is the one exception. Its ProteinMPNN files predate the convention
(`mpnn_retention_frozen_2026-08-18.csv`), so the empty variant maps to those
explicitly rather than pretending they fit the pattern.
"""
from __future__ import annotations

import argparse
from pathlib import Path

SCORE_SETS = ("dataset", "controls", "secretory")

# Which manifest each retention run was produced from. The manifest supplies
# occupancy_status, so a retention table cannot be interpreted without it.
RETENTION_TAGS = (
    ("scoring_manifest", "results/manifests/scoring_manifest.csv"),
    ("secretory", "results/manifests/manifest_matched_secretory.csv"),
)

# ProteinMPNN's retention files predate the naming convention.
LEGACY_RETENTION = {
    "scoring_manifest": "results/designs/mpnn_retention_frozen_2026-08-18.csv",
    "secretory": "results/designs/mpnn_retention_secretory.csv",
}


def suffix(variant: "str | None") -> str:
    """`_esm_if` for a named variant, empty string for the legacy run."""
    variant = (variant or "").strip().strip("_")
    return f"_{variant}" if variant else ""


def scores(set_tag: str, variant: "str | None" = None) -> Path:
    if set_tag not in SCORE_SETS:
        raise ValueError(f"unknown score set {set_tag!r}; expected one of {SCORE_SETS}")
    return Path(f"results/scores/scores_{set_tag}{suffix(variant)}.csv")


def retention_sources(variant: "str | None" = None) -> "list[tuple[Path, Path]]":
    """`(retention table, manifest)` pairs for this variant."""
    tag_suffix = suffix(variant)
    pairs = []
    for tag, manifest in RETENTION_TAGS:
        if tag_suffix:
            retention = Path(f"results/designs/retention_{tag}{tag_suffix}.csv")
        else:
            retention = Path(LEGACY_RETENTION[tag])
        pairs.append((retention, Path(manifest)))
    return pairs


def retention_all_classes(variant: "str | None" = None) -> Path:
    return Path(f"results/designs/retention_all_classes{suffix(variant)}.csv")


def analysis_json(label: str, variant: "str | None" = None) -> Path:
    return Path(f"results/analysis/analysis_{label}{suffix(variant)}.json")


def contrasts(label: str, variant: "str | None" = None) -> Path:
    return Path(f"results/analysis/contrasts_{label}{suffix(variant)}.csv")


def retention_paired(label: str, variant: "str | None" = None) -> Path:
    return Path(f"results/analysis/retention_paired_{label}{suffix(variant)}.csv")


def retention_paired_json(variant: "str | None" = None) -> Path:
    return Path(f"results/analysis/retention_paired{suffix(variant)}.json")


def significance(variant: "str | None" = None, extension: str = "csv") -> Path:
    return Path(f"results/analysis/significance{suffix(variant)}.{extension}")


def add_variant_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--variant", default="",
        help="which run's numbers to read and write, e.g. alphabet_corrected, "
             "esm_if, esmc_single. Default: the original ProteinMPNN filenames.")


def describe_source(path: Path) -> str:
    """`model / conditioning / n_orders` as recorded in a score file.

    Read from the data rather than restated in the analysis, so a stage can never
    label one model's numbers with another's provenance. Stage 09 previously
    hard-coded "ProteinMPNN v_48_020", which would have mislabelled every ESM-IF
    and ESMC result as ProteinMPNN.
    """
    import pandas as pd

    if not Path(path).exists():
        return "unknown"
    frame = pd.read_csv(path, low_memory=False, nrows=1)
    parts = [str(frame[col].iloc[0]) for col in ("model", "conditioning", "n_orders")
             if col in frame.columns]
    return " / ".join(parts) if parts else "unknown"
