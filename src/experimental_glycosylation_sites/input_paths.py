"""Where the read-only bulk inputs live, resolved by configuration.

This repository holds code, tests and the tables it generates. The inputs it
reads but never writes -- the UniProt release, the cached structures -- are
large and live outside it. Stages used to reach them with literal relative
paths like `../../data/raw/uniprot/...`, which encodes one checkout layout into
the source and breaks the moment the analysis is run from anywhere else.

Locations come from the environment, falling back to the layouts that already
existed so no current invocation changes behaviour:

    GCA_DATA_ROOTS        os.pathsep-separated roots searched in order
    GCA_STRUCTURE_DIRS    os.pathsep-separated structure caches, searched first
    GCA_DATASETS_DIR      directory holding the upstream dataset tables

A missing input raises `MissingInput` naming every location searched, because
the failure this module exists to prevent is a path that silently resolves to
nothing.
"""
from __future__ import annotations

import os
from pathlib import Path

ENV_DATA_ROOTS = "GCA_DATA_ROOTS"
ENV_STRUCTURE_DIRS = "GCA_STRUCTURE_DIRS"
ENV_DATASETS_DIR = "GCA_DATASETS_DIR"

# The in-repo layout first, then the location the analysis directory used when
# it lived inside the SugarFix tree.
DEFAULT_DATA_ROOTS = ("data", "../../data")
DEFAULT_STRUCTURE_DIRS = (
    "data/cache/pdb",
    "../ortholog_sequon_conservation/results/database_current/structures/pdb",
)


class MissingInput(FileNotFoundError):
    """An input file was not found under any configured root."""


def _from_env(name: str) -> "list[Path]":
    raw = os.environ.get(name, "")
    return [Path(part).expanduser() for part in raw.split(os.pathsep) if part.strip()]


def data_roots() -> "list[Path]":
    """Configured roots first, then the historical defaults."""
    return _from_env(ENV_DATA_ROOTS) + [Path(p) for p in DEFAULT_DATA_ROOTS]


def structure_dirs(extra: "tuple[str, ...] | None" = None) -> "list[Path]":
    """Structure caches, most specific first.

    `extra` carries a stage's own --structure-dir flags, which outrank both the
    environment and the defaults.
    """
    return ([Path(p) for p in (extra or ())]
            + _from_env(ENV_STRUCTURE_DIRS)
            + [Path(p) for p in DEFAULT_STRUCTURE_DIRS])


def resolve_input(relative: str) -> Path:
    """The first existing `relative` under a configured root.

    Raises MissingInput naming every location searched, so a stage that cannot
    find its input says where it looked rather than failing obscurely later.
    """
    roots = data_roots()
    for root in roots:
        candidate = root / relative
        if candidate.exists():
            return candidate
    searched = "\n  ".join(str(root / relative) for root in roots)
    raise MissingInput(
        f"{relative} not found. Searched:\n  {searched}\n"
        f"Set {ENV_DATA_ROOTS} to the directory holding it "
        f"({os.pathsep}-separated for more than one)."
    )


def datasets_dir() -> Path:
    """Where the upstream dataset tables are read from.

    Separate from the output directory on purpose: a stage may read the
    evidence tables produced by an earlier analysis tree while writing its own
    results into this repository.
    """
    configured = os.environ.get(ENV_DATASETS_DIR, "").strip()
    return Path(configured).expanduser() if configured else Path("results/datasets")


def resolve_optional_input(relative: str) -> "Path | None":
    """`resolve_input` for a file that legitimately may not exist.

    Callers must still decide what a missing file means. Skipping one silently
    is how the context manifest came out with no sequence context for its
    largest population and no error to say so.
    """
    try:
        return resolve_input(relative)
    except MissingInput:
        return None
