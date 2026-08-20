"""Shared plumbing for the model-dependent pipeline stages.

Stages 05, 07 and 08 all have to answer the same three questions — which model,
on which device, and where is this PDB id cached — and they used to answer the
third by repeating the same glob in each file. Centralising it means a new
structure directory is added once rather than three times, and it gives the
`--model` flag one definition instead of three.

Defaults are chosen so that every command in the README behaves exactly as it
did before this module existed: `--model proteinmpnn`, `--device cpu`.
"""
from __future__ import annotations

import argparse
from pathlib import Path

# Searched in order; the first directory holding a given PDB id wins.
STRUCTURE_DIRS = (
    "data/cache/pdb",
    "../ortholog_sequon_conservation/results/database_current/structures/pdb",
)


def structure_paths(extra_dirs: "tuple[str, ...] | None" = None) -> "dict[str, Path]":
    """PDB id (upper case) -> cached structure file."""
    paths: "dict[str, Path]" = {}
    for directory in tuple(extra_dirs or ()) + STRUCTURE_DIRS:
        base = Path(directory)
        if not base.is_dir():
            continue
        for path in list(base.glob("*.pdb")) + list(base.glob("*.cif")):
            paths.setdefault(path.stem.upper(), path)
    return paths


def parse_args(argv, default_manifest: str, default_out: str, *,
               description: str = "") -> argparse.Namespace:
    """Positional manifest and output, plus the model/device flags.

    Positional-with-defaults keeps every existing invocation working unchanged
    while letting a second model be selected without editing the file.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("manifest", nargs="?", default=default_manifest)
    parser.add_argument("out", nargs="?", default=default_out)
    parser.add_argument("--model", default="proteinmpnn",
                        help="registered adapter name (default: proteinmpnn)")
    parser.add_argument("--device", default="cpu",
                        help="torch device, e.g. cpu, cuda, mps (default: cpu)")
    parser.add_argument("--structure-dir", action="append", default=[],
                        help="extra directory of cached structures; repeatable")
    return parser.parse_args(argv)


# Where ProteinMPNN's checkout might be. `../../ProteinMPNN` is correct inside
# SugarFix, where this module sits two levels below the repo root -- but wrong
# anywhere the module is the root, which is exactly the standalone/Colab layout.
# Hardcoding it meant the ProteinMPNN runs failed off-laptop while ESM-IF ran
# fine, so probe instead, and let the environment override.
PROTEINMPNN_CANDIDATES = ("../../ProteinMPNN", "/content/ProteinMPNN",
                          "../ProteinMPNN", "ProteinMPNN", "~/ProteinMPNN")


def proteinmpnn_dir() -> Path:
    """Locate ProteinMPNN, honouring $PROTEINMPNN_DIR first.

    Identified by `protein_mpnn_utils.py` rather than by the directory merely
    existing, so an empty or half-cloned checkout is rejected here with a clear
    message instead of failing later inside an import.
    """
    import os

    candidates = []
    override = os.environ.get("PROTEINMPNN_DIR")
    if override:
        candidates.append(override)
    candidates.extend(PROTEINMPNN_CANDIDATES)

    for candidate in candidates:
        path = Path(candidate).expanduser()
        if (path / "protein_mpnn_utils.py").exists():
            return path

    raise FileNotFoundError(
        "ProteinMPNN not found. Looked for protein_mpnn_utils.py in: "
        + ", ".join(str(Path(c).expanduser()) for c in candidates)
        + ". Clone it (https://github.com/dauparas/ProteinMPNN) or set "
          "PROTEINMPNN_DIR to its checkout.")


def build_adapter(name: str, device: str = "cpu"):
    """Construct a registered adapter, passing the device it understands."""
    from . import adapters

    if name == "proteinmpnn":
        return adapters.load(name, device=device,
                             proteinmpnn_dir=proteinmpnn_dir())
    return adapters.load(name, device=device)


def resolve_device(requested: str) -> str:
    """Fall back to CPU rather than dying when the requested device is absent.

    A Colab session that loses its accelerator should finish the sweep slowly,
    not lose the run: every stage here is resumable, so a downgraded device
    costs time and nothing else.
    """
    import torch

    if requested in ("auto", None):
        if torch.cuda.is_available():
            return "cuda"
        return "cpu"
    if requested.startswith("cuda") and not torch.cuda.is_available():
        print(f"requested {requested} but CUDA is unavailable; using cpu", flush=True)
        return "cpu"
    if requested == "mps" and not torch.backends.mps.is_available():
        print("requested mps but it is unavailable; using cpu", flush=True)
        return "cpu"
    return requested
