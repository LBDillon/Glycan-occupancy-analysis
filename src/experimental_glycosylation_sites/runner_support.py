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


def build_adapter(name: str, device: str = "cpu"):
    """Construct a registered adapter, passing the device it understands."""
    from . import adapters

    if name == "proteinmpnn":
        return adapters.load(name, device=device,
                             proteinmpnn_dir=Path("../../ProteinMPNN"))
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
