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

from .input_paths import structure_dirs

# Searched in order; the first directory holding a given PDB id wins.
STRUCTURE_DIRS = (
    "data/cache/pdb",
    "../ortholog_sequon_conservation/results/database_current/structures/pdb",
)


def structure_paths(extra_dirs: "tuple[str, ...] | None" = None) -> "dict[str, Path]":
    """PDB id (upper case) -> cached structure file."""
    paths: "dict[str, Path]" = {}
    for base in structure_dirs(tuple(extra_dirs or ())):
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
    parser.add_argument("--mask-mode", default=None,
                        help="ESMC only: 'single' (default) or 'joint'")
    parser.add_argument("--max-batch", type=int, default=None,
                        help="cap designs decoded at once. Default: chosen from "
                             "chain length so memory stays bounded.")
    parser.add_argument("--shard", default=None, metavar="K/N",
                        help="process only chain group K of N (0-based), for "
                             "SLURM job arrays. Each shard must write its own "
                             "output file; merge them afterwards.")
    return parser.parse_args(argv)


def apply_shard(groups: list, shard: "str | None") -> list:
    """Keep chain groups belonging to this shard.

    Sharding is by chain group rather than by site so that a chain's expensive
    per-chain work happens exactly once, in exactly one task. Interleaving
    (`index % n == k`) rather than contiguous blocks, because chains are ordered
    by PDB id and length correlates with neither -- contiguous blocks would give
    one task all the long chains.
    """
    if not shard:
        return groups
    try:
        k, n = (int(part) for part in str(shard).split("/"))
    except ValueError as exc:
        raise SystemExit(f"--shard must look like K/N, got {shard!r}") from exc
    if not 0 <= k < n:
        raise SystemExit(f"--shard K must satisfy 0 <= K < N, got {shard!r}")
    selected = [g for i, g in enumerate(groups) if i % n == k]
    print(f"shard {k}/{n}: {len(selected)} of {len(groups)} chain groups", flush=True)
    return selected


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


def build_adapter(name: str, device: str = "cpu", **options):
    """Construct a registered adapter, passing what it understands.

    Model-specific options travel as keywords rather than as attributes the
    runners know about, so a new model's knobs never require editing a stage.
    Options that are None are dropped, so an unset flag means "the adapter's
    default" rather than an explicit None the adapter has to interpret.
    """
    from . import adapters

    options = {k: v for k, v in options.items() if v is not None}
    if name == "proteinmpnn":
        return adapters.load(name, device=device,
                             proteinmpnn_dir=proteinmpnn_dir(), **options)
    return adapters.load(name, device=device, **options)


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
