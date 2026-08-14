from __future__ import annotations

import datetime as dt
import hashlib
import subprocess
from pathlib import Path

from .config import Config, OUTPUT_PATH_KEYS


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_state() -> dict:
    def run(*args: str) -> str:
        try:
            return subprocess.run(
                args, capture_output=True, text=True, timeout=10, check=False
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""

    return {
        "commit": run("git", "rev-parse", "HEAD"),
        "dirty": bool(run("git", "status", "--porcelain")),
    }


def build_manifest(config: Config, counts: dict, extra: dict) -> dict:
    """Everything needed to explain how a result set was produced."""
    inputs = {}
    for key, path in sorted(config.paths.items()):
        if key in OUTPUT_PATH_KEYS or not path.is_file():
            continue
        stat = path.stat()
        inputs[key] = {
            "path": str(path),
            "sha256": hash_file(path),
            "size_bytes": stat.st_size,
            "modified": dt.datetime.fromtimestamp(stat.st_mtime, dt.UTC).isoformat(),
        }

    return {
        "generated_at": dt.datetime.now(dt.UTC).isoformat(),
        "inputs": inputs,
        "config": {
            "path": str(config.source_path),
            "hash": config.config_hash,
            "layers": config.layers,
            "policy": config.policy,
        },
        "git": _git_state(),
        "counts": counts,
        "extra": extra,
    }
