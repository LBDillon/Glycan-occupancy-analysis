from __future__ import annotations

import hashlib
import json
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

OUTPUT_PATH_KEYS = {"cache_dir", "results_dir"}

# Declared for provenance and possible future use, but read by no module today.
# Validating them would hard-fail every command — and silently skip the frozen
# regression test — on files nothing actually consumes.
UNUSED_PATH_KEYS = {"proteins_master", "existing_structural_context"}


@dataclass(frozen=True)
class Config:
    paths: dict[str, Path]
    layers: dict[str, bool]
    policy: dict
    api: dict
    source_path: Path
    config_hash: str
    raw: dict = field(repr=False, default_factory=dict)

    def validate_inputs(self) -> list[str]:
        """Return one message per configured input path that does not exist."""
        errors = []
        for key, path in sorted(self.paths.items()):
            if key in OUTPUT_PATH_KEYS or key in UNUSED_PATH_KEYS:
                continue
            if not path.exists():
                errors.append(f"config key '{key}' points at a missing path: {path}")
        return errors


def load_config(path: Path) -> Config:
    path = Path(path).resolve()
    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    base = path.parent
    paths = {
        key: (base / value).resolve()
        for key, value in raw.get("paths", {}).items()
    }
    payload = json.dumps(raw, sort_keys=True).encode("utf-8")
    return Config(
        paths=paths,
        layers=raw.get("layers", {}),
        policy=raw.get("policy", {}),
        api=raw.get("api", {}),
        source_path=path,
        config_hash=hashlib.sha256(payload).hexdigest(),
        raw=raw,
    )
