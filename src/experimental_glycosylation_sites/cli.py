from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import load_config
from .glyconnect import fetch_details as fetch_glyconnect
from .glygen import fetch_details as fetch_glygen
from .orthologs import build_candidate_sites
from .pipeline import run_full

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "default.toml"


def _accessions(config) -> list[str]:
    import pandas as pd

    pairs = pd.read_csv(config.paths["pairs_master"], low_memory=False)
    sites = build_candidate_sites(
        pairs, bool(config.policy.get("require_analysis_ready", True))
    )
    return sorted(set(sites["accession"]))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="experimental_glycosylation_sites")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "command",
        choices=["validate", "fetch-glygen", "fetch-glyconnect", "run"],
    )
    parser.add_argument(
        "--fetch", action="store_true",
        help="refresh API caches before running (run command only)",
    )
    args = parser.parse_args(argv)

    config = load_config(args.config)
    errors = config.validate_inputs()
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    if args.command == "validate":
        print("All configured inputs resolve.")
        return 0
    if args.command == "fetch-glygen":
        print("cache:", fetch_glygen(_accessions(config), config))
        return 0
    if args.command == "fetch-glyconnect":
        print("cache:", fetch_glyconnect(_accessions(config), config))
        return 0

    print(json.dumps(run_full(config, fetch=args.fetch), indent=2))
    return 0
