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


def _accessions(config, xref_column: str | None = None) -> list[str]:
    """Candidate accessions, optionally restricted to those a source indexes.

    GlyGen answers HTTP 500 — not 404 — for an accession it has no entry for, so
    requesting every candidate spends most of the run on responses that can never
    succeed. Only 1,714 of 2,878 candidates carry a GlyGen cross-reference.
    """
    import csv
    import gzip

    import pandas as pd

    pairs = pd.read_csv(config.paths["pairs_master"], low_memory=False)
    sites = build_candidate_sites(
        pairs, bool(config.policy.get("require_analysis_ready", True))
    )
    accessions = set(sites["accession"])
    if xref_column is None:
        return sorted(accessions)

    tsv_path = config.paths["uniprot_tsv"]
    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    indexed = set()
    with opener(tsv_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("Entry") or "").strip()
            if accession in accessions and (row.get(xref_column) or "").strip():
                indexed.add(accession)
    return sorted(indexed)


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
        print("cache:", fetch_glygen(_accessions(config, "GlyGen"), config))
        return 0
    if args.command == "fetch-glyconnect":
        print("cache:", fetch_glyconnect(_accessions(config, "GlyConnect"), config))
        return 0

    print(json.dumps(run_full(config, fetch=args.fetch), indent=2))
    return 0
