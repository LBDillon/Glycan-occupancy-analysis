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
from .uniprot import accessions_with_xref

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "default.toml"


def _accessions(config, xref_column: str | None = None) -> list[str]:
    """Candidate accessions, optionally restricted to those a source indexes.

    Only 1,714 of 2,878 candidates carry a GlyGen cross-reference; see
    `uniprot.accessions_with_xref` for why the rest are not worth requesting.
    """
    import pandas as pd

    pairs = pd.read_csv(config.paths["pairs_master"], low_memory=False)
    sites = build_candidate_sites(
        pairs, bool(config.policy.get("require_analysis_ready", True))
    )
    accessions = set(sites["accession"])
    if xref_column is None:
        return sorted(accessions)
    return accessions_with_xref(config.paths["uniprot_tsv"], accessions, xref_column)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="experimental_glycosylation_sites")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "command",
        choices=[
            "validate", "fetch-glygen", "fetch-glyconnect", "fetch-structures",
            "fetch-controls", "run",
        ],
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
    if args.command == "fetch-controls":
        import pandas as pd

        from .controls import (
            CONTROL_SETS, build_control_sites, fetch_control_proteins, summarise,
        )

        frames = []
        for name in CONTROL_SETS:
            print(f"fetching {name}...", flush=True)
            frames.append(fetch_control_proteins(name))
            print(f"  {len(frames[-1])} proteins")
        proteins = pd.concat(frames, ignore_index=True)
        sites = build_control_sites(proteins)

        results = Path(config.paths["results_dir"])
        results.mkdir(parents=True, exist_ok=True)
        sites.to_csv(results / "negative_control_sites.csv", index=False)
        (results / "negative_control_summary.json").write_text(
            json.dumps(summarise(sites, proteins), indent=2)
        )
        # Sequences are needed by the feature stage but are bulky; keep them in
        # the cache rather than the results tables.
        cache = Path(config.paths["cache_dir"])
        cache.mkdir(parents=True, exist_ok=True)
        proteins.to_csv(cache / "negative_control_proteins.csv.gz", index=False)

        print(f"\n{len(sites)} control sequons written")
        print(sites.groupby("control_set").size().to_string())
        return 0

    if args.command == "fetch-structures":
        from .structures import fetch_structures, load_manifest

        manifest = load_manifest(
            config.paths["structure_manifest"], config.paths.get("structure_dir")
        )
        wanted: dict[str, set[str]] = {}
        for accession in _accessions(config):
            row = manifest.get(accession)
            if row is None:
                continue
            ids = {x.strip() for x in str(row.get("all_pdb_ids") or "").split(";") if x.strip()}
            ids.discard(str(row.get("pdb_id", "")).strip())
            if ids:
                wanted[accession] = ids

        cap = int(config.api.get("structures_per_accession", 20))
        print(
            f"{len(wanted)} proteins have structures beyond the one already cached; "
            f"fetching up to {cap} each"
        )
        stats = fetch_structures(
            wanted,
            Path(config.paths["cache_dir"]) / "pdb",
            delay=float(config.api.get("delay_seconds", 0.34)),
            timeout=int(config.api.get("timeout_seconds", 60)),
            per_accession_cap=cap,
        )
        print(json.dumps(stats, indent=2))
        return 0

    if args.command == "fetch-glyconnect":
        print("cache:", fetch_glyconnect(_accessions(config, "GlyConnect"), config))
        return 0

    print(json.dumps(run_full(config, fetch=args.fetch), indent=2))
    return 0
