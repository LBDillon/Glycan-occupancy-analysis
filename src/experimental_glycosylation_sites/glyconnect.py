from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from .config import Config

_SITE_LABEL = re.compile(r"^Asn-(\d+)$", re.IGNORECASE)


def parse_site_label(label: str | None) -> int | None:
    """Parse a GlyConnect site string such as "Asn-80". Non-asparagine returns None."""
    if not label:
        return None
    match = _SITE_LABEL.match(str(label).strip())
    return int(match.group(1)) if match else None


def extract_positions(detail: dict) -> dict[int, dict]:
    """Collapse GlyConnect structures to one record per N-linked asparagine position."""
    found: dict[int, dict] = {}
    for structure in detail.get("structures", []) or []:
        if str(structure.get("type", "")).lower() != "n-linked":
            continue
        composition = str(structure.get("composition_string") or "")
        for site in structure.get("sites", []) or []:
            position = parse_site_label(site.get("site"))
            if position is None:
                continue
            record = found.setdefault(position, {"glycan_count": 0, "compositions": set()})
            record["glycan_count"] += 1
            if composition:
                record["compositions"].add(composition)
    return {
        position: {
            "glycan_count": record["glycan_count"],
            "compositions": "|".join(sorted(record["compositions"])),
        }
        for position, record in found.items()
    }


def _get_json(url: str, timeout: int, retries: int, delay: float):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return json.load(response), ""
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            time.sleep(delay * (attempt + 1))
    return None, error


def fetch_details(accessions: list[str], config: Config) -> Path:
    """Two-step GlyConnect retrieval into a resumable JSONL cache."""
    cache_path = Path(config.paths["cache_dir"]) / "glyconnect_protein_detail.jsonl"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    seen = set(load_cache(cache_path))

    api = config.api
    search = api.get("glyconnect_search_url", "https://glyconnect.expasy.org/api/proteins?uniprot={accession}")
    detail_url = api.get("glyconnect_detail_url", "https://glyconnect.expasy.org/api/proteins/{protein_id}")
    delay = float(api.get("delay_seconds", 0.4))
    timeout = int(api.get("timeout_seconds", 60))
    retries = int(api.get("max_retries", 3))

    with cache_path.open("a", encoding="utf-8") as handle:
        for accession in accessions:
            if accession in seen:
                continue
            hits, error = _get_json(search.format(accession=accession), timeout, retries, delay)
            detail = None
            if hits:
                protein_id = hits[0].get("id")
                if protein_id is not None:
                    time.sleep(delay)
                    detail, error = _get_json(
                        detail_url.format(protein_id=protein_id), timeout, retries, delay
                    )
            handle.write(json.dumps(
                {"accession": accession, "detail": detail, "error": error}
            ) + "\n")
            handle.flush()
            time.sleep(delay)
    return cache_path


def load_cache(path: Path) -> dict[str, dict]:
    path = Path(path)
    if not path.exists():
        return {}
    cache: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("detail"):
                cache[record["accession"]] = record["detail"]
    return cache


def build_site_evidence(candidates: pd.DataFrame, cache: dict[str, dict]) -> pd.DataFrame:
    """One row per candidate site. Corroboration only; never a primary label."""
    per_accession = {
        accession: extract_positions(detail) for accession, detail in cache.items()
    }
    rows = []
    for accession, position in zip(candidates["accession"], candidates["position"]):
        accession, position = str(accession), int(position)
        record = per_accession.get(accession, {}).get(position)
        rows.append({
            "accession": accession,
            "position": position,
            "glyconnect_supported": record is not None,
            "glyconnect_glycan_count": record["glycan_count"] if record else 0,
            "glyconnect_compositions": record["compositions"] if record else "",
        })
    return pd.DataFrame(rows).sort_values(["accession", "position"]).reset_index(drop=True)
