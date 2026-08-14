from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from .config import Config

GLYGEN_TIER_ORDER: tuple[str, ...] = (
    "glygen_reported_with_glycan",
    "glygen_reported_independent",
    "glygen_reported_uniprot_derived",
    "glygen_predicted",
    "glygen_unclassified",
)

# GlyGen "predicted" entries are UniProt's sequence-model annotation re-exported.
# Counting them would readmit predictions into a set whose value is that it
# contains none. "reported" citing only UniProtKB is circular for the same
# reason and is kept separate so it never corroborates the UniProt call.
_UNIPROT_DATABASES = {"UniProtKB", "UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"}


def _databases(entry: dict) -> set[str]:
    return {
        str(item.get("database") or "")
        for item in entry.get("evidence", []) or []
        if item.get("database")
    }


def classify_glygen_entry(entry: dict) -> str:
    categories = set((entry.get("site_category_dict") or {}).keys())
    if not categories and entry.get("site_category"):
        categories = {entry["site_category"]}
    databases = _databases(entry)

    if "reported_with_glycan" in categories:
        return "glygen_reported_with_glycan"
    if "reported" in categories:
        independent = databases - _UNIPROT_DATABASES
        return "glygen_reported_independent" if independent else "glygen_reported_uniprot_derived"
    if "predicted" in categories or "predicted_with_glycan" in categories:
        return "glygen_predicted"
    return "glygen_unclassified"


def _ids(entry: dict, database: str) -> list[str]:
    return [
        str(item.get("id"))
        for item in entry.get("evidence", []) or []
        if item.get("database") == database and item.get("id")
    ]


def extract_sites(detail: dict) -> dict[int, dict]:
    """Collapse a GlyGen protein detail payload to one record per N-linked position."""
    rank = {tier: index for index, tier in enumerate(GLYGEN_TIER_ORDER)}
    sites: dict[int, dict] = {}

    for entry in detail.get("glycosylation", []) or []:
        if str(entry.get("type", "")).lower() != "n-linked":
            continue
        start, end = entry.get("start_pos"), entry.get("end_pos")
        if start is None or start != end:
            continue

        position = int(start)
        tier = classify_glygen_entry(entry)
        record = sites.setdefault(position, {
            "tier": tier, "categories": set(), "evidence_databases": set(),
            "pubmed_ids": set(), "glytoucan_ids": set(),
        })
        if rank[tier] < rank[record["tier"]]:
            record["tier"] = tier
        record["categories"].update((entry.get("site_category_dict") or {}).keys())
        record["evidence_databases"].update(_databases(entry))
        record["pubmed_ids"].update(_ids(entry, "PubMed"))
        if entry.get("glytoucan_ac"):
            record["glytoucan_ids"].add(str(entry["glytoucan_ac"]))

    return {
        position: {
            "tier": record["tier"],
            "categories": "|".join(sorted(record["categories"])),
            "evidence_databases": "|".join(sorted(record["evidence_databases"])),
            "pubmed_ids": "|".join(sorted(record["pubmed_ids"])),
            "glytoucan_ids": "|".join(sorted(record["glytoucan_ids"])),
        }
        for position, record in sites.items()
    }


def fetch_details(accessions: list[str], config: Config) -> Path:
    """Fetch GlyGen protein details into a resumable JSONL cache."""
    cache_path = Path(config.paths["cache_dir"]) / "glygen_protein_detail.jsonl"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    already = set(load_cache(cache_path))

    api = config.api
    url_template = api.get("glygen_url", "https://api.glygen.org/protein/detail/{accession}/")
    delay = float(api.get("delay_seconds", 0.4))
    timeout = int(api.get("timeout_seconds", 60))
    retries = int(api.get("max_retries", 3))

    with cache_path.open("a", encoding="utf-8") as handle:
        for accession in accessions:
            if accession in already:
                continue
            payload, error = None, ""
            for attempt in range(retries):
                try:
                    with urllib.request.urlopen(
                        url_template.format(accession=accession), timeout=timeout
                    ) as response:
                        payload = json.load(response)
                    break
                except urllib.error.HTTPError as exc:
                    # An HTTP status from this API means "no entry for this
                    # accession" (it answers 500, not 404). That is permanent, so
                    # retrying triples the cost of every miss for no benefit.
                    # HTTPError subclasses URLError, so it must be caught first.
                    error = f"HTTPError: {exc}"
                    break
                except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                    error = f"{type(exc).__name__}: {exc}"
                    time.sleep(delay * (attempt + 1))
            record = {"accession": accession, "detail": payload, "error": error}
            handle.write(json.dumps(record) + "\n")
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
    """One row per candidate site with its GlyGen tier (empty when unsupported)."""
    per_accession = {
        accession: extract_sites(detail) for accession, detail in cache.items()
    }
    blank = {"tier": "", "categories": "", "evidence_databases": "",
             "pubmed_ids": "", "glytoucan_ids": ""}

    rows = []
    for accession, position in zip(candidates["accession"], candidates["position"]):
        accession, position = str(accession), int(position)
        record = per_accession.get(accession, {}).get(position, blank)
        rows.append({
            "accession": accession,
            "position": position,
            "glygen_tier": record["tier"],
            "glygen_categories": record["categories"],
            "glygen_evidence_databases": record["evidence_databases"],
            "glygen_pubmed_ids": record["pubmed_ids"],
            "glygen_glytoucan_ids": record["glytoucan_ids"],
        })
    return pd.DataFrame(rows).sort_values(["accession", "position"]).reset_index(drop=True)
