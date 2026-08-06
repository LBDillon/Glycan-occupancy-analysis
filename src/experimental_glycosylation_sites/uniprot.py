from __future__ import annotations

import gzip
import re
from pathlib import Path

from .models import UniProtFeature

_ENTRY_SPLIT = re.compile(r"(?=CARBOHYD\s)")
_EXACT_POSITION = re.compile(r"^CARBOHYD\s+(\d+)\s*;")
_ANY_POSITION = re.compile(r"^CARBOHYD\s+(\S+)")
_NOTE = re.compile(r'/note="([^"]*)"')
_ECO = re.compile(r"ECO:\d{7}")


def _glyco_type(note: str) -> str:
    note = note.lower()
    if "n-linked" in note:
        if "glycation" in note:
            return "glycation"
        return "N-linked" if "asparagine" in note else "N-linked-other"
    for prefix, label in (("o-linked", "O-linked"), ("c-linked", "C-linked"), ("s-linked", "S-linked")):
        if prefix in note:
            return label
    return "unknown"


def parse_glycosylation_column(text: str | None, accession: str) -> list[UniProtFeature]:
    """Parse a UniProt Glycosylation (ft_carbohyd) column into features.

    Only single-residue positions yield a usable position. Ranges ("10..12")
    and uncertain positions ("?") are returned with position=None and
    parse_status="uncertain_or_range_position" so they can be excluded with a
    reason rather than silently truncated to their start coordinate.
    """
    if not text or not str(text).strip():
        return []

    features: list[UniProtFeature] = []
    for entry in _ENTRY_SPLIT.split(str(text)):
        entry = entry.strip()
        if not entry.startswith("CARBOHYD"):
            continue

        note_match = _NOTE.search(entry)
        note = note_match.group(1) if note_match else ""
        codes = frozenset(_ECO.findall(entry))

        exact = _EXACT_POSITION.match(entry)
        if exact:
            position, status = int(exact.group(1)), "ok"
        elif _ANY_POSITION.match(entry):
            position, status = None, "uncertain_or_range_position"
        else:
            position, status = None, "malformed_feature"

        features.append(UniProtFeature(
            accession=accession,
            position=position,
            glyco_type=_glyco_type(note),
            evidence_codes=codes,
            raw_note=note,
            parse_status=status,
        ))
    return features


def load_uniprot_features(
    tsv_path: Path, accessions: set[str]
) -> tuple[list[UniProtFeature], set[str]]:
    """Parse CARBOHYD features for the requested accessions.

    Returns (features, accessions_absent_from_snapshot). Streams the gzipped
    TSV rather than loading it whole.
    """
    import csv

    opener = gzip.open if Path(tsv_path).suffix == ".gz" else open
    seen: set[str] = set()
    features: list[UniProtFeature] = []

    with opener(tsv_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = (row.get("Entry") or "").strip()
            if accession not in accessions:
                continue
            seen.add(accession)
            features.extend(
                parse_glycosylation_column(row.get("Glycosylation"), accession)
            )
    return features, accessions - seen
