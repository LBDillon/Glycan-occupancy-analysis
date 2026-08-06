from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SiteKey:
    accession: str
    position: int


@dataclass(frozen=True)
class UniProtFeature:
    accession: str
    position: int | None
    glyco_type: str
    evidence_codes: frozenset[str]
    raw_note: str
    parse_status: str
