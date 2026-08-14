from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from experimental_glycosylation_sites.cli import _accessions
from experimental_glycosylation_sites.config import load_config


def _write_config(tmp_path: Path) -> Path:
    """A tiny synthetic config + pairs CSV + gzipped UniProt TSV.

    Deliberately built in tmp_path rather than pointed at the real 40MB
    snapshot, so these tests stay hermetic and fast.
    """
    (tmp_path / "config").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sibling").mkdir(parents=True, exist_ok=True)

    pairs = pd.DataFrame([
        # P00001: analysis-ready, has both cross-references.
        {"pos_accession": "P00001", "best_loss_positive_position": 10, "analysis_ready": True},
        # P00002: analysis-ready, GlyGen column is blank, GlyConnect populated.
        {"pos_accession": "P00002", "best_loss_positive_position": 20, "analysis_ready": True},
        # P00003: analysis-ready, GlyGen populated, GlyConnect blank.
        {"pos_accession": "P00003", "best_loss_positive_position": 30, "analysis_ready": True},
        # P00004: NOT analysis-ready, must be excluded from every result below,
        # including the unrestricted (xref_column=None) case.
        {"pos_accession": "P00004", "best_loss_positive_position": 40, "analysis_ready": False},
    ])
    pairs.to_csv(tmp_path / "sibling" / "pairs.csv", index=False)

    tsv_path = tmp_path / "sibling" / "uniprot.tsv.gz"
    with gzip.open(tsv_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("Entry\tGlyGen\tGlyConnect\n")
        handle.write("P00001\tG12345\tC12345\n")
        handle.write("P00002\t\tC99999\n")
        handle.write("P00003\tG99999\t\n")
        # P00004 intentionally absent: it must never surface since it fails
        # the analysis_ready filter before the TSV is even consulted.

    cfg = tmp_path / "config" / "test.toml"
    cfg.write_text(
        "[paths]\n"
        'pairs_master = "../sibling/pairs.csv"\n'
        'uniprot_tsv = "../sibling/uniprot.tsv.gz"\n'
        'cache_dir = "cache"\n'
        "[layers]\n"
        "uniprot = true\n"
        "[policy]\n"
        "require_analysis_ready = true\n"
        "[api]\n"
        "delay_seconds = 0.0\n"
    )
    return cfg


def test_accessions_without_xref_column_returns_all_candidates(tmp_path):
    config = load_config(_write_config(tmp_path))
    assert _accessions(config, None) == ["P00001", "P00002", "P00003"]


def test_accessions_restricted_to_glygen_xref(tmp_path):
    config = load_config(_write_config(tmp_path))
    # P00002's GlyGen column is blank, so it must be excluded even though it
    # is an analysis-ready candidate. GlyGen answers HTTP 500 (not 404) for
    # accessions it has no entry for, so this restriction is what keeps
    # fetch-glygen from burning its whole run on unanswerable requests.
    assert _accessions(config, "GlyGen") == ["P00001", "P00003"]


def test_accessions_restricted_to_glyconnect_xref(tmp_path):
    config = load_config(_write_config(tmp_path))
    assert _accessions(config, "GlyConnect") == ["P00001", "P00002"]
