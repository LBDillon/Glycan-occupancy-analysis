from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from experimental_glycosylation_sites.structures import (
    GlycanLink,
    assess_site,
    build_site_evidence,
    parse_link_records,
)
from pdb_lines import link_line

FIXTURE = Path(__file__).parent / "fixtures" / "mini_link.pdb"


def write_link(tmp_path: Path, name: str, line: str) -> Path:
    path = tmp_path / name
    path.write_text(line + "\n")
    return path


def test_parses_asn_first_link_order():
    links = parse_link_records(FIXTURE)
    assert links == [GlycanLink(chain_id="A", resseq=2, icode="", glycan_resname="NAG")]


def test_parses_glycan_first_link_order(tmp_path):
    line = link_line("NAG", "A", 101, "ASN", "A", 2, name1="C1", name2="ND2")
    assert parse_link_records(write_link(tmp_path, "flipped.pdb", line)) == [
        GlycanLink(chain_id="A", resseq=2, icode="", glycan_resname="NAG")
    ]


def test_parses_insertion_code_on_the_asparagine(tmp_path):
    line = link_line("ASN", "N", 86, "NAG", "N", 477, icode1="A")
    links = parse_link_records(write_link(tmp_path, "icode.pdb", line))
    assert links[0].icode == "A"
    assert links[0].resseq == 86


def test_ignores_links_without_an_asparagine(tmp_path):
    line = link_line("SER", "A", 10, "MAN", "A", 200, name1="OG")
    assert parse_link_records(write_link(tmp_path, "serine.pdb", line)) == []


def test_site_with_linked_glycan_is_tier_linked():
    links = parse_link_records(FIXTURE)
    result = assess_site("MNKTA", 2, FIXTURE, "MINI", links)
    assert result["tier"] == "structure_linked_glycan"
    assert result["chain_id"] == "A"
    assert result["resseq"] == 2
    assert result["observed_residue"] == "N"


def test_resolved_site_without_glycan_is_not_called_unmodified():
    result = assess_site("MNKTA", 4, FIXTURE, "MINI", parse_link_records(FIXTURE))
    assert result["tier"] == "structure_residue_resolved"
    assert "unmodified" not in result["tier"]


def test_position_outside_the_model_is_unresolved():
    result = assess_site("MNKTAWWWWW", 9, FIXTURE, "MINI", parse_link_records(FIXTURE))
    assert result["tier"] == "structure_residue_unresolved"


def test_missing_structure_file_is_not_assessed(tmp_path):
    result = assess_site("MNKTA", 2, tmp_path / "absent.pdb", "NONE", [])
    assert result["tier"] == "structure_not_assessed"
    assert result["detail"] == "structure_file_missing"


def test_mmcif_linkage_is_recorded_as_unsupported(tmp_path):
    path = tmp_path / "entry.cif"
    path.write_text("data_TEST\n")
    result = assess_site("MNKTA", 2, path, "TEST", [])
    assert result["tier"] == "structure_not_assessed"
    assert result["detail"] == "mmcif_linkage_unsupported"


def test_accession_without_a_cached_structure_is_not_assessed():
    candidates = pd.DataFrame([{"accession": "P99999", "position": 5}])
    frame = build_site_evidence(candidates, {"P99999": "MNKTA"}, {})
    assert frame.iloc[0].structure_tier == "structure_not_assessed"
    assert frame.iloc[0].structure_detail == "no_cached_structure"


def test_build_site_evidence_uses_the_manifest(tmp_path):
    candidates = pd.DataFrame([{"accession": "P1", "position": 2}])
    manifest = {"P1": {"accession": "P1", "pdb_id": "MINI",
                       "output_path": str(FIXTURE), "status": "already_present"}}
    frame = build_site_evidence(candidates, {"P1": "MNKTA"}, manifest)
    assert frame.iloc[0].structure_tier == "structure_linked_glycan"
