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
from pdb_lines import LONG_RESIDUES, LONG_SEQUENCE, chain_lines, link_line

FIXTURE = Path(__file__).parent / "fixtures" / "mini_link.pdb"

# A short chain that matches only the query's first two residues. Local
# alignment maps the queried position onto it at identity 1.0, so it is exactly
# the "short perfect sub-block on an unrelated chain" that the credibility floor
# (MIN_CHAIN_COVERAGE / MIN_ALIGNED_RESIDUES) exists to reject.
SHORT_DECOY_RESIDUES = ["MET", "ASN", "TRP", "TRP", "TRP"]


def write_link(tmp_path: Path, name: str, line: str) -> Path:
    path = tmp_path / name
    path.write_text(line + "\n")
    return path


def build_structure(
    chains: list[tuple[str, list[str]]],
    links: list[tuple] = (),
) -> str:
    """A PDB file from (chain_id, residue names) blocks plus optional LINK records."""
    lines = [link_line(*link) for link in links]
    serial = 1
    for chain_id, resnames in chains:
        block, serial = chain_lines(chain_id, resnames, start_serial=serial)
        lines.extend(block)
    lines.append("END")
    return "\n".join(lines) + "\n"


def substituted(positions: dict[int, str]) -> list[str]:
    """LONG_RESIDUES with 1-indexed positions replaced, for building homologues
    of the query at a controlled identity. Position 2 is never touched, so the
    asparagine every test targets stays put."""
    residues = list(LONG_RESIDUES)
    for position, resname in positions.items():
        assert position != 2, "position 2 must stay ASN"
        residues[position - 1] = resname
    return residues


def homodimer_structure() -> str:
    """Two copies of the same protein (an 1AVD-style homodimer). Chain A is
    missing its last eight C-terminal residues (unresolved) but carries a glycan
    LINK on its Asn. Chain B is fully resolved (higher coverage) but has no
    LINK. Both chains clear the credibility floor and align at identity 1.0, so
    both are the same protein; a coverage-sensitive tie would incorrectly drop
    chain A even though it is the copy actually carrying the glycan.
    """
    return build_structure(
        [("A", LONG_RESIDUES[:32]), ("B", LONG_RESIDUES)],
        [("ASN", "A", 2, "NAG", "A", 101)],
    )


def linked_but_unrelated_chain_structure() -> str:
    """Chain A: a full-length but divergent chain (eight substitutions, identity
    0.8) that carries a synthetic glycan LINK on its own Asn 2. Chain B: an
    exact match to the query, with no LINK. Chain A clears the credibility floor
    on its own — full coverage over 40 aligned residues — so only the identity
    gate can reject it, which is what this pins: a glycan link is never trusted
    just because it exists somewhere in the file.
    """
    divergent = substituted({
        5: "GLY", 9: "ALA", 13: "GLY", 17: "ALA",
        23: "GLY", 27: "ALA", 31: "GLY", 35: "ALA",
    })
    return build_structure(
        [("A", divergent), ("B", LONG_RESIDUES)],
        [("ASN", "A", 2, "NAG", "A", 101)],
    )


def two_chain_structure() -> str:
    """Chain A: the short decoy, sharing only the query's first two residues.
    Chain B: an exact match to the full query sequence. Chain A is written
    first, so file-order-wins chain selection would pick it over the chain that
    actually matches the query.
    """
    return build_structure([("A", SHORT_DECOY_RESIDUES), ("B", LONG_RESIDUES)])


def short_linked_chain_structure() -> str:
    """The short decoy as the only chain, carrying a glycan LINK on its Asn 2.
    Nothing in the file is credibly the queried protein, so the linkage must not
    be credited to the site however perfectly that two-residue block aligns.
    """
    return build_structure(
        [("A", SHORT_DECOY_RESIDUES)], [("ASN", "A", 2, "NAG", "A", 101)]
    )


def point_mutant_structure() -> str:
    """A single full-length chain one substitution away from the query (identity
    0.975), carrying a glycan LINK on its Asn 2. Real structures routinely
    differ from the canonical UniProt sequence by a point mutation or an
    isoform difference, and must still count as the same protein.
    """
    return build_structure(
        [("A", substituted({20: "LYS"}))], [("ASN", "A", 2, "NAG", "A", 101)]
    )


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
    result = assess_site(LONG_SEQUENCE, 2, FIXTURE, "MINI", links)
    assert result["tier"] == "structure_linked_glycan"
    assert result["chain_id"] == "A"
    assert result["resseq"] == 2
    assert result["observed_residue"] == "N"
    assert result["detail"] == ""


def test_resolved_site_without_glycan_is_not_called_unmodified():
    result = assess_site(LONG_SEQUENCE, 4, FIXTURE, "MINI", parse_link_records(FIXTURE))
    assert result["tier"] == "structure_residue_resolved"
    assert "unmodified" not in result["tier"]
    assert result["observed_residue"] == "T"


def test_position_outside_the_model_is_unresolved():
    # The query runs five residues past the end of the modelled chain.
    result = assess_site(
        LONG_SEQUENCE + "KKKKK", 43, FIXTURE, "MINI", parse_link_records(FIXTURE)
    )
    assert result["tier"] == "structure_residue_unresolved"


def test_missing_structure_file_is_not_assessed(tmp_path):
    result = assess_site(LONG_SEQUENCE, 2, tmp_path / "absent.pdb", "NONE", [])
    assert result["tier"] == "structure_not_assessed"
    assert result["detail"] == "structure_file_missing"


def test_mmcif_glycan_bonds_are_parsed_from_struct_conn():
    """mmCIF records connectivity in _struct_conn, not LINK lines.

    Recent and large depositions exist only in mmCIF, so refusing them would
    silently drop glycan evidence for exactly the newest structures.
    """
    links = parse_link_records(Path(__file__).parent / "fixtures" / "mini_conn.cif")
    assert GlycanLink("A", 2, "", "NAG") in links       # ASN listed first
    assert GlycanLink("B", 47, "A", "NAG") in links     # glycan first, insertion code kept


def test_mmcif_ignores_bonds_that_are_not_asn_glycan():
    """A disulfide is covalent too, and an O-linked serine bond is not N-linked."""
    links = parse_link_records(Path(__file__).parent / "fixtures" / "mini_conn.cif")
    assert not any(link.resseq in (5, 9) for link in links)   # CYS-CYS disulfide
    assert not any(link.resseq == 12 for link in links)       # SER-MAN, O-linked

def test_accession_without_a_cached_structure_is_not_assessed():
    candidates = pd.DataFrame([{"accession": "P99999", "position": 5}])
    frame = build_site_evidence(candidates, {"P99999": LONG_SEQUENCE}, {})
    assert frame.iloc[0].structure_tier == "structure_not_assessed"
    assert frame.iloc[0].structure_detail == "no_cached_structure"


def test_build_site_evidence_uses_the_manifest(tmp_path):
    candidates = pd.DataFrame([{"accession": "P1", "position": 2}])
    manifest = {"P1": {"accession": "P1", "pdb_id": "MINI",
                       "output_path": str(FIXTURE), "status": "already_present"}}
    frame = build_site_evidence(candidates, {"P1": LONG_SEQUENCE}, manifest)
    assert frame.iloc[0].structure_tier == "structure_linked_glycan"


def test_chain_selection_prefers_alignment_quality_over_file_order(tmp_path):
    path = tmp_path / "two_chain.pdb"
    path.write_text(two_chain_structure())
    result = assess_site(LONG_SEQUENCE, 2, path, "MULTI", [])
    assert result["chain_id"] == "B"
    assert result["resseq"] == 2
    assert result["observed_residue"] == "N"
    assert result["detail"] == ""


def test_glycan_link_on_lower_coverage_homodimer_copy_still_wins(tmp_path):
    path = tmp_path / "homodimer.pdb"
    path.write_text(homodimer_structure())
    links = parse_link_records(path)
    result = assess_site(LONG_SEQUENCE, 2, path, "DIMER", links)
    assert result["tier"] == "structure_linked_glycan"
    assert result["chain_id"] == "A"
    assert result["resseq"] == 2


def test_identity_gate_excludes_linked_but_unrelated_chain(tmp_path):
    path = tmp_path / "linked_unrelated.pdb"
    path.write_text(linked_but_unrelated_chain_structure())
    links = parse_link_records(path)
    result = assess_site(LONG_SEQUENCE, 2, path, "MULTI2", links)
    assert result["tier"] == "structure_residue_resolved"
    assert result["chain_id"] == "B"


def test_short_block_match_on_a_linked_chain_is_low_confidence(tmp_path):
    path = tmp_path / "short_linked.pdb"
    path.write_text(short_linked_chain_structure())
    links = parse_link_records(path)
    result = assess_site(LONG_SEQUENCE, 2, path, "SHORT", links)
    assert result["tier"] == "structure_residue_resolved"
    assert result["detail"] == "low_confidence_chain_match"


def test_full_length_point_mutant_is_still_the_same_protein(tmp_path):
    path = tmp_path / "point_mutant.pdb"
    path.write_text(point_mutant_structure())
    links = parse_link_records(path)
    result = assess_site(LONG_SEQUENCE, 2, path, "MUT", links)
    assert result["tier"] == "structure_linked_glycan"
    assert result["chain_id"] == "A"
    assert result["resseq"] == 2
    assert result["detail"] == ""


# --- manifest path resolution -------------------------------------------------
# The manifest records absolute paths from the machine that downloaded the
# structures. If those cannot be resolved elsewhere the structural layer silently
# disappears, so resolution is pinned in all three directions.

def _manifest_csv(tmp_path, output_path: str):
    path = tmp_path / "manifest.csv"
    path.write_text(
        "accession,pdb_id,output_path,status\n"
        f"P1,MINI,{output_path},already_present\n"
    )
    return path


def test_manifest_uses_the_recorded_path_when_it_exists(tmp_path):
    from experimental_glycosylation_sites.structures import load_manifest

    manifest = load_manifest(_manifest_csv(tmp_path, str(FIXTURE)), None)
    assert manifest["P1"]["output_path"] == str(FIXTURE)


def test_manifest_recovers_via_structure_dir_when_recorded_path_is_absent(tmp_path):
    from experimental_glycosylation_sites.structures import load_manifest

    stale = "/nonexistent/machine/structures/pdb/mini_link.pdb"
    manifest = load_manifest(_manifest_csv(tmp_path, stale), FIXTURE.parent)
    assert manifest["P1"]["output_path"] == str(FIXTURE)


def test_manifest_raises_rather_than_silently_returning_nothing(tmp_path):
    from experimental_glycosylation_sites.structures import (
        StructureManifestError,
        load_manifest,
    )

    stale = "/nonexistent/machine/structures/pdb/mini_link.pdb"
    with pytest.raises(StructureManifestError) as excinfo:
        load_manifest(_manifest_csv(tmp_path, stale), tmp_path / "also-missing")
    assert "structure_dir" in str(excinfo.value)
