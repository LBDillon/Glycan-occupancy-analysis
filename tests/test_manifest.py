from __future__ import annotations

from pathlib import Path

import pandas as pd

from experimental_glycosylation_sites.manifest import (
    build_manifest,
    map_sequon,
    within_structure_pairs,
)
from pdb_lines import LONG_RESIDUES, LONG_SEQUENCE, chain_lines

FIXTURE = Path(__file__).parent / "fixtures" / "mini_link.pdb"


def write_structure(tmp_path, name, chains):
    lines, serial = [], 1
    for chain_id, resnames in chains:
        block, serial = chain_lines(chain_id, resnames, start_serial=serial)
        lines.extend(block)
    lines.append("END")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def test_maps_all_three_sequon_residues():
    """LONG_SEQUENCE begins MNKT..., so position 2 opens an N-K-T sequon."""
    result = map_sequon(LONG_SEQUENCE, 2, FIXTURE, "MINI", "A")
    assert result["ok"] is True
    assert result["triplet"] == "NKT"
    assert result["subtype"] == "NXT"
    assert result["n_aa"] == "N" and result["plus2_aa"] == "T"
    # model indices are 0-based ordinals in the chain the model reads
    assert result["plus1_model_index"] == result["n_model_index"] + 1
    assert result["plus2_model_index"] == result["n_model_index"] + 2


def test_missing_plus2_is_excluded_not_approximated(tmp_path):
    """The N resolving is not enough: a two-residue observation is not a sequon."""
    truncated = write_structure(tmp_path, "trunc.pdb", [("A", LONG_RESIDUES[:3])])
    result = map_sequon(LONG_SEQUENCE, 2, truncated, "TRUNC", "A")
    assert result["ok"] is False
    assert result["reason"] == "plus2_not_mapped"


def test_absent_chain_is_reported():
    result = map_sequon(LONG_SEQUENCE, 2, FIXTURE, "MINI", "Z")
    assert result["ok"] is False and result["reason"] == "chain_absent"


def test_triplet_must_match_the_sequon_pattern(tmp_path):
    """A structure whose residues do not spell N-X-S/T is not scoreable here."""
    residues = list(LONG_RESIDUES)
    residues[3] = "ALA"  # position 4 -> the +2 of the sequon at 2
    altered = write_structure(tmp_path, "alt.pdb", [("A", residues)])
    sequence = LONG_SEQUENCE[:3] + "A" + LONG_SEQUENCE[4:]
    result = map_sequon(sequence, 2, altered, "ALT", "A")
    assert result["ok"] is False
    assert result["reason"] == "triplet_not_a_sequon"
    assert result["triplet"] == "NKA"


def test_build_manifest_partitions_every_input_row():
    rows = pd.DataFrame([
        {"accession": "P1", "position": 2, "structure_pdb_id": "MINI",
         "structure_chain_id": "A", "occupancy_status": "occupied_supported"},
        {"accession": "P2", "position": 2, "structure_pdb_id": "MINI",
         "structure_chain_id": "Z", "occupancy_status": "occupied_supported"},
        {"accession": "P3", "position": 2, "structure_pdb_id": "ABSENT",
         "structure_chain_id": "A", "occupancy_status": "occupied_supported"},
    ])
    manifest, exclusions = build_manifest(
        rows, {"P1": LONG_SEQUENCE, "P2": LONG_SEQUENCE, "P3": LONG_SEQUENCE},
        {"MINI": FIXTURE},
    )
    assert len(manifest) + len(exclusions) == len(rows)
    assert set(exclusions.exclusion_reason) == {"chain_absent", "structure_not_cached"}


def test_missing_sequence_is_excluded_with_a_reason():
    rows = pd.DataFrame([{"accession": "PX", "position": 2, "structure_pdb_id": "MINI",
                          "structure_chain_id": "A", "occupancy_status": "occupied_supported"}])
    manifest, exclusions = build_manifest(rows, {}, {"MINI": FIXTURE})
    assert manifest.empty
    assert exclusions.iloc[0].exclusion_reason == "no_uniprot_sequence"


def test_within_structure_pairs_finds_shared_proteins():
    manifest = pd.DataFrame([
        {"accession": "P1", "position": 2, "structure_pdb_id": "AAAA",
         "occupancy_status": "occupied_supported"},
        {"accession": "P1", "position": 40, "structure_pdb_id": "AAAA",
         "occupancy_status": "observed_unmodified"},
        {"accession": "P9", "position": 5, "structure_pdb_id": "BBBB",
         "occupancy_status": "occupied_supported"},
    ])
    pairs = within_structure_pairs(manifest)
    assert len(pairs) == 1, "one distinct pair, not one row per level"
    row = pairs.iloc[0]
    assert bool(row.same_protein) and bool(row.same_structure)
    assert row.occupied_position == 2 and row.unmodified_position == 40


def test_same_structure_different_protein_is_kept_and_flagged():
    """Different proteins in one entry still share the experimental context."""
    manifest = pd.DataFrame([
        {"accession": "P1", "position": 2, "structure_pdb_id": "AAAA",
         "occupancy_status": "occupied_supported"},
        {"accession": "P2", "position": 7, "structure_pdb_id": "AAAA",
         "occupancy_status": "observed_unmodified"},
    ])
    pairs = within_structure_pairs(manifest)
    assert len(pairs) == 1
    assert bool(pairs.iloc[0].same_protein) is False
    assert bool(pairs.iloc[0].same_structure) is True
