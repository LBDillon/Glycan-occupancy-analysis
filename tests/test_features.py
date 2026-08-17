from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio.PDB.Polypeptide import is_aa

from experimental_glycosylation_sites.features import (
    _clean_icode,
    _model_with_sasa,
    build_features,
    residue_features,
)

FIXTURE = Path(__file__).parent / "fixtures" / "mini_link.pdb"


def test_sasa_model_excludes_glycans():
    """The decisive one: an occupied site's own sugar must not shade it.

    mini_link.pdb carries a NAG. If it survived into the accessibility
    calculation, every occupied site would look more buried than an unmodified
    one purely because it is glycosylated — burial would encode the label.
    """
    _model_with_sasa.__globals__["_MODEL_CACHE"].clear()
    model = _model_with_sasa(FIXTURE)
    resnames = {r.get_resname() for chain in model for r in chain}
    assert "NAG" not in resnames
    assert all(is_aa(r, standard=False) for chain in model for r in chain)


def test_clean_icode_handles_missing_values():
    """Insertion codes arrive as NaN from CSV; str(nan) matches nothing."""
    assert _clean_icode(float("nan")) == ""
    assert _clean_icode(None) == ""
    assert _clean_icode("nan") == ""
    assert _clean_icode("") == ""
    assert _clean_icode(" A ") == "A"


def test_residue_features_for_a_known_residue():
    features = residue_features(FIXTURE, "A", 2)
    assert features is not None
    assert features["observed_residue"] == "N"
    assert features["rsa"] is not None
    assert features["rsa_bin"] in {"buried", "intermediate", "exposed"}
    assert features["n_neighbours_8a"] >= 0


def test_absent_residue_and_chain_return_none():
    assert residue_features(FIXTURE, "A", 9999) is None
    assert residue_features(FIXTURE, "Z", 2) is None


def test_build_features_keeps_every_site():
    """Sites without coordinates stay in the output, flagged, not filtered away."""
    sites = pd.DataFrame([
        {"accession": "P1", "position": 2, "occupancy_status": "occupied_supported",
         "structure_pdb_id": "MINI", "structure_chain_id": "A",
         "structure_resseq": 2, "structure_icode": float("nan")},
        {"accession": "P2", "position": 5, "occupancy_status": "unknown",
         "structure_pdb_id": "", "structure_chain_id": "",
         "structure_resseq": None, "structure_icode": float("nan")},
    ])
    frame = build_features(sites, {"MINI": FIXTURE})
    assert len(frame) == 2
    by_accession = frame.set_index("accession")
    assert bool(by_accession.loc["P1", "features_available"]) is True
    assert bool(by_accession.loc["P2", "features_available"]) is False
    assert by_accession.loc["P1", "observed_residue"] == "N"


def test_carried_columns_survive_the_output_resort():
    """The regression that mislabelled 43.5% of control rows.

    build_features sorts its output. A caller that assigns a provenance column
    afterwards with .values pairs labels from input order against rows in sorted
    order, silently scrambling them. Here the input order is deliberately the
    reverse of the sorted order, so a positional assignment would swap the labels.
    """
    sites = pd.DataFrame([
        {"accession": "Z9", "position": 2, "structure_pdb_id": "MINI",
         "structure_chain_id": "A", "structure_resseq": 2,
         "structure_icode": float("nan"), "control_set": "set_last"},
        {"accession": "A1", "position": 2, "structure_pdb_id": "MINI",
         "structure_chain_id": "A", "structure_resseq": 2,
         "structure_icode": float("nan"), "control_set": "set_first"},
    ])
    frame = build_features(sites, {"MINI": FIXTURE}, carry_columns=("control_set",))
    by_accession = frame.set_index("accession")["control_set"]
    assert by_accession.loc["Z9"] == "set_last"
    assert by_accession.loc["A1"] == "set_first"
