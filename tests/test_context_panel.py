"""The frozen biological panel: definitions that must mean what they say.

Two failures this guards against. A feature named for one measurement while
computing another -- `sidechain_contacts_5a` counted residues, not contacts --
and a run length reported as though observed when its boundary was never
resolved, which biases loop lengths downward exactly where density is poor.
"""
from __future__ import annotations

import numpy as np
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from experimental_glycosylation_sites.context_features import (_loop_run,
                                                               _ss_key,
                                                               sequon_context)


def _res(resseq, icode=" ", name="ALA"):
    residue = Residue((" ", resseq, icode), name, "")
    residue.add(Atom("CA", np.array([0.0, 0.0, 0.0]), 0, 1, " ", "CA", 1, "C"))
    return residue


def _dssp(resolved, codes):
    return {_ss_key(r.id[1], r.id[2]): (code, 0.5)
            for r, code in zip(resolved, codes)}


def test_loop_run_counts_the_coil_residues_containing_the_asn():
    """H H T T T H H -- the Asn sits in a three-residue loop."""
    resolved = [_res(n) for n in range(1, 8)]
    dssp = _dssp(resolved, ["H", "H", "T", "S", "-", "H", "H"])
    length, censored = _loop_run(resolved, 3, dssp)
    assert length == 3
    assert censored is False


def test_loop_run_is_censored_at_a_chain_end():
    """A loop running off the end of the model has no measured boundary."""
    resolved = [_res(n) for n in range(1, 5)]
    dssp = _dssp(resolved, ["-", "-", "-", "H"])
    length, censored = _loop_run(resolved, 0, dssp)
    assert length == 3
    assert censored is True


def test_loop_run_is_censored_at_a_numbering_gap():
    """The loop may continue through residues the deposition never resolved."""
    resolved = [_res(10), _res(11), _res(12), _res(40)]
    dssp = _dssp(resolved, ["H", "-", "-", "-"])
    length, censored = _loop_run(resolved, 1, dssp)
    assert length == 2, "residue 40 is not part of this run"
    assert censored is True


def test_loop_run_is_none_when_the_asn_is_not_in_a_loop():
    resolved = [_res(n) for n in range(1, 4)]
    dssp = _dssp(resolved, ["H", "H", "H"])
    assert _loop_run(resolved, 1, dssp) == (None, False)


def test_loop_run_is_none_without_secondary_structure():
    resolved = [_res(n) for n in range(1, 4)]
    assert _loop_run(resolved, 1, {}) == (None, False)


ASN_SIDE_CHAIN = [
    # ND2 at the origin, so neighbour distances read directly off the coordinates
    ("ND2", 0.0, 0.0, 0.0, "N"), ("CG", 0.0, 1.4, 0.0, "C"),
    ("OD1", 1.2, 2.0, 0.0, "O"), ("CB", -1.0, 2.2, 0.0, "C"),
    ("CA", -0.8, 3.7, 0.0, "C"), ("N", -2.2, 4.1, 0.0, "N"),
    ("C", 0.2, 4.6, 0.0, "C"), ("O", 0.3, 5.8, 0.0, "O"),
]


def _atom_record(serial, name, resname, chain, resseq, x, y, z, element):
    return (f"ATOM  {serial:5d} {name:^4s} {resname:>3s} {chain:1s}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}")


def _neighbourhood_fixture(tmp_path, include_nd2=True):
    lines, serial = [], 1
    for name, x, y, z, element in ASN_SIDE_CHAIN:
        if name == "ND2" and not include_nd2:
            continue
        lines.append(_atom_record(serial, name, "ASN", "A", 10, x, y, z, element))
        serial += 1
    for name, resname, resseq, x, y, z, element in [
        ("CA", "ALA", 11, 5.0, 0.0, 0.0, "C"),      # 5.0 A from ND2
        ("CA", "PHE", 12, 6.0, 0.0, 0.0, "C"),      # 6.0 A, backbone
        ("CZ", "PHE", 12, 7.0, 0.0, 0.0, "C"),      # 7.0 A, aromatic side chain
        ("CA", "ALA", 13, 20.0, 0.0, 0.0, "C"),     # far outside
    ]:
        lines.append(_atom_record(serial, name, resname, "A", resseq, x, y, z, element))
        serial += 1
    lines.append(_atom_record(serial, "CA", "ALA", "B", 1, 0.0, 7.0, 0.0, "C"))
    path = tmp_path / "neighbourhood.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_nd2_neighbourhood_counts_atoms_and_residues_separately(tmp_path):
    """The spec asks for both: atoms measure crowding, residues measure how
    many distinct partners contribute it. One is not a substitute for the other.
    """
    out = sequon_context(_neighbourhood_fixture(tmp_path), "A", 10)
    assert out is not None
    assert out["nd2_atoms_8a_same_chain"] == 3, "ALA11 CA, PHE12 CA, PHE12 CZ"
    assert out["nd2_residues_8a_same_chain"] == 2, "ALA11 and PHE12"


def test_nd2_neighbourhood_separates_other_chains(tmp_path):
    """Oligomer interfaces and crystal contacts must not be silently mixed into
    local sequence context."""
    out = sequon_context(_neighbourhood_fixture(tmp_path), "A", 10)
    assert out["nd2_atoms_8a_other_chain"] == 1
    assert out["nd2_residues_8a_other_chain"] == 1


def test_nd2_features_are_absent_when_the_side_chain_is_unresolved(tmp_path):
    """No ND2 means no ND2-centred measurement -- not a zero, which would read
    as an uncrowded site."""
    out = sequon_context(_neighbourhood_fixture(tmp_path, include_nd2=False), "A", 10)
    assert out is not None
    assert out["has_nd2"] is False
    assert out["nd2_atoms_8a_same_chain"] is None
    assert out["nd2_residues_8a_same_chain"] is None


def test_nearest_aromatic_is_measured_to_the_side_chain(tmp_path):
    """A CA distance describes where the backbone is, not where the ring is."""
    out = sequon_context(_neighbourhood_fixture(tmp_path), "A", 10)
    assert out["nearest_aromatic_sidechain_nd2"] == pytest.approx(7.0, abs=0.01)


def test_residue_contact_count_is_named_for_what_it_counts(tmp_path):
    """`sidechain_contacts_5a` counted residues while reading as atom contacts.

    One phenylalanine contributing three heavy atoms inside the shell is one
    residue, not three. A name that does not distinguish the two is how the
    original 8 A ND2 atom count and this 5 A residue count were conflated.
    """
    lines, serial = [], 1
    for name, x, y, z, element in ASN_SIDE_CHAIN:
        lines.append(_atom_record(serial, name, "ASN", "A", 10, x, y, z, element))
        serial += 1
    for name, x, y, z in [("CA", 4.0, 0.0, 0.0), ("CZ", 4.5, 0.0, 0.0),
                          ("CE1", 4.2, 0.5, 0.0)]:
        lines.append(_atom_record(serial, name, "PHE", "A", 11, x, y, z, "C"))
        serial += 1
    path = tmp_path / "one_residue_many_atoms.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")

    out = sequon_context(path, "A", 10)
    assert "sidechain_contacts_5a" not in out
    assert out["sidechain_neighbour_residues_5a"] == 1, "one PHE, three atoms"


def _gap_fixture(tmp_path):
    """ASN 10, SER 11, then THR 13 -- the +2 was never resolved."""
    lines, serial = [], 1
    for name, x, y, z, element in ASN_SIDE_CHAIN:
        lines.append(_atom_record(serial, name, "ASN", "A", 10, x, y, z, element))
        serial += 1
    for resname, resseq, x in [("SER", 11, 4.0), ("THR", 13, 8.0)]:
        lines.append(_atom_record(serial, "CA", resname, "A", resseq, x, 0.0, 0.0, "C"))
        serial += 1
    path = tmp_path / "gap.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_qc_records_residue_numbers_and_insertion_codes(tmp_path):
    """Without these the row cannot be traced back to the residues measured."""
    out = sequon_context(_gap_fixture(tmp_path), "A", 10)
    assert out["n_resseq"] == 10
    assert out["n_icode"] == ""
    assert out["plus1_resseq"] == 11
    assert out["plus2_resseq"] is None, "the +2 was never resolved"


def test_mapping_continuity_is_reported(tmp_path):
    """The triplet check only catches a gap when the residue identity differs;
    continuity catches it regardless."""
    out = sequon_context(_gap_fixture(tmp_path), "A", 10)
    assert out["mapping_continuous"] is False


def test_dssp_availability_is_reported_per_position(tmp_path):
    """Chain-level dssp_ok overstates coverage: it says DSSP ran, not that all
    three sequon positions have a secondary-structure call."""
    out = sequon_context(_gap_fixture(tmp_path), "A", 10)
    for prefix in ("n", "plus1", "plus2"):
        assert f"{prefix}_dssp_ok" in out
    assert out["plus2_dssp_ok"] is False, "there is no +2 to assign"


def test_loop_run_appears_in_the_feature_row(tmp_path):
    out = sequon_context(_gap_fixture(tmp_path), "A", 10)
    assert "loop_run_length" in out
    assert "loop_run_censored" in out
