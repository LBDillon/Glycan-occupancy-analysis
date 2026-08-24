"""Structure provenance: how the coordinates were obtained, not what they say.

Resolution and method belong with the QC fields. A 3.5 A structure and a 1.2 A
structure do not support the same claim about a side chain's environment, and
without the field an analysis cannot tell which it is looking at.
"""
from __future__ import annotations

from glyco_context.context_features import structure_metadata

PDB_HEADER = """HEADER    HYDROLASE                               01-JAN-00   1ABC
EXPDTA    X-RAY DIFFRACTION
REMARK   2 RESOLUTION.    1.85 ANGSTROMS.
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C
END
"""

CIF_HEADER = """data_TEST
_exptl.method 'ELECTRON MICROSCOPY'
_refine.ls_d_res_high 3.40
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 10.000 10.000 10.000 1.00 20.00 1 A 1
#
"""


def test_reads_resolution_and_method_from_pdb(tmp_path):
    path = tmp_path / "1abc.pdb"
    path.write_text(PDB_HEADER)
    meta = structure_metadata(path)
    assert meta["structure_resolution"] == 1.85
    assert "X-RAY" in meta["structure_method"].upper()


def test_reads_resolution_and_method_from_mmcif(tmp_path):
    path = tmp_path / "test.cif"
    path.write_text(CIF_HEADER)
    meta = structure_metadata(path)
    assert meta["structure_resolution"] == 3.40
    assert "ELECTRON" in meta["structure_method"].upper()


def test_missing_metadata_is_none_not_zero(tmp_path):
    """An absent resolution must not read as a perfect one."""
    path = tmp_path / "bare.pdb"
    path.write_text("ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C\nEND\n")
    meta = structure_metadata(path)
    assert meta["structure_resolution"] is None
