"""Sequence-to-structure mapping in the v2 context extractor.

These cover the four ways a sequon's +1/+2 environment can be attributed to the
wrong residue: numbering gaps, insertion codes, author-number arithmetic, and
chain identifiers DSSP cannot round-trip through legacy PDB.
"""
from __future__ import annotations

import pytest

from experimental_glycosylation_sites.context_features import sequon_context
from pdb_lines import atom_line


def _write(tmp_path, lines, name="fixture.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _residue(serial, resname, chain, resseq, icode=" ", x=10.0):
    """CA-only residue; sequon_context only requires a CA to accept a residue."""
    line = atom_line(serial, "CA", resname, chain, resseq, x, 10.0, 10.0, "C")
    # atom_line has no insertion-code parameter; column 27 (0-indexed 26) is iCode.
    return line[:26] + icode + line[27:]


def test_plus_two_is_unresolved_when_a_residue_is_missing(tmp_path):
    """A numbering gap must not be walked over.

    Chain has ASN 10, SER 11, then THR 13 -- residue 12 is absent from the
    model. Taking "the next resolved residue" makes THR 13 the +2 and reports a
    complete NST triplet, so a site whose real +2 was never observed is scored
    as though it had been. The gap has to surface as unresolved.
    """
    lines = [
        _residue(1, "ASN", "A", 10),
        _residue(2, "SER", "A", 11),
        _residue(3, "THR", "A", 13),
    ]
    out = sequon_context(_write(tmp_path, lines), "A", 10)
    assert out is not None
    assert out["plus1_residue"] == "S"
    assert out["plus2_residue"] is None, "residue 13 is not the sequon's +2"
    assert out["triplet_observed"] == "NS?"


def test_insertion_coded_residues_are_contiguous(tmp_path):
    """36, 36A, 36B, 36C is a run, not a gap.

    Chymotrypsin-numbered proteases put the sequon inside an insertion block.
    Walking author numbers alone would treat 36 -> 36A as a repeat; treating the
    block as discontinuous would wrongly blank the triplet.
    """
    lines = [
        _residue(1, "LEU", "A", 36),
        _residue(2, "LYS", "A", 36, "A"),
        _residue(3, "ASN", "A", 36, "B"),
        _residue(4, "ASP", "A", 36, "C"),
        _residue(5, "THR", "A", 37),
    ]
    out = sequon_context(_write(tmp_path, lines), "A", 36, icode="B")
    assert out is not None
    assert out["triplet_observed"] == "NDT"


def test_terminal_distances_count_residues_not_author_numbers(tmp_path):
    """`*_resolved` distances must be commensurate with chain_length_resolved.

    Author numbering skips; a residue count does not. Deriving one from
    max/min author numbers makes n + c + 1 != chain_length wherever the
    deposition has a gap, which is most of the PDB.
    """
    lines = [
        _residue(1, "ALA", "A", 10),
        _residue(2, "ASN", "A", 11),
        _residue(3, "SER", "A", 12),
        _residue(4, "THR", "A", 90),  # large numbering jump, still the 4th residue
    ]
    out = sequon_context(_write(tmp_path, lines), "A", 11)
    assert out is not None
    assert out["chain_length_resolved"] == 4
    assert out["distance_to_n_terminus_resolved"] == 1
    assert out["distance_to_c_terminus_resolved"] == 2
    assert (out["distance_to_n_terminus_resolved"]
            + out["distance_to_c_terminus_resolved"] + 1
            == out["chain_length_resolved"])


def _backbone(serial, resname, chain, resseq, origin, icode=" "):
    """N, CA, C, O for one residue, spaced so DSSP sees a connected backbone."""
    atoms = [("N", origin, 10.0), ("CA", origin + 1.45, 10.0),
             ("C", origin + 2.9, 10.0), ("O", origin + 3.1, 11.2)]
    lines = []
    for offset, (name, x, y) in enumerate(atoms):
        line = atom_line(serial + offset, name, resname, chain, resseq, x, y, 10.0,
                         name[0])
        lines.append(line[:26] + icode + line[27:])
    return lines, serial + len(atoms)


def _backbone_chain(chain_id, spec):
    """spec: list of (resname, resseq, icode). Returns PDB lines."""
    lines, serial, origin = [], 1, 0.0
    for resname, resseq, icode in spec:
        block, serial = _backbone(serial, resname, chain_id, resseq, origin, icode)
        lines.extend(block)
        origin += 3.8
    return lines


def _cif_chain(chain_id, spec):
    """Minimal mmCIF for one chain. Multi-character chain ids only reach the
    pipeline through mmCIF -- the legacy PDB chain field is one column wide --
    so reproducing that failure needs this format, not a .pdb fixture."""
    header = ["data_TEST", "#", "loop_", "_atom_site.group_PDB",
              "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
              "_atom_site.label_alt_id", "_atom_site.label_comp_id",
              "_atom_site.label_asym_id", "_atom_site.label_entity_id",
              "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code",
              "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
              "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
              "_atom_site.auth_seq_id", "_atom_site.auth_asym_id",
              "_atom_site.pdbx_PDB_model_num"]
    rows, serial, origin = [], 1, 0.0
    for seq, (resname, resseq, icode) in enumerate(spec, start=1):
        for name, dx, dy in (("N", 0.0, 10.0), ("CA", 1.45, 10.0),
                             ("C", 2.9, 10.0), ("O", 3.1, 11.2)):
            rows.append(
                f"ATOM {serial} {name[0]} {name} . {resname} {chain_id} 1 {seq} "
                f"{icode.strip() or '?'} {origin + dx:.3f} {dy:.3f} 10.000 "
                f"1.00 20.00 {resseq} {chain_id} 1")
            serial += 1
        origin += 3.8
    return "\n".join(header + rows + ["#", ""])


def test_dihedrals_are_not_computed_across_a_gap(tmp_path):
    """Psi needs the next residue's N; across a gap that atom belongs to
    another part of the chain entirely, and the angle is meaningless."""
    spec = [("ASN", 10, " "), ("SER", 11, " "), ("THR", 13, " ")]
    out = sequon_context(_write(tmp_path, _backbone_chain("A", spec)), "A", 10)
    assert out is not None
    assert out["n_psi"] is not None, "10 -> 11 is adjacent, psi is well defined"
    assert out["plus1_psi"] is None, "11 -> 13 spans a gap"


def test_dssp_runs_on_multi_character_chain_ids(tmp_path):
    """Chain 'AB' must not be lost to the legacy PDB single-column chain field.

    Every site in the run whose chain id was longer than one character failed
    DSSP, so secondary structure was missing from exactly the large assemblies
    the control arms draw most heavily on.
    """
    from experimental_glycosylation_sites.context_features import (_DSSP_CACHE,
                                                                   dssp_for_chain)
    _DSSP_CACHE.clear()
    spec = [("ALA", n, " ") for n in range(1, 9)]
    path = tmp_path / "multichain.cif"
    path.write_text(_cif_chain("AB", spec))
    result, reason = dssp_for_chain(path, "AB")
    assert "exceeds PDB format limit" not in reason
    assert result, f"no DSSP entries returned (reason: {reason!r})"


def test_dssp_entries_keep_their_insertion_code(tmp_path):
    """36 and 36A are different residues and can hold different structure.

    Keying the DSSP table on residue number alone collapses an insertion block
    onto one entry, so whichever residue is parsed last supplies the secondary
    structure for all of them.
    """
    from experimental_glycosylation_sites.context_features import (_DSSP_CACHE,
                                                                   dssp_for_chain)
    _DSSP_CACHE.clear()
    spec = [("ALA", 1, " "), ("ALA", 2, " "), ("ALA", 3, " "), ("ALA", 3, "A"),
            ("ALA", 3, "B"), ("ALA", 4, " "), ("ALA", 5, " "), ("ALA", 6, " ")]
    path = _write(tmp_path, _backbone_chain("A", spec), "icodes.pdb")
    result, reason = dssp_for_chain(path, "A")
    assert result, f"no DSSP entries returned (reason: {reason!r})"
    assert (3, "A") in result or ("3", "A") in result, (
        f"insertion-coded residues collapsed onto plain numbers: {sorted(result)}")
