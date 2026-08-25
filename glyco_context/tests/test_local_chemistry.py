"""The local chemical environment panel.

This is what fixed-backbone design can actually change, so it is the outcome of
the context-retention test. Positions come from the structure and never move;
only the identities at those positions differ between wild type and design.

The sequon itself is excluded from every count. Including it would put the three
residues we deliberately held fixed into the measurement, guaranteeing partial
agreement and flattering the result.
"""
from __future__ import annotations

import pytest

from glyco_context.local_chemistry import (CLASSES, chemistry_panel,
                                           flank_indices)


def test_flank_window_excludes_the_sequon():
    """+-2 around an Asn at index 10 is 8,9 and 13,14 -- never 10,11,12."""
    assert flank_indices(10, length=40, window=2) == [8, 9, 13, 14]


def test_flank_window_is_clipped_at_the_chain_start():
    assert flank_indices(1, length=40, window=3) == [0, 4, 5, 6]


def test_flank_window_is_clipped_at_the_chain_end():
    assert flank_indices(36, length=40, window=3) == [33, 34, 35, 39]


def test_panel_reports_a_fraction_for_every_class():
    sequence = "A" * 40
    panel = chemistry_panel(sequence, n_index=10, shell_indices=[20, 21])
    for name in CLASSES:
        assert f"flank_{name}_fraction" in panel
        assert f"shell_{name}_fraction" in panel
    assert "shell_net_charge" in panel


def test_hydrophobic_flank_is_measured():
    sequence = list("A" * 40)
    for i in (5, 6, 7, 8, 9, 13, 14, 15, 16, 17):
        sequence[i] = "L"
    panel = chemistry_panel("".join(sequence), n_index=10, shell_indices=[])
    assert panel["flank_hydrophobic_fraction"] == pytest.approx(1.0)


def test_shell_net_charge_counts_acids_and_bases_with_histidine_neutral():
    sequence = list("A" * 40)
    sequence[20], sequence[21], sequence[22], sequence[23] = "D", "E", "K", "H"
    panel = chemistry_panel("".join(sequence), n_index=10,
                            shell_indices=[20, 21, 22, 23])
    assert panel["shell_net_charge"] == -1     # -1 -1 +1, histidine neutral


def test_shell_excludes_the_sequon_positions():
    """A shell that happens to contain the sequon must not count it."""
    sequence = list("A" * 40)
    sequence[10], sequence[11], sequence[12] = "N", "G", "T"
    sequence[20] = "W"
    panel = chemistry_panel("".join(sequence), n_index=10,
                            shell_indices=[10, 11, 12, 20])
    assert panel["shell_aromatic_fraction"] == pytest.approx(1.0), "only W counts"


def test_empty_shell_gives_none_rather_than_zero():
    """No neighbours is not the same as neighbours that are all non-aromatic."""
    panel = chemistry_panel("A" * 40, n_index=10, shell_indices=[])
    assert panel["shell_aromatic_fraction"] is None
    assert panel["shell_net_charge"] is None


def _write(tmp_path, lines):
    path = tmp_path / "chain.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _atom(serial, name, resname, resseq, x, element="C"):
    return (f"ATOM  {serial:5d} {name:^4s} {resname:>3s} A{resseq:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}")


def test_shell_returns_none_when_the_sequence_does_not_align(tmp_path):
    """A length or identity mismatch means every index is meaningless."""
    from glyco_context.local_chemistry import shell_indices_from_structure
    lines = [_atom(1, "CA", "ALA", 1, 0.0), _atom(2, "CA", "ALA", 2, 4.0)]
    # structure says AA, caller claims a three-residue chain
    assert shell_indices_from_structure(_write(tmp_path, lines), "A", 1, "",
                                        "AAA", 0) is None


def test_shell_returns_none_when_the_asn_is_not_where_claimed(tmp_path):
    from glyco_context.local_chemistry import shell_indices_from_structure
    lines = [_atom(1, "CA", "ALA", 1, 0.0), _atom(2, "CA", "ALA", 2, 4.0)]
    assert shell_indices_from_structure(_write(tmp_path, lines), "A", 99, "",
                                        "AA", 0) is None
