"""Fixed-column PDB record formatters, used by fixtures and tests alike.

LINK column layout (1-indexed, per the PDB format spec):
  13-16 name1  18-20 resName1  22 chainID1  23-26 resSeq1  27 iCode1
  43-46 name2  48-50 resName2  52 chainID2  53-56 resSeq2  57 iCode2
"""
from __future__ import annotations


def link_line(
    res1: str, ch1: str, seq1: int, res2: str, ch2: str, seq2: int,
    name1: str = "ND2", name2: str = "C1",
    icode1: str = " ", icode2: str = " ", distance: float = 1.45,
) -> str:
    head = f"LINK        {name1:^4s} {res1:>3s} {ch1:1s}{seq1:4d}{icode1:1s}"
    tail = f"{name2:^4s} {res2:>3s} {ch2:1s}{seq2:4d}{icode2:1s}  1555   1555  {distance:5.2f}"
    return head + " " * 15 + tail


def atom_line(
    serial: int, name: str, resname: str, chain: str, resseq: int,
    x: float, y: float, z: float, element: str, record: str = "ATOM  ",
) -> str:
    return (
        f"{record}{serial:5d} {name:^4s} {resname:>3s} {chain:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
    )


def mini_structure() -> str:
    """Chain A holding MNKTA, with a NAG linked to the asparagine at resseq 2."""
    lines = [link_line("ASN", "A", 2, "NAG", "A", 101)]
    serial = 1
    for offset, (resname, element_seq) in enumerate([
        ("MET", "NCC"), ("ASN", "NCC"), ("LYS", "NCC"), ("THR", "NCC"), ("ALA", "NCC"),
    ]):
        for name, element in zip(("N", "CA", "C"), element_seq):
            lines.append(atom_line(
                serial, name, resname, "A", offset + 1,
                10.0 + serial, 10.0, 10.0, element,
            ))
            serial += 1
    lines.append(atom_line(serial, "C1", "NAG", "A", 101, 15.0, 12.0, 10.0, "C", "HETATM"))
    lines.append("END")
    return "\n".join(lines) + "\n"
