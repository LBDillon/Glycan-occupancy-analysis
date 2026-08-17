"""Fixed-column PDB record formatters, used by fixtures and tests alike.

NOTE: the coordinates these produce are collinear placeholders, adequate for
residue identity, numbering and linkage parsing but meaningless as geometry.
Never feed a fixture built here to a structure-based model — ProteinMPNN on this
backbone returns saturated nonsense, because no real protein looks like a
straight line.

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


# A 40-residue, repeat-free amino acid sequence (all 20 canonical residues,
# then a second differently-ordered pass through the same 20) long enough to
# clear structures.MIN_ALIGNED_RESIDUES (30) on a full self-match. Used
# wherever a fixture needs a chain that's credibly "the same protein" as the
# query, not just a short coincidental local-alignment match. One-letter form
# verified against Bio.SeqUtils.seq1 rather than transcribed by hand:
# MNKTALWDGYSPVHCEQFIRGSTVLIPFYHQEDKMWNACR
LONG_RESIDUES = [
    "MET", "ASN", "LYS", "THR", "ALA", "LEU", "TRP", "ASP", "GLY", "TYR",
    "SER", "PRO", "VAL", "HIS", "CYS", "GLU", "GLN", "PHE", "ILE", "ARG",
    "GLY", "SER", "THR", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "HIS",
    "GLN", "GLU", "ASP", "LYS", "MET", "TRP", "ASN", "ALA", "CYS", "ARG",
]
LONG_SEQUENCE = "MNKTALWDGYSPVHCEQFIRGSTVLIPFYHQEDKMWNACR"


def chain_lines(
    chain_id: str, resnames: list[str], start_serial: int = 1, start_resseq: int = 1,
    x0: float = 10.0, y: float = 10.0, z: float = 10.0,
) -> tuple[list[str], int]:
    """CA-only ATOM records for a chain (sufficient for _parse_chains, which
    only requires a CA atom per residue). Returns the lines and the next free
    atom serial number, so callers can continue numbering across chains."""
    lines = []
    serial = start_serial
    for offset, resname in enumerate(resnames):
        lines.append(atom_line(
            serial, "CA", resname, chain_id, start_resseq + offset, x0 + serial, y, z, "C",
        ))
        serial += 1
    return lines, serial


def mini_structure() -> str:
    """Chain A holding the 40-residue LONG_RESIDUES sequence, with a NAG
    linked to the asparagine at resseq 2."""
    lines = [link_line("ASN", "A", 2, "NAG", "A", 101)]
    chain, serial = chain_lines("A", LONG_RESIDUES)
    lines.extend(chain)
    lines.append(atom_line(serial, "C1", "NAG", "A", 101, 15.0, 12.0, 10.0, "C", "HETATM"))
    lines.append("END")
    return "\n".join(lines) + "\n"
