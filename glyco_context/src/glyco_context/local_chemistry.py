"""The local chemical environment of a sequon, as fixed-backbone design sees it.

Fixed-backbone design cannot move an atom: solvent accessibility, secondary
structure and distance to the termini are inherited from the input structure. It
can only change which residue sits at each position. This module measures that —
the composition of the sequence window around the Asn, and of the
three-dimensional shell around ND2 — so a wild-type site and a design of the same
backbone are directly comparable.

The sequon itself is excluded everywhere. Those three positions are the ones held
fixed during design, so counting them would guarantee partial agreement and
flatter the result.
"""
from __future__ import annotations

CLASSES = ("hydrophobic", "aromatic", "charged", "polar", "glycine", "proline",
           "cysteine")

AA_CLASS = {
    **{a: "hydrophobic" for a in "AVLIM"},
    **{a: "aromatic" for a in "FWY"},
    **{a: "charged" for a in "DEKR"},
    **{a: "polar" for a in "STNQH"},
    "G": "glycine", "P": "proline", "C": "cysteine",
}
CHARGE = {"D": -1, "E": -1, "K": +1, "R": +1}   # histidine neutral


def flank_indices(n_index: int, length: int, window: int = 5) -> "list[int]":
    """Sequence positions around the sequon, excluding the sequon itself."""
    sequon = {n_index, n_index + 1, n_index + 2}
    lower = range(max(0, n_index - window), n_index)
    upper = range(n_index + 3, min(length, n_index + 3 + window))
    return [i for i in list(lower) + list(upper) if i not in sequon]


def _composition(residues: "list[str]", prefix: str) -> dict:
    """Class fractions over a set of residues, or None where there are none.

    An empty set gives None rather than zero: no neighbours is a different fact
    from neighbours that happen to contain no aromatics, and averaging the two
    together would understate aromatic content wherever density is poor.
    """
    out = {f"{prefix}_{name}_fraction": None for name in CLASSES}
    if not residues:
        return out
    total = len(residues)
    for name in CLASSES:
        out[f"{prefix}_{name}_fraction"] = round(
            sum(1 for r in residues if AA_CLASS.get(r) == name) / total, 4)
    return out


def chemistry_panel(sequence: str, n_index: int,
                    shell_indices) -> dict:
    """The fifteen-feature panel for one site on one sequence.

    `shell_indices` are chain positions whose residues sit within the 8 A shell
    around ND2. They come from the structure and are the same for the wild type
    and every design of it, so any difference in the panel is a difference in
    sequence alone.
    """
    length = len(sequence)
    flank = [sequence[i] for i in flank_indices(n_index, length)]
    sequon = {n_index, n_index + 1, n_index + 2}
    shell = [sequence[i] for i in shell_indices
             if 0 <= i < length and i not in sequon]

    panel = {}
    panel.update(_composition(flank, "flank"))
    panel.update(_composition(shell, "shell"))
    panel["shell_net_charge"] = (sum(CHARGE.get(r, 0) for r in shell)
                                 if shell else None)
    return panel


def shell_indices_from_structure(path, chain_id: str, resseq: int, icode: str,
                                 sequence: str, n_index: int,
                                 radius: float = 8.0):
    """Chain indices within `radius` of ND2, or None if the mapping is unsafe.

    The design sequence is indexed by ProteinMPNN's parse; the geometry comes
    from BioPython. Those two must describe the same residues in the same order
    or every index is meaningless, so the alignment is verified rather than
    assumed: same length, same one-letter sequence, and an asparagine where the
    manifest says the Asn is. Any mismatch returns None and the site is dropped
    and counted, which is the behaviour that would have caught this class of bug
    the first time.
    """
    import numpy as np
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import is_aa
    from Bio.SeqUtils import seq1

    parser = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    try:
        chain = parser.get_structure("x", str(path))[0][str(chain_id)]
    except Exception:
        return None
    resolved = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    if len(resolved) != len(sequence):
        return None
    for residue, letter in zip(resolved, sequence):
        if seq1(residue.get_resname(), undef_code="X") != letter:
            return None
    if not (0 <= n_index < len(resolved)):
        return None
    asn = resolved[n_index]
    from experimental_glycosylation_sites.features import _clean_icode
    if (int(asn.id[1]) != int(resseq)
            or _clean_icode(asn.id[2]) != _clean_icode(icode)):
        return None
    if "ND2" not in asn:
        return None

    origin = asn["ND2"].coord
    indices = []
    for index, residue in enumerate(resolved):
        if residue is asn:
            continue
        for atom in residue:
            if atom.element == "H":
                continue
            if float(np.linalg.norm(atom.coord - origin)) <= radius:
                indices.append(index)
                break
    return indices


def disulfide_indices(path, chain_id: str, sequence: str, radius: float = 2.5):
    """Chain indices of cysteines whose SG is bonded to another chain's SG.

    Held fixed during design because a redesign that removes half a disulfide
    produces a structure the backbone no longer supports, and the resulting
    sequence is not a fair test of anything. Partners on other chains count: the
    bond constrains this residue either way.
    """
    import numpy as np
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import is_aa

    parser = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    try:
        model = parser.get_structure("x", str(path))[0]
        chain = model[str(chain_id)]
    except Exception:
        return []
    resolved = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    if len(resolved) != len(sequence):
        return []

    partners = [r for other in model for r in other
                if is_aa(r, standard=False) and "SG" in r]
    indices = []
    for index, residue in enumerate(resolved):
        if "SG" not in residue:
            continue
        origin = residue["SG"].coord
        for partner in partners:
            if partner is residue:
                continue
            if float(np.linalg.norm(partner["SG"].coord - origin)) <= radius:
                indices.append(index)
                break
    return indices
