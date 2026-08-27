"""CARBonAra's alphabet, structure preparation, index mapping and score.

Every test here is hermetic: no CARBonAra checkout, no weights, no gemmi, no
network, no CUDA. The upstream parser is replaced by a stub that reproduces what
`src/structure.clean_structure` does to a PDB file, and the model by a fake whose
output we choose. Nothing here feeds the collinear geometry from `pdb_lines` to a
real structure model, which that module warns against.

The mapping tests carry real numbering gaps and insertion codes rather than only
gapless chains, because the 2026-08-25 ProteinMPNN correction was a gapless-only
test suite failing to notice a 25.3% misindexing.
"""
from __future__ import annotations

import numpy as np
import pytest

from experimental_glycosylation_sites import carbonara_scoring as cs
from pdb_lines import atom_line

# --------------------------------------------------------------------------
# Fixtures: PDB text with the awkward features real depositions have.
# --------------------------------------------------------------------------

BACKBONE = ("N", "CA", "C", "O")


def residue_lines(serial, resname, chain, resseq, icode=" ", atoms=BACKBONE,
                  record="ATOM  ", element=None, x0=0.0):
    """One residue's atoms. Coordinates are placeholders, not geometry."""
    lines = []
    for i, name in enumerate(atoms):
        el = element or ("SE" if name == "SE" else
                         "N" if name.startswith("N") else
                         "O" if name.startswith("O") else
                         "S" if name.startswith("S") else "C")
        line = atom_line(serial, name, resname, chain, resseq,
                         x0 + i * 1.4, 10.0, 10.0, el, record=record)
        # atom_line has no insertion-code slot; splice one in at column 27.
        if icode != " ":
            line = line[:26] + icode + line[27:]
        lines.append(line)
        serial += 1
    return lines, serial


def build_pdb(residues, chain="A"):
    """residues: list of (resname, resseq, icode, atoms, record)."""
    lines, serial = [], 1
    for resname, resseq, icode, atoms, record in residues:
        new, serial = residue_lines(serial, resname, chain, resseq, icode,
                                    atoms, record)
        lines += new
    return "\n".join(lines) + "\nTER\nEND\n"


def simple_chain(resnames, start=1, chain="A"):
    return [(rn, start + i, " ", BACKBONE, "ATOM  ") for i, rn in enumerate(resnames)]


# A short chain carrying an N-X-S/T sequon at manifest indices 2,3,4.
SEQUON_CHAIN = ["ALA", "GLY", "ASN", "LYS", "SER", "VAL", "LEU"]
SEQUON_SEQ = "AGNKSVL"


def write(tmp_path, text, name="test.pdb"):
    path = tmp_path / name
    path.write_text(text)
    return path


# --------------------------------------------------------------------------
# A stub of CARBonAra's own parse.
# --------------------------------------------------------------------------

def stub_structure(pdb_text, rm_hetatm=True, rm_wat=True):
    """What `read_pdb` + `clean_structure` produce, for our generated files.

    Reproduces the three behaviours the mapping depends on: hydrogens and water
    are dropped, HETATM records are dropped when `rm_hetatm`, and residues are
    renumbered consecutively from 1 in file order with the insertion code folded
    into that renumbering and then discarded.
    """
    names, elements, resnames, resids, het, xyz = [], [], [], [], [], []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip()
        resname = line[17:20].strip()
        flag = "H" if line.startswith("HETATM") else "A"
        if element in ("H", "D"):
            continue
        if rm_wat and resname in ("HOH", "DOD"):
            continue
        if rm_hetatm and flag != "A":
            continue
        names.append(line[12:16].strip())
        elements.append(element)
        resnames.append(resname)
        resids.append((line[21], int(line[22:26]), line[26]))
        het.append(flag)
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    renumbered, seen, current = [], None, 0
    for key in resids:
        if key != seen:
            current += 1
            seen = key
        renumbered.append(current)

    return {"name": np.array(names), "element": np.array(elements),
            "resname": np.array(resnames), "resid": np.array(renumbered),
            "het_flag": np.array(het), "xyz": np.array(xyz, dtype=np.float32)}


@pytest.fixture
def stub_parser(monkeypatch):
    monkeypatch.setattr(cs, "_load_structure",
                        lambda pdb_text, carbonara_dir=None: stub_structure(pdb_text))


# --------------------------------------------------------------------------
# 1. Alphabet — derived from upstream, never assumed.
# --------------------------------------------------------------------------

def test_alphabet_is_carbonaras_abundance_order_not_alphabetical():
    """`std_aminoacids` in src/data_encoding.py is sorted by abundance.

    Reading it as alphabetical would put N at 11 instead of 13 and score
    P(Pro) as P(Asn) — the shape of the 2026-08-20 defect.
    """
    assert cs.ALPHABET == "LERKVIFDYATSQNPGHWMC"
    assert cs.CARBONARA_RESNAMES[0] == "LEU"
    assert cs.CARBONARA_RESNAMES[-1] == "CYS"
    assert len(cs.CARBONARA_RESNAMES) == 20


def test_the_three_positions_the_score_reads():
    assert cs.AA_INDEX["N"] == 13
    assert cs.AA_INDEX["S"] == 11
    assert cs.AA_INDEX["T"] == 10
    assert cs.AA_INDEX["P"] == 14


def test_alphabet_round_trips_three_letter_to_one_letter():
    assert "".join(cs.RES3TO1[r] for r in cs.CARBONARA_RESNAMES) == cs.ALPHABET


def test_verify_alphabet_accepts_the_upstream_order():
    cs.verify_alphabet(list(cs.CARBONARA_RESNAMES))


def test_verify_alphabet_rejects_a_reordered_upstream():
    """If upstream ever reorders, scoring must stop rather than shift silently."""
    shuffled = list(cs.CARBONARA_RESNAMES)
    shuffled[13], shuffled[11] = shuffled[11], shuffled[13]
    with pytest.raises(cs.AlphabetMismatchError):
        cs.verify_alphabet(shuffled)


def test_verify_alphabet_rejects_a_different_length():
    with pytest.raises(cs.AlphabetMismatchError):
        cs.verify_alphabet(list(cs.CARBONARA_RESNAMES)[:19])


# --------------------------------------------------------------------------
# 3-6. Protein-only structure preparation.
# --------------------------------------------------------------------------

def test_preparation_removes_glycan_water_ligand_ion_and_other_chains(tmp_path):
    """The baseline arm must be protein-only, or it is not ProteinMPNN's input."""
    residues = simple_chain(SEQUON_CHAIN)
    residues += [("NAG", 500, " ", ("C1", "C2", "O5"), "HETATM"),
                 ("HOH", 600, " ", ("O",), "HETATM"),
                 ("SO4", 601, " ", ("S", "O1"), "HETATM"),
                 ("ZN", 602, " ", ("ZN",), "HETATM")]
    text = build_pdb(residues, chain="A")
    other, _ = residue_lines(900, "LEU", "B", 1)
    text += "\n".join(other) + "\nEND\n"

    path = write(tmp_path, text)
    prepared, residue_ids, native = cs.protein_only_pdb(path, "A")

    assert native == SEQUON_SEQ
    assert len(residue_ids) == len(SEQUON_CHAIN)
    for absent in ("NAG", "HOH", "SO4", " ZN"):
        assert absent not in prepared
    assert "HETATM" not in prepared
    # chain B contributed nothing
    assert all(line[21] == "A" for line in prepared.splitlines()
               if line.startswith("ATOM  "))


def test_mse_becomes_canonical_met_rather_than_disappearing(tmp_path):
    """CARBonAra's `clean_structure` has no MSE special case.

    Dropped as a HETATM it would remove a backbone position and shift every
    index after it; kept as MSE it is not in `std_aminoacids`, so it would be
    treated as a ligand and excluded from the amino-acid residue list. Either
    way the sequon moves. Converting it here is the only safe option.
    """
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("MSE", 2, " ", BACKBONE, "HETATM"),
                ("ASN", 3, " ", BACKBONE, "ATOM  "),
                ("LYS", 4, " ", BACKBONE, "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    prepared, residue_ids, native = cs.protein_only_pdb(path, "A")

    assert native == "AMNKS"
    assert "MSE" not in prepared
    assert "MET" in prepared
    assert "HETATM" not in prepared
    assert len(residue_ids) == 5


def test_selenium_atom_of_an_mse_is_not_emitted_as_a_met_atom(tmp_path):
    """SE is not a MET atom name; leaving it in would feed a stray element."""
    residues = [("MSE", 1, " ", ("N", "CA", "C", "O", "SE"), "HETATM")]
    path = write(tmp_path, build_pdb(residues))
    prepared, _, _ = cs.protein_only_pdb(path, "A")
    assert " SE " not in prepared


def test_hydrogens_are_removed(tmp_path):
    residues = [("ALA", 1, " ", ("N", "CA", "C", "O"), "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    text = path.read_text().replace(
        "END\n",
        atom_line(90, "H", "ALA", "A", 1, 1.0, 1.0, 1.0, "H") + "\nEND\n")
    path.write_text(text)
    prepared, _, native = cs.protein_only_pdb(path, "A")
    assert native == "A"
    assert not any(line[76:78].strip() == "H"
                   for line in prepared.splitlines() if line.startswith("ATOM  "))


def test_author_numbering_and_insertion_codes_are_kept_as_metadata(tmp_path):
    residues = [("ALA", 36, " ", BACKBONE, "ATOM  "),
                ("GLY", 36, "A", BACKBONE, "ATOM  "),
                ("ASN", 40, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    _, residue_ids, native = cs.protein_only_pdb(path, "A")
    assert native == "AGN"
    assert residue_ids == [(36, ""), (36, "A"), (40, "")]


def test_prepared_residues_are_renumbered_consecutively(tmp_path):
    """Author numbering can exceed the PDB column width and carry icodes.

    CARBonAra renumbers internally regardless, so emitting 1..N loses nothing
    and removes both hazards. The author numbering survives in the metadata.
    """
    residues = [("ALA", 36, " ", BACKBONE, "ATOM  "),
                ("GLY", 36, "A", BACKBONE, "ATOM  "),
                ("ASN", 4000, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    prepared, _, _ = cs.protein_only_pdb(path, "A")
    numbers = [int(line[22:26]) for line in prepared.splitlines()
               if line.startswith("ATOM  ") and line[12:16].strip() == "CA"]
    assert numbers == [1, 2, 3]
    assert all(line[26] == " " for line in prepared.splitlines()
               if line.startswith("ATOM  "))


def test_a_missing_chain_is_reported_not_guessed(tmp_path):
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    with pytest.raises(cs.ChainUnreadableError):
        cs.protein_only_pdb(path, "Z")


def test_pdb_and_mmcif_of_the_same_chain_agree(tmp_path):
    """The cache holds both suffixes; CARBonAra's loader reads only PDB.

    Preparation must therefore accept either and produce the same target chain,
    or the model silently scores a different set of sites depending on which
    format a structure happened to be cached in.
    """
    pytest.importorskip("Bio.PDB.MMCIFParser")
    pdb_path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    from Bio.PDB import MMCIFIO, PDBParser

    structure = PDBParser(QUIET=True).get_structure("t", str(pdb_path))
    cif_path = tmp_path / "test.cif"
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))

    from_pdb = cs.protein_only_pdb(pdb_path, "A")
    from_cif = cs.protein_only_pdb(cif_path, "A")
    assert from_pdb[2] == from_cif[2] == SEQUON_SEQ
    assert from_pdb[1] == from_cif[1]


# --------------------------------------------------------------------------
# 6-8. Index mapping.
# --------------------------------------------------------------------------

def test_contiguous_numbering_maps_directly():
    assert cs.build_index_map("ACD", "ACD") == {0: 0, 1: 1, 2: 2}


def test_mapping_is_rejected_when_lengths_disagree():
    """A dropped residue would shift the sequon; refuse rather than shift."""
    assert cs.build_index_map("ACD", "AC") == {}
    assert cs.build_index_map("ACD", "ACDE") == {}


def test_mapping_is_rejected_when_a_residue_disagrees():
    assert cs.build_index_map("ACD", "AWD") == {}


def test_noncanonical_residues_map_to_carbonaras_non_amino_acid_slot():
    """A residue outside `std_aminoacids` is a ligand to CARBonAra.

    It keeps its slot, so nothing shifts, but it has no amino-acid identity to
    imprint or read. Recorded as X on both sides and left unimprinted.
    """
    assert cs.build_index_map("ASXD", "ASXD") == {0: 0, 1: 1, 2: 2, 3: 3}


def test_a_noncanonical_residue_read_as_canonical_is_rejected():
    assert cs.build_index_map("ASXD", "ASND") == {}


def test_empty_inputs_give_an_empty_map():
    assert cs.build_index_map("", "") == {}


def test_numbering_gaps_do_not_shift_the_sequon(tmp_path, stub_parser):
    """The failure that inverted the ProteinMPNN headline, in CARBonAra's terms.

    ProteinMPNN inserts a placeholder per absent author number; CARBonAra
    renumbers observed residues consecutively instead, so a gap must NOT shift
    anything. Asserting the fixed point matters as much as asserting the shift.
    """
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("GLY", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 90, " ", BACKBONE, "ATOM  "),
                ("LYS", 91, " ", BACKBONE, "ATOM  "),
                ("SER", 92, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    mapping = cs.chain_mapping(path, "A")
    assert mapping.model_seq == "AGNKS"
    assert mapping.map_indices((2, 3, 4)) == [2, 3, 4]
    mapping.check_triplet([2, 3, 4], "NKS")


def test_insertion_codes_occupy_their_own_slots(tmp_path, stub_parser):
    residues = [("ASN", 36, " ", BACKBONE, "ATOM  "),
                ("LYS", 36, "A", BACKBONE, "ATOM  "),
                ("SER", 36, "B", BACKBONE, "ATOM  "),
                ("VAL", 37, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    mapping = cs.chain_mapping(path, "A")
    assert mapping.model_seq == "NKSV"
    assert mapping.map_indices((0, 1, 2)) == [0, 1, 2]
    mapping.check_triplet([0, 1, 2], "NKS")


def test_a_mapping_inconsistent_with_the_parsed_sequence_is_rejected(
        tmp_path, monkeypatch):
    """The identity mapping is a hypothesis about CARBonAra's parse.

    It is checked against what CARBonAra actually read, and a chain that
    disagrees is dropped rather than scored at an index we cannot justify.
    """
    def wrong(pdb_text, carbonara_dir=None):
        structure = stub_structure(pdb_text)
        structure["resname"] = np.array(["TRP"] * len(structure["resname"]))
        return structure

    monkeypatch.setattr(cs, "_load_structure", wrong)
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    with pytest.raises(cs.ChainUnreadableError):
        cs.chain_mapping(path, "A")


def test_an_expected_triplet_mismatch_is_rejected(tmp_path, stub_parser):
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    mapping = cs.chain_mapping(path, "A")
    mapping.check_triplet([2, 3, 4], "NKS")
    with pytest.raises(cs.SequonMismatchError):
        mapping.check_triplet([2, 3, 4], "NAS")


def test_an_out_of_range_index_is_rejected(tmp_path, stub_parser):
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    mapping = cs.chain_mapping(path, "A")
    with pytest.raises(cs.IncompleteBackboneError):
        mapping.map_indices((0, 1, 99))


# --------------------------------------------------------------------------
# 10. Backbone completeness and scoreability.
# --------------------------------------------------------------------------

def test_an_incomplete_backbone_makes_a_position_unscoreable(tmp_path, stub_parser):
    """CARBonAra needs N, CA and C to place its virtual C-beta, and O for the
    backbone encoding. A residue missing any of them is an exclusion."""
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("GLY", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 3, " ", ("N", "CA", "C"), "ATOM  "),   # no O
                ("LYS", 4, " ", BACKBONE, "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    mapping = cs.chain_mapping(path, "A")
    assert not bool(mapping.backbone_ok[2])
    with pytest.raises(cs.IncompleteBackboneError):
        mapping.map_indices((2, 3, 4))


def test_decodable_positions_are_returned_in_manifest_space(tmp_path, stub_parser):
    """Stage 05 indexes the array with `n_model_index`, a manifest ordinal."""
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("GLY", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 3, " ", ("N", "CA", "C"), "ATOM  "),
                ("LYS", 4, " ", BACKBONE, "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    decodable = cs.decodable_positions(path, "A")
    assert decodable.dtype == bool
    assert len(decodable) == 5
    assert list(decodable) == [True, True, False, True, True]


def test_an_unreadable_chain_is_entirely_undecodable_not_an_exception(tmp_path):
    """One bad structure must not stop a sweep of thousands."""
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    assert not cs.decodable_positions(path, "Z").any()


def test_decodable_positions_needs_no_model(tmp_path, stub_parser):
    """Scoreability precedes matching, so it cannot depend on a forward pass."""
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    assert cs.decodable_positions(path, "A").all()


# --------------------------------------------------------------------------
# 11-12. The conditional, and calibration.
# --------------------------------------------------------------------------

class FakeCARBonAra:
    """Records the imprint it was given and returns a chosen raw confidence.

    Mirrors the three upstream entry points `conditional_probabilities` uses:
    `process_structure`, `apply_model` and `conf`.
    """

    def __init__(self, n_residues=None, raw=None):
        # `n` is re-derived from whatever structure it is handed, so one fake
        # serves chains of different lengths.
        self.n = n_residues
        self.calls = []
        self.raw = raw

    def process_structure(self, structure):
        import torch as pt

        resids = np.unique(structure["resid"])
        n = self.n = len(resids)
        resnames = [structure["resname"][structure["resid"] == rid][0]
                    for rid in resids]
        y = pt.zeros((n, 20))
        mr_aa = pt.zeros(n, dtype=pt.bool)
        for i, rn in enumerate(resnames):
            if rn in cs.CARBONARA_RESNAMES:
                y[i, cs.CARBONARA_RESNAMES.index(rn)] = 1.0
                mr_aa[i] = True
        X = pt.zeros((n, 3))
        qe = pt.zeros((n, 33))
        Mr = pt.eye(n)
        return X, qe, None, None, Mr, None, y, mr_aa

    def apply_model(self, X, qe, Mr, yt=None):
        import torch as pt

        self.calls.append(None if yt is None else yt.clone())
        if self.raw is not None:
            return pt.as_tensor(self.raw, dtype=pt.float32)
        # Deterministic, non-uniform, and NOT normalised: raw sigmoid output.
        rows = pt.arange(self.n, dtype=pt.float32).unsqueeze(1)
        cols = pt.arange(20, dtype=pt.float32).unsqueeze(0)
        return pt.sigmoid(pt.sin(rows + cols))

    @staticmethod
    def conf(p):
        """Stands in for the empirical CDF map: clip, then normalise rows."""
        c = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
        return c / c.sum(axis=1, keepdims=True)


def test_the_scored_position_is_not_in_its_own_imprint(tmp_path, stub_parser):
    """The whole conditional depends on this.

    Leaving the scored residue in `yt` hands the model the answer, and P(N)
    would measure the imprint rather than the structure.
    """
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    mapping = cs.chain_mapping(path, "A")
    model = FakeCARBonAra(len(SEQUON_CHAIN))
    cs.conditional_probabilities(mapping, model, [2])

    assert len(model.calls) == 1
    yt = model.calls[0].numpy()
    assert yt[2].sum() == 0.0, "the scored position was imprinted"
    for other in (0, 1, 3, 4, 5, 6):
        assert yt[other].sum() == 1.0, f"position {other} was not imprinted"
    # and the imprint carries the native identity, not an arbitrary residue
    assert yt[4, cs.AA_INDEX["S"]] == 1.0


def test_each_position_gets_its_own_pass(tmp_path, stub_parser):
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    mapping = cs.chain_mapping(path, "A")
    model = FakeCARBonAra(len(SEQUON_CHAIN))
    probabilities = cs.conditional_probabilities(mapping, model, [2, 3, 4])
    assert len(model.calls) == 3
    assert set(probabilities) == {2, 3, 4}
    for index, call in zip((2, 3, 4), model.calls):
        assert call.numpy()[index].sum() == 0.0


def test_a_noncanonical_residue_is_never_imprinted(tmp_path, stub_parser):
    """It has no amino-acid identity, so there is nothing honest to imprint."""
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("SEP", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 3, " ", BACKBONE, "ATOM  "),
                ("LYS", 4, " ", BACKBONE, "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    mapping = cs.chain_mapping(path, "A")
    model = FakeCARBonAra(5)
    cs.conditional_probabilities(mapping, model, [2])
    yt = model.calls[0].numpy()
    assert yt[1].sum() == 0.0


def test_raw_confidences_are_calibrated_and_normalised(tmp_path, stub_parser):
    """`apply_model` returns independent sigmoids that do not sum to one.

    Scoring them directly would put a logit on something that is not a
    probability, so `conf` runs first and the result is checked.
    """
    path = write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))
    mapping = cs.chain_mapping(path, "A")
    model = FakeCARBonAra(len(SEQUON_CHAIN))

    raw = model.apply_model(None, None, None, yt=None).numpy()
    assert abs(raw[2].sum() - 1.0) > 0.1, "fixture is not exercising the point"

    probabilities = cs.conditional_probabilities(mapping, model, [2])
    assert probabilities[2].shape == (20,)
    assert abs(float(probabilities[2].sum()) - 1.0) < cs.PROBABILITY_SUM_TOLERANCE
    model.calls.clear()


# --------------------------------------------------------------------------
# 13. Malformed vectors are refused.
# --------------------------------------------------------------------------

def uniform(value=1.0 / 20):
    return np.full(20, value)


def test_a_wrong_length_vector_is_rejected():
    with pytest.raises(cs.InvalidProbabilityVector):
        cs.check_scoreable({0: np.full(21, 1.0 / 21)}, (0, 0, 0))


def test_an_unnormalised_vector_is_rejected():
    with pytest.raises(cs.InvalidProbabilityVector):
        cs.check_scoreable({0: uniform(0.5)}, (0, 0, 0))


def test_a_non_finite_vector_is_rejected():
    bad = uniform()
    bad[3] = np.nan
    with pytest.raises(cs.InvalidProbabilityVector):
        cs.check_scoreable({0: bad}, (0, 0, 0))


def test_a_negative_entry_is_rejected():
    bad = uniform()
    bad[3], bad[4] = -0.1, bad[4] + 0.1
    with pytest.raises(cs.InvalidProbabilityVector):
        cs.check_scoreable({0: bad}, (0, 0, 0))


def test_a_missing_position_is_rejected_rather_than_defaulted():
    with pytest.raises(cs.IncompleteBackboneError):
        cs.check_scoreable({0: uniform()}, (0, 1, 2))


def test_a_uniform_vector_is_accepted():
    cs.check_scoreable({0: uniform(), 1: uniform(), 2: uniform()}, (0, 1, 2))


# --------------------------------------------------------------------------
# 14. The score, against a hand calculation.
# --------------------------------------------------------------------------

def test_the_score_matches_a_hand_calculation():
    """0.5 * [logit P(N at i) + logit (P(S) + P(T)) at i+2].

    Chosen so the arithmetic is checkable by hand: P(N)=0.5 gives logit 0, and
    P(S)+P(T)=0.75 gives log 3. The mean is therefore log(3)/2.
    """
    n_row = np.full(20, 0.5 / 19)
    n_row[cs.AA_INDEX["N"]] = 0.5

    plus2 = np.zeros(20)
    plus2[cs.AA_INDEX["S"]] = 0.5
    plus2[cs.AA_INDEX["T"]] = 0.25
    plus2[cs.AA_INDEX["A"]] = 0.25

    plus1 = np.full(20, 0.9 / 19)
    plus1[cs.AA_INDEX["P"]] = 0.1

    score = cs.sequon_score({0: n_row, 1: plus1, 2: plus2}, 0, 1, 2)

    assert score["p_asn_at_n"] == pytest.approx(0.5)
    assert score["p_ser_at_plus2"] == pytest.approx(0.5)
    assert score["p_thr_at_plus2"] == pytest.approx(0.25)
    assert score["p_ser_or_thr_at_plus2"] == pytest.approx(0.75)
    assert score["p_pro_at_plus1"] == pytest.approx(0.1)
    assert score["logit_p_asn"] == pytest.approx(0.0)
    assert score["logit_p_ser_or_thr"] == pytest.approx(np.log(3.0))
    assert score["conditional_sequon_score"] == pytest.approx(np.log(3.0) / 2)


def test_the_score_reports_a_single_deterministic_conditional():
    row = uniform()
    score = cs.sequon_score({0: row, 1: row, 2: row}, 0, 1, 2)
    assert score["conditional_sequon_score_sd"] == 0.0
    assert score["n_decoding_orders"] == 1


def test_the_score_carries_full_twenty_entry_vectors():
    row = uniform()
    score = cs.sequon_score({0: row, 1: row, 2: row}, 0, 1, 2)
    for key in ("probs_n", "probs_plus1", "probs_plus2"):
        assert len(score[key]) == 20
        assert isinstance(score[key], list)


def test_the_middle_residue_does_not_enter_the_score():
    """Any residue but proline permits a sequon, so a preference there is not a
    preference for the motif. Proline is kept as a diagnostic only."""
    n_row = np.full(20, 0.5 / 19)
    n_row[cs.AA_INDEX["N"]] = 0.5
    plus2 = np.zeros(20)
    plus2[cs.AA_INDEX["S"]] = 0.75
    plus2[cs.AA_INDEX["A"]] = 0.25

    first = cs.sequon_score({0: n_row, 1: uniform(), 2: plus2}, 0, 1, 2)
    other = np.zeros(20)
    other[cs.AA_INDEX["P"]] = 1.0
    second = cs.sequon_score({0: n_row, 1: other, 2: plus2}, 0, 1, 2)

    assert (first["conditional_sequon_score"]
            == pytest.approx(second["conditional_sequon_score"]))
    assert first["p_pro_at_plus1"] != second["p_pro_at_plus1"]


def test_a_saturated_probability_gives_a_large_finite_score():
    """Clamping keeps an infinity out of every average downstream."""
    n_row = np.zeros(20)
    n_row[cs.AA_INDEX["N"]] = 1.0
    plus2 = np.zeros(20)
    plus2[cs.AA_INDEX["S"]] = 1.0
    score = cs.sequon_score({0: n_row, 1: uniform(), 2: plus2}, 0, 1, 2)
    assert np.isfinite(score["conditional_sequon_score"])
    assert score["conditional_sequon_score"] > 10
