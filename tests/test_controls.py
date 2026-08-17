from __future__ import annotations

import pandas as pd

from experimental_glycosylation_sites.controls import (
    CONTROL_SETS,
    OST_BEARING_TAXA,
    build_control_sites,
    composition,
    find_sequons,
    summarise,
)


def test_finds_overlapping_sequons():
    """NNSS opens a sequon at both asparagines; a regex scan would drop the second."""
    assert find_sequons("NNSS") == [1, 2]


def test_proline_blocks_the_sequon():
    assert find_sequons("NPS") == []
    assert find_sequons("NPT") == []


def test_third_position_must_be_serine_or_threonine():
    assert find_sequons("NAT") == [1]
    assert find_sequons("NAS") == [1]
    assert find_sequons("NAC") == []
    assert find_sequons("NAA") == []


def test_truncated_motifs_at_the_terminus_are_not_counted():
    assert find_sequons("NA") == []
    assert find_sequons("N") == []


def test_every_ost_bearing_clade_is_excluded_from_the_bacterial_query():
    """Excluding by machinery, not by annotation - the whole point of the set."""
    query = CONTROL_SETS["bacterial_extracytoplasmic"]["query"]
    for taxon in OST_BEARING_TAXA:
        assert f"NOT taxonomy_id:{taxon}" in query
    # archaea in particular: a genuine OST, easy to overlook as "prokaryote"
    assert "NOT taxonomy_id:2157" in query


def test_control_queries_require_a_structure_and_exclude_glycoproteins():
    for name, spec in CONTROL_SETS.items():
        assert "database:pdb" in spec["query"], name
        assert "NOT keyword:KW-0325" in spec["query"], name


def test_cytosolic_query_excludes_signal_peptide_and_transmembrane():
    """A cytosolic label alone is not enough; secretory-pathway markers must be absent."""
    query = CONTROL_SETS["cytosolic_eukaryotic"]["query"]
    assert "NOT keyword:KW-0732" in query  # signal peptide
    assert "NOT keyword:KW-0812" in query  # transmembrane


def _proteins() -> pd.DataFrame:
    return pd.DataFrame([
        {"Entry": "P1", "Organism": "Escherichia coli", "Organism (ID)": "83333",
         "Sequence": "MNKTAAAANAS", "PDB": "1ABC;2DEF",
         "control_set": "bacterial_extracytoplasmic"},
        {"Entry": "P2", "Organism": "Homo sapiens", "Organism (ID)": "9606",
         "Sequence": "MAAAAAAA", "PDB": "",
         "control_set": "cytosolic_eukaryotic"},
    ])


def test_build_control_sites_labels_provenance_and_records_sequon():
    sites = build_control_sites(_proteins())
    assert set(sites.control_set) == {"bacterial_extracytoplasmic"}
    assert list(sites.position) == [2, 9]
    assert list(sites.sequon) == ["NKT", "NAS"]
    assert int(sites.iloc[0].n_pdb_entries) == 2


def test_protein_without_a_sequon_contributes_nothing():
    sites = build_control_sites(_proteins())
    assert "P2" not in set(sites.accession)


def test_composition_sums_to_one():
    result = composition(["AAAC", "CCGG"])
    assert abs(sum(result.values()) - 1.0) < 1e-9
    assert result["C"] == 0.375


def test_summary_reports_density_and_excluded_taxa():
    proteins = _proteins()
    summary = summarise(build_control_sites(proteins), proteins)
    entry = summary["bacterial_extracytoplasmic"]
    assert entry["sequons"] == 2
    assert entry["sequons_per_1000_residues"] is not None
    assert "composition" in entry
    # the exclusions travel with the data, so a reader can audit them
    assert summary["_excluded_taxa"] == OST_BEARING_TAXA


def test_repeated_protein_does_not_duplicate_its_sequons():
    """Paginated UniProt responses can repeat an entry across page boundaries.

    A repeat would emit that protein's sequons twice and then square in any join
    keyed on (accession, position) — 2 copies became 4 rows in the first build.
    """
    doubled = pd.concat([_proteins(), _proteins()], ignore_index=True)
    sites = build_control_sites(doubled)
    assert not sites.duplicated(["control_set", "accession", "position"]).any()
    assert len(sites) == len(build_control_sites(_proteins()))


def test_cytosolic_query_is_restricted_to_eukaryota():
    """GO:0005829 is domain-agnostic.

    Without an explicit Eukaryota restriction the set filled with bacteria (E.
    coli was its third most common organism) and 302 archaeal proteins — and
    archaea glycosylate via AglB, so those are potential true glycoproteins in a
    set defined as unable to be glycosylated.
    """
    assert "taxonomy_id:2759" in CONTROL_SETS["cytosolic_eukaryotic"]["query"]


def test_control_sets_target_disjoint_domains():
    """A protein must not be able to satisfy both queries.

    Bacterial proteins annotated both cytosolic and periplasmic previously
    appeared in both sets, which makes the two comparisons non-independent.
    """
    cyto = CONTROL_SETS["cytosolic_eukaryotic"]["query"]
    bact = CONTROL_SETS["bacterial_extracytoplasmic"]["query"]
    assert "taxonomy_id:2759" in cyto and "taxonomy_id:2" in bact
    assert "taxonomy_id:2759" not in bact
