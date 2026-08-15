"""Negative control sets for model-scoring experiments.

These are sequons that cannot be N-glycosylated, for reasons of biology rather
than absence of annotation. They are NOT part of the resource: they never join
the candidate universe, never count toward a site total, and always carry a
`control_set` label recording their provenance.

Two sets, chosen so their confounds do not overlap. Cytosolic eukaryotic
proteins never meet oligosaccharyltransferase because they never enter the
secretory pathway; they match the positives on taxonomy but differ in
compartment. Bacterial periplasmic and outer-membrane proteins cross a membrane
and fold in an oxidising compartment much like the ER lumen, but their clades
have no OST; they match on compartment but differ in taxonomy. A result that
survives both, and the 32 structural internal controls, cannot be explained by
either confound alone.

See docs/negative_controls.md for the full rationale.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

SEQUON = re.compile(r"N[^P][ST]")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Bacterial N-glycosylation is rare but real. Excluding it by annotation alone
# would repeat the absence-of-evidence error this project exists to avoid, so
# clades are excluded by known machinery. Archaea glycosylate via AglB — a
# genuine OST — and are excluded wholesale.
OST_BEARING_TAXA = {
    2157: "Archaea (AglB oligosaccharyltransferase)",
    194: "Campylobacter (PglB, the best-characterised bacterial OST)",
    209: "Helicobacter (PglB-family machinery)",
    724: "Haemophilus (HMW1C-type N-glycosyltransferase, acts on N-X-S/T)",
    713: "Actinobacillus (HMW1C homologue)",
    629: "Yersinia (HMW1C homologue)",
    32257: "Kingella (HMW1C homologue)",
}

_EXCLUDE_TAXA = " ".join(f"AND NOT taxonomy_id:{t}" for t in sorted(OST_BEARING_TAXA))

CONTROL_SETS: dict[str, dict] = {
    "cytosolic_eukaryotic": {
        "query": (
            "reviewed:true AND go:0005829 AND database:pdb "
            "AND NOT keyword:KW-0325 AND NOT keyword:KW-0732 "
            "AND NOT keyword:KW-0812"
        ),
        "rationale": (
            "Cytosolic, no signal peptide, no transmembrane region: never enters "
            "the secretory pathway, so its sequons never meet OST. Matches the "
            "positives on taxonomy; differs in compartment."
        ),
    },
    "bacterial_extracytoplasmic": {
        "query": (
            "reviewed:true AND taxonomy_id:2 AND database:pdb "
            "AND (keyword:KW-0574 OR keyword:KW-0998 OR keyword:KW-0964) "
            "AND NOT keyword:KW-0325 " + _EXCLUDE_TAXA
        ),
        "rationale": (
            "Periplasmic, outer-membrane or secreted bacterial proteins from clades "
            "with no known N-glycosylation machinery. Crosses a membrane and folds "
            "in an oxidising compartment like the ER lumen; differs in taxonomy."
        ),
    },
}

FIELDS = "accession,id,organism_name,organism_id,length,sequence,xref_pdb"


def fetch_control_proteins(
    name: str, limit: int | None = None, delay: float = 0.2, timeout: int = 120
) -> pd.DataFrame:
    """Download one control set from UniProt, following pagination."""
    if name not in CONTROL_SETS:
        raise KeyError(f"unknown control set {name!r}; expected {sorted(CONTROL_SETS)}")

    params = {
        "query": CONTROL_SETS[name]["query"],
        "fields": FIELDS,
        "format": "tsv",
        "size": "500",
    }
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"

    frames, header, rows = [], None, []
    while url:
        request = urllib.request.Request(url, headers={"Accept": "text/plain"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            text = response.read().decode("utf-8")
            link = response.headers.get("Link", "")

        lines = text.splitlines()
        if not lines:
            break
        if header is None:
            header = lines[0].split("\t")
        rows.extend(line.split("\t") for line in lines[1:] if line.strip())

        if limit is not None and len(rows) >= limit:
            rows = rows[:limit]
            break

        match = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = match.group(1) if match else None
        time.sleep(delay)

    frame = pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header or [])
    if "Entry" in frame.columns:
        # Paginated responses can repeat an entry across page boundaries. A repeat
        # would emit that protein's sequons twice and then square in any join.
        frame = frame.drop_duplicates("Entry")
    frame["control_set"] = name
    return frame


def find_sequons(sequence: str) -> list[int]:
    """1-indexed asparagine positions of every N-X-S/T motif (X is not proline).

    Uses a lookahead-free scan with manual stepping so overlapping motifs are not
    lost: in NNSS the asparagine at 1 and at 2 both open a valid sequon.
    """
    positions = []
    for index in range(len(sequence) - 2):
        if sequence[index] != "N":
            continue
        if sequence[index + 1] == "P":
            continue
        if sequence[index + 2] in ("S", "T"):
            positions.append(index + 1)
    return positions


def build_control_sites(proteins: pd.DataFrame) -> pd.DataFrame:
    """One row per sequon in the control proteins."""
    rows = []
    for record in proteins.itertuples(index=False):
        sequence = str(getattr(record, "Sequence", "") or "")
        if not sequence:
            continue
        pdb_ids = [
            x.strip() for x in str(getattr(record, "PDB", "") or "").split(";") if x.strip()
        ]
        for position in find_sequons(sequence):
            rows.append({
                "control_set": record.control_set,
                "accession": getattr(record, "Entry", ""),
                "position": position,
                "organism": getattr(record, "Organism", ""),
                "organism_id": getattr(record, "Organism__ID_", ""),
                "protein_length": len(sequence),
                "sequon": sequence[position - 1: position + 2],
                "n_pdb_entries": len(pdb_ids),
                "pdb_ids": ";".join(pdb_ids[:20]),
            })

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return (
        frame.drop_duplicates(["control_set", "accession", "position"])
        .sort_values(["control_set", "accession", "position"])
        .reset_index(drop=True)
    )


def composition(sequences: list[str]) -> dict[str, float]:
    """Amino acid usage, so sequon density differences can be attributed properly.

    Bacterial and eukaryotic proteomes differ in composition, which changes how
    often N-X-S/T arises by chance alone — nothing to do with glycosylation.
    """
    counts: dict[str, int] = {}
    total = 0
    for sequence in sequences:
        for residue in sequence:
            counts[residue] = counts.get(residue, 0) + 1
            total += 1
    return {k: round(v / total, 5) for k, v in sorted(counts.items())} if total else {}


def summarise(sites: pd.DataFrame, proteins: pd.DataFrame) -> dict:
    """Per-set counts, sequon density and composition, for the provenance record."""
    out: dict[str, dict] = {}
    for name, group in sites.groupby("control_set"):
        source = proteins[proteins.control_set == name]
        sequences = [str(s) for s in source.get("Sequence", pd.Series(dtype=str)) if s]
        residues = sum(len(s) for s in sequences)
        out[name] = {
            "rationale": CONTROL_SETS[name]["rationale"],
            "proteins": int(source.shape[0]),
            "sequons": int(group.shape[0]),
            "sequons_per_protein": round(group.shape[0] / max(source.shape[0], 1), 3),
            "sequons_per_1000_residues": round(1000 * group.shape[0] / residues, 3) if residues else None,
            "proteins_with_structure": int((group.groupby("accession").n_pdb_entries.first() > 0).sum()),
            "composition": composition(sequences),
        }
    out["_excluded_taxa"] = OST_BEARING_TAXA
    return out


def control_structure_targets(
    sites: pd.DataFrame,
    proteins_per_set: dict[str, int | None] | None = None,
    seed: int = 0,
) -> dict[str, set[str]]:
    """One deposited structure per control protein, optionally subsampled.

    Feature extraction needs coordinates, and one entry per protein is enough:
    the residue either resolves in it or it does not. Fetching every
    cross-referenced entry for 7,499 proteins would cost hours and gigabytes for
    no gain in matching power, since a few thousand controls already dwarf the
    332 occupied sites they are matched against.

    Subsampling is by protein and seeded, so the selection is reproducible and
    does not depend on which sequons a protein happens to carry.
    """
    proteins_per_set = proteins_per_set or {}
    wanted: dict[str, set[str]] = {}

    for name, group in sites.groupby("control_set"):
        first_entry = (
            group[group.pdb_ids.fillna("") != ""]
            .groupby("accession")
            .pdb_ids.first()
            .str.split(";")
            .str[0]
        )
        cap = proteins_per_set.get(name)
        if cap is not None and len(first_entry) > cap:
            first_entry = first_entry.sample(n=cap, random_state=seed)
        for accession, pdb_id in first_entry.items():
            if pdb_id:
                wanted[accession] = {pdb_id.strip()}

    return wanted
