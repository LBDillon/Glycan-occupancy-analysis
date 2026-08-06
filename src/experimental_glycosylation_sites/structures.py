from __future__ import annotations

import csv
import warnings
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import seq1

GLYCAN_RESNAMES = {"NAG", "NDG", "BGC", "GLC", "MAN", "BMA", "FUC", "GAL", "XYS"}
MMCIF_SUFFIXES = {".cif", ".mmcif"}


@dataclass(frozen=True)
class GlycanLink:
    chain_id: str
    resseq: int
    icode: str
    glycan_resname: str


@dataclass(frozen=True)
class ChainData:
    chain_id: str
    sequence: str
    residue_ids: list[tuple[int, str]]


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner


def _alignment_pairs(seq_a: str, seq_b: str) -> list[tuple[int, int]]:
    """1-indexed (position_in_a, position_in_b) pairs from the best local alignment."""
    if not seq_a or not seq_b:
        return []
    try:
        alignment = _aligner().align(seq_a, seq_b)[0]
    except (ValueError, IndexError):
        return []
    pairs = []
    for (a_start, a_end), (b_start, b_end) in zip(alignment.aligned[0], alignment.aligned[1]):
        for offset in range(a_end - a_start):
            pairs.append((a_start + offset + 1, b_start + offset + 1))
    return pairs


def _one_letter(resname: str) -> str:
    try:
        return seq1(resname, custom_map={"MSE": "M", "SEC": "U", "PYL": "O"}, undef_code="X")
    except (TypeError, KeyError):
        return "X"


def parse_link_records(path: Path) -> list[GlycanLink]:
    """Asparagine-side coordinates of every ASN-glycan LINK record.

    A LINK between an asparagine and a glycan residue is direct physical
    evidence that the site carried a glycan in that structure. The absence of
    such a record is NOT evidence that the site was unmodified.
    """
    path = Path(path)
    if not path.exists() or path.suffix.lower() in MMCIF_SUFFIXES:
        return []

    links: list[GlycanLink] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("LINK"):
                continue
            line = line.rstrip("\n").ljust(60)
            first = (line[17:20].strip(), line[21].strip(), line[22:26].strip(), line[26].strip())
            second = (line[47:50].strip(), line[51].strip(), line[52:56].strip(), line[56].strip())

            for asn, glycan in ((first, second), (second, first)):
                if asn[0] != "ASN" or glycan[0] not in GLYCAN_RESNAMES:
                    continue
                if not asn[2].lstrip("-").isdigit():
                    continue
                links.append(GlycanLink(
                    chain_id=asn[1], resseq=int(asn[2]), icode=asn[3], glycan_resname=glycan[0],
                ))
                break
    return links


def _parse_chains(path: Path, pdb_id: str) -> list[ChainData]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = PDBParser(QUIET=True).get_structure(pdb_id, str(path))

    chains = []
    for chain in next(iter(structure)):
        sequence, residue_ids = [], []
        for residue in chain:
            if not is_aa(residue, standard=False) or "CA" not in residue:
                continue
            _, resseq, icode = residue.id
            sequence.append(_one_letter(residue.get_resname()))
            residue_ids.append((int(resseq), str(icode).strip()))
        if sequence:
            chains.append(ChainData(chain.id, "".join(sequence), residue_ids))
    return chains


def assess_site(
    uniprot_sequence: str,
    position: int,
    structure_path: Path,
    pdb_id: str,
    links: list[GlycanLink],
) -> dict:
    """Place one site on the structural resolution ladder.

    structure_residue_resolved means the residue is present in the model with no
    glycan linkage. It must never be read as "observed unmodified": glycans are
    routinely removed before crystallisation, expressed in bacterial systems, or
    left unmodelled through disorder.
    """
    blank = {"tier": "structure_not_assessed", "pdb_id": pdb_id, "chain_id": "",
             "resseq": None, "icode": "", "observed_residue": "", "detail": ""}
    structure_path = Path(structure_path)

    if not structure_path.exists():
        return {**blank, "detail": "structure_file_missing"}
    if structure_path.suffix.lower() in MMCIF_SUFFIXES:
        return {**blank, "detail": "mmcif_linkage_unsupported"}

    try:
        chains = _parse_chains(structure_path, pdb_id)
    except Exception as exc:  # malformed structure files are data, not bugs
        return {**blank, "detail": f"structure_unreadable: {type(exc).__name__}"}
    if not chains:
        return {**blank, "detail": "no_protein_chain"}

    linked = {(link.chain_id, link.resseq, link.icode) for link in links}
    best = None
    for chain in chains:
        pairs = _alignment_pairs(uniprot_sequence, chain.sequence)
        mapped = {u: c for u, c in pairs}
        index = mapped.get(position)
        if index is None or not 1 <= index <= len(chain.residue_ids):
            continue
        resseq, icode = chain.residue_ids[index - 1]
        candidate = {
            "tier": "structure_residue_resolved", "pdb_id": pdb_id,
            "chain_id": chain.chain_id, "resseq": resseq, "icode": icode,
            "observed_residue": chain.sequence[index - 1], "detail": "",
        }
        if (chain.chain_id, resseq, icode) in linked:
            candidate["tier"] = "structure_linked_glycan"
            return candidate
        best = best or candidate

    if best is not None:
        return best
    return {**blank, "tier": "structure_residue_unresolved", "detail": "position_not_in_model"}


def load_manifest(path: Path) -> dict[str, dict]:
    """Accession to manifest row, keeping only rows whose file exists."""
    manifest: dict[str, dict] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("accession") or "").strip()
            output_path = (row.get("output_path") or "").strip()
            if not accession or not output_path or (row.get("status") or "").strip() == "failed":
                continue
            candidate = Path(output_path)
            if candidate.exists() and candidate.stat().st_size > 0:
                manifest[accession] = row
    return manifest


def build_site_evidence(
    candidates: pd.DataFrame,
    sequences: dict[str, str],
    manifest: dict[str, dict],
) -> pd.DataFrame:
    """One row per candidate site on the structural resolution ladder."""
    link_cache: dict[str, list[GlycanLink]] = {}
    rows = []

    for accession, position in zip(candidates["accession"], candidates["position"]):
        accession, position = str(accession), int(position)
        entry = manifest.get(accession)
        sequence = sequences.get(accession, "")

        if entry is None:
            result = {"tier": "structure_not_assessed", "pdb_id": "", "chain_id": "",
                      "resseq": None, "icode": "", "detail": "no_cached_structure"}
        elif not sequence:
            result = {"tier": "structure_not_assessed", "pdb_id": entry.get("pdb_id", ""),
                      "chain_id": "", "resseq": None, "icode": "",
                      "detail": "no_uniprot_sequence"}
        else:
            path = Path(entry["output_path"])
            key = str(path)
            if key not in link_cache:
                link_cache[key] = parse_link_records(path)
            result = assess_site(
                sequence, position, path, entry.get("pdb_id", ""), link_cache[key]
            )

        rows.append({
            "accession": accession,
            "position": position,
            "structure_tier": result["tier"],
            "structure_pdb_id": result.get("pdb_id", ""),
            "structure_chain_id": result.get("chain_id", ""),
            "structure_resseq": result.get("resseq"),
            "structure_icode": result.get("icode", ""),
            "structure_detail": result.get("detail", ""),
        })

    return pd.DataFrame(rows).sort_values(["accession", "position"]).reset_index(drop=True)
