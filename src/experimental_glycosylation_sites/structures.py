from __future__ import annotations

import csv
import time
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

# A chain counts as "the same protein" when its alignment identity is within this
# much of the best-matching chain's. Copies in a homodimer align at essentially
# identical identity while differing in coverage, so identity — not coverage — is
# what distinguishes a genuine copy from an unrelated chain that local alignment
# happened to match.
SAME_PROTEIN_IDENTITY_TOLERANCE = 0.02

# Identity alone is not enough: it is measured over the aligned block, so an
# unrelated chain sharing a short motif scores identity 1.0 on a two-residue block
# and would pass the gate. A chain that genuinely is this protein aligns across
# most of itself, so require that too before trusting a glycan linkage from it.
MIN_CHAIN_COVERAGE = 0.5
MIN_ALIGNED_RESIDUES = 30


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
    scored: list[tuple[tuple[float, float, float, int], bool, dict]] = []

    for chain in chains:
        pairs = _alignment_pairs(uniprot_sequence, chain.sequence)
        if not pairs:
            continue
        index = {u: c for u, c in pairs}.get(position)
        if index is None or not 1 <= index <= len(chain.residue_ids):
            continue

        matches = sum(
            1 for u, c in pairs if uniprot_sequence[u - 1] == chain.sequence[c - 1]
        )
        identity = matches / len(pairs)
        coverage_chain = len({c for _, c in pairs}) / len(chain.sequence)
        coverage_uniprot = len({u for u, _ in pairs}) / len(uniprot_sequence)

        resseq, icode = chain.residue_ids[index - 1]
        has_link = (chain.chain_id, resseq, icode) in linked
        scored.append((
            (identity, coverage_chain, coverage_uniprot, len(pairs)),
            has_link,
            {
                "tier": "structure_linked_glycan" if has_link else "structure_residue_resolved",
                "pdb_id": pdb_id, "chain_id": chain.chain_id,
                "resseq": resseq, "icode": icode,
                "observed_residue": chain.sequence[index - 1],
                # An N-linked site must land on an asparagine. Anything else is a
                # mapping failure — often an engineered N->Q sequon knockout — and
                # must not be read as "examined and found bare", which is exactly
                # what a blank detail would imply.
                "detail": (
                    ""
                    if chain.sequence[index - 1] == "N"
                    else f"residue_mismatch:{chain.sequence[index - 1]}"
                ),
            },
        ))

    if not scored:
        return {**blank, "tier": "structure_residue_unresolved", "detail": "position_not_in_model"}

    # Identity is the "is this the same protein" test; coverage is not. So gate on
    # identity, then prefer a linked chain, then rank the rest by coverage. But
    # admit only chains that align across enough of themselves to be credibly this
    # protein at all — otherwise a short perfect sub-block on an unrelated chain
    # passes at identity 1.0 and its glycan is wrongly credited.
    credible = [
        item for item in scored
        if item[0][1] >= MIN_CHAIN_COVERAGE and item[0][3] >= MIN_ALIGNED_RESIDUES
    ]

    if not credible:
        # Nothing aligns well enough to identify. Report the best positional guess,
        # but never assert a glycan from a chain we cannot credibly call this protein.
        scored.sort(key=lambda item: item[0], reverse=True)
        fallback = dict(scored[0][2])
        fallback["tier"] = "structure_residue_resolved"
        fallback["detail"] = "low_confidence_chain_match"
        return fallback

    best_identity = max(item[0][0] for item in credible)
    same_protein = [
        item for item in credible
        if item[0][0] >= best_identity - SAME_PROTEIN_IDENTITY_TOLERANCE
    ]
    same_protein.sort(key=lambda item: item[0], reverse=True)
    for _, has_link, candidate in same_protein:
        if has_link:
            return candidate
    return same_protein[0][2]


class StructureManifestError(RuntimeError):
    """The manifest lists structures but none of them could be located."""


def load_manifest(path: Path, structure_dir: Path | None = None) -> dict[str, dict]:
    """Accession to manifest row, keeping only rows whose file can be located.

    The manifest stores absolute `output_path` values recorded on the machine that
    downloaded the structures. On any other checkout those paths do not exist, so
    each row is also tried as a basename under `structure_dir`. Without that
    fallback the whole structural layer silently disappears — every lookup misses,
    the run completes, and it reports a smaller site count as though that were the
    answer.

    A manifest with usable rows but zero resolvable files is a broken setup, not an
    empty result, and raises rather than returning {}.
    """
    manifest: dict[str, dict] = {}
    considered = 0
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("accession") or "").strip()
            output_path = (row.get("output_path") or "").strip()
            if not accession or not output_path or (row.get("status") or "").strip() == "failed":
                continue
            considered += 1

            located = None
            recorded = Path(output_path)
            if recorded.exists() and recorded.stat().st_size > 0:
                located = recorded
            elif structure_dir is not None:
                fallback = Path(structure_dir) / recorded.name
                if fallback.exists() and fallback.stat().st_size > 0:
                    located = fallback

            if located is not None:
                row = dict(row)
                row["output_path"] = str(located)
                manifest[accession] = row

    if considered and not manifest:
        raise StructureManifestError(
            f"{path}: {considered} usable rows but none could be located. "
            f"Recorded paths are absolute and this is not the machine that wrote them; "
            f"set [paths] structure_dir to the directory holding the structure files "
            f"(tried: {structure_dir})."
        )
    return manifest


# Best-answer precedence when a site is examined in several structures. A glycan
# seen in any structure is an existence claim and outranks every silence; "not
# assessed" is the weakest because it records only that we did not look.
TIER_RANK = {
    "structure_linked_glycan": 3,
    "structure_residue_resolved": 2,
    "structure_residue_unresolved": 1,
    "structure_not_assessed": 0,
}


def extra_structure_paths(accession: str, pdb_ids: set[str], cache_dir: Path) -> list[Path]:
    """Locally cached structures for an accession beyond the manifest's one entry."""
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return []
    found = []
    for pdb_id in sorted(pdb_ids):
        path = cache_dir / f"{pdb_id.upper()}.pdb"
        if path.exists() and path.stat().st_size > 0:
            found.append(path)
    return found


def fetch_structures(
    wanted: dict[str, set[str]],
    cache_dir: Path,
    delay: float = 0.34,
    timeout: int = 60,
    per_accession_cap: int = 20,
) -> dict[str, int]:
    """Download PDB entries from RCSB into a module-local cache.

    `wanted` maps accession to the PDB ids worth fetching. Resumable: anything
    already cached is skipped. Never writes outside `cache_dir`, so the ortholog
    database's own structure cache stays read-only.
    """
    import urllib.error
    import urllib.request

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    stats = {"requested": 0, "downloaded": 0, "cached": 0, "failed": 0}

    for accession in sorted(wanted):
        for pdb_id in sorted(wanted[accession])[:per_accession_cap]:
            target = cache_dir / f"{pdb_id.upper()}.pdb"
            stats["requested"] += 1
            if target.exists() and target.stat().st_size > 0:
                stats["cached"] += 1
                continue
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            try:
                with urllib.request.urlopen(url, timeout=timeout) as response:
                    payload = response.read()
            except (urllib.error.URLError, TimeoutError, OSError):
                # Large entries are mmCIF-only and have no .pdb form; that is data,
                # not a bug, so record and continue rather than aborting the sweep.
                stats["failed"] += 1
                continue
            target.write_bytes(payload)
            stats["downloaded"] += 1
            time.sleep(delay)
    return stats


def build_site_evidence(
    candidates: pd.DataFrame,
    sequences: dict[str, str],
    manifest: dict[str, dict],
    extra_cache_dir: Path | None = None,
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
            paths = [(Path(entry["output_path"]), str(entry.get("pdb_id", "")))]
            if extra_cache_dir is not None:
                ids = {
                    x.strip()
                    for x in str(entry.get("all_pdb_ids") or "").split(";")
                    if x.strip()
                }
                ids.discard(str(entry.get("pdb_id", "")).strip())
                paths += [
                    (p, p.stem)
                    for p in extra_structure_paths(accession, ids, extra_cache_dir)
                ]

            # Examine every structure available for this protein and keep the
            # strongest answer. A glycan modelled in any one of them is a real
            # observation; seeing none in a particular crystal is not evidence.
            result = None
            for path, pdb_id in paths:
                key = str(path)
                if key not in link_cache:
                    link_cache[key] = parse_link_records(path)
                candidate = assess_site(sequence, position, path, pdb_id, link_cache[key])
                if result is None or TIER_RANK[candidate["tier"]] > TIER_RANK[result["tier"]]:
                    result = candidate
                if result["tier"] == "structure_linked_glycan":
                    break
            result["n_structures_examined"] = len(paths)

        rows.append({
            "accession": accession,
            "position": position,
            "structure_tier": result["tier"],
            "structure_pdb_id": result.get("pdb_id", ""),
            "structure_chain_id": result.get("chain_id", ""),
            "structure_resseq": result.get("resseq"),
            "structure_icode": result.get("icode", ""),
            # Carried through so a consumer can detect a mis-mapped residue: an
            # N-linked site must map to an asparagine, and anything else is a
            # mapping failure that would otherwise be invisible.
            "structure_observed_residue": result.get("observed_residue", ""),
            "structure_n_examined": result.get("n_structures_examined", 0),
            "structure_detail": result.get("detail", ""),
        })

    return pd.DataFrame(rows).sort_values(["accession", "position"]).reset_index(drop=True)
