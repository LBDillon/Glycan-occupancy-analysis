"""The scoring manifest: one definitive row per site to be scored.

A model is asked about three specific residues, so all three must be located
exactly in the structure the model will see. Anything less is excluded and
reported rather than approximated: a site whose asparagine is resolved but whose
+1 or +2 residue is missing cannot be scored as a sequon, and silently scoring it
would put a two-residue observation into a three-residue benchmark.

The manifest also records the model index of each residue — its ordinal position
in the chain as the model reads it — because that, not the UniProt position or
the author residue number, is what a scoring call has to address.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .structures import _alignment_pairs, _parse_chains

SEQUON_PATTERN = re.compile(r"^N[^P][ST]$")

FAILURE_REASONS = (
    "structure_unreadable",
    "chain_absent",
    "n_not_mapped",
    "plus1_not_mapped",
    "plus2_not_mapped",
    "triplet_not_a_sequon",
)


def map_sequon(
    uniprot_sequence: str,
    position: int,
    structure_path: Path,
    pdb_id: str,
    chain_id: str,
) -> dict:
    """Locate all three sequon residues in one chain, or explain why not.

    The chain is fixed by the caller rather than chosen here: it must be the
    chain the site's matching features were computed from, or the manifest would
    describe a different residue environment than the one the sites were matched
    on.
    """
    try:
        chains = _parse_chains(Path(structure_path), pdb_id)
    except Exception:
        return {"ok": False, "reason": "structure_unreadable"}

    chain = next((c for c in chains if c.chain_id == chain_id), None)
    if chain is None:
        return {"ok": False, "reason": "chain_absent"}

    mapped = {u: c for u, c in _alignment_pairs(uniprot_sequence, chain.sequence)}

    residues = []
    for offset, label in enumerate(("n", "plus1", "plus2")):
        index = mapped.get(position + offset)
        if index is None or not 1 <= index <= len(chain.residue_ids):
            return {"ok": False, "reason": f"{label}_not_mapped"}
        resseq, icode = chain.residue_ids[index - 1]
        residues.append({
            "model_index": index - 1,          # 0-based, as a model reads the chain
            "resseq": int(resseq),
            "icode": icode,
            "aa": chain.sequence[index - 1],
        })

    triplet = "".join(r["aa"] for r in residues)
    if not SEQUON_PATTERN.match(triplet):
        return {"ok": False, "reason": "triplet_not_a_sequon", "triplet": triplet}

    return {
        "ok": True,
        "triplet": triplet,
        "subtype": f"NX{triplet[2]}",
        "chain_length_model": len(chain.sequence),
        **{
            f"{label}_{field}": residues[i][field]
            for i, label in enumerate(("n", "plus1", "plus2"))
            for field in ("model_index", "resseq", "icode", "aa")
        },
    }


def build_manifest(
    rows: pd.DataFrame,
    sequences: dict[str, str],
    structure_paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (manifest, exclusions).

    `rows` needs accession, position, structure_pdb_id, structure_chain_id and
    whatever provenance columns should travel through. Every input row appears in
    exactly one of the two outputs, so the pair is a complete account of the input.
    """
    kept, dropped = [], []

    for row in rows.itertuples(index=False):
        record = row._asdict()
        accession = str(record.get("accession", ""))
        pdb_id = str(record.get("structure_pdb_id") or "").upper()
        chain_id = str(record.get("structure_chain_id") or "")
        sequence = sequences.get(accession, "")
        path = structure_paths.get(pdb_id)

        if not sequence:
            dropped.append({**record, "exclusion_reason": "no_uniprot_sequence"})
            continue
        if path is None:
            dropped.append({**record, "exclusion_reason": "structure_not_cached"})
            continue

        result = map_sequon(sequence, int(record["position"]), path, pdb_id, chain_id)
        if not result.pop("ok"):
            dropped.append({
                **record,
                "exclusion_reason": result["reason"],
                "observed_triplet": result.get("triplet", ""),
            })
            continue
        kept.append({**record, **result})

    manifest = pd.DataFrame(kept)
    exclusions = pd.DataFrame(dropped)
    if not manifest.empty:
        manifest = manifest.sort_values(["accession", "position"]).reset_index(drop=True)
    return manifest, exclusions


def within_structure_pairs(manifest: pd.DataFrame) -> pd.DataFrame:
    """Occupied and observed-unmodified sites sharing a protein or a structure.

    The strongest available sensitivity subset: protein identity, expression
    system, crystallisation conditions and refinement are all held constant, so a
    score difference cannot be attributed to any of them.
    """
    occupied = manifest[manifest.occupancy_status == "occupied_supported"]
    unmodified = manifest[manifest.occupancy_status == "observed_unmodified"]
    if occupied.empty or unmodified.empty:
        return pd.DataFrame(columns=[
            "accession", "occupied_position", "unmodified_position",
            "occupied_pdb", "unmodified_pdb", "same_protein", "same_structure",
        ])

    # A pair sharing a protein almost always shares its structure too, so listing
    # per level would count the same pair twice. Emit one row per distinct pair and
    # flag which controls hold, which is what a sensitivity analysis needs to know.
    merged = occupied.merge(
        unmodified, on="accession", suffixes=("_occ", "_unm"), how="inner"
    )
    rows = [
        {
            "accession": r.accession,
            "occupied_position": r.position_occ,
            "unmodified_position": r.position_unm,
            "occupied_pdb": r.structure_pdb_id_occ,
            "unmodified_pdb": r.structure_pdb_id_unm,
            "same_protein": True,
            "same_structure": r.structure_pdb_id_occ == r.structure_pdb_id_unm,
        }
        for r in merged.itertuples(index=False)
    ]

    # Different proteins that happen to be solved in the same entry: experimental
    # context is shared even though protein identity is not.
    cross = occupied.merge(
        unmodified, on="structure_pdb_id", suffixes=("_occ", "_unm"), how="inner"
    )
    rows += [
        {
            "accession": f"{r.accession_occ}|{r.accession_unm}",
            "occupied_position": r.position_occ,
            "unmodified_position": r.position_unm,
            "occupied_pdb": r.structure_pdb_id,
            "unmodified_pdb": r.structure_pdb_id,
            "same_protein": False,
            "same_structure": True,
        }
        for r in cross.itertuples(index=False)
        if r.accession_occ != r.accession_unm
    ]

    frame = pd.DataFrame(rows)
    return frame.drop_duplicates() if not frame.empty else frame
