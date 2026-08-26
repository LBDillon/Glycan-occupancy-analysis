"""Are occupied sequons lost more often than the chain around them?

Unconstrained redesign: nothing is held fixed, so the sequon is free to
disappear. Three quantities from the same designs, which is the point — they are
only comparable if they come from one run:

    sequon        does the N-X-S/T pattern survive at its own position?
    background    what fraction of all positions change at all?
    control       does an arbitrary three-residue window survive intact?

The control matters because a sequon is three residues, and three consecutive
residues surviving redesign is unlikely for reasons that have nothing to do with
glycosylation. Comparing sequon loss against the per-residue mutation rate would
overstate the case; comparing it against other triplets does not.

Both an exact and a pattern reading are reported. Retention as the benchmark
defines it is a pattern — asparagine, not proline, then serine or threonine — so
a sequon can survive while all three residues change. An arbitrary triplet has
no such latitude, so the exact comparison is the like-for-like one.

Usage:
    56_sequon_retention_rate.py [--proteins 250] [--designs 32]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from experimental_glycosylation_sites.mpnn_scoring import (chain_mapping, load_model,
                                                           to_manifest_space)
from experimental_glycosylation_sites.retention import design_sequences
from experimental_glycosylation_sites.runner_support import (proteinmpnn_dir,
                                                             structure_paths)
from experimental_glycosylation_sites.structures import _parse_chains
from glyco_context.fixed_design import verify_sequon_index
from experimental_glycosylation_sites.provenance import _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--indices", default="results/manifests/candidate_manifest_dataset.csv")
parser.add_argument("--proteins", type=int, default=250)
parser.add_argument("--designs", type=int, default=32)
parser.add_argument("--temperature", type=float, default=0.1)
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--out", default="glyco_context/results/analysis/sequon_retention_rate.csv")
args = parser.parse_args()

manifest = pd.read_csv(args.indices, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
manifest = manifest.dropna(subset=["n_model_index"])
groups = list(manifest.groupby(["structure_pdb_id", "structure_chain_id"]))[: args.proteins]
print(f"{len(manifest)} occupied sites in {len(groups)} chains")

paths = structure_paths(())
model = load_model(proteinmpnn_dir(), device="cpu")
rows, dropped, t0 = [], [], time.time()

for n, ((pdb_id, chain_id), group) in enumerate(groups, start=1):
    path = paths.get(str(pdb_id).upper())
    if path is None:
        dropped.append({"pdb": pdb_id, "reason": "structure_not_cached"}); continue
    try:
        mapping, _ = chain_mapping(path, str(chain_id), pdb_id)
        native = next((c for c in _parse_chains(path, str(pdb_id))
                       if c.chain_id == str(chain_id)), None)
        if not mapping or native is None:
            raise KeyError("chain cannot be mapped onto ProteinMPNN's parse")
        wild = native.sequence
    except Exception as exc:
        dropped.append({"pdb": pdb_id, "reason": f"{type(exc).__name__}"}); continue

    sites = []
    for site in group.itertuples():
        n_index = int(site.n_model_index)
        if verify_sequon_index(wild, n_index, str(site.triplet)):
            sites.append((site.accession, int(site.position), n_index))
    if not sites:
        continue

    # Nothing fixed: this is the design condition the retention benchmark uses.
    designs = [to_manifest_space(d, mapping) for d in
               design_sequences(path, str(chain_id), model, n_designs=args.designs,
                                temperature=args.temperature, seed=args.seed)]

    sequon_positions = {i for _, _, k in sites for i in (k, k + 1, k + 2)}
    # Every three-residue window that does not touch a sequon: the baseline for
    # "three consecutive residues survive", measured on the same designs.
    windows = [i for i in range(len(wild) - 2)
               if not sequon_positions & {i, i + 1, i + 2}]

    for design in designs:
        if len(design) != len(wild):
            continue
        identical = [a == b for a, b in zip(wild, design)]
        background = 1.0 - float(np.mean(identical))
        control_exact = float(np.mean([identical[i] and identical[i + 1] and identical[i + 2]
                                       for i in windows])) if windows else np.nan
        for accession, position, k in sites:
            trip = design[k:k + 3]
            rows.append({
                "accession": accession, "position": position,
                "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                "wild_triplet": wild[k:k + 3], "design_triplet": trip,
                "sequon_exact": trip == wild[k:k + 3],
                "sequon_pattern": (len(trip) == 3 and trip[0] == "N"
                                   and trip[2] in ("S", "T") and trip[1] != "P"),
                "background_mutation_rate": background,
                "control_triplet_exact": control_exact,
                "n_control_windows": len(windows),
            })
    if n % 25 == 0:
        print(f"  {n}/{len(groups)} chains, {len(rows)} rows, "
              f"{(time.time()-t0)/60:.1f} min", flush=True)

frame = pd.DataFrame(rows)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
frame.to_csv(args.out, index=False)
print(f"\n{len(frame)} design-site rows over "
      f"{frame.groupby(['accession','position']).ngroups} sites -> {Path(args.out).resolve()}")
if dropped:
    print(f"dropped chains: {len(dropped)}")
