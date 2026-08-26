"""Is the composition shift local to the sequon, or global to the chain?

The context-retention result says designs put more proline and glycine and fewer
aromatics near a protected sequon. That is only a statement about sequons if
ProteinMPNN is not doing the same thing everywhere. It has its own amino-acid
preferences, and the random control cannot separate them: the control draws
replacements from the wild-type chain's own frequencies, so it changes which
residues sit where without changing the overall mix.

This compares the same classes in two regions of the same designed chains:

    flank   the +-5 residues around each sequon, sequon itself excluded
    rest    every other designable position in the chain

Both are measured on the same sequences from the same run, so anything global to
ProteinMPNN appears in both and only a difference between them is local.

Usage:  55_composition_control.py
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.local_chemistry import AA_CLASS, CLASSES, flank_indices
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sequences",
                    default="glyco_context/results/analysis/fixed_sequon_sequences.csv")
parser.add_argument("--panels",
                    default="glyco_context/results/analysis/fixed_sequon_panels.csv")
parser.add_argument("--indices", default="results/manifests/candidate_manifest_dataset.csv")
parser.add_argument("--boot", type=int, default=2000)
parser.add_argument("--out", default="glyco_context/results/analysis")
args = parser.parse_args()

sequences = pd.read_csv(args.sequences, low_memory=False)
panels = pd.read_csv(args.panels, low_memory=False)
indices = pd.read_csv(args.indices, low_memory=False)
sites = panels[panels.variant == "wild_type"][
    ["accession", "position", "structure_pdb_id", "structure_chain_id"]].merge(
    indices[["accession", "position", "n_model_index"]].drop_duplicates(
        ["accession", "position"]), on=["accession", "position"], how="inner")
print(f"{len(sites)} sites on {sites.groupby(['structure_pdb_id','structure_chain_id']).ngroups} chains")


def composition(sequence, positions):
    residues = [sequence[i] for i in positions if 0 <= i < len(sequence)]
    if not residues:
        return {}
    n = len(residues)
    return {c: sum(1 for r in residues if AA_CLASS.get(r) == c) / n for c in CLASSES}


rows = []
for (pdb, chain), group in sites.groupby(["structure_pdb_id", "structure_chain_id"]):
    chain_seqs = sequences[(sequences.structure_pdb_id == pdb)
                           & (sequences.structure_chain_id == chain)]
    wt_rows = chain_seqs[chain_seqs.variant == "wild_type"]
    designs = chain_seqs[chain_seqs.variant == "design"].sequence.tolist()
    if not len(wt_rows) or not designs:
        continue
    wild = wt_rows.sequence.iloc[0]

    sequon, flank = set(), set()
    for site in group.itertuples():
        n_index = int(site.n_model_index)
        sequon.update({n_index, n_index + 1, n_index + 2})
        flank.update(flank_indices(n_index, len(wild)))
    flank -= sequon
    rest = set(range(len(wild))) - sequon - flank
    if len(flank) < 4 or len(rest) < 20:
        continue

    wild_flank, wild_rest = composition(wild, flank), composition(wild, rest)
    for design in designs:
        d_flank, d_rest = composition(design, flank), composition(design, rest)
        for c in CLASSES:
            rows.append({"accession": group.accession.iloc[0], "pdb": pdb, "chain": chain,
                         "amino_acid_class": c,
                         "flank_shift": d_flank.get(c, np.nan) - wild_flank.get(c, np.nan),
                         "rest_shift": d_rest.get(c, np.nan) - wild_rest.get(c, np.nan)})

shifts = pd.DataFrame(rows)
# Designs of one chain are replicates of one draw, so average within chain first.
per_chain = shifts.groupby(["accession", "pdb", "chain", "amino_acid_class"])[
    ["flank_shift", "rest_shift"]].mean().reset_index()
per_chain["local_excess"] = per_chain.flank_shift - per_chain.rest_shift


def bootstrap(values, clusters, n_boot, seed=0):
    values = values.dropna()
    if len(values) < 3:
        return np.nan, np.nan, np.nan, np.nan
    groups = [g.to_numpy(float) for _, g in values.groupby(clusters.loc[values.index])]
    rng = np.random.default_rng(seed)
    draws = np.array([np.concatenate([groups[i] for i in
                                      rng.integers(0, len(groups), len(groups))]).mean()
                      for _ in range(n_boot)])
    p = min(1.0, max(2 * min((draws <= 0).mean(), (draws >= 0).mean()), 1 / n_boot))
    return values.mean(), np.percentile(draws, 2.5), np.percentile(draws, 97.5), p


print(f"\n{'class':14}{'near sequon':>13}{'rest of chain':>15}{'local excess':>15}"
      f"{'95% CI':>20}{'p':>8}")
out = []
for c in CLASSES:
    sub = per_chain[per_chain.amino_acid_class == c]
    if len(sub) < 3:
        continue
    proteins = sub.accession
    flank_mean = sub.flank_shift.mean()
    rest_mean = sub.rest_shift.mean()
    mean, low, high, p = bootstrap(sub.local_excess, proteins, args.boot)
    star = "*" if p < 0.05 else " "
    print(f"{c:14}{flank_mean:+13.4f}{rest_mean:+15.4f}{mean:+15.4f}"
          f"  [{low:+7.4f},{high:+7.4f}]{p:8.4f}{star}")
    out.append({"amino_acid_class": c, "flank_shift": flank_mean, "rest_shift": rest_mean,
                "local_excess": mean, "ci_low": low, "ci_high": high, "p": p,
                "chains": int(len(sub))})

table = pd.DataFrame(out)
outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
table.to_csv(outdir / "composition_control.csv", index=False)
(outdir / "composition_control.json").write_text(json.dumps(
    {"n_boot": args.boot, "chains": int(per_chain.pdb.nunique()), "git": _git_state()},
    indent=2, default=str))
print(f"\nwrote {(outdir / 'composition_control.csv').resolve()}")
