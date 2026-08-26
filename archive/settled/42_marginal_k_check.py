"""How many samples does ESM-IF's marginalised joint score need?

`marginalised_probabilities` estimates a 400-term sum by sampling the hidden
residues from the model's own belief. The sample count is a free parameter, and
picking one without checking is how an arbitrary number ends up in a methods
section. This measures three things and prints them:

  * **convergence** — does the estimate stop moving as K grows?
  * **seed-to-seed spread** — two seeds at the same K should agree better as K
    rises; if they do not, the estimator is not doing what it claims.
  * **the artefact check** — substituting `<mask>` moved 93% of the probability
    at +2 onto aromatics against 0.3% natively. Marginalisation should stay near
    the native value. If it does not, it has not escaped the problem it exists
    to solve and the approach is wrong.

Usage:  42_marginal_k_check.py [--sites 8] [--k 4,8,16,32] [--device cpu]
"""
import argparse, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites import esmif_scoring as E
from experimental_glycosylation_sites.runner_support import (resolve_device,
                                                             structure_paths)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sites", type=int, default=8)
parser.add_argument("--k", default="4,8,16,32")
parser.add_argument("--device", default="cpu")
parser.add_argument("--manifest", default="results/manifests/candidate_manifest_dataset.csv")
parser.add_argument("--out", default="results/analysis/marginal_k_check.csv")
args = parser.parse_args()

Ks = [int(k) for k in args.k.split(",")]
device = resolve_device(args.device)
model, alphabet = E.load_model(device)
paths = structure_paths()

manifest = pd.read_csv(args.manifest, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
manifest = manifest.drop_duplicates(["accession", "position"]).head(args.sites)

i_s, i_t = alphabet.get_idx("S"), alphabet.get_idx("T")
aromatic = [alphabet.get_idx(a) for a in "FYWH"]

rows, t0 = [], time.time()
for r in manifest.itertuples(index=False):
    path = paths.get(str(r.structure_pdb_id).upper())
    if path is None:
        continue
    try:
        mapping = E.chain_mapping(path, r.structure_chain_id, str(r.structure_pdb_id))
        idx = mapping.map_indices((int(r.n_model_index), int(r.plus1_model_index),
                                   int(r.plus2_model_index)))
    except Exception:
        continue
    base = E.conditional_probabilities(mapping, model, alphabet, device)
    record = {"site": f"{r.accession}:{r.position}",
              "chain_length": len(mapping.esm_seq),
              "native_pst": float(base[idx[2], i_s] + base[idx[2], i_t]),
              "native_aromatic": float(base[idx[2], aromatic].sum())}
    for k in Ks:
        draws = []
        for seed in (0, 1):
            probs, _ = E.marginalised_probabilities(
                mapping, model, alphabet, *idx, device=device, n_samples=k, seed=seed)
            draws.append(float(probs[idx[2], i_s] + probs[idx[2], i_t]))
            if k == max(Ks) and seed == 0:
                record["aromatic_marginal"] = float(probs[idx[2], aromatic].sum())
        record[f"k{k}"] = float(np.mean(draws))
        record[f"k{k}_spread"] = float(abs(draws[0] - draws[1]))
    rows.append(record)
    print(f"  {record['site']:22} len {record['chain_length']:4d}  "
          f"({time.time()-t0:.0f}s)", flush=True)

frame = pd.DataFrame(rows)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
frame.to_csv(args.out, index=False)
largest = max(Ks)

print(f"\n{len(frame)} sites in {(time.time()-t0)/60:.1f} min  ->  {args.out}\n")
print("P(S/T) at +2, hidden residues marginalised")
print(f"  {'native':>10}" + "".join(f"{f'K={k}':>10}" for k in Ks))
print(f"  {frame.native_pst.mean():>10.4f}"
      + "".join(f"{frame[f'k{k}'].mean():>10.4f}" for k in Ks))

print("\nseed-to-seed spread at the same K (should fall as K rises)")
for k in Ks:
    print(f"  K={k:<4} mean {frame[f'k{k}_spread'].mean():.4f}   max {frame[f'k{k}_spread'].max():.4f}")

print(f"\nconvergence, |K - K{largest}| relative to K{largest}")
for k in Ks[:-1]:
    rel = (frame[f"k{k}"] - frame[f"k{largest}"]).abs() / frame[f"k{largest}"].clip(lower=1e-6)
    print(f"  K={k:<4} mean {rel.mean():.3f}   max {rel.max():.3f}")

print("\nartefact check, aromatic mass at +2")
print(f"  native {frame.native_aromatic.mean():.4f}   "
      f"marginalised {frame.aromatic_marginal.mean():.4f}   "
      f"(<mask> substitution gave 0.926)")
if frame.aromatic_marginal.mean() > 5 * max(frame.native_aromatic.mean(), 1e-4):
    print("  WARNING: still elevated -- marginalisation has not escaped the artefact")
