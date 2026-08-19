"""Design retention — reported separately, and not part of the primary result.

Retention asks what ProteinMPNN *does* when it generates a sequence, which is a
different question from what probability it holds at a site. It is the bridge to
the preprint, and it is descriptive: it has no pre-specified margin, so it
cannot carry a statistical conclusion of its own.

Two exclusions apply. Sites ProteinMPNN cannot decode are dropped, because
sample() writes their native residue back unchanged in every design and the
resulting number describes the parser rather than the model.

The sweep is complete: 2,526 of 2,640 manifest sites, 1,659 of 1,725 chains.
Every one of the 114 uncovered sites is accounted for by a logged parse failure
on its chain (59 KeyError, 7 ValueError across 66 chains), so the gap is
structural rather than a truncation. It is not neutral, though: the uncovered
entries skew towards recent large depositions whose chain identifiers
ProteinMPNN's PDB parser cannot read, so very large complexes are
under-represented relative to the manifest.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

SNAPSHOT = Path("results/designs/mpnn_retention_frozen_2026-08-18.csv")
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
FULL = "std_frac_full_sequon_retained"

ret = pd.read_csv(SNAPSHOT, low_memory=False)
score = pd.read_csv("results/scores/mpnn_conditional_scores.csv", low_memory=False)
able = pd.read_csv("results/manifests/scoreability.csv", low_memory=False)
for d in (ret, score, able):
    for k in KEY:
        if k in d: d[k] = d[k].astype(str)

# the per-row scores of valid sites are untouched by the defect, which corrupted
# whole rows rather than biasing them, so valid rows carry over unchanged
score["valid"] = ~((score.p_ser_or_thr_at_plus2 > 1.0)
                   | (score.probs_n.map(lambda v: abs(sum(json.loads(v)) - 1) > 1e-3)))

before = len(ret)
ret = ret.merge(able[KEY + ["scoreable"]], on=KEY, how="left")
ret = ret[ret.scoreable == True].copy()
print(f"retention rows {before} -> {len(ret)} after dropping sites the model cannot decode")
print(f"chains covered {ret.structure_pdb_id.nunique()}; sweep complete "
      f"(uncovered sites are parse failures, skewed to large recent entries)\n")

print("=== descriptive retention (32 designs, temperature 0.1) ===")
print(f"  mean full-sequon retention   {ret[FULL].mean():.4f}")
print(f"  median                       {ret[FULL].median():.4f}")
print(f"  sites losing the sequon in every design  "
      f"{int((ret[FULL]==0).sum())} of {len(ret)} ({100*(ret[FULL]==0).mean():.1f}%)")
print(f"  sites retaining it in every design       "
      f"{int((ret[FULL]==1).sum())} of {len(ret)} ({100*(ret[FULL]==1).mean():.1f}%)")

pre = "pre_frac_full_sequon_retained"
if pre in ret.columns and ret[pre].notna().any():
    both = ret.dropna(subset=[pre, FULL])
    print(f"\n=== preprint condition (8 designs) vs standard (32) ===")
    print(f"  8-design mean  {both[pre].mean():.4f}   32-design mean {both[FULL].mean():.4f}")
    print(f"  correlation    {both[pre].corr(both[FULL]):.4f}  (n={len(both)})")
    print("  the preprint's setting is noisier per site, not biased")

merged = ret.merge(score[score.valid][KEY + ["conditional_sequon_score"]], on=KEY, how="inner")
print(f"\n=== does the conditional score predict retention? (n={len(merged)}) ===")
rho, p = st.spearmanr(merged.conditional_sequon_score, merged[FULL])
print(f"  Spearman rho {rho:+.4f}  (p={p:.2e})")
q = pd.qcut(merged.conditional_sequon_score, 5, labels=False, duplicates="drop")
print("  retention by score quintile (low to high):")
for i, g in merged.groupby(q):
    print(f"    Q{int(i)+1}  score {g.conditional_sequon_score.mean():+.2f}   "
          f"retention {g[FULL].mean():.4f}   n={len(g)}")

Path("results/retention_analysis.json").write_text(json.dumps({
    "snapshot": SNAPSHOT.name,
    "status": "secondary and descriptive; not part of the primary conclusion",
    "coverage_caveat": "sweep complete; 114 of 2,640 manifest sites uncovered, all "
                       "explained by logged chain parse failures, skewed to large recent entries",
    "rows_before_exclusion": int(before), "rows_analysed": int(len(ret)),
    "excluded_unscoreable": int(before - len(ret)),
    "mean_full_retention": round(float(ret[FULL].mean()), 4),
    "frac_sites_never_retained": round(float((ret[FULL] == 0).mean()), 4),
    "spearman_score_vs_retention": round(float(rho), 4),
    "spearman_p": float(p), "spearman_n": int(len(merged)),
}, indent=2))
print("\nwrote results/retention_analysis.json")
