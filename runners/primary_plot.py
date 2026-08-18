"""The primary figure: every matched pair, and the difference they average to.

A slope graph rather than two box plots. With sixteen pairs the individual
pairings are the evidence, and a summary that hides them would imply more
precision than sixteen contrasts can carry. The right-hand panel shows the mean
paired difference against the pre-specified equivalence margin, which is what
the verdict is read from.
"""
import json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LABEL = sys.argv[1] if len(sys.argv) > 1 else "primary"
contrasts = pd.read_csv(f"results/contrasts_{LABEL}.csv")
stats = json.loads(Path(f"results/analysis_{LABEL}.json").read_text())

INK, OCC, CTL = "#22252b", "#1f6f8b", "#b05c3b"
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 5.6),
                             gridspec_kw={"width_ratios": [1.45, 1]})

# ---- left: every pair -------------------------------------------------------
order = contrasts.sort_values("contrast").reset_index(drop=True)
for i, r in order.iterrows():
    up = r.contrast > 0
    ax.plot([0, 1], [r.control_mean_score, r.case_score],
            color=(OCC if up else CTL), alpha=0.75, lw=1.4, zorder=2)
ax.scatter(np.zeros(len(order)), order.control_mean_score, s=34, color=CTL,
           zorder=3, label="no modelled glycan\n(internal control)")
ax.scatter(np.ones(len(order)), order.case_score, s=34, color=OCC,
           zorder=3, label="experimentally occupied")
ax.set_xlim(-0.42, 1.42); ax.set_xticks([0, 1])
ax.set_xticklabels(["internal\ncontrol", "occupied"], fontsize=10)
ax.set_ylabel("conditional sequon score (log-odds)", fontsize=10)
ax.set_title(f"{len(order)} matched pairs", fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
for side in ("top", "right"): ax.spines[side].set_visible(False)

# ---- right: the difference, against the margin ------------------------------
lo, hi = stats["ci95_raw"]; mean = stats["mean_difference_raw"]
margin = stats["margin_raw"]
bx.axhspan(-margin, margin, color="#9aa3ad", alpha=0.16, zorder=1)
bx.axhline(0, color=INK, lw=0.9, zorder=2)
bx.scatter(np.random.default_rng(0).normal(0.62, 0.035, len(order)), order.contrast,
           s=26, color="#9aa3ad", alpha=0.85, zorder=3)
bx.errorbar([0.18], [mean], yerr=[[mean - lo], [hi - mean]], fmt="o", color=INK,
            capsize=5, lw=1.8, markersize=7, zorder=4)
bx.set_xlim(-0.1, 0.95); bx.set_xticks([0.18, 0.62])
bx.set_xticklabels(["mean\n(95% CI)", "per pair"], fontsize=9.5)
bx.set_ylabel("occupied − control (log-odds)", fontsize=10)
bx.set_title(f"{stats['verdict']}", fontsize=11, loc="left")
bx.annotate(f"±0.2 SD margin\n(±{margin:.2f} log-odds)", xy=(0.88, margin),
            fontsize=8, color="#5c6570", ha="right", va="bottom")
for side in ("top", "right"): bx.spines[side].set_visible(False)

fig.suptitle("ProteinMPNN conditional sequon score: occupied vs matched internal controls",
             fontsize=12.5, x=0.02, ha="left", y=0.98)
rb = stats["robustness"]
fig.text(0.02, 0.055,
         f"ProteinMPNN v_48_020, 8 decoding orders. Matched on relative accessibility, "
         f"neighbour count and hydrophobic fraction, with NXS/NXT required to be identical.",
         fontsize=8, color="#5c6570")
fig.text(0.02, 0.030,
         f"Mean {mean:+.3f} log-odds ({stats['mean_difference_sd']:+.3f} SD); 95% CI "
         f"[{stats['ci95_sd'][0]:+.3f}, {stats['ci95_sd'][1]:+.3f}] SD from a cluster "
         f"bootstrap over {stats['n_resampling_units']} units of shared ortholog clusters "
         f"and control proteins.",
         fontsize=8, color="#5c6570")
fig.text(0.02, 0.005,
         f"Direction is inconsistent: occupied scores higher in {rb['n_positive']} of "
         f"{rb['n_contrasts']} pairs (sign test p={rb['sign_test_p']:.2f}); the effect is "
         f"carried by magnitude, not count (Wilcoxon p={rb['wilcoxon_p']:.3f}). "
         f"Provisional — 16 pairs.",
         fontsize=8, color="#8a5a3b")
fig.tight_layout(rect=[0, 0.085, 1, 0.955])
out = f"results/figures/primary_{LABEL}.png"
Path("results/figures").mkdir(exist_ok=True)
fig.savefig(out, dpi=200); print("wrote", out)
