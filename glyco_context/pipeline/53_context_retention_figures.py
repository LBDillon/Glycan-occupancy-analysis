"""Figures for the fixed-sequon context-retention test, from the scored table.

Two figures.

  fig4  the result: where each site's wild type sits relative to natural
        occupied context, where its designs sit, and how that compares to
        changing the same number of residues at random
  fig5  which features move, and by how much

The paired panel shows every site rather than a summary. With fifty of them the
individual pairings are the evidence, and the random control is the reason the
shift is readable at all -- a wild type sits inside the reference by
construction, so any perturbation moves outward.

Usage:  53_context_retention_figures.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5,
                     "axes.titlelocation": "left", "axes.spines.top": False,
                     "axes.spines.right": False})
INK, OCC, CTL, DIM = "#22252b", "#1f6f8b", "#b05c3b", "#9aa0a6"

ANALYSIS = Path("glyco_context/results/analysis")
OUT = Path("glyco_context/results/figures")
OUT.mkdir(parents=True, exist_ok=True)

scored = pd.read_csv(ANALYSIS / "context_retention_distances.csv")
site = scored.pivot_table(index=["accession", "position"], columns="variant",
                          values="distance", aggfunc="mean").dropna(subset=["wild_type"])
site = site.sort_values("wild_type").reset_index()

# =============================================================================
# fig4 — the result
# =============================================================================
fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.4, 5.8),
                             gridspec_kw={"width_ratios": [1.5, 1]})

for i, r in site.iterrows():
    ax.plot([0, 1], [r.wild_type, r.design], color=(OCC if r.design > r.wild_type else CTL),
            alpha=0.55, lw=1.2, zorder=2)
ax.scatter(np.zeros(len(site)), site.wild_type, s=30, color=INK, zorder=3,
           label="wild type")
ax.scatter(np.ones(len(site)), site.design, s=30, color=OCC, zorder=3,
           label="design (sequon held fixed)")
ax.set_xlim(-0.35, 1.35); ax.set_xticks([0, 1])
ax.set_xticklabels(["wild type", "design"])
ax.set_ylabel("distance from natural occupied context  (median |z| over the panel)")
moved = int((site.design > site.wild_type).sum())
ax.set_title(f"Every design keeps the sequon; {moved} of {len(site)} still move\n"
             "further from natural occupied context", fontsize=12)
ax.legend(frameon=False, fontsize=9, loc="upper left")

deltas = {"design\n- wild type": site.design - site.wild_type,
          "random control\n- wild type": site.random - site.wild_type,
          "design\n- random": site.design - site.random}
positions = np.arange(len(deltas))
for i, (name, values) in enumerate(deltas.items()):
    bx.scatter(np.random.default_rng(i).normal(i, 0.055, len(values)), values,
               s=16, color=DIM, alpha=0.6, zorder=2)
    mean = values.mean()
    bx.plot([i - 0.22, i + 0.22], [mean, mean], color=OCC, lw=2.6, zorder=4)
    bx.text(i + 0.27, mean, f"{mean:+.3f}", fontsize=9.5, color=INK, va="center")
bx.axhline(0, color=INK, lw=1.0, zorder=1)
bx.set_xticks(positions); bx.set_xticklabels(list(deltas), fontsize=9)
bx.set_ylabel("change in distance from natural context")
bx.set_title("Designs drift; arbitrary change of the\nsame size drifts about half as far",
             fontsize=12)
bx.text(0.0, -0.22, "Each point is one site. Bars are means. A wild type sits inside the "
        "reference by\nconstruction, so any change moves outward — which is what the "
        "random control measures.",
        transform=bx.transAxes, fontsize=8.6, color=INK, va="top")
fig.tight_layout()
fig.savefig(OUT / "fig4_context_retention.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig4_context_retention.png")

# =============================================================================
# fig5 — which features move
# =============================================================================
path = ANALYSIS / "context_retention_features.csv"
if path.exists():
    features = pd.read_csv(path)
    col = "shift" if "shift" in features.columns else "mean"
    features = features.reindex(features[col].abs().sort_values().index)
    pretty = (features.feature.str.replace("_fraction", "", regex=False)
              .str.replace("flank_", "flank: ", regex=False)
              .str.replace("shell_", "ND2 shell: ", regex=False)
              .str.replace("_", " ", regex=False))
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    y = np.arange(len(features))
    for i, r in enumerate(features.itertuples()):
        ax.plot([r.ci_low, r.ci_high], [i, i], color=DIM, lw=1.6, zorder=2)
    ax.scatter(features[col], y, s=46,
               color=[OCC if q < 0.10 else DIM for q in features.q], zorder=3)
    ax.axvline(0, color=INK, lw=1.0, zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(pretty, fontsize=9.5)
    ax.set_xlabel("shift, design minus wild type (reference standard deviations)")
    ax.set_title("No single feature carries the drift\n"
                 "(shaded: q < 0.10; none survives correction at 0.05)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_context_features.png", dpi=200, bbox_inches="tight")
    print("wrote", OUT / "fig5_context_features.png")
