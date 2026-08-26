"""Figures for the fixed-sequon context-retention test.

Two plain figures, each answering one question.

  fig4  how far designs sit from natural occupied context, and by how much that
        differs from the wild type and from arbitrary change of the same size
  fig5  which features move

Deliberately spare. An earlier version drew all fifty pairings as crossing lines
and the differences as jittered clouds; both were honest and neither was
readable. The numbers here are small and few, so the figures should be too.

Usage:  53_context_retention_figures.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12.5,
                     "axes.titlelocation": "left", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.spines.left": False})
INK, OCC, CTL, DIM = "#22252b", "#1f6f8b", "#b05c3b", "#9aa0a6"

ANALYSIS = Path("glyco_context/results/analysis")
OUT = Path("glyco_context/results/figures")
OUT.mkdir(parents=True, exist_ok=True)

scored = pd.read_csv(ANALYSIS / "context_retention_distances.csv")
site = scored.pivot_table(index=["accession", "position"], columns="variant",
                          values="distance", aggfunc="mean").dropna(subset=["wild_type"])

def ci(values, n_boot=2000, seed=0):
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, (n_boot, len(values)), replace=True).mean(axis=1)
    return values.mean(), np.percentile(draws, 2.5), np.percentile(draws, 97.5)

# =============================================================================
# fig4 — the result, as two simple rows of dots
# =============================================================================
fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.6, 6.2),
                             gridspec_kw={"height_ratios": [1, 1.15]})

# where each group sits
groups = [("wild type", site.wild_type, INK),
          ("random control", site.random, DIM),
          ("design", site.design, OCC)]
for i, (name, values, colour) in enumerate(groups):
    mean, low, high = ci(values)
    ax.plot([low, high], [i, i], color=colour, lw=3, solid_capstyle="round")
    ax.plot(mean, i, "o", color=colour, ms=10)
    ax.text(high + 0.012, i, f"{mean:.3f}", va="center", fontsize=10.5, color=INK)
ax.set_yticks(range(3)); ax.set_yticklabels([g[0] for g in groups])
ax.set_ylim(-0.6, 2.6); ax.invert_yaxis()
ax.set_xlabel("distance from natural occupied context")
ax.set_title("Designs sit further from natural context than the wild type", pad=10)
ax.tick_params(left=False)

# the paired change
changes = [("design\nminus wild type", site.design - site.wild_type, OCC),
           ("random control\nminus wild type", site.random - site.wild_type, DIM)]
for i, (name, values, colour) in enumerate(changes):
    mean, low, high = ci(values)
    bx.plot([low, high], [i, i], color=colour, lw=3, solid_capstyle="round")
    bx.plot(mean, i, "o", color=colour, ms=10)
    bx.text(high + 0.004, i, f"+{mean:.3f}", va="center", fontsize=10.5, color=INK)
bx.axvline(0, color=INK, lw=1.0, zorder=1)
bx.set_yticks(range(len(changes))); bx.set_yticklabels([c[0] for c in changes])
bx.set_ylim(-0.6, len(changes) - 0.4); bx.invert_yaxis()
bx.set_xlabel("change in distance  (right = further from natural context)")
bx.set_title("The design moves about twice as far as arbitrary change\n"
             "of the same size — but the gap between them is not significant", pad=10)
bx.tick_params(left=False)
bx.text(0, -0.30, f"{len(site)} sites, 38 proteins. Bars are 95% intervals. Every design "
        "keeps the sequon.\nA wild type sits inside the reference by construction, so the "
        "random control\nshows how much of the drift is simply that.",
        transform=bx.transAxes, fontsize=9.5, color=INK, va="top")
fig.tight_layout()
fig.savefig(OUT / "fig4_context_retention.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig4_context_retention.png")

# =============================================================================
# fig5 — only the features that move at all
# =============================================================================
path = ANALYSIS / "context_retention_features.csv"
if path.exists():
    features = pd.read_csv(path)
    col = "shift" if "shift" in features.columns else "mean"
    notable = features[features.q < 0.25].copy()
    notable = notable.reindex(notable[col].abs().sort_values().index)
    name = (notable.feature.str.replace("_fraction", "", regex=False)
            .str.replace("flank_", "flanking ", regex=False)
            .str.replace("shell_", "near ND2: ", regex=False))
    fig, ax = plt.subplots(figsize=(9.0, 3.4))
    y = np.arange(len(notable))
    for i, r in enumerate(notable.itertuples()):
        colour = OCC if getattr(r, col) > 0 else CTL
        # Weaker evidence drawn fainter, so strength is visible without a legend.
        alpha = 1.0 if r.q < 0.10 else 0.42
        ax.plot([r.ci_low, r.ci_high], [i, i], color=colour, lw=3, alpha=alpha,
                solid_capstyle="round")
        ax.plot(getattr(r, col), i, "o", color=colour, ms=9, alpha=alpha)
        ax.text(0.995, i, f"q={r.q:.2f}", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=9, color=DIM)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(y); ax.set_yticklabels(name)
    ax.set_xlabel("shift after design  (reference standard deviations)")
    ax.set_title("Five of fifteen features move at all\nthe strongest: more glycine, fewer aromatics near the site", pad=10)
    ax.tick_params(left=False)
    ax.text(0, -0.42, "Faint bars are weaker evidence. None survives multiple-testing "
            "correction at\nq < 0.05, so this is a lead rather than a result.",
            transform=ax.transAxes, fontsize=9.5, color=INK, va="top")
    fig.tight_layout()
    fig.savefig(OUT / "fig5_context_features.png", dpi=200, bbox_inches="tight")
    print("wrote", OUT / "fig5_context_features.png")
