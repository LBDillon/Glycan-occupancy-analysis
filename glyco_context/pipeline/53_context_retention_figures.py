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
    """Mean, percentile interval, and the two-sided bootstrap p against zero."""
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, (n_boot, len(values)), replace=True).mean(axis=1)
    p = min(1.0, max(2 * min((draws <= 0).mean(), (draws >= 0).mean()), 1 / n_boot))
    return values.mean(), np.percentile(draws, 2.5), np.percentile(draws, 97.5), p

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
    mean, low, high, _ = ci(values)
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
    mean, low, high, p = ci(values)
    bx.plot([low, high], [i, i], color=colour, lw=3, solid_capstyle="round")
    bx.plot(mean, i, "o", color=colour, ms=10)
    star = "*" if p < 0.05 else ""
    bx.text(high + 0.004, i, f"+{mean:.3f}   p={p:.3f}{star}", va="center",
            fontsize=10.5, color=INK)
bx.axvline(0, color=INK, lw=1.0, zorder=1)
bx.set_yticks(range(len(changes))); bx.set_yticklabels([c[0] for c in changes])
bx.set_ylim(-0.6, len(changes) - 0.4); bx.invert_yaxis()
bx.set_xlabel("change in distance  (right = further from natural context)")
bx.set_title("The design moves about twice as far as arbitrary change\n"
             "of the same size — but the gap between them is not significant", pad=10)
bx.tick_params(left=False)
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
    fig, ax = plt.subplots(figsize=(9.4, 3.8))
    y = np.arange(len(notable))
    for i, r in enumerate(notable.itertuples()):
        colour = OCC if getattr(r, col) > 0 else CTL
        # Weaker evidence drawn fainter, so strength is visible without a legend.
        alpha = 1.0 if r.q < 0.10 else 0.42
        ax.plot([r.ci_low, r.ci_high], [i, i], color=colour, lw=3, alpha=alpha,
                solid_capstyle="round")
        ax.plot(getattr(r, col), i, "o", color=colour, ms=9, alpha=alpha)
        label = "q<0.01" if r.q < 0.01 else f"q={r.q:.2f}"
        ax.text(1.02, i, label, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=9, color=DIM)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(y); ax.set_yticklabels(name)
    ax.set_xlabel("shift after design  (reference standard deviations)")
    top = notable.reindex(notable[col].abs().sort_values(ascending=False).index).head(3)
    rising = sorted({n.split()[-1] for n, v in zip(
        top.feature.str.replace("_fraction", "", regex=False).str.replace("_", " "),
        top[col]) if v > 0})
    falling = sorted({n.split()[-1] for n, v in zip(
        top.feature.str.replace("_fraction", "", regex=False).str.replace("_", " "),
        top[col]) if v < 0})
    movement = " and ".join(filter(None, [
        "more " + ", ".join(rising) if rising else "",
        "fewer " + ", ".join(falling) if falling else ""]))
    ax.set_title(f"{len(notable)} of 15 features move at all\n"
                 f"the strongest: {movement} near the site", pad=10)
    ax.tick_params(left=False)
    fig.savefig(OUT / "fig5_context_features.png", dpi=200, bbox_inches="tight")
    print("wrote", OUT / "fig5_context_features.png")
