"""Feature distributions: natural occupied sites against protected-sequon designs.

One small panel per feature, natural sites filled and designs outlined, on shared
bins. Histograms rather than smoothed densities because several of these features
are heavily zero-inflated -- 91% of natural sites have no cysteine in the ND2
shell -- and a kernel density would spread mass across values that never occur.

Text on the plot is held to feature names, axis labels and a legend. Everything
else belongs in the caption.

Usage:  54_context_distributions.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "glyco_context/src")
from glyco_context.local_chemistry import CLASSES

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10,
                     "axes.titlelocation": "left", "axes.spines.top": False,
                     "axes.spines.right": False})
INK, OCC, NAT = "#22252b", "#1f6f8b", "#c7ccd1"

ANALYSIS = Path("glyco_context/results/analysis")
OUT = Path("glyco_context/results/figures")
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ([f"flank_{c}_fraction" for c in CLASSES]
         + [f"shell_{c}_fraction" for c in CLASSES] + ["shell_net_charge"])
LABEL = {f"flank_{c}_fraction": f"flanking {c}" for c in CLASSES}
LABEL.update({f"shell_{c}_fraction": f"ND2 shell {c}" for c in CLASSES})
LABEL["shell_net_charge"] = "ND2 shell net charge"

natural = pd.read_csv(ANALYSIS / "natural_reference_panels.csv")
natural = natural[natural.variant == "wild_type"]
designs = pd.read_csv(ANALYSIS / "fixed_sequon_panels.csv")
designs = designs[designs.variant == "design"]

rows, cols = 5, 3
fig, axes = plt.subplots(rows, cols, figsize=(11.0, 12.4))
for ax, feature in zip(axes.ravel(), PANEL):
    a = natural[feature].dropna().to_numpy(float)
    b = designs[feature].dropna().to_numpy(float)
    if not len(a) or not len(b):
        ax.axis("off"); continue
    # A modest fixed number of bins. Aligning bins to the data's own levels was
    # tried and is worse: the flanking window is clipped at chain ends, so the
    # denominator varies between sites and the "levels" are not shared.
    lo, hi = float(min(a.min(), b.min())), float(max(a.max(), b.max()))
    bins = np.linspace(lo, hi, 15) if hi > lo else np.linspace(lo - 0.5, hi + 0.5, 3)

    ax.hist(a, bins=bins, density=True, color=NAT, label="natural occupied")
    ax.hist(b, bins=bins, density=True, histtype="step", lw=1.8, color=OCC,
            label="design")
    ax.set_title(LABEL[feature])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)

for ax in axes.ravel()[len(PANEL):]:
    ax.axis("off")
handles, labels = axes.ravel()[0].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=10, ncol=2,
           loc="lower center", bbox_to_anchor=(0.5, 0.008))
fig.text(0.5, 0.052, "value", ha="center", fontsize=10)
fig.text(0.02, 0.5, "density", va="center", rotation="vertical", fontsize=10)
fig.tight_layout(rect=(0.03, 0.085, 1, 1))
fig.savefig(OUT / "fig6_feature_distributions.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig6_feature_distributions.png")
print(f"natural n={len(natural)}   designs n={len(designs)}")
