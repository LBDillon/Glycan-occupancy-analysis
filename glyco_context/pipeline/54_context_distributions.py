"""Feature distributions as empirical CDFs: natural, wild type, design.

Overlaid histograms were tried first and are the wrong tool here. Binning hides
small shifts, the choice of bin width changes the picture, and several of these
features are mostly zero -- 91% of natural sites have no cysteine in the ND2
shell -- so the zero bar dominates and everything else is invisible.

An empirical CDF has no binning choice, uses every observation exactly, and
reads at these sample sizes as a smooth curve without any density being
invented. A shift between two distributions appears as horizontal separation; a
mass of identical values appears as a vertical jump, which is what those
zero-inflated features honestly are.

Three curves per panel, because two different questions are being asked at once:
where the tested sites sat to begin with (wild type against natural), and where
design moved them (design against wild type). The q value marks the second.

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
INK, OCC, NAT, DIM = "#22252b", "#1f6f8b", "#8d949b", "#9aa0a6"

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
panels = pd.read_csv(ANALYSIS / "fixed_sequon_panels.csv")
wild = panels[panels.variant == "wild_type"]
design = panels[panels.variant == "design"]

qs = {}
features_path = ANALYSIS / "context_retention_features.csv"
if features_path.exists():
    table = pd.read_csv(features_path)
    qs = dict(zip(table.feature, table.q))


def ecdf(ax, values, colour, lw, label, zorder=2):
    x = np.sort(np.asarray(values, float))
    if not len(x):
        return
    y = np.arange(1, len(x) + 1) / len(x)
    ax.step(np.concatenate([x[:1], x]), np.concatenate([[0.0], y]), where="post",
            color=colour, lw=lw, label=label, zorder=zorder)


fig, axes = plt.subplots(5, 3, figsize=(11.0, 12.6))
for ax, feature in zip(axes.ravel(), PANEL):
    ecdf(ax, natural[feature].dropna(), NAT, 2.6, "natural occupied", 2)
    ecdf(ax, wild[feature].dropna(), INK, 1.3, "wild type of tested sites (subset of natural)", 3)
    ecdf(ax, design[feature].dropna(), OCC, 1.8, "design", 4)
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.5, 1])
    q = qs.get(feature)
    star = "*" if q is not None and q < 0.05 else ""
    title = LABEL[feature] + (f"   q={q:.2f}{star}" if q is not None else "")
    ax.set_title(title)

for ax in axes.ravel()[len(PANEL):]:
    ax.axis("off")
handles, labels = axes.ravel()[0].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=10, ncol=3,
           loc="lower center", bbox_to_anchor=(0.5, 0.008))
fig.text(0.5, 0.052, "value", ha="center", fontsize=10)
fig.text(0.02, 0.5, "cumulative fraction of sites", va="center",
         rotation="vertical", fontsize=10)
fig.tight_layout(rect=(0.03, 0.085, 1, 1))
fig.savefig(OUT / "fig6_feature_distributions.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig6_feature_distributions.png")
print(f"natural n={len(natural)}  wild type n={len(wild)}  design n={len(design)}")
