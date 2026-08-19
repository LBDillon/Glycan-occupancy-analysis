"""The two measurements, all five classes, side by side.

Left: what the model *writes* — retention for every class on one axis, so the
occupied sites can be read against each control set directly.

Right: what the model *believes* — the conditional-score difference against each
control set, with the pre-specified equivalence margin.

They are separate panels rather than one because they are different quantities
measured on different scales. The point of putting them together is that they
agree: neither shows the model distinguishing occupied sequons from unoccupied
ones by any margin the design can resolve.
"""
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OCC, CTL, INK, MUTE, NEW = "#1f6f8b", "#b05c3b", "#22252b", "#9aa3ad", "#3f7d5a"
OUT = Path("results/figures"); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5, "axes.titlelocation": "left",
                     "axes.spines.top": False, "axes.spines.right": False})

byclass = json.loads(Path("results/analysis/retention_by_class.json").read_text())["classes"]
ORDER = ["occupied_supported", "observed_unmodified",
         "control_secretory_eukaryotic_unannotated",
         "control_bacterial_extracytoplasmic", "control_cytosolic_eukaryotic"]
COLOUR = {"occupied_supported": OCC, "observed_unmodified": CTL,
          "control_secretory_eukaryotic_unannotated": NEW,
          "control_bacterial_extracytoplasmic": MUTE,
          "control_cytosolic_eukaryotic": MUTE}
SHORT = {"occupied_supported": "occupied\n(experimental\nglycan)",
         "observed_unmodified": "internal\ncontrol",
         "control_secretory_eukaryotic_unannotated": "eukaryotic\nsecretory",
         "control_bacterial_extracytoplasmic": "bacterial",
         "control_cytosolic_eukaryotic": "cytosolic"}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.4, 4.4),
                             gridspec_kw={"width_ratios": [1.15, 1]})

present = [k for k in ORDER if k in byclass]
labels = []
for i, key in enumerate(present):
    entry = byclass[key]
    mean = entry["mean_retention"]; low, high = entry["ci95"]
    ax.bar(i, mean, color=COLOUR[key], width=0.6)
    if np.isfinite(low):
        ax.errorbar(i, mean, yerr=[[mean - low], [high - mean]], fmt="none",
                    ecolor=INK, capsize=5, lw=1.4)
    ax.text(i, (high if np.isfinite(high) else mean) + 0.004, f"{mean:.3f}",
            ha="center", fontsize=9.2)
    labels.append(f"{SHORT[key]}\nn={entry['n_sites']:,}")
ax.set_xticks(range(len(present))); ax.set_xticklabels(labels, fontsize=8.4)
ax.set_ylabel("fraction of 32 designs keeping the sequon")
ax.set_title("A. What the model writes (design retention)")
ax.set_ylim(0, max(byclass[k]["ci95"][1] for k in present) * 1.28)

COMPARISONS = [("optimal", "internal control", CTL),
               ("secretory", "eukaryotic secretory", NEW),
               ("bacterial", "bacterial", MUTE),
               ("cytosolic", "cytosolic", MUTE)]
margin = None
for i, (key, label, colour) in enumerate(COMPARISONS):
    path = Path(f"results/analysis/analysis_{key}.json")
    if not path.exists():
        continue
    d = json.loads(path.read_text())
    margin = d["margin_standardised"]
    mean, (low, high) = d["mean_difference_sd"], d["ci95_sd"]
    bx.errorbar([mean], [-i], xerr=[[mean - low], [high - mean]], fmt="o",
                color=colour, capsize=4, lw=1.7, markersize=7, zorder=3)
    bx.text(high + 0.05, -i, f"{mean:+.3f}  (n={d['n_contrasts']})",
            va="center", fontsize=8.6, color=colour)
if margin:
    bx.axvspan(-margin, margin, color=MUTE, alpha=0.18, zorder=0)
bx.axvline(0, color=INK, lw=0.9, zorder=1)
bx.set_yticks([-i for i in range(len(COMPARISONS))])
bx.set_yticklabels([c[1] for c in COMPARISONS], fontsize=9)
bx.set_xlabel("occupied − control (SD units)")
bx.set_xlim(-0.75, 1.7); bx.set_ylim(-len(COMPARISONS) + 0.45, 0.55)
bx.set_title("B. What the model believes (conditional score)")

caption = (
    "A. Every class on one axis, intervals bootstrapped over proteins rather than sites "
    "because sequons on a chain share one set of designs. Occupied sequons are destroyed "
    "as readily as controls that cannot be glycosylated at all.  "
    "B. Occupied minus control in units of the reference SD; shaded band is the "
    "pre-specified equivalence margin. The internal-control comparison has the strongest "
    "negative label and only 16 pairs; the eukaryotic secretory comparison matches on both "
    "taxonomy and compartment and supplies 262.  "
    "The two panels are different quantities and agree: nothing here shows the model "
    "distinguishing occupied sequons by a margin this design can resolve."
)
lines = textwrap.wrap(caption, width=148)
caption_inches = 0.155 * len(lines) + 0.2
w, h = fig.get_size_inches()
fig.set_size_inches(w, h + caption_inches)
fig.tight_layout(rect=[0, caption_inches / (h + caption_inches), 1, 1])
fig.text(0.008, 0.012, "\n".join(lines), fontsize=7.7, color="#5c6570",
         va="bottom", linespacing=1.5)
fig.savefig(OUT / "fig9_all_classes.png", dpi=200)
print("wrote", OUT / "fig9_all_classes.png")
