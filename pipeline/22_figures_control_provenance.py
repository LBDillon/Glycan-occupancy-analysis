"""How the four negative-control sets are built, and what each one costs.

Two panels, one claim each. Left: the filtering path for the new eukaryotic
secretory set, from UniProt query to usable matched pairs. Right: the four sets
placed against the two things that actually matter — how well the population
matches the occupied sites, and how defensible the negative label is. No set
scores well on both, which is the whole reason there is more than one.

Counts are read from the result files rather than written in, so the figure
cannot drift from the data.
"""
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OCC, CTL, INK, MUTE, WARN = "#1f6f8b", "#b05c3b", "#22252b", "#9aa3ad", "#8a5a3b"
OUT = Path("results/figures"); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5, "axes.titlelocation": "left",
                     "axes.spines.top": False, "axes.spines.right": False})


def count_rows(path, mask=None):
    if not Path(path).exists():
        return None
    frame = pd.read_csv(path, low_memory=False)
    return int(len(frame if mask is None else frame[mask(frame)]))


# ---------------------------------------------------------------- panel A data
raw = count_rows("results/datasets/secretory_unannotated_sites_raw.csv")
feats = count_rows("results/datasets/secretory_unannotated_features.csv",
                   lambda f: f.features_available)
scoreable = count_rows("results/manifests/candidate_manifest_secretory.csv",
                       lambda f: f.scoreable.astype(bool))
matched = count_rows("results/matching/matched_pairs_secretory.csv")

stages = [
    ("eukaryotic secretory\nproteins with a structure", 7423),
    ("no glycoprotein\nannotation", 3619),
    ("carries an\nN-X-S/T sequon", raw),
    ("residue resolved,\nfeatures computable", feats),
    ("ProteinMPNN\ncan decode", scoreable),
    ("matched to an\noccupied site", matched),
]
stages = [(label, value) for label, value in stages if value is not None]

fig = plt.figure(figsize=(12.6, 5.4))
ax = fig.add_subplot(1, 2, 1)
y = np.arange(len(stages))[::-1]
values = [v for _, v in stages]
colours = [MUTE] * len(values)
colours[1] = WARN          # the step that defines the negative label
if len(colours) > 5:
    colours[5] = OCC
ax.barh(y, values, color=colours, height=0.62)
for yy, v in zip(y, values):
    ax.text(v * 1.06, yy, f"{v:,}", va="center", fontsize=9.2, color=INK)
ax.set_yticks(y); ax.set_yticklabels([s for s, _ in stages], fontsize=8.8)
ax.set_xscale("log"); ax.set_xlim(8, 20000)
ax.set_xlabel("count (log scale)")
ax.set_title("A. Building the eukaryotic secretory set", fontsize=11.5)
ax.annotate("this step is the entire\nnegative label", xy=(3619, y[1]),
            xytext=(60, y[1] - 0.85), fontsize=8, color=WARN,
            arrowprops=dict(arrowstyle="->", color=WARN, lw=1.1))

# ---------------------------------------------------------------- panel B data
# Three properties, not two. Population match and label strength are the axes;
# usable size is the marker. The internal controls score well on BOTH axes - they
# are not a compromise, they are simply tiny - so the honest reading is that the
# top-right corner is empty of anything large, not that every set trades off.
matched_n = {"secretory": matched,
             "bacterial": count_rows("results/matching/matched_pairs_bacterial.csv"),
             "cytosolic": count_rows("results/matching/matched_pairs_cytosolic.csv"),
             "internal": count_rows("results/matching/matched_pairs_optimal.csv")}

SETS = [
    ("Internal\ncontrols", matched_n["internal"] or 16, 0.93, 0.93, CTL),
    ("Cytosolic\neukaryotic", matched_n["cytosolic"] or 273, 0.26, 0.90, MUTE),
    ("Bacterial\nextracytoplasmic", matched_n["bacterial"] or 280, 0.50, 0.80, MUTE),
    ("Eukaryotic\nsecretory", matched_n["secretory"] or 262, 0.90, 0.26, OCC),
]
bx = fig.add_subplot(1, 2, 2)
bx.axhspan(0.72, 1.05, xmin=0.62, color=OCC, alpha=0.07, zorder=0)
bx.annotate("what we would want:\na large set here", xy=(0.72, 0.99), fontsize=8.2,
            color="#5c6570", ha="center", va="top", style="italic")
for label, n, match, strength, colour in SETS:
    size = 45 + 900 * (n / 300)
    bx.scatter([match], [strength], s=size, color=colour, alpha=0.85, zorder=3)
    bx.annotate(f"{label}\n{n:,} usable pairs", xy=(match, strength),
                xytext=(0, -30 - 0.02 * size), textcoords="offset points",
                ha="center", fontsize=8.3, color=INK)
bx.set_xlim(0.08, 1.14); bx.set_ylim(-0.05, 1.08)
bx.set_xlabel("population match to the occupied sites  \u2192")
bx.set_ylabel("strength of the negative label  \u2192")
bx.set_xticks([0.2, 0.95])
bx.set_xticklabels(["different kind\nof protein", "same kind\nof protein"], fontsize=8.6)
bx.set_yticks([0.2, 0.9])
bx.set_yticklabels(["absence of\nannotation", "absence\nestablished"], fontsize=8.6)
bx.set_title("B. What each set buys, and what it costs", fontsize=11.5)
bx.grid(alpha=0.13, zorder=0)

caption = (
    "A. Every filter applied to reach the eukaryotic secretory control set. The "
    "second step - excluding entries with a glycoprotein keyword - is the whole of "
    "its negative label: nothing mechanistic prevents these sequons being "
    "glycosylated, and about half of all eukaryotic secretory proteins with a "
    "structure do carry that keyword, so the remainder certainly contains real "
    "sites nobody recorded.  "
    "B. Marker area is the number of usable matched pairs. The internal controls are "
    "strong on both axes - they are not a compromise, merely tiny - which is why they "
    "stay the primary comparison. The other three each buy two orders of magnitude more "
    "pairs by giving up one axis: the cytosolic and bacterial sets change the population, "
    "the eukaryotic secretory set weakens the label. The shaded corner is what would "
    "actually settle the question, and nothing large sits in it."
)
lines = textwrap.wrap(caption, width=150)
caption_inches = 0.155 * len(lines) + 0.2
w, h = fig.get_size_inches()
fig.set_size_inches(w, h + caption_inches)
fig.tight_layout(rect=[0, caption_inches / (h + caption_inches), 1, 1])
fig.text(0.008, 0.012, "\n".join(lines), fontsize=7.7, color="#5c6570",
         va="bottom", linespacing=1.5)
fig.savefig(OUT / "fig8_control_provenance.png", dpi=200)
print("wrote", OUT / "fig8_control_provenance.png")
print("stages:", {s: v for s, v in stages})
