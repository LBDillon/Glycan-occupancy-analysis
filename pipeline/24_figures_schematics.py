"""Draft schematics: the mechanism, the two measurements, and the spectrum.

Deliberately plain. These explain what the pipeline does and what is at stake,
rather than reporting a result — the only one carrying data is the dataset
summary. Drafts: the shapes and wording are meant to be argued with.
"""
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

OCC, CTL, INK, MUTE, NEW = "#1f6f8b", "#b05c3b", "#22252b", "#9aa3ad", "#3f7d5a"
GLY = "#c98b2e"
OUT = Path("results/figures"); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5,
                     "axes.titlelocation": "left"})


def caption(fig, text, width=140):
    lines = textwrap.wrap(text, width=width)
    inches = 0.155 * len(lines) + 0.2
    w, h = fig.get_size_inches()
    fig.set_size_inches(w, h + inches)
    fig.subplots_adjust(bottom=inches / (h + inches) + 0.02)
    fig.text(0.01, 0.012, "\n".join(lines), fontsize=7.7, color="#5c6570",
             va="bottom", linespacing=1.5)


def blob(ax, cx, cy, w, h, colour, alpha=0.18):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor=colour, edgecolor=colour,
                                alpha=alpha, linewidth=1.4))


def arrow(fig, ax, x0, y0, x1, y1, colour=INK, style="-|>", lw=1.5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=14, color=colour, lw=lw,
                                 shrinkA=2, shrinkB=2))


# ===========================================================================
# fig10 — why this matters: the model cannot see the glycan
# ===========================================================================
fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
for ax in axes:
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

def surface(ax, y=0.42):
    xs = np.linspace(0.08, 0.92, 200)
    ax.plot(xs, y + 0.035 * np.sin(xs * 11), color=INK, lw=2.2, zorder=3)

# --- A: what is actually there ---
ax = axes[0]
ax.set_title("A. What is really there", fontsize=11.5)
blob(ax, 0.5, 0.30, 0.84, 0.34, OCC)
surface(ax)
ax.add_patch(Circle((0.5, 0.455), 0.045, facecolor=OCC, edgecolor="white",
                    lw=1.4, zorder=4))
ax.text(0.5, 0.455, "N", ha="center", va="center", fontsize=10,
        color="white", zorder=5, fontweight="bold")
for dx, dy in [(0, 0.13), (-0.06, 0.215), (0.06, 0.215), (-0.115, 0.30)]:
    ax.add_patch(Circle((0.5 + dx, 0.455 + dy), 0.034, facecolor=GLY,
                        edgecolor="white", lw=1.1, zorder=4))
ax.plot([0.5, 0.5], [0.500, 0.585], color=GLY, lw=2.2, zorder=3)
ax.text(0.5, 0.86, "a glycan, attached at\nthe exposed asparagine",
        ha="center", fontsize=9, color=GLY)
ax.text(0.5, 0.08, "solvent-exposed patch\nis exposed FOR A REASON",
        ha="center", fontsize=8.6, color=INK)

# --- B: what the model sees ---
ax = axes[1]
ax.set_title("B. What ProteinMPNN sees", fontsize=11.5)
blob(ax, 0.5, 0.30, 0.84, 0.34, MUTE)
surface(ax)
ax.add_patch(Circle((0.5, 0.455), 0.045, facecolor=MUTE, edgecolor="white",
                    lw=1.4, zorder=4))
ax.text(0.5, 0.455, "N", ha="center", va="center", fontsize=10,
        color="white", zorder=5, fontweight="bold")
ax.text(0.5, 0.86, "glycans are HETATM records.\nthe parser never reads them.",
        ha="center", fontsize=9, color=INK)
ax.text(0.5, 0.08, "an ordinary exposed patch,\nindistinguishable from any other",
        ha="center", fontsize=8.6, color=INK)

# --- C: what it designs ---
ax = axes[2]
ax.set_title("C. What it designs", fontsize=11.5)
blob(ax, 0.5, 0.30, 0.84, 0.34, CTL)
surface(ax)
ax.add_patch(Circle((0.5, 0.455), 0.045, facecolor=CTL, edgecolor="white",
                    lw=1.4, zorder=4))
ax.text(0.5, 0.455, "L", ha="center", va="center", fontsize=10,
        color="white", zorder=5, fontweight="bold")
for dx in (-0.115, 0.115):
    ax.add_patch(Circle((0.5 + dx, 0.448), 0.030, facecolor=CTL, alpha=0.7,
                        edgecolor="white", lw=1.0, zorder=4))
ax.text(0.5, 0.86, "a hydrophobic residue packs\nbetter. the sequon is gone.",
        ha="center", fontsize=9, color=CTL)
ax.text(0.5, 0.08, "usually an improvement.\nhere it aglycosylates the protein.",
        ha="center", fontsize=8.6, color=CTL)

for x in (0.335, 0.665):
    fig.add_artist(FancyArrowPatch((x - 0.018, 0.50), (x + 0.018, 0.50),
                                   transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=22, color=INK, lw=2.0))
fig.suptitle("The problem: an exposed patch can be exposed on purpose",
             fontsize=13, x=0.012, ha="left", y=0.98)
caption(fig,
        "ProteinMPNN parses only ATOM records, so a glycosylated asparagine and a bare one "
        "present it with identical information. Removing a solvent-exposed polar residue in "
        "favour of tighter hydrophobic packing is usually a sound design move, and for most "
        "proteins it is. For a glycoprotein it can delete a required modification site "
        "without anything in the model's input indicating a cost. This is what 'lack of "
        "glycosylation relevance awareness' means concretely, and what the rest of the "
        "project tries to measure. SCHEMATIC - not a real structure.")
fig.savefig(OUT / "fig10_mechanism.png", dpi=200)
plt.close(fig); print("wrote", OUT / "fig10_mechanism.png")


# ===========================================================================
# fig11 — what is actually in each dataset (the only data-carrying schematic)
# ===========================================================================
feat = pd.read_csv("results/datasets/site_structural_features.csv", low_memory=False)
sec = pd.read_csv("results/datasets/secretory_unannotated_sites_raw.csv", low_memory=False)
neg = pd.read_csv("results/datasets/negative_control_sites.csv", low_memory=False)
occ_info = pd.read_csv("results/datasets/occupied_protein_info.csv", low_memory=False)
secf = pd.read_csv("results/datasets/secretory_unannotated_features.csv", low_memory=False)
negf = pd.read_csv("results/datasets/negative_control_features.csv", low_memory=False)

SETS = [("occupied", OCC), ("internal\ncontrol", CTL),
        ("eukaryotic\nsecretory", NEW), ("bacterial", MUTE), ("cytosolic", MUTE)]

fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.2))

# --- A: sites at each stage --------------------------------------------------
ax = axes[0]
sites = [922, 32, len(sec),
         int((neg.control_set == "bacterial_extracytoplasmic").sum()),
         int((neg.control_set == "cytosolic_eukaryotic").sum())]
usable = [16, 16, 262, 280, 273]
x = np.arange(len(SETS))
ax.bar(x - 0.19, sites, width=0.36, color=[c for _, c in SETS], alpha=0.45,
       label="sequons in the set")
ax.bar(x + 0.19, usable, width=0.36, color=[c for _, c in SETS],
       label="usable matched pairs")
for i, (a, b) in enumerate(zip(sites, usable)):
    ax.text(i - 0.19, a * 1.12, f"{a:,}", ha="center", fontsize=7.6)
    ax.text(i + 0.19, b * 1.12, f"{b:,}", ha="center", fontsize=7.6)
ax.set_yscale("log"); ax.set_ylim(6, 60000)
ax.set_xticks(x); ax.set_xticklabels([s for s, _ in SETS], fontsize=8)
ax.set_ylabel("sites (log scale)")
ax.set_title("A. Size, before and after matching", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper left")
for side in ("top", "right"): ax.spines[side].set_visible(False)

# --- B: taxonomic spread -----------------------------------------------------
ax = axes[1]
def top_taxa(series, k=4):
    clean = series.dropna().astype(str).str.split(" (", regex=False).str[0]
    return clean.value_counts().head(k)

panels = [("occupied", occ_info["Organism"] if "Organism" in occ_info else pd.Series(dtype=str)),
          ("eukaryotic secretory", sec.organism),
          ("cytosolic", neg[neg.control_set == "cytosolic_eukaryotic"].organism)]
offset = 0
yticks, ylabels = [], []
for name, series in panels:
    counts = top_taxa(series)
    total = max(len(series.dropna()), 1)
    for organism, n in counts.items():
        yticks.append(offset); ylabels.append(organism[:26])
        ax.barh(offset, 100 * n / total, color=MUTE, height=0.62)
        ax.text(100 * n / total + 1, offset, f"{100*n/total:.0f}%",
                va="center", fontsize=7.4, color="#5c6570")
        offset -= 1
    ax.text(0, offset + len(counts) + 0.55, name, fontsize=8.8, color=INK,
            fontweight="bold", va="bottom")
    offset -= 1.6
ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=7.6)
ax.set_xlabel("% of the set"); ax.set_xlim(0, 78)
ax.set_ylim(offset + 0.6, 1.4)
ax.set_title("B. Taxonomic spread (top 4 each)", fontsize=11)
for side in ("top", "right", "left"): ax.spines[side].set_visible(False)

# --- C: burial ---------------------------------------------------------------
ax = axes[2]
dists = [
    ("occupied", feat[feat.features_available & (feat.occupancy_status == "occupied_supported")].rsa, OCC),
    ("internal ctrl", feat[feat.features_available & (feat.occupancy_status == "observed_unmodified")].rsa, CTL),
    ("euk. secretory", secf[secf.features_available].rsa, NEW),
    ("bacterial", negf[negf.features_available & (negf.control_set == "bacterial_extracytoplasmic")].rsa, MUTE),
    ("cytosolic", negf[negf.features_available & (negf.control_set == "cytosolic_eukaryotic")].rsa, MUTE),
]
parts = ax.boxplot([d.dropna().to_numpy() for _, d, _ in dists], vert=True,
                   widths=0.55, patch_artist=True, showfliers=False)
for patch, (_, _, colour) in zip(parts["boxes"], dists):
    patch.set_facecolor(colour); patch.set_alpha(0.65); patch.set_edgecolor(INK)
for element in ("medians", "whiskers", "caps"):
    for item in parts[element]: item.set_color(INK)
ax.set_xticklabels([n for n, _, _ in dists], fontsize=8, rotation=18, ha="right")
ax.set_ylabel("relative solvent accessibility")
ax.set_title("C. How exposed the sequons are", fontsize=11)
ax.axhline(feat[feat.features_available &
                (feat.occupancy_status == "occupied_supported")].rsa.median(),
           color=OCC, ls=":", lw=1.3)
for side in ("top", "right"): ax.spines[side].set_visible(False)

fig.suptitle("What is in each dataset", fontsize=13, x=0.012, ha="left", y=0.98)
caption(fig,
        "A. Log scale. Set size and usable pairs are different quantities: the cytosolic set "
        "is by far the largest and still yields fewer pairs than the eukaryotic secretory "
        "one, because pairs are limited by the caliper rather than by supply.  "
        "B. A caveat rather than a reassurance. The eukaryotic secretory set matches the "
        "occupied sites on kingdom and compartment but NOT on species: it is 28% human "
        "against the occupied set's 43%, with far more yeast and plant. The cytosolic set is "
        "more human-dominated (54%) than either. Species composition is not matched anywhere "
        "in this study, and this panel is how you would see that mattering.  "
        "C. Dotted line is the occupied median. Occupied sequons are the most exposed - an "
        "oligosaccharyltransferase has to reach them - which is why matching on accessibility "
        "precedes any score comparison. DRAFT.", width=150)
fig.savefig(OUT / "fig11_dataset_summary.png", dpi=200)
plt.close(fig); print("wrote", OUT / "fig11_dataset_summary.png")


# ===========================================================================
# fig12 / fig13 — what each measurement actually does
# ===========================================================================
def flow(ax, steps, colour, y=0.56, box_h=0.42, gap=0.035):
    """Left-to-right boxes with arrows between them.

    Box width is derived from the gap rather than fixed, so the arrows always
    have room; hard-coding a width leaves a sliver the arrowhead cannot fit in
    and the marker renders backwards.
    """
    n = len(steps)
    box_w = (1.0 - 0.04 - gap * (n - 1)) / n
    centres = [0.02 + box_w * (i + 0.5) + gap * i for i in range(n)]
    for cx, (title, body) in zip(centres, steps):
        ax.add_patch(FancyBboxPatch((cx - box_w / 2, y - box_h / 2), box_w, box_h,
                                    boxstyle="round,pad=0.012,rounding_size=0.03",
                                    facecolor=colour, alpha=0.16,
                                    edgecolor=colour, lw=1.5))
        ax.text(cx, y + box_h / 2 - 0.06, title, ha="center", va="top",
                fontsize=9.4, fontweight="bold", color=INK)
        ax.text(cx, y + box_h / 2 - 0.145, body, ha="center", va="top",
                fontsize=7.9, color="#3d434b", linespacing=1.55)
    for a, b in zip(centres, centres[1:]):
        ax.add_patch(FancyArrowPatch((a + box_w / 2 + 0.004, y),
                                     (b - box_w / 2 - 0.004, y),
                                     arrowstyle="-|>", mutation_scale=13,
                                     color=INK, lw=1.4, shrinkA=0, shrinkB=0))
    return centres


# --- fig12: scoring ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.6, 3.4))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
centres = flow(ax, [
    ("structure", "one deposited chain\nbackbone atoms only\nN, CA, C, O"),
    ("check first", "can the model decode\nall three residues?\nif not, EXCLUDE"),
    ("model reads", "backbone + the rest of\nthe native sequence\nnothing is altered"),
    ("probabilities", "20 amino acids at\neach of the three\nsequon positions"),
    ("log-odds", "log(p / 1-p) at N,\nand at S+T\naveraged over 8 orders"),
    ("site score", "one number\nper site\n-5 to +1"),
], OCC)
ax.text(centres[0], 0.20, "glycans are HETATM —\nnever read", ha="center",
        fontsize=7.6, color=GLY, style="italic")
ax.text(centres[3], 0.20, "middle residue excluded:\nanything but proline works",
        ha="center", fontsize=7.6, color="#5c6570", style="italic")
fig.suptitle("Measurement 1 — the conditional sequon score: what the model believes",
             fontsize=12.5, x=0.012, ha="left", y=0.97)
caption(fig,
        "Nothing is generated and the sequon is never altered; the model is asked what it "
        "would put at those three positions given everything else. Scoreability is checked "
        "BEFORE matching, from coordinates alone — a residue with an incomplete backbone is "
        "returned as a row of zeros that exponentiates to twenty-one ones and scores +13.8, "
        "which is the defect that inverted the first version of this study. Run by "
        "pipeline/07_score.py. DRAFT.", width=150)
fig.savefig(OUT / "fig12_scoring_process.png", dpi=200)
plt.close(fig); print("wrote", OUT / "fig12_scoring_process.png")

# --- fig13: design ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.6, 3.4))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
centres = flow(ax, [
    ("structure", "same chain,\nsame backbone\nas measurement 1"),
    ("generate", "32 complete sequences\nunconstrained, T=0.1\nnothing held fixed"),
    ("read back", "look at the three\noriginal sequon\npositions in each"),
    ("classify", "N present?\nS or T present?\nproline at the middle?"),
    ("retention", "fraction of 32 designs\nkeeping the whole motif\n0.00 to 1.00"),
], NEW)
ax.text(centres[1], 0.20, "the sequon is free to vanish —\nthat is the point",
        ha="center", fontsize=7.6, color=NEW, style="italic")
ax.text(centres[4], 0.20, "median across all sites: 0.00",
        ha="center", fontsize=7.6, color="#5c6570", style="italic")
fig.suptitle("Measurement 2 — design retention: what the model actually writes",
             fontsize=12.5, x=0.012, ha="left", y=0.97)
caption(fig,
        "No positions are fixed and no residues biased, so a sequon survives only if the "
        "model independently chooses to keep it. 32 designs rather than 8 because 8 gives a "
        "per-site standard error near 0.18 — too coarse for a site-level analysis. Designs "
        "are generated per chain, so every sequon on a chain is read from the same 32. "
        "Run by pipeline/08_design.py. The two measurements track each other at Spearman "
        "+0.547, which is what licenses treating them as two views of one quantity. DRAFT.",
        width=150)
fig.savefig(OUT / "fig13_design_process.png", dpi=200)
plt.close(fig); print("wrote", OUT / "fig13_design_process.png")


# ===========================================================================
# fig14 — the glycosylation-awareness spectrum (DRAFT: the axis is a proposal)
# ===========================================================================
fig, ax = plt.subplots(figsize=(13.2, 4.4))
ax.set_xlim(-0.04, 1.06); ax.set_ylim(0, 1); ax.axis("off")

ax.annotate("", xy=(1.02, 0.62), xytext=(0.0, 0.62),
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=2.2))
ax.text(0.0, 0.665, "glycan-blind", fontsize=9.5, color=MUTE, fontweight="bold")
ax.text(1.02, 0.665, "glycan-aware", fontsize=9.5, color=MUTE,
        fontweight="bold", ha="right")

STATIONS = [
    (0.03, "0", "cannot see glycans",
     "glycans are not in the input\nat all — HETATM, stripped"),
    (0.30, "1", "sees structure, not function",
     "treats an occupied sequon\nexactly like any other motif"),
    (0.62, "2", "implicit awareness",
     "treats occupied sequons\ndifferently without being told,\nfrom structural context alone"),
    (0.92, "3", "explicit representation",
     "the glycan is modelled;\ndesign preserves the site\nbecause it knows it is there"),
]
for x, num, title, body in STATIONS:
    ax.plot([x], [0.62], "o", markersize=13, color=INK, zorder=4)
    ax.text(x, 0.62, num, ha="center", va="center", fontsize=8,
            color="white", zorder=5, fontweight="bold")
    ax.text(x, 0.50, title, ha="center", va="top", fontsize=9.2,
            fontweight="bold", color=INK)
    ax.text(x, 0.435, body, ha="center", va="top", fontsize=7.8,
            color="#3d434b", linespacing=1.5)

# where the evidence puts ProteinMPNN
ax.add_patch(FancyBboxPatch((0.005, 0.775), 0.47, 0.165,
                            boxstyle="round,pad=0.014,rounding_size=0.03",
                            facecolor=OCC, alpha=0.15, edgecolor=OCC, lw=1.6))
ax.text(0.24, 0.915, "ProteinMPNN sits between 0 and 1", ha="center",
        fontsize=10.2, fontweight="bold", color=OCC)
ax.text(0.24, 0.868, "structurally at 0: it never receives a glycan.\n"
                     "empirically no better than 1: no comparison, on either\n"
                     "measurement, distinguishes occupied sequons.",
        ha="center", va="top", fontsize=8.1, color="#3d434b", linespacing=1.5)
ax.plot([0.03, 0.03], [0.775, 0.645], color=OCC, lw=1.5, ls="--")
ax.plot([0.30, 0.30], [0.775, 0.645], color=OCC, lw=1.5, ls="--")
ax.annotate("", xy=(0.30, 0.72), xytext=(0.03, 0.72),
            arrowprops=dict(arrowstyle="<->", color=OCC, lw=1.5))

ax.text(0.62, 0.24, "what would move a model to 2", ha="center",
        fontsize=9, fontweight="bold", color=NEW)
ax.text(0.62, 0.19, "a positive, pre-registered difference on a matched\n"
                    "comparison large enough to exclude the margin",
        ha="center", va="top", fontsize=7.9, color="#3d434b", linespacing=1.5)
ax.text(0.92, 0.24, "what would move it to 3", ha="center",
        fontsize=9, fontweight="bold", color=GLY)
ax.text(0.92, 0.19, "a different model — glycans in the input,\n"
                    "not a better readout of this one",
        ha="center", va="top", fontsize=7.9, color="#3d434b", linespacing=1.5)

fig.suptitle("Where a model can sit on glycosylation awareness  —  DRAFT, the axis is a proposal",
             fontsize=12.5, x=0.012, ha="left", y=0.97)
caption(fig,
        "The stations are a proposed scale, not an established one, and the gap between 1 and "
        "2 is the only part this project can currently measure. Note that 0 and 1 are "
        "genuinely different claims: 0 is about the input representation and is simply true "
        "of ProteinMPNN, while 1 is about behaviour and is what the experiments test. A model "
        "at 0 could in principle still reach 2, if occupied sequons sat in backbone "
        "environments distinctive enough to be learned indirectly — that is exactly the "
        "hypothesis under test, and the answer so far is no detectable difference. Worth "
        "arguing about before this is used in a write-up. DRAFT.", width=150)
fig.savefig(OUT / "fig14_awareness_spectrum.png", dpi=200)
plt.close(fig); print("wrote", OUT / "fig14_awareness_spectrum.png")
