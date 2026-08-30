"""The model comparison, as separate figures rather than four panels of one.

Split because they answer different questions and are read at different times:

  fig_scoring               what each model scores, and the spread behind it
  fig_masking               what survives hiding the motif
  fig_retention             what each model keeps when it redesigns the chain
  fig_score_distributions   the raw log-odds, occupied against no-annotation

Everything except the last is SD-standardised. Raw magnitudes are not comparable
between models -- eight averaged decoding orders, one causal pass, one masked
position, one shot -- so each model's contrasts are divided by its own reference
SD, which is the quantity the analysis rests on anyway.

Deliberately spare: titles, axes, keys and significance markers only. Exact
numbers go to a companion values file and to docs/figures_and_captions.md.

A model with no data for a figure is omitted and named, and one that *cannot*
have data is distinguished from one that has simply not been run -- saying the
wrong reason is worse than saying none.

Usage:  25_figures_model_comparison.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 12,
                     "axes.titlelocation": "left", "axes.spines.top": False,
                     "axes.spines.right": False})

INK, OCC, CTL, DIM = "#22252b", "#1f6f8b", "#b05c3b", "#9aa0a6"
ANALYSIS = Path("results/analysis")
OUT = Path("results/figures")
OUT.mkdir(parents=True, exist_ok=True)

# Ordered so ESM3's two arms sit adjacent: the difference between them is the
# only within-model measurement of what the structure track contributes.
#
#   (name, motif-visible variant, motif-hidden variant, retention variant, colour)
#
# A None hidden variant means the model HAS no masked condition; a named one
# that has not been run means only that. The two are shown differently.
MODELS = [
    ("ProteinMPNN", "proteinmpnn_index_corrected",
     "proteinmpnn_joint_index_corrected", "proteinmpnn_index_corrected", OCC),
    ("ESM-IF", "esm_if_index_corrected",
     "esm_if_joint_index_corrected", "esm_if", CTL),
    ("CARBonAra", "carbonara", None, "carbonara", "#8a6fa8"),
    ("ESM3 structure", "esm3_struct_single", "esm3_struct_joint", "esm3", "#1f6f8b"),
    ("ESM3 sequence", "esm3_seq_single", "esm3_seq_joint", None, "#7fb3c8"),
    ("ESMC", "esmc_single", "esmc_joint", None, "#5b8c5a"),
    ("ProGen2", "progen2_direction1", "progen2_joint_direction1", None, "#c58a3d"),
]

# Why a model can never have a retention row, as opposed to not having one yet.
# Why a model has no retention row. "ESM3 structure" was listed here as
# conditioning on sequence not backbone, which is ESMC's reason and false for
# it: the structure track IS a backbone conditioning, and it can redesign a
# chain. It is absent from this figure only until its design arm has been run.
NO_DESIGNER = {
    "ESMC": "masked LM, conditions on sequence not backbone",
    "ESM3 sequence": "structure track withheld, so no backbone to redesign",
    "ProGen2": "generation is unconditioned by any backbone",
}
SECRETORY = "eukaryotic secretory"
values = {}


def scoring(variant):
    """The contrast for one variant, or None if it has not been run."""
    if variant is None:
        return None
    path = ANALYSIS / f"analysis_secretory_{variant}.json"
    return json.loads(path.read_text()) if path.exists() else None


def retention(variant):
    if variant is None:
        return None
    path = ANALYSIS / f"retention_paired_{variant}.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["comparisons"].get(SECRETORY)


def masking(visible, hidden):
    if hidden is None:
        return None
    path = ANALYSIS / f"masking_secretory_{visible}_vs_{hidden}.json"
    return json.loads(path.read_text()) if path.exists() else None


def note(ax, start, texts):
    """Italic reason on a row with no data, so a gap is never unexplained."""
    for offset, text in enumerate(texts):
        ax.text(0.02, start + offset, text, transform=ax.get_yaxis_transform(),
                fontsize=9.5, color=INK, alpha=0.65, style="italic", va="center")


# =============================================================================
# fig_scoring
# =============================================================================
scored = [(n, scoring(v), c) for n, v, _, _, c in MODELS if scoring(v)]
unscored = [n for n, v, _, _, _ in MODELS if not scoring(v)]

fig, (ax, dx) = plt.subplots(1, 2, figsize=(12.8, 5.6),
                             gridspec_kw={"width_ratios": [1, 1.05]})
for i, (name, stat, colour) in enumerate(scored):
    low, high = stat["ci95_sd"]
    ax.plot([low, high], [i, i], color=colour, lw=3.4, solid_capstyle="round")
    ax.plot(stat["mean_difference_sd"], i, "o", color=colour, ms=10, zorder=3)
    if low > 0:
        ax.text(high + 0.03, i, "*", va="center", fontsize=15, color=INK)
note(ax, len(scored), ["not yet run"] * len(unscored))
ax.axvline(0, color=INK, lw=1.0, zorder=1)
ax.set_yticks(range(len(scored) + len(unscored)))
ax.set_yticklabels([s[0] for s in scored] + unscored)
ax.set_ylim(-0.7, len(scored) + len(unscored) - 0.3)
ax.invert_yaxis()
ax.tick_params(left=False)
ax.set_xlabel("occupied - matched control  (reference SD)")
ax.set_title("A   Occupancy-associated score shifts across model conditionals",
             pad=10)

for name, variant, _, _, colour in MODELS:
    stat = scoring(variant)
    path = ANALYSIS / f"contrasts_secretory_{variant}.csv"
    if stat is None or not path.exists():
        continue
    standardised = np.sort(pd.read_csv(path).contrast.values / stat["reference_sd"])
    ecdf = np.arange(1, len(standardised) + 1) / len(standardised)
    dx.step(standardised, ecdf, where="post", color=colour, lw=2.0, label=name)
dx.axvline(0, color=INK, lw=1.0, zorder=1)
dx.set_xlim(-4, 4)
dx.set_ylim(0, 1)
dx.set_xlabel("per-site contrast  (reference SD)")
dx.set_ylabel("cumulative proportion of sites")
dx.set_title("B   The averages come from broad, overlapping distributions", pad=10)
dx.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(pad=1.8)
fig.savefig(OUT / "fig_scoring.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig_scoring.png")

values["scoring"] = {
    n: {"mean_sd": s["mean_difference_sd"], "ci95_sd": s["ci95_sd"],
        "n_contrasts": s["n_contrasts"], "verdict": s["verdict"]}
    for n, s, _ in scored}
values["scoring_not_run"] = unscored

# =============================================================================
# fig_masking
# =============================================================================
maskable = [m for m in MODELS if scoring(m[1]) and scoring(m[2])]
absent = {m[0]: ("no masked condition - one-shot, so no motif-visible arm"
                 if m[2] is None else "masked arm not yet run")
          for m in MODELS if scoring(m[1]) and not scoring(m[2])}

height = 0.62 * (len(maskable) + len(absent)) + 1.9
fig, bx = plt.subplots(figsize=(9.4, height))
for i, (name, vis, hid, _, colour) in enumerate(maskable):
    v, h = scoring(vis)["mean_difference_sd"], scoring(hid)["mean_difference_sd"]
    bx.annotate("", xy=(h, i), xytext=(v, i),
                arrowprops=dict(arrowstyle="-|>", color=colour, lw=3.0,
                                shrinkA=0, shrinkB=0, mutation_scale=22))
    bx.plot(v, i, "o", color=colour, ms=10, zorder=3,
            label="motif visible" if i == 0 else None)
    bx.plot(h, i, "o", color="white", ms=10, zorder=3, markeredgecolor=colour,
            markeredgewidth=2.4, label="motif hidden" if i == 0 else None)
    change = masking(vis, hid)
    # Above the arrow's midpoint: the marker refers to the CHANGE, so it has to
    # read as attached to the arrow rather than to either endpoint.
    if change and change["p"] < 0.05:
        bx.text((v + h) / 2, i - 0.28, "*", va="center", ha="center",
                fontsize=15, color=INK)
note(bx, len(maskable), list(absent.values()))
bx.axvline(0, color=INK, lw=1.0, zorder=1)
bx.margins(x=0.10)
rows = len(maskable) + len(absent)
bx.set_yticks(range(rows))
bx.set_yticklabels([m[0] for m in maskable] + list(absent))
bx.set_ylim(-0.8, rows - 0.3)
bx.invert_yaxis()
bx.tick_params(left=False)
bx.set_xlabel("occupied - matched control  (reference SD)")
bx.set_title("Effect of hiding native sequon identity", pad=10)
bx.legend(frameon=False, fontsize=9.5, loc="upper left")
fig.tight_layout(pad=1.6)
fig.savefig(OUT / "fig_masking.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig_masking.png")

values["masking"] = {}
for name, vis, hid, _, _ in maskable:
    change = masking(vis, hid)
    values["masking"][name] = {
        "visible_sd": scoring(vis)["mean_difference_sd"],
        "hidden_sd": scoring(hid)["mean_difference_sd"],
        "change_logodds": change["mean"] if change else None,
        "ci95_logodds": [change["ci_low"], change["ci_high"]] if change else None,
        "p": change["p"] if change else None}
values["masking_absent"] = absent

# =============================================================================
# fig_retention
# =============================================================================
designers = [(n, retention(r), c) for n, _, _, r, c in MODELS if retention(r)]
pending = [n for n, _, _, r, _ in MODELS if r is not None and not retention(r)]

height = 0.62 * (len(designers) + len(pending)) + 1.9
fig, cx = plt.subplots(figsize=(9.4, height))
for i, (name, r, colour) in enumerate(designers):
    cx.plot([r["control_mean"], r["occupied_mean"]], [i, i], color=colour,
            lw=2.0, alpha=0.55, zorder=1)
    cx.plot(r["control_mean"], i, "o", color="white", ms=10, zorder=3,
            markeredgecolor=DIM, markeredgewidth=2.2,
            label="matched control" if i == 0 else None)
    cx.plot(r["occupied_mean"], i, "o", color=colour, ms=10, zorder=3,
            label="occupied" if i == 0 else None)
    if r["excludes_zero"]:
        cx.text(r["occupied_mean"] + 0.006, i, "*", va="center", fontsize=15,
                color=INK)
note(cx, len(designers), ["not yet run"] * len(pending))
cx.set_yticks(range(len(designers) + len(pending)))
cx.set_yticklabels([d[0] for d in designers] + pending)
cx.set_ylim(-0.8, len(designers) + len(pending) - 0.3)
cx.set_xlim(left=0)
cx.invert_yaxis()
cx.tick_params(left=False)
without = [n for n in NO_DESIGNER if n in {m[0] for m in MODELS}]
cx.set_xlabel("sequon retained in redesign  (proportion of sites)"
              + ("\n" + ", ".join(without) + ": no designer" if without else ""),
              fontsize=10)
cx.set_title("Design: the sequon is usually not kept", pad=10)
# Upper right, not lower: with ESM3 added the bottom row's occupied dot lands
# beside the legend markers, and the two read as data on the same row.
cx.legend(frameon=False, fontsize=9.5, loc="upper right")
fig.tight_layout(pad=1.6)
fig.savefig(OUT / "fig_retention.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig_retention.png")

values["retention"] = {
    n: {"occupied": r["occupied_mean"], "control": r["control_mean"],
        "paired_difference": r["paired_difference"], "ci95": r["ci95"],
        "n_pairs": r["n_pairs"], "n_tied": r["n_tied"],
        "wilcoxon_p": r["wilcoxon_p"], "ci_excludes_zero": r["excludes_zero"]}
    for n, r, _ in designers}
values["retention_no_designer"] = NO_DESIGNER
values["retention_not_run"] = pending

# =============================================================================
# fig_score_distributions
# =============================================================================
# NOT standardised, because this is what each model actually emits. The panels
# therefore cannot share an axis and must not be read against one another.
drawn = [m for m in MODELS if scoring(m[1])
         and (ANALYSIS / f"contrasts_secretory_{m[1]}.csv").exists()]
# Laid out as a grid rather than one row: at seven models a single row is
# thirty inches wide and cannot be printed or read. The panels never shared an
# axis anyway -- each model emits its own scale, which is the whole point of
# this figure -- so wrapping them costs no comparability.
COLS = 4
n_panels = len(drawn)
cols = min(COLS, n_panels)
rows = -(-n_panels // cols)
fig, axes = plt.subplots(rows, cols, figsize=(4.1 * cols, 3.6 * rows),
                         squeeze=False)
flat = axes.ravel()
for spare in flat[n_panels:]:
    spare.set_visible(False)
for i, (panel, (name, variant, _, _, colour)) in enumerate(zip(flat, drawn)):
    pairs = pd.read_csv(ANALYSIS / f"contrasts_secretory_{variant}.csv")
    occupied, unknown = pairs.case_score.values, pairs.control_mean_score.values
    edges = np.histogram_bin_edges(np.concatenate([occupied, unknown]), bins=26)
    panel.hist(unknown, bins=edges, color=DIM, alpha=0.55, density=True,
               label="no annotated glycan")
    panel.hist(occupied, bins=edges, color=colour, alpha=0.62, density=True,
               label="experimentally occupied")
    panel.axvline(np.median(unknown), color=DIM, lw=1.8, ls="--")
    panel.axvline(np.median(occupied), color=colour, lw=1.8, ls="--")
    if scoring(variant)["ci95_sd"][0] > 0:
        panel.text(0.97, 0.97, "*", transform=panel.transAxes, ha="right",
                   va="top", fontsize=17, color=INK)
    panel.set_title(name, pad=8)
    # Axis labels only on the edges of the grid: an x label under every panel
    # repeats the same string seven times and crowds out the ticks.
    if i + cols >= n_panels:
        panel.set_xlabel("conditional sequon score  (log odds)")
    if i % cols == 0:
        panel.set_ylabel("density of sites")
    panel.set_yticks([])
    panel.spines["left"].set_visible(False)
    values.setdefault("raw_distributions", {})[name] = {
        "n_pairs": int(len(pairs)),
        "occupied_median": float(np.median(occupied)),
        "unknown_median": float(np.median(unknown)),
        "occupied_below_unknown_median":
            float((occupied < np.median(unknown)).mean())}
flat[0].legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(pad=1.6)
fig.savefig(OUT / "fig_score_distributions.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig_score_distributions.png")

(OUT / "fig_model_comparison_values.json").write_text(json.dumps(values, indent=2))
print("wrote", OUT / "fig_model_comparison_values.json")
if unscored:
    print(f"\nnot yet run, omitted: {unscored}")
