"""Three models side by side: what each one scores, and what each one designs.

Four panels, one claim each.

  A  the scoring effect, occupied against matched control
  B  what survives hiding the motif -- the sharpest result in the benchmark
  C  the design result, sequon retention, for the two models that generate
  D  the per-site contrasts panel A averages, so the spread is visible

Everything is SD-standardised. Raw log-odds are NOT comparable between models --
ProteinMPNN averages eight decoding orders, ESM-IF runs one causal pass, ESMC
never sees the backbone -- so each model's contrasts are divided by its own
reference SD, which is the quantity the analysis rests on anyway.

Deliberately spare: titles, axes, keys and significance markers only. The exact
numbers go to a companion values file rather than onto the plot.

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

# Current variants only. The superseded ones (alphabet_corrected, the
# pre-index-correction esm_if) are deliberately absent: a figure that quietly
# mixed corrected and uncorrected numbers is how the last two corrections
# stayed invisible for a month.
MODELS = [
    ("ProteinMPNN", "proteinmpnn_index_corrected",
     "proteinmpnn_joint_index_corrected", "proteinmpnn_index_corrected", OCC),
    ("ESM-IF", "esm_if_index_corrected",
     "esm_if_joint_index_corrected", "esm_if", CTL),
    ("ESMC", "esmc_single", "esmc_joint", None, "#5b8c5a"),
    # Scorer only, and one-shot: no motif-visible condition, so no within-model
    # masking contrast and no retention. Panel B necessarily omits it.
    ("CARBonAra", "carbonara", "carbonara", None, "#8a6fa8"),
]
SECRETORY = "eukaryotic secretory"


def scoring(variant):
    return json.loads((ANALYSIS / f"analysis_secretory_{variant}.json").read_text())


def retention(variant):
    if variant is None:
        return None
    path = ANALYSIS / f"retention_paired_{variant}.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["comparisons"][SECRETORY]


def masking(visible, hidden):
    path = ANALYSIS / f"masking_secretory_{visible}_vs_{hidden}.json"
    return json.loads(path.read_text()) if path.exists() else None


fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4))
(ax, bx), (cx, dx) = axes
values = {}

# =============================================================================
# A — the scoring effect, motif visible
# =============================================================================
rows = [(name, scoring(vis), colour) for name, vis, _, _, colour in MODELS]
for i, (name, stat, colour) in enumerate(rows):
    low, high = stat["ci95_sd"]
    mean = stat["mean_difference_sd"]
    ax.plot([low, high], [i, i], color=colour, lw=3.4, solid_capstyle="round")
    ax.plot(mean, i, "o", color=colour, ms=10, zorder=3)
    if low > 0:
        ax.text(high + 0.03, i, "*", va="center", fontsize=15, color=INK)
ax.axvline(0, color=INK, lw=1.0, zorder=1)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows])
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.invert_yaxis()
ax.tick_params(left=False)
ax.set_xlabel("occupied − matched control  (reference SD)")
ax.set_title("A   Scoring: all three models separate occupied sequons", pad=10)

values["A_scoring_visible"] = {
    name: {"mean_sd": s["mean_difference_sd"], "ci95_sd": s["ci95_sd"],
           "n_contrasts": s["n_contrasts"], "wilcoxon_p": s["robustness"]["wilcoxon_p"]}
    for name, s, _ in rows}

# =============================================================================
# B — what survives hiding the motif
# =============================================================================
for i, (name, vis, hid, _, colour) in enumerate(MODELS):
    v, h = scoring(vis)["mean_difference_sd"], scoring(hid)["mean_difference_sd"]
    bx.annotate("", xy=(h, i), xytext=(v, i),
                arrowprops=dict(arrowstyle="-|>", color=colour, lw=3.0,
                                shrinkA=0, shrinkB=0, mutation_scale=22))
    bx.plot(v, i, "o", color=colour, ms=10, zorder=3,
            label="motif visible" if i == 0 else None)
    bx.plot(h, i, "o", color="white", ms=10, zorder=3,
            markeredgecolor=colour, markeredgewidth=2.4,
            label="motif hidden" if i == 0 else None)
    change = masking(vis, hid)
    # Above the arrow's midpoint: the marker refers to the CHANGE, so it has to
    # read as attached to the arrow rather than to either endpoint.
    if change and change["p"] < 0.05:
        bx.text((v + h) / 2, i - 0.26, "*", va="center", ha="center",
                fontsize=15, color=INK)
bx.axvline(0, color=INK, lw=1.0, zorder=1)
bx.margins(x=0.09)
bx.set_yticks(range(len(MODELS)))
bx.set_yticklabels([m[0] for m in MODELS])
bx.set_ylim(-0.85, len(MODELS) - 0.4)
bx.invert_yaxis()
bx.tick_params(left=False)
bx.set_xlabel("occupied − matched control  (reference SD)")
bx.set_title("B   Hiding the motif collapses only the sequence-only model", pad=10)
bx.legend(frameon=False, fontsize=9.5, loc="upper left")

values["B_masking"] = {}
for name, vis, hid, _, _ in MODELS:
    change = masking(vis, hid)
    values["B_masking"][name] = {
        "visible_sd": scoring(vis)["mean_difference_sd"],
        "hidden_sd": scoring(hid)["mean_difference_sd"],
        "change_logodds": None if not change else change["mean"],
        "ci95_logodds": None if not change else [change["ci_low"], change["ci_high"]],
        "p": None if not change else change["p"]}

# =============================================================================
# C — the design result: sequon retention
# =============================================================================
designers = [(name, retention(ret), colour)
             for name, _, _, ret, colour in MODELS if retention(ret) is not None]
absent = [name for name, _, _, ret, _ in MODELS if retention(ret) is None]

for i, (name, r, colour) in enumerate(designers):
    cx.plot([r["control_mean"], r["occupied_mean"]], [i, i],
            color=colour, lw=2.0, alpha=0.55, zorder=1)
    cx.plot(r["control_mean"], i, "o", color="white", ms=10, zorder=3,
            markeredgecolor=DIM, markeredgewidth=2.2,
            label="matched control" if i == 0 else None)
    cx.plot(r["occupied_mean"], i, "o", color=colour, ms=10, zorder=3,
            label="occupied" if i == 0 else None)
    if r["excludes_zero"]:
        cx.text(r["occupied_mean"] + 0.012, i, "*", va="center",
                fontsize=15, color=INK)
for j, name in enumerate(absent):
    i = len(designers) + j
    cx.text(0.004, i, "no designer — scorer only", va="center", fontsize=9.5,
            color=DIM, style="italic")

cx.set_yticks(range(len(designers) + len(absent)))
cx.set_yticklabels([d[0] for d in designers] + absent)
cx.set_ylim(-0.6, len(designers) + len(absent) - 0.4)
cx.set_xlim(left=-0.005)
cx.invert_yaxis()
cx.tick_params(left=False)
cx.set_xlabel("sequon retained in redesign  (proportion of sites)")
cx.set_title("C   Design: the sequon is usually not kept, by either model", pad=10)
cx.legend(frameon=False, fontsize=9.5, loc="upper right")

values["C_retention"] = {
    name: {"occupied": r["occupied_mean"], "control": r["control_mean"],
           "paired_difference": r["paired_difference"], "ci95": r["ci95"],
           "n_pairs": r["n_pairs"], "n_tied": r["n_tied"],
           "wilcoxon_p": r["wilcoxon_p"], "ci_excludes_zero": r["excludes_zero"]}
    for name, r, _ in designers}
values["C_retention_absent"] = absent

# =============================================================================
# D — the per-site contrasts panel A averages
# =============================================================================
for name, vis, _, _, colour in MODELS:
    stat = scoring(vis)
    contrasts = pd.read_csv(ANALYSIS / f"contrasts_secretory_{vis}.csv")
    standardised = np.sort(contrasts.contrast.values / stat["reference_sd"])
    ecdf = np.arange(1, len(standardised) + 1) / len(standardised)
    dx.step(standardised, ecdf, where="post", color=colour, lw=2.0, label=name)
dx.axvline(0, color=INK, lw=1.0, zorder=1)
dx.set_xlim(-4, 4)
dx.set_ylim(0, 1)
dx.set_xlabel("per-site contrast  (reference SD)")
dx.set_ylabel("cumulative proportion of sites")
dx.set_title("D   The averages come from broad, overlapping distributions", pad=10)
dx.legend(frameon=False, fontsize=9.5, loc="upper left")

values["D_contrasts"] = {}
for name, vis, _, _, _ in MODELS:
    stat = scoring(vis)
    c = pd.read_csv(ANALYSIS / f"contrasts_secretory_{vis}.csv")
    s = c.contrast.values / stat["reference_sd"]
    values["D_contrasts"][name] = {
        "n": int(len(s)), "reference_sd": stat["reference_sd"],
        "median_sd": float(np.median(s)),
        "fraction_positive": float((s > 0).mean())}

fig.tight_layout(pad=2.0)
path = OUT / "fig_model_comparison.png"
fig.savefig(path, dpi=200, bbox_inches="tight")
print("wrote", path)

values_path = OUT / "fig_model_comparison_values.json"
values_path.write_text(json.dumps(values, indent=2))
print("wrote", values_path)


# =============================================================================
# Second figure — the raw score distributions, one panel per model
# =============================================================================
# Raw log-odds, NOT standardised, because that is what each model actually
# emits. The consequence is that the three panels cannot share an axis: eight
# averaged decoding orders, one causal pass and a sequence-only masked position
# produce numbers on genuinely different scales. Each panel is therefore its own
# comparison, and the panels must not be read against each other.
#
# The matched pairs rather than the whole populations. Every case here has
# exactly one control, so nothing is averaged and no spread is compressed --
# and matching is what stops this being the confounded comparison the archived
# population analysis showed collapsing once composition was controlled.
fig2, axes2 = plt.subplots(1, 3, figsize=(14.4, 4.6))

for panel, (name, vis, _, _, colour) in zip(axes2, MODELS):
    stat = scoring(vis)
    pairs = pd.read_csv(ANALYSIS / f"contrasts_secretory_{vis}.csv")
    occupied = pairs.case_score.values
    unknown = pairs.control_mean_score.values

    edges = np.histogram_bin_edges(np.concatenate([occupied, unknown]), bins=26)
    panel.hist(unknown, bins=edges, color=DIM, alpha=0.55, density=True,
               label="no annotated glycan")
    panel.hist(occupied, bins=edges, color=colour, alpha=0.62, density=True,
               label="experimentally occupied")
    panel.axvline(np.median(unknown), color=DIM, lw=1.8, ls="--")
    panel.axvline(np.median(occupied), color=colour, lw=1.8, ls="--")

    if stat["ci95_sd"][0] > 0:
        panel.text(0.97, 0.97, "*", transform=panel.transAxes, ha="right",
                   va="top", fontsize=17, color=INK)
    panel.set_title(f"{name}", pad=8)
    panel.set_xlabel("conditional sequon score  (log odds)")
    panel.set_yticks([])
    panel.spines["left"].set_visible(False)

    values.setdefault("E_raw_distributions", {})[name] = {
        "n_pairs": int(len(pairs)),
        "occupied_median": float(np.median(occupied)),
        "unknown_median": float(np.median(unknown)),
        "occupied_mean": float(occupied.mean()),
        "unknown_mean": float(unknown.mean()),
        "overlap_fraction_occupied_below_unknown_median":
            float((occupied < np.median(unknown)).mean()),
    }

axes2[0].set_ylabel("density of sites")
axes2[0].legend(frameon=False, fontsize=9.5, loc="upper left")
fig2.tight_layout(pad=1.8)
path2 = OUT / "fig_score_distributions.png"
fig2.savefig(path2, dpi=200, bbox_inches="tight")
print("wrote", path2)

values_path.write_text(json.dumps(values, indent=2))
