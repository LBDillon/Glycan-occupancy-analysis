"""ARCHIVED figures for the comparative context analysis.

Three figures, each carrying one claim.

  fig1  the effects collapse once composition is matched, and why they were
        there in the first place
  fig2  what an occupied site actually looks like
  fig3  the single feature that survives every framing, and the one that does not

Deliberately not a gallery. The result is largely negative, and a figure that
made it look otherwise would be the most misleading thing in the repository.

Usage:  run from the repository root. See ../README.md for why this is archived.
"""
import json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5,
                     "axes.titlelocation": "left", "axes.spines.top": False,
                     "axes.spines.right": False})

INK, OCC, CTL, DIM = "#22252b", "#1f6f8b", "#b05c3b", "#9aa0a6"
DATA = Path("glyco_context/results/datasets")
ANALYSIS = Path("glyco_context/archive/comparative_analysis/results")
OUT = Path("glyco_context/archive/comparative_analysis/results")
OUT.mkdir(parents=True, exist_ok=True)

core = pd.read_csv(DATA / "context_triplet_core.csv", low_memory=False)
population = pd.read_csv(ANALYSIS / "context_differences.csv", low_memory=False)
paired = pd.read_csv(ANALYSIS / "context_paired_differences.csv", low_memory=False)
occ = core[core.population == "occupied"]

PRETTY = {
    "n_rsa": "exposure at Asn", "plus1_rsa": "exposure at +1",
    "plus2_rsa": "exposure at +2", "loop_run_length": "loop length",
    "n_neighbours_8a": "neighbours within 8 Å",
    "nd2_atoms_8a_same_chain": "atoms near ND2 (same chain)",
    "nd2_residues_8a_same_chain": "residues near ND2",
    "nd2_atoms_8a_other_chain": "atoms near ND2 (other chain)",
    "sidechain_neighbour_residues_5a": "side-chain neighbours (5 Å)",
    "neighbour_net_charge_8a": "local net charge",
    "neighbour_hydrophobic_fraction_8a": "hydrophobic fraction",
    "neighbour_aromatic_fraction_8a": "aromatic fraction",
    "nearest_aromatic_sidechain_nd2": "nearest aromatic side chain",
    "uniprot_residues_after_asn": "residues after Asn",
    "uniprot_residues_after_sequon": "residues after sequon",
    "distance_to_n_terminus_resolved": "distance to N-terminus",
    "distance_to_c_terminus_resolved": "distance to C-terminus",
    "n_ss_coarse==helix": "helix at Asn", "n_ss_coarse==sheet": "sheet at Asn",
    "n_ss_coarse==loop": "loop at Asn",
    "plus1_ss_coarse==helix": "helix at +1", "plus1_ss_coarse==sheet": "sheet at +1",
    "plus1_ss_coarse==loop": "loop at +1",
    "plus2_ss_coarse==helix": "helix at +2", "plus2_ss_coarse==sheet": "sheet at +2",
    "plus2_ss_coarse==loop": "loop at +2",
    "aromatic_within_8a==True": "aromatic within 8 Å",
}
MATCHED_AWAY = {"n_rsa", "n_neighbours_8a", "neighbour_hydrophobic_fraction_8a"}
label = lambda k: PRETTY.get(k, k)

# =============================================================================
# fig1 — the collapse, and the confound that produced it
# =============================================================================
pop = population[(population.comparison == "secretory_unannotated")].set_index("feature")
pair = paired[paired.comparison == "secretory"].set_index("feature")
rows = []
for f in pop.index:
    rows.append({"feature": f, "before": pop.loc[f, "estimate"],
                 "after": pair.loc[f, "standardised"] if f in pair.index else np.nan,
                 "q_after": pair.loc[f, "q"] if f in pair.index else np.nan,
                 "matched_away": f in MATCHED_AWAY})
comparison = pd.DataFrame(rows)
comparison = comparison.sort_values("before").reset_index(drop=True)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.8, 7.2),
                             gridspec_kw={"width_ratios": [2.0, 1]})
y = np.arange(len(comparison))
for i, r in enumerate(comparison.itertuples()):
    if np.isfinite(r.after):
        ax.plot([r.before, r.after], [i, i], color=DIM, lw=1.3, zorder=1)
ax.scatter(comparison.before, y, s=48, color=CTL, zorder=3,
           label="populations compared directly")
testable = comparison[np.isfinite(comparison.after)]
ax.scatter(testable.after, testable.index, s=48, color=OCC, zorder=3,
           label="matched pairs (composition controlled)")
away = comparison[comparison.matched_away]
if len(away):
    ax.scatter(away.before, away.index, s=150, facecolor="none", edgecolor=INK,
               lw=1.4, zorder=4, label="matched away — untestable here")
survivor = comparison[(comparison.q_after < 0.05)]
if len(survivor):
    ax.scatter(survivor.after, survivor.index, s=185, facecolor="none",
               edgecolor=OCC, lw=2.2, zorder=4, label="survives correction")
ax.axvline(0, color=INK, lw=1.0, zorder=2)
ax.set_yticks(y)
ax.set_yticklabels([label(f) for f in comparison.feature], fontsize=9.5)
ax.set_ylim(-1.0, len(comparison) - 0.3)
ax.set_xlabel("standardised difference (occupied − comparison)")
ax.set_title("Apparent context differences collapse\nwhen composition is matched",
             fontsize=12.5, pad=12)
ax.legend(frameon=False, fontsize=9, loc="lower right", ncol=1,
          borderpad=0.9)

# right: the confound, as a share so the two sets are comparable
sec = core[core.population == "secretory_unannotated"]
ic = core[core.population == "internal_control"]
bars = [("internal-control\ncomparison",
         len(set(ic.accession) & set(occ.accession)), ic.accession.nunique()),
        ("secretory-unannotated\ncomparison",
         len(set(sec.accession) & set(occ.accession)), sec.accession.nunique())]
for i, (name, sharedn, total) in enumerate(bars):
    share = 100 * sharedn / total
    bx.barh(i, 100, color="#eceff1", height=0.42)
    bx.barh(i, share, color=OCC, height=0.42)
    bx.text(2.5, i + 0.33, f"{sharedn} of {total} proteins  ({share:.0f}%)",
            fontsize=10, va="center", color=INK)
bx.set_yticks(range(len(bars)))
bx.set_yticklabels([b[0] for b in bars], fontsize=10)
bx.set_ylim(-0.6, len(bars) - 0.15)
bx.set_xlim(0, 100)
bx.set_xlabel("% of comparison proteins that also carry an occupied site")
bx.set_title("Why: occupancy was confounded\nwith protein identity", fontsize=12.5, pad=12)
fig.subplots_adjust(wspace=0.42)
fig.savefig(OUT / "fig1_effects_collapse.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig1_effects_collapse.png")

# =============================================================================
# fig3 — the same feature under four framings
# =============================================================================
FRAMINGS = [("populations, vs secretory", population, "secretory_unannotated", "estimate"),
            ("populations, vs internal", population, "internal_control", "estimate"),
            ("matched pairs, secretory", paired, "secretory", "standardised"),
            ("matched pairs, internal", paired, "internal", "standardised")]
FEATURES = ["plus2_ss_coarse==sheet", "nd2_atoms_8a_same_chain",
            "nearest_aromatic_sidechain_nd2", "n_ss_coarse==helix",
            "distance_to_n_terminus_resolved"]
fig, ax = plt.subplots(figsize=(10.4, 5.2))
width = 0.19
for j, (name, table, comparison_name, column) in enumerate(FRAMINGS):
    sub = table[table.comparison == comparison_name].set_index("feature")
    values = [sub.loc[f, column] if f in sub.index else np.nan for f in FEATURES]
    qs = [sub.loc[f, "q"] if f in sub.index else np.nan for f in FEATURES]
    xs = np.arange(len(FEATURES)) + (j - 1.5) * width
    colour = [OCC, "#5b8c5a", CTL, "#8a6fa8"][j]
    ax.bar(xs, values, width=width * 0.92, color=colour, label=name,
           alpha=1.0 if j >= 2 else 0.55,
           edgecolor=INK if j >= 2 else "none", lw=0.8)
    for x, v, q in zip(xs, values, qs):
        if np.isfinite(q) and q < 0.05:
            ax.text(x, v + (0.02 if v >= 0 else -0.05), "*", ha="center",
                    fontsize=12, color=INK)
ax.axhline(0, color=INK, lw=1.0)
ax.set_xticks(np.arange(len(FEATURES)))
ax.set_xticklabels([label(f).replace(" (same chain)", "\n(same chain)")
                    .replace("nearest aromatic side chain", "nearest aromatic\nside chain")
                    .replace("distance to N-terminus", "distance to\nN-terminus")
                    for f in FEATURES], fontsize=9.5)
ax.set_ylabel("standardised difference")
ax.set_title("One feature holds its direction under every framing; the large ones do not",
             fontsize=12)
ax.legend(frameon=False, fontsize=9, ncol=2)
fig.tight_layout()
fig.savefig(OUT / "fig3_framings.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig3_framings.png")
