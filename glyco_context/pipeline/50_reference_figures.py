"""The natural reference picture: what an occupied site looks like.

One figure. This is the distribution the fixed-sequon context-retention test
scores designs against, so it is kept with the reference rather than with any
particular comparison.

Usage:  50_reference_figures.py
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
ANALYSIS = Path("glyco_context/results/analysis")
OUT = Path("glyco_context/results/figures")
OUT.mkdir(parents=True, exist_ok=True)

core = pd.read_csv(DATA / "context_triplet_core.csv", low_memory=False)
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
# fig2 — what an occupied site looks like
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3))
ax = axes[0]
for name, column, colour in (("Asn", "n_rsa", OCC), ("+1", "plus1_rsa", "#5b8c5a"),
                             ("+2", "plus2_rsa", CTL)):
    ax.hist(occ[column].dropna(), bins=np.linspace(0, 1.2, 31), histtype="step",
            lw=1.9, color=colour, label=name)
ax.set_xlabel("relative solvent accessibility")
ax.set_ylabel("occupied sites")
ax.set_title(f"Exposure across the sequon (n={len(occ)})")
ax.legend(frameon=False, fontsize=9)

bx = axes[1]
ss = pd.DataFrame({p: occ[f"{p}_ss_coarse"].value_counts(normalize=True)
                   for p in ("n", "plus1", "plus2")}).fillna(0).T
order = [c for c in ("loop", "sheet", "helix", "unknown") if c in ss.columns]
bottom = np.zeros(len(ss))
for column, colour in zip(order, [OCC, CTL, "#5b8c5a", DIM]):
    bx.bar(range(len(ss)), ss[column] * 100, bottom=bottom, color=colour,
           label=column, width=0.62)
    bottom += ss[column].to_numpy() * 100
bx.set_xticks(range(len(ss)))
bx.set_xticklabels(["Asn", "+1", "+2"])
bx.set_ylabel("% of sites")
bx.set_title("Secondary structure")
bx.legend(frameon=False, fontsize=9, ncol=2)

cx = axes[2]
region = occ.assign(r=[None if pd.isna(a) or pd.isna(b) else
                       ("alpha_L" if a > 0 and -60 <= b <= 90 else
                        "other" if a > 0 else
                        "alpha_R" if -120 <= b <= 50 else "beta")
                       for a, b in zip(occ.n_phi, occ.n_psi)])
counts = region.r.value_counts()
cx.bar(range(len(counts)), counts.to_numpy(), color=[OCC, CTL, "#5b8c5a", DIM][:len(counts)],
       width=0.62)
cx.set_xticks(range(len(counts)))
cx.set_xticklabels([{"alpha_R": "α-right", "beta": "β", "alpha_L": "α-left",
                     "other": "other"}.get(i, i) for i in counts.index])
cx.set_ylabel("occupied sites")
cx.set_title("Backbone region at the Asn")
fig.savefig(OUT / "fig2_occupied_context.png", dpi=200, bbox_inches="tight")
print("wrote", OUT / "fig2_occupied_context.png")

