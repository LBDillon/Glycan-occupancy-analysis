"""What each model predicts at the three sequon positions, as a heatmap.

Two figures.

  fig_sequon_heatmap        the mean predicted distribution at each of the three
                            sequon positions, per model
  fig_sequon_heatmap_diff   the PAIRED difference, occupied minus its matched
                            partner, so it shows which amino acids the occupied
                            and unannotated arms are actually separated by

The score files already hold a full distribution over amino acids at the N, X and
S/T positions -- that is what `probs_n`, `probs_plus1` and `probs_plus2` are. This
draws them, which explains what the score measures better than the prose does.

The difference figure is PAIRED, using `matched_pairs_<label>.csv`. Taking two
population means instead would be the confounded comparison the archived context
analysis showed collapsing once composition was controlled.

Each model's own alphabet is used to find the columns, and the constant is
checked against the columns the scorer already wrote rather than trusted.

Usage:  26_figures_sequon_heatmap.py [label]      (default: secretory)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, "src")

LABEL = sys.argv[1] if len(sys.argv) > 1 else "secretory"
# Grouped by side-chain chemistry rather than alphabetically, so neighbouring
# columns are residues that could plausibly substitute for one another. The
# groups are the conventional ones; that S, T and N end up adjacent is a
# property of the polar-uncharged group, not an arrangement chosen to flatter
# the result.
GROUPS = [("hydrophobic", "AVLIM"), ("aromatic", "FWY"),
          ("polar", "STNQ"), ("special", "CGP"),
          ("basic", "HKR"), ("acidic", "DE")]
ORDER = "".join(letters for _, letters in GROUPS)
SEQUON = set("NST")                 # the residues the score actually reads
assert len(ORDER) == 20 and len(set(ORDER)) == 20
SCORES, MATCHING = Path("results/scores"), Path("results/matching")
OUT = Path("results/figures")
OUT.mkdir(parents=True, exist_ok=True)

COLUMNS = [("probs_n", "1 (Asn)"), ("probs_plus1", "2 (X)"),
           ("probs_plus2", "3 (Ser/Thr)")]
KEY = ["accession", "position"]


def carbonara_alphabet():
    from experimental_glycosylation_sites.carbonara_scoring import ALPHABET
    return list(ALPHABET)


def mpnn_alphabet():
    from experimental_glycosylation_sites.mpnn_scoring import ALPHABET
    return list(ALPHABET)


def esmif_alphabet():
    from esm.data import Alphabet
    # The token LIST, not a joined string: joining would let "N".index() match a
    # letter inside a multi-character token such as <null_1>.
    return list(Alphabet.from_architecture("invariant_gvp").all_toks)


# EvolutionaryScale's `esm` cannot be installed beside `fair-esm` -- both claim
# the import name -- so ESMC's tokens are recorded here rather than read live.
# Dumped from EsmSequenceTokenizer and verified: N lands at 17, S at 8, T at 11,
# each reproducing the p_asn_at_n / p_ser_at_plus2 / p_thr_at_plus2 columns the
# scorer wrote. The stored vectors are 64 wide because ESM-C's output head is
# padded; columns 33-63 are identically zero. The assertion in `load` re-checks
# this on every run, so a drift upstream stops the figure rather than shifting it.
ESMC_TOKENS = ["<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E",
               "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y", "M", "H", "W",
               "C", "X", "B", "U", "Z", "O", ".", "-", "|", "<mask>"]


def esmc_alphabet():
    return list(ESMC_TOKENS)


MODELS = [("ProteinMPNN", "proteinmpnn_index_corrected", mpnn_alphabet),
          ("ESM-IF", "esm_if_index_corrected", esmif_alphabet),
          ("ESMC", "esmc_single", esmc_alphabet),
          ("CARBonAra", "carbonara", carbonara_alphabet)]


def load(variant, alphabet):
    """(accession, position) -> [3, 20] of probabilities in ORDER."""
    path = SCORES / f"scores_{{}}_{variant}.csv"
    frames = []
    for tag in ("dataset", LABEL):
        p = Path(str(path).format(tag))
        if p.exists():
            frames.append(pd.read_csv(p, low_memory=False))
    if not frames:
        return None
    d = pd.concat(frames, ignore_index=True).drop_duplicates(KEY)

    take = [alphabet.index(aa) for aa in ORDER]
    grids = {}
    for _, r in d.iterrows():
        rows = [np.asarray(json.loads(r[col]), dtype=float)[take]
                for col, _ in COLUMNS]
        grids[(r.accession, int(r.position))] = np.vstack(rows)
    # the alphabet is a claim about column order; check it against what the
    # scorer already wrote for a column it names explicitly
    first = d.iloc[0]
    assert abs(np.asarray(json.loads(first.probs_n))[alphabet.index("N")]
               - first.p_asn_at_n) < 1e-5, f"{variant}: probs_n[N] != p_asn_at_n"
    return grids


def cluster_ci(values, units, n_boot=2000, seed=0):
    """Bootstrap over resampling units, not sites. Returns (mean, low, high)."""
    values, units = np.asarray(values), np.asarray(units)
    groups = [values[units == u] for u in np.unique(units)]
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for b in range(n_boot):
        picked = rng.integers(0, len(groups), len(groups))
        draws[b] = np.concatenate([groups[i] for i in picked]).mean()
    return values.mean(), np.percentile(draws, 2.5), np.percentile(draws, 97.5)


pairs = pd.read_csv(MATCHING / f"matched_pairs_{LABEL}.csv", low_memory=False)
contrasts = pd.read_csv(f"results/analysis/contrasts_{LABEL}_esm_if_index_corrected.csv")
unit_of = {(r.case_accession, int(r.case_position)): r.resample_unit
           for r in contrasts.itertuples()}

panels, difference, notes = [], [], []
for name, variant, get_alphabet in MODELS:
    try:
        alphabet = get_alphabet()
    except Exception as exc:
        notes.append(f"{name}: alphabet unavailable ({type(exc).__name__})")
        continue
    grids = load(variant, alphabet)
    if grids is None:
        notes.append(f"{name}: no score files")
        continue

    occupied_arm, control_arm, deltas, units = [], [], [], []
    for r in pairs.itertuples():
        case = grids.get((r.case_accession, int(r.case_position)))
        control = grids.get((r.control_accession, int(r.control_position)))
        if case is None or control is None:
            continue
        occupied_arm.append(case)
        control_arm.append(control)
        deltas.append(case - control)
        units.append(unit_of.get((r.case_accession, int(r.case_position)),
                                 r.case_accession))
    if not deltas:
        notes.append(f"{name}: no pair had both arms scored")
        continue

    deltas = np.stack(deltas)
    mean_delta = deltas.mean(axis=0)
    significant = np.zeros_like(mean_delta, dtype=bool)
    for i in range(mean_delta.shape[0]):
        for j in range(mean_delta.shape[1]):
            _, low, high = cluster_ci(deltas[:, i, j], units)
            significant[i, j] = (low > 0) or (high < 0)

    panels.append((name, np.stack(occupied_arm).mean(axis=0),
                   np.stack(control_arm).mean(axis=0), significant, len(deltas)))
    difference.append((name, mean_delta, significant, len(deltas)))

for n in notes:
    print(n)


ARMS = ("occupied", "no annotated glycan")
# Boundaries between chemistry groups, in category-index space.
BOUNDARIES = np.cumsum([len(letters) for _, letters in GROUPS])[:-1] - 0.5
TICKS = [f"<b>{a}</b>" if a in SEQUON else a for a in ORDER]


def draw(panels, filename):
    """One row per model, the two arms side by side.

    Amino acids across and positions down, which is the transpose of the usual
    inverse-folding heatmap. That layout puts hundreds of residue positions on
    the x axis; here there are three, so the panel would be a tall thin sliver
    the wrong way round.
    """
    rows = len(panels)
    # The arm is named once, at the top of its column, rather than on every row.
    titles = list(ARMS) + [""] * (2 * (rows - 1))
    fig = make_subplots(rows=rows, cols=2, shared_xaxes=True, shared_yaxes=True,
                        subplot_titles=titles,
                        horizontal_spacing=0.055, vertical_spacing=0.075)
    limit = max(max(occ.max(), ctl.max()) for _, occ, ctl, _, _ in panels)
    ylabels = [c[1] for c in COLUMNS]

    for r, (name, occupied, control, significant, n) in enumerate(panels, start=1):
        for c, grid in ((1, occupied), (2, control)):
            fig.add_trace(go.Heatmap(
                z=grid, x=list(ORDER), y=ylabels, zmin=0, zmax=limit,
                colorscale="Blues", showscale=(r == 1 and c == 2),
                colorbar=dict(title="probability", thickness=14, len=0.5, y=0.8),
                hovertemplate="%{x} at position %{y}<br>%{z:.4f}<extra></extra>"),
                row=r, col=c)
            for boundary in BOUNDARIES:
                fig.add_vline(x=boundary, line=dict(color="white", width=2),
                              row=r, col=c)
        ys, xs = np.where(significant)
        if len(ys):
            fig.add_trace(go.Scatter(
                x=[ORDER[i] for i in xs], y=[ylabels[i] for i in ys],
                mode="markers",
                marker=dict(symbol="circle-open", size=9, line=dict(width=2),
                            color="rgba(20,20,20,0.9)"),
                showlegend=False, hoverinfo="skip"), row=r, col=1)
        for c in (1, 2):
            fig.update_yaxes(autorange="reversed", row=r, col=c)
            fig.update_xaxes(tickmode="array", tickvals=list(ORDER),
                             ticktext=TICKS, row=r, col=c)
        # the model is named once too, down the left-hand side
        fig.update_yaxes(title_text=f"<b>{name}</b>", row=r, col=1)

    fig.update_layout(
        template="simple_white",
        title=dict(text="Predicted amino-acid distribution at the sequon",
                   x=0.02, xanchor="left"),
        height=210 * rows + 130, width=1220,
        margin=dict(t=105, l=110, r=90, b=70))
    for annotation in fig.layout.annotations[:2]:
        annotation.update(font=dict(size=15))
    fig.write_image(OUT / f"{filename}.png", scale=2)
    fig.write_html(OUT / f"{filename}.html", include_plotlyjs="cdn")
    print("wrote", OUT / f"{filename}.png", "and .html")


if panels:
    n = panels[0][4]
    # Everything explanatory lives in docs/figures_and_captions.md, so the
    # figure carries only its title, axes, key and significance markers.
    draw(panels, "fig_sequon_heatmap")
    print(f"{n} matched pairs; caption in docs/figures_and_captions.md")

summary = {name: {"n_pairs": int(n),
                  "delta_asn_at_1": float(delta[0, ORDER.index("N")]),
                  "delta_ser_at_3": float(delta[2, ORDER.index("S")]),
                  "delta_thr_at_3": float(delta[2, ORDER.index("T")]),
                  "significant_cells": int(sig.sum())}
           for name, delta, sig, n in difference}
(OUT / "fig_sequon_heatmap_values.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
