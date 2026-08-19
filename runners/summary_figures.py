"""Six simple figures, one claim each.

Deliberately plain. Each answers a single question a reader would ask, and
nothing is combined that could be read separately. Regenerate with:

    python runners/summary_figures.py
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

def save(fig, name, caption, width=112):
    """Grow the canvas by exactly the caption's height and reserve that strip.

    Working in inches rather than figure fractions keeps the plot area the size
    it was designed at, however many lines the caption runs to.
    """
    lines = textwrap.wrap(caption, width=width)
    caption_inches = 0.155 * len(lines) + 0.18
    width_in, height_in = fig.get_size_inches()
    fig.set_size_inches(width_in, height_in + caption_inches)
    fig.tight_layout(rect=[0, caption_inches / (height_in + caption_inches), 1, 1])
    fig.text(0.012, 0.012, "\n".join(lines), fontsize=7.6, color="#5c6570",
             va="bottom", linespacing=1.5)
    fig.savefig(OUT / name, dpi=200)
    plt.close(fig); print("wrote", OUT / name)

feat = pd.read_csv("results/site_structural_features.csv", low_memory=False)
man = pd.read_csv("results/candidate_manifest_dataset.csv", low_memory=False)
opt = json.loads(Path("results/analysis_optimal.json").read_text())

# ----------------------------------------------------------------- 1. funnel
fig, ax = plt.subplots(figsize=(8.2, 4.0))
occ = [922, int((feat.features_available & (feat.occupancy_status == "occupied_supported")).sum()),
       int((man.scoreable & (man.occupancy_status == "occupied_supported")).sum()), 16]
ctl = [32, int((feat.features_available & (feat.occupancy_status == "observed_unmodified")).sum()),
       int((man.scoreable & (man.occupancy_status == "observed_unmodified")).sum()), 16]
stages = ["evidence\nlayer", "mapped to\na structure", "ProteinMPNN\ncan decode", "survives\nmatching"]
y = np.arange(len(stages))[::-1]
ax.barh(y + 0.19, occ, height=0.34, color=OCC, label="occupied")
ax.barh(y - 0.19, ctl, height=0.34, color=CTL, label="internal control")
for yy, a, b in zip(y, occ, ctl):
    ax.text(a * 1.15, yy + 0.19, f"{a:,}", va="center", fontsize=9, color=OCC)
    ax.text(b * 1.15, yy - 0.19, f"{b:,}", va="center", fontsize=9, color=CTL)
ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(stages)
ax.set_xlim(8, 4000); ax.set_xlabel("sites (log scale)")
ax.set_title("Where the sites go: 32 internal controls become 16 usable pairs")
ax.legend(frameon=False, loc="lower right", fontsize=9)
save(fig, "fig1_attrition.png",
     "Log scale. The occupied class is never limiting. The study is bounded by the internal "
     "controls: 32 exist, 28 have a backbone ProteinMPNN will decode, and 16 find a partner "
     "inside the matching caliper.")

# ------------------------------------------------------- 2. three comparisons
fig, ax = plt.subplots(figsize=(8.6, 3.4))
rows = [("Internal control\n(primary)", "optimal", OCC),
        ("Eukaryotic secretory\n(parallel)", "secretory", "#3f7d5a"),
        ("Bacterial\n(diagnostic)", "bacterial", MUTE),
        ("Cytosolic\n(diagnostic)", "cytosolic", MUTE)]
margin = opt["margin_standardised"]
ax.axvspan(-margin, margin, color=MUTE, alpha=0.18, zorder=0)
ax.axvline(0, color=INK, lw=0.9, zorder=1)
for i, (label, key, colour) in enumerate(rows):
    d = json.loads(Path(f"results/analysis_{key}.json").read_text())
    m, (lo, hi) = d["mean_difference_sd"], d["ci95_sd"]
    ax.errorbar([m], [-i], xerr=[[m - lo], [hi - m]], fmt="o", color=colour,
                capsize=4, lw=1.7, markersize=7, zorder=3)
    ax.text(hi + 0.07, -i, f"{m:+.3f}  (n={d['n_contrasts']})", va="center",
            fontsize=8.8, color=colour)
ax.set_yticks([-i for i in range(len(rows))])
ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel("occupied − control (SD units)"); ax.set_xlim(-0.85, 1.65)
ax.set_title("Four comparisons: the best-powered one sits on zero")
ax.set_ylim(-len(rows) + 0.4, 0.6)
save(fig, "fig2_three_comparisons.png",
     "Shaded band is the pre-specified equivalence margin. The internal-control comparison "
     "has the strongest negative label but only 16 pairs, so its interval is enormous. The "
     "eukaryotic secretory set matches the occupied sites on both taxonomy and compartment "
     "and supplies 262 pairs, at the cost of a weaker label; it lands almost exactly on zero "
     "and is the best-powered estimate available. The two diagnostics are confounded by "
     "construction. An earlier reading in which all comparisons were negative and shrank as "
     "matching improved was an artefact of a scoring defect and is withdrawn.")

# ------------------------------------------------- 3. matching sensitivity
sweep = pd.read_csv("results/matching_seed_sweep.csv")
sens = json.loads(Path("results/matching_sensitivity.json").read_text())
sd = sens["reference_sd"]
fig, ax = plt.subplots(figsize=(8.2, 3.5))
ax.hist(sweep["mean"] / sd, bins=26, color=MUTE, edgecolor="white", linewidth=0.6)
ax.axvline(0, color=INK, lw=1.2)
ax.axvline(opt["mean_difference_sd"], color=OCC, lw=2.0,
           label=f"deterministic optimum ({opt['mean_difference_sd']:+.3f})")
ax.set_xlabel("mean occupied − control (SD units), one value per greedy seed")
ax.set_ylabel("seeds")
ax.set_title("Every seed gives a positive estimate; only 38% exclude zero")
ax.legend(frameon=False, fontsize=8.8, loc="upper left")
save(fig, "fig3_matching_sensitivity.png",
     f"200 greedy matchings, one per seed. All 200 point estimates are positive, spanning "
     f"{sweep['mean'].min()/sd:+.3f} to {sweep['mean'].max()/sd:+.3f} SD, but the confidence "
     f"interval excludes zero in only {int(sweep.excludes_zero.sum())} of 200. Whether the "
     "result 'reached significance' was being decided by the seed, which is why matching is "
     "now deterministic.")

# ------------------------------------------------------- 4. retention bridge
ret = pd.read_csv("results/mpnn_retention_frozen_2026-08-18.csv", low_memory=False)
able = pd.read_csv("results/scoreability.csv", low_memory=False)
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
score = pd.concat([pd.read_csv("results/mpnn_conditional_scores.csv", low_memory=False),
                   pd.read_csv("results/mpnn_conditional_scores_unmatched.csv",
                               low_memory=False)],
                  ignore_index=True).drop_duplicates(KEY)
for d in (ret, able, score):
    for k in KEY: d[k] = d[k].astype(str)
score["valid"] = ~((score.p_ser_or_thr_at_plus2 > 1.0)
                   | (score.probs_n.map(lambda v: abs(sum(json.loads(v)) - 1) > 1e-3)))
ret = ret.merge(able[KEY + ["scoreable"]], on=KEY, how="left")
ret = ret[ret.scoreable == True]
merged = ret.merge(score[score.valid][KEY + ["conditional_sequon_score"]], on=KEY, how="inner")
FULL = "std_frac_full_sequon_retained"
q = pd.qcut(merged.conditional_sequon_score, 5, labels=False, duplicates="drop")
grp = merged.groupby(q).agg(score=("conditional_sequon_score", "mean"),
                            ret=(FULL, "mean"), n=(FULL, "size"))
fig, ax = plt.subplots(figsize=(7.4, 3.6))
labels = [f"Q{i+1}\nmean score {v:+.1f}" for i, v in zip(grp.index, grp.score)]
bars = ax.bar(labels, grp.ret, color=OCC, width=0.62)
for b, r in zip(bars, grp.ret):
    ax.text(b.get_x() + b.get_width()/2, r + 0.008, f"{r:.3f}", ha="center", fontsize=9)
ax.set_ylabel("fraction of designs keeping the sequon")
ax.set_xlabel("conditional sequon score, low to high quintile")
ax.set_ylim(0, 0.36)
ax.set_title("The conditional score predicts what the model actually does")
save(fig, "fig4_retention_bridge.png",
     f"n = {len(merged):,} scoreable sites, 32 unconstrained designs each. Spearman "
     f"+0.547. This is the bridge between the site-level probability the model holds and "
     "the sequences it generates: the two describe one underlying quantity.")

# -------------------------------------------------- 5. retention distribution
fig, ax = plt.subplots(figsize=(7.4, 3.4))
ax.hist(ret[FULL], bins=33, color=CTL, edgecolor="white", linewidth=0.5)
ax.set_yscale("log")
ax.set_xlabel("fraction of 32 designs keeping the sequon (per site)")
ax.set_ylabel("sites (log scale)")
ax.set_title(f"{100*(ret[FULL]==0).mean():.1f}% of sequons are lost in every single design")
save(fig, "fig5_retention_distribution.png",
     f"n = {len(ret):,} sites. Mean retention {ret[FULL].mean():.3f}; "
     f"{int((ret[FULL]==1).sum())} sites ({100*(ret[FULL]==1).mean():.1f}%) keep the sequon in "
     "all 32. This replicates the preprint's finding at site level: ProteinMPNN removes "
     "natural glycosylation sequons as a matter of course.")

# ------------------------------------------------------------ 6. the defect
fig, ax = plt.subplots(figsize=(7.8, 3.4))
ax.hist(score[score.valid].conditional_sequon_score, bins=50, color=OCC,
        edgecolor="white", linewidth=0.4, label=f"valid ({int(score.valid.sum()):,})")
ax.hist(score[~score.valid].conditional_sequon_score, bins=8, color=WARN,
        edgecolor="white", linewidth=0.4, label=f"invalid ({int((~score.valid).sum())})")
ax.set_yscale("log"); ax.set_xlabel("conditional sequon score (log-odds)")
ax.set_ylabel("sites (log scale)")
ax.set_title(f"The scoring defect: {int((~score.valid).sum())} sites returned "
             f"an impossible +13.8")
ax.legend(frameon=False, fontsize=8.8)
save(fig, "fig6_scorer_defect.png",
     "ProteinMPNN never decodes a residue whose backbone is incomplete; it returns a row of "
     "zeros that exponentiates to twenty-one ones, giving P(N)=1 and P(S)+P(T)=2. Only 4% of "
     "sites were affected, but they inflated the reference SD from 1.33 to 2.62 — every "
     "standardised effect was divided by that number, and the first result had the wrong sign.")

# ------------------------------------------------- 7. retention by class
# The best-powered comparison in the project, and it had no figure. Uncertainty
# is bootstrapped over PROTEINS, not sites: several sequons on one protein share
# a structure and a set of designs, so treating them as independent would make
# every interval far too narrow.
man_all = pd.read_csv("results/scoring_manifest.csv", low_memory=False).drop_duplicates(KEY)
for k in KEY:
    man_all[k] = man_all[k].astype(str)
by_class = (ret.merge(man_all[KEY + ["occupancy_status"]], on=KEY, how="left")
               .dropna(subset=["occupancy_status"]))

ORDER = [("occupied_supported", "occupied\n(experimental glycan)", OCC),
         ("control_bacterial_extracytoplasmic", "bacterial\ncontrol", MUTE),
         ("control_cytosolic_eukaryotic", "cytosolic\ncontrol", MUTE),
         ("observed_unmodified", "internal\ncontrol", CTL)]

def protein_bootstrap(frame, n_boot=4000, seed=11):
    """Resample whole proteins; a protein's sequons share one set of designs."""
    rng = np.random.default_rng(seed)
    groups = [g[FULL].to_numpy() for _, g in frame.groupby("accession")]
    if len(groups) < 2:
        return np.nan, np.nan
    draws = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        draws[i] = np.concatenate([groups[j] for j in pick]).mean()
    return np.percentile(draws, [2.5, 97.5])

fig, ax = plt.subplots(figsize=(7.6, 3.8))
tick_labels = []
for i, (key, label, colour) in enumerate(ORDER):
    grp = by_class[by_class.occupancy_status == key]
    if grp.empty:
        continue
    mean = grp[FULL].mean()
    lo, hi = protein_bootstrap(grp)
    ax.bar(i, mean, color=colour, width=0.6)
    if np.isfinite(lo):
        ax.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="none",
                    ecolor=INK, capsize=5, lw=1.4)
    ax.text(i, (hi if np.isfinite(hi) else mean) + 0.004,
            f"{mean:.3f}", ha="center", fontsize=9.5)
    tick_labels.append(f"{label}\nn={len(grp):,}")

ax.set_xticks(range(len(tick_labels))); ax.set_xticklabels(tick_labels, fontsize=9)
ax.set_ylabel("fraction of 32 designs keeping the sequon")
ax.set_ylim(0, max(by_class.groupby("occupancy_status")[FULL].mean()) * 1.9)
ax.set_title("Occupied sequons are destroyed at the same rate as unoccupied ones")
save(fig, "fig7_retention_by_class.png",
     "Mean per-site retention by class, with 95% intervals bootstrapped over proteins "
     "rather than sites. Occupied sites and the two large control sets are "
     "indistinguishable: whether a sequon actually carries a glycan makes no difference "
     "to how often ProteinMPNN removes it. This is the best-powered comparison in the "
     "project - roughly a thousand sites per control group against sixteen matched pairs "
     "for the primary analysis - and it points the same way. The internal-control bar is "
     "lower but rests on 21 sites.")
