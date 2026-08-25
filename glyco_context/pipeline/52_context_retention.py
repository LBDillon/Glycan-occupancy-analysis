"""Does fixing the motif protect the biology around it, or only the three letters?

Scores each design's local chemistry against the distribution of natural
occupied sites, and compares that to where its own wild type sits. The paired
quantity is what matters:

    dD = D(design) - D(wild type)

positive meaning the design has moved away from natural occupied context.

The reference for a site excludes its own protein, so nothing is scored partly
against itself. Designs of one chain are replicates of a single draw, so they
are averaged within site before proteins are resampled.

Usage:
    52_context_retention.py [--panels ...] [--reference ...] [--boot 2000]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.context_distance import context_distance, reference_moments
from glyco_context.context_stats import benjamini_hochberg
from glyco_context.local_chemistry import CLASSES
from experimental_glycosylation_sites.provenance import hash_file, _git_state

PANEL = ([f"flank_{c}_fraction" for c in CLASSES]
         + [f"shell_{c}_fraction" for c in CLASSES] + ["shell_net_charge"])

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--panels", default="glyco_context/results/analysis/fixed_sequon_panels.csv")
parser.add_argument("--reference", default="glyco_context/results/analysis/natural_reference_panels.csv")
parser.add_argument("--boot", type=int, default=2000)
parser.add_argument("--out", default="glyco_context/results/analysis")
args = parser.parse_args()

panels = pd.read_csv(args.panels, low_memory=False)
reference = pd.read_csv(args.reference, low_memory=False)
reference = reference[reference.variant == "wild_type"]
print(f"reference: {len(reference)} natural occupied sites, "
      f"{reference.accession.nunique()} proteins")
print(f"test set : {panels[panels.variant=='wild_type'].shape[0]} sites, "
      f"{panels.accession.nunique()} proteins, variants {sorted(panels.variant.unique())}")

# --- distance from natural context, protein held out ------------------------
rows = []
for (accession, position), group in panels.groupby(["accession", "position"]):
    mu, sigma = reference_moments(reference, PANEL, exclude_accession=accession)
    for row in group.itertuples():
        rows.append({"accession": accession, "position": position,
                     "variant": row.variant, "replicate": row.replicate,
                     "distance": context_distance(row._asdict(), mu, sigma, PANEL)})
scored = pd.DataFrame(rows)

# Designs of one chain are replicates of a single draw, not separate sites.
per_site = scored.pivot_table(index=["accession", "position"], columns="variant",
                              values="distance", aggfunc="mean")
per_site = per_site.dropna(subset=["wild_type"])
print(f"\nsites scored: {len(per_site)}")

def bootstrap(values: pd.Series, proteins: pd.Series, n_boot: int, seed: int = 0):
    values = values.dropna()
    if len(values) < 3:
        return {"n": int(len(values)), "mean": np.nan, "ci_low": np.nan,
                "ci_high": np.nan, "p": np.nan}
    groups = [g.to_numpy(float) for _, g in values.groupby(proteins.loc[values.index])]
    rng = np.random.default_rng(seed)
    draws = np.array([np.concatenate([groups[i] for i in
                                      rng.integers(0, len(groups), len(groups))]).mean()
                      for _ in range(n_boot)])
    below, above = float((draws <= 0).mean()), float((draws >= 0).mean())
    return {"n": int(len(values)), "mean": float(values.mean()),
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)),
            "p": float(min(1.0, max(2 * min(below, above), 1.0 / n_boot)))}

proteins = pd.Series(per_site.index.get_level_values("accession"), index=per_site.index)
summary = {}
print(f"\n{'quantity':34}{'n':>5}{'mean':>9}{'95% CI':>20}{'p':>9}")
for name, series in [
        ("D(wild type)", per_site.get("wild_type")),
        ("D(design)", per_site.get("design")),
        ("D(random control)", per_site.get("random")),
        ("dD = design - wild type", per_site.get("design", pd.Series(dtype=float)) - per_site.wild_type),
        ("dD = random - wild type", per_site.get("random", pd.Series(dtype=float)) - per_site.wild_type),
        ("design - random (same n mutations)",
         per_site.get("design", pd.Series(dtype=float)) - per_site.get("random", pd.Series(dtype=float)))]:
    if series is None or not len(series.dropna()):
        continue
    stats = bootstrap(series, proteins, args.boot)
    summary[name] = stats
    print(f"{name:34}{stats['n']:5d}{stats['mean']:9.3f}"
          f"  [{stats['ci_low']:7.3f},{stats['ci_high']:7.3f}]{stats['p']:9.4f}")

# --- which features move ----------------------------------------------------
print(f"\nper-feature shift, design minus wild type (standardised on the reference):")
feature_rows = []
wt = panels[panels.variant == "wild_type"].set_index(["accession", "position"])
des = panels[panels.variant == "design"].groupby(["accession", "position"])[PANEL].mean()
for feature in PANEL:
    if feature not in wt.columns:
        continue
    mu, sigma = reference_moments(reference, [feature], exclude_accession=None)
    spread = sigma.get(feature)
    if not spread or not np.isfinite(spread) or spread <= 0:
        continue
    common = wt.index.intersection(des.index)
    shift = ((des.loc[common, feature] - wt.loc[common, feature]) / spread).dropna()
    if len(shift) < 3:
        continue
    stats = bootstrap(shift, pd.Series(common.get_level_values("accession"), index=common),
                      args.boot)
    feature_rows.append({"feature": feature, **stats})
features = pd.DataFrame(feature_rows)
if len(features):
    features["q"] = benjamini_hochberg(features.p.tolist())
    # 'mean' does not survive itertuples as an attribute; name it for what it is.
    features = features.rename(columns={"mean": "shift"})
    for row in features.reindex(features["shift"].abs().sort_values(ascending=False).index).itertuples():
        flag = "*" if row.q < 0.05 else " "
        print(f"  {row.feature:34}{row.shift:+8.3f}  [{row.ci_low:+7.3f},{row.ci_high:+7.3f}]"
              f"  q={row.q:.4f}{flag}")

outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
scored.to_csv(outdir / "context_retention_distances.csv", index=False)
if len(features):
    features.to_csv(outdir / "context_retention_features.csv", index=False)
(outdir / "context_retention_summary.json").write_text(json.dumps({
    "prespecification": "glyco_context/docs/prespecification_fixed_sequon_context_retention.md",
    "summary": summary, "n_boot": args.boot,
    "reference_sites": int(len(reference)), "git": _git_state()}, indent=2, default=str))
print(f"\nwrote {(outdir / 'context_retention_distances.csv').resolve()}")
