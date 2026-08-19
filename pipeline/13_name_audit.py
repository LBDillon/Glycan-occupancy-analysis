"""Are any 'unannotated' control proteins obviously glycoproteins by name?

The eukaryotic secretory control set rests entirely on the absence of a UniProt
glycoprotein keyword. That is a database fact, not a biochemical one, so it is
worth asking whether the names give the game away — a protein whose orthologue
is a known glycoprotein, or one from a family where N-glycosylation is the rule,
is a likely false negative arising from patchy curation rather than biology.

Three checks, then a sensitivity analysis with the flagged sites removed.

Lectins are deliberately NOT treated as suspect. Galectins bind glycans and are
cytosolic and bare, so flagging them would inflate the suspect count with
proteins that are correctly in the control set.
"""
import re, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.contrasts import (
    build_contrasts, classify, cluster_bootstrap)

KEY = ["accession", "position"]
STRICT_FAMILIES = (
    r"\b(immunoglobulin|antibody|Fc receptor|mucin|collagen|integrin|cadherin"
    r"|complement C[0-9]|coagulation factor|fibrinogen|von Willebrand"
    r"|histocompatibility|growth factor receptor)\b")


def gene_tokens(value):
    if not isinstance(value, str) or not value.strip():
        return set()
    return {g.upper() for g in re.split(r"[\s;]+", value.strip()) if g}


occ = pd.read_csv("results/datasets/occupied_protein_info.csv", low_memory=False)
ctl = pd.read_csv("results/datasets/secretory_unannotated_protein_info.csv", low_memory=False)
sites = pd.read_csv("results/datasets/secretory_unannotated_sites_raw.csv", low_memory=False)
ctl = ctl[ctl.Entry.astype(str).isin(set(sites.accession.astype(str)))].copy()

# Check 0 — did the keyword exclusion also remove automated annotation?
carb = [c for c in ctl.columns if "Glycosylation" in c or "Carbohyd" in c]
if carb:
    has = ctl[carb[0]].notna() & ctl[carb[0]].astype(str).str.strip().ne("") \
          & ctl[carb[0]].astype(str).ne("nan")
    print(f"control proteins with ANY CARBOHYD feature (manual or automated): "
          f"{int(has.sum())} of {len(ctl)}")

# Check 1 — same gene symbol as a protein in the occupied set
occ_genes = set().union(*(gene_tokens(v) for v in occ["Gene Names"].fillna("")))
ctl["shared_gene"] = ctl["Gene Names"].fillna("").map(lambda v: bool(gene_tokens(v) & occ_genes))

# Check 2 — family where the protein itself is normally N-glycosylated
ctl["suspect_family"] = ctl["Protein names"].fillna("").astype(str).str.contains(
    STRICT_FAMILIES, case=False, regex=True, na=False)

ctl["suspect"] = ctl.shared_gene | ctl.suspect_family
ctl[["Entry", "Protein names", "shared_gene", "suspect_family", "suspect"]].to_csv(
    "results/datasets/secretory_unannotated_name_audit.csv", index=False)
print(f"suspect by shared gene symbol : {int(ctl.shared_gene.sum())}")
print(f"suspect by protein family     : {int(ctl.suspect_family.sum())}")
print(f"suspect overall               : {int(ctl.suspect.sum())} of {len(ctl)}")

# Check 3 — how many actually reach the analysis, and does removing them matter?
pairs = pd.read_csv("results/matching/matched_pairs_secretory.csv")
suspects = set(ctl.loc[ctl.suspect, "Entry"].astype(str))
pairs["suspect"] = pairs.control_accession.astype(str).isin(suspects)
pairs.to_csv("results/matching/matched_pairs_secretory_audited.csv", index=False)
print(f"\nmatched control sites from suspect proteins: "
      f"{int(pairs.suspect.sum())} of {len(pairs)} ({100*pairs.suspect.mean():.1f}%)")


def load(manifest_path, score_path):
    manifest = pd.read_csv(manifest_path, low_memory=False)
    manifest = manifest[manifest.scoreable.astype(bool)].copy()
    scores = pd.read_csv(score_path, low_memory=False)
    for frame in (manifest, scores):
        frame["accession"] = frame.accession.astype(str)
        frame["position"] = frame.position.astype(int)
    return manifest.merge(scores[KEY + ["conditional_sequon_score"]], on=KEY, how="inner")


dataset = load("results/manifests/candidate_manifest_dataset.csv", "results/scores/scores_dataset.csv")
controls = load("results/manifests/manifest_matched_secretory.csv", "results/scores/scores_secretory.csv")
site = pd.concat([dataset, controls], ignore_index=True)
if "ortholog_clusters" not in site.columns:
    site["ortholog_clusters"] = pd.NA
sd = float(dataset.conditional_sequon_score.std(ddof=1))
margin = 0.2 * sd

print("\n=== sensitivity: does removing the suspects move the result? ===")
for label, subset in (("all pairs", pairs),
                      ("suspects removed", pairs[~pairs.suspect])):
    contrasts = build_contrasts(subset, site)
    draws = cluster_bootstrap(contrasts, 10000, 20260818)
    low, high = np.percentile(draws, [2.5, 97.5])
    mean = contrasts.contrast.mean()
    print(f"  {label:18s} n={len(contrasts):3d}  {mean/sd:+.3f} SD  "
          f"[{low/sd:+.3f}, {high/sd:+.3f}]  {classify(low, high, margin)}")
