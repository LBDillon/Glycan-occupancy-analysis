"""One row per biological site, for the natural-context atlas.

The occupancy benchmark's manifests are keyed by (site, structure) and exist to
be *scored*. This one is keyed by (accession, position) and exists to be
*described*: it carries the evidence behind the site, the sequence around it,
which structure was used and why, and which reference tier it belongs to.

Two things it deliberately records rather than resolves.

**Why a structure was chosen.** The evidence mapper picks the structure carrying
the strongest occupancy evidence, which is usually a glycan-bearing one. That is
right for gathering evidence and not label-blind for comparing contexts. But only
14% of occupied sites (45 of 332) had more than one structure examined at all —
for the other 287 there was nothing to choose between, so no selection bias is
possible. `n_structures_examined` and `structure_choice` mark the 45 that need a
sensitivity check, instead of re-running selection for everyone.

**Tier membership, not a single reference set.** Evidence strength trades off
against sample size, and which trade is right depends on the question. Every site
carries its tier flags so an analysis can pick, and so "consistent across tiers"
can be checked rather than assumed. The tiers nest, so their overlap is reported.

Usage:  41_context_manifest.py [--window 11]
"""
import argparse, csv, gzip, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.table_io import write_table

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--window", type=int, default=11,
                    help="residues of UniProt sequence either side of the Asn")
parser.add_argument("--out", default="results/datasets/context_manifest.csv")
args = parser.parse_args()

KEY = ["accession", "position"]
D = Path("results/datasets")

AA_CLASS = {**{a: "acidic" for a in "DE"}, **{a: "basic" for a in "KRH"},
            **{a: "polar" for a in "STNQ"}, **{a: "hydrophobic" for a in "AVLIM"},
            **{a: "aromatic" for a in "FWY"},
            "G": "glycine", "P": "proline", "C": "cysteine"}


def load(name):
    path = D / name
    if not path.exists():
        raise SystemExit(f"missing input: {path}")
    frame = pd.read_csv(path, low_memory=False)
    if "accession" in frame:
        frame["accession"] = frame.accession.astype(str)
    if "position" in frame:
        frame["position"] = pd.to_numeric(frame.position, errors="coerce")
        frame = frame.dropna(subset=["position"])
        frame["position"] = frame.position.astype(int)
    return frame


sites = load("experimental_sites_all.csv")
features = load("site_structural_features.csv")
evidence = load("structure_site_evidence.csv")
assoc = load("site_pair_associations.csv")

# --- UniProt sequences, for the window and the +1 residue -------------------
sequences = {}
with gzip.open("../../data/raw/uniprot/uniprot_reviewed_glycoproteins_2026-04-27.tsv.gz",
               "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        sequences[row["Entry"]] = row.get("Sequence", "")

rows = []
for r in sites.itertuples(index=False):
    seq = sequences.get(r.accession, "")
    pos = int(r.position)                      # 1-indexed
    i = pos - 1
    triplet = seq[i:i + 3] if 0 <= i and i + 3 <= len(seq) else ""
    plus1 = triplet[1] if len(triplet) == 3 else ""
    plus2 = triplet[2] if len(triplet) == 3 else ""
    lo, hi = max(0, i - args.window), min(len(seq), i + 3 + args.window)
    rows.append({
        "accession": r.accession,
        "position": pos,
        # --- evidence -----------------------------------------------------
        "occupancy_status": getattr(r, "occupancy_status", ""),
        "support_sources": getattr(r, "support_sources", ""),
        "support_count": int(getattr(r, "support_count", 0) or 0),
        "uniprot_tier": getattr(r, "uniprot_tier", ""),
        "glygen_tier": getattr(r, "glygen_tier", ""),
        "structure_tier_evidence": getattr(r, "structure_tier", ""),
        # --- sequence context ---------------------------------------------
        "sequon_triplet": triplet,
        "subtype": f"NX{plus2}" if plus2 else "",
        "plus1_residue": plus1,
        "plus1_class": AA_CLASS.get(plus1, "other" if plus1 else ""),
        "sequence_window": seq[lo:hi],
        "window_offset": i - lo,
        "window_full": (i - lo == args.window) and (hi - (i + 3) == args.window),
        "uniprot_length": len(seq),
        "distance_to_n_terminus": pos,
        "distance_to_c_terminus": len(seq) - (pos + 2) if seq else None,
    })

manifest = pd.DataFrame(rows)

# --- ortholog cluster: the resampling unit, never an independence claim ----
clusters = (assoc.groupby(KEY).cluster_id
            .agg(lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x)})[:3]))
            .rename("ortholog_clusters").reset_index())
manifest = manifest.merge(clusters, on=KEY, how="left")
manifest["n_ortholog_associations"] = manifest[KEY].merge(
    assoc.groupby(KEY).size().rename("n").reset_index(), on=KEY, how="left")["n"].fillna(0).astype(int)

# --- structure, and why it was the one used --------------------------------
ev = evidence.drop_duplicates(KEY)
manifest = manifest.merge(
    ev[KEY + ["structure_pdb_id", "structure_chain_id", "structure_resseq",
              "structure_tier", "structure_n_examined", "structure_glycans_elsewhere",
              "structure_expression_system"]], on=KEY, how="left")
manifest["n_structures_examined"] = pd.to_numeric(
    manifest.structure_n_examined, errors="coerce").fillna(0).astype(int)
manifest["glycan_modelled_at_site"] = manifest.structure_tier.eq("structure_linked_glycan")
manifest["structure_choice"] = manifest.n_structures_examined.map(
    lambda n: "no_alternative" if n <= 1 else "selected_from_alternatives")

# --- structural features, and whether they are usable ----------------------
feat_cols = [c for c in ("features_available", "rsa", "rsa_bin", "sasa",
                         "n_neighbours_8a", "hydrophobic_fraction_8a",
                         "charged_fraction_8a", "chain_length_resolved",
                         "distance_to_chain_terminus", "observed_residue")
             if c in features.columns]
manifest = manifest.merge(features.drop_duplicates(KEY)[KEY + feat_cols], on=KEY, how="left")
manifest["features_available"] = manifest.get("features_available", False).fillna(False).astype(bool)

# --- reference tiers, nested ------------------------------------------------
occupied = manifest.occupancy_status.eq("occupied_supported")
manifest["tier_all_occupied"] = occupied
manifest["tier_features"] = occupied & manifest.features_available
manifest["tier_two_layer"] = manifest.tier_features & manifest.support_count.ge(2)
manifest["tier_linked_glycan"] = manifest.tier_features & manifest.glycan_modelled_at_site
manifest["tier_three_layer"] = manifest.tier_features & manifest.support_count.ge(3)

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
write_table(manifest, out)

tiers = ["tier_all_occupied", "tier_features", "tier_two_layer",
         "tier_linked_glycan", "tier_three_layer"]
print(f"wrote {out}  ({len(manifest)} sites, {manifest.accession.nunique()} proteins)\n")
print(f"{'tier':22}{'sites':>7}{'proteins':>10}{'clusters':>10}")
for t in tiers:
    sub = manifest[manifest[t]]
    print(f"{t:22}{len(sub):>7}{sub.accession.nunique():>10}"
          f"{sub.ortholog_clusters.fillna('').nunique():>10}")

print("\ntier overlap (rows shared with the column tier):")
core = [t for t in tiers if t != "tier_all_occupied"]
print(f"{'':22}" + "".join(f"{t.replace('tier_',''):>14}" for t in core))
for a in core:
    print(f"{a:22}" + "".join(
        f"{int((manifest[a] & manifest[b]).sum()):>14}" for b in core))

sel = manifest[manifest.tier_features]
print(f"\nstructure choice among {len(sel)} occupied sites with features:")
print("  " + str(sel.structure_choice.value_counts().to_dict()))
print(f"  glycan modelled at the site: {int(sel.glycan_modelled_at_site.sum())}")
print(f"\nsequence window: {int(manifest.window_full.sum())} of {len(manifest)} full "
      f"(+/-{args.window}); subtype {manifest[occupied].subtype.value_counts().to_dict()}")
