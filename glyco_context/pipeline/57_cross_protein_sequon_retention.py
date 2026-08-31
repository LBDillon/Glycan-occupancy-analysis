"""Occupied sequons against motif-only sequons in the SAME chains.

The matched-secretory comparison pairs an occupied site in one protein with a
control site in another, so protein identity is never fully controlled. This
does not: an occupied sequon and a motif-only sequon on one chain share the
chain, the fold, the depositor and — because designs are generated per chain —
the same 32 designs. The contrast is therefore within-protein by construction.

Four populations, all read off one design run so they are comparable:

    occupied sequon    N-X-S/T whose Asn is an `occupied_supported` manifest site
    motif-only sequon  N-X-S/T in the same chain with no occupancy support
    non-sequon Asn     any other asparagine in the same chain
    control triplet    any 3-residue window touching no sequon

Three readings, per the vocabulary in the brief:

    n_retained        design[k] == "N"
    pattern_retained  N, not-P, S/T  -- `classify_retention`'s full_sequon_retained
    exact_retained    design[k:k+3] == wild[k:k+3]

Pattern is the biologically meaningful reading at a sequon, but an arbitrary
three-residue window has no equivalent latitude, so **every sequon-against-
control-triplet comparison here is exact against exact**. Pattern-against-exact
overstates the case and appears nowhere.

Reads the sequences persisted by `08_design.py --save-sequences`, so it is
offline analysis: no model is run, and adding a model means pointing `--sequences`
at that model's file.

Usage:
    57_cross_protein_sequon_retention.py --sequences <designs>_sequences.csv
"""
import argparse, json, re, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from experimental_glycosylation_sites.retention import classify_retention
from experimental_glycosylation_sites.runner_support import structure_paths
from experimental_glycosylation_sites.structures import _parse_chains
from experimental_glycosylation_sites.provenance import _git_state
from glyco_context.fixed_design import verify_sequon_index

# Lookahead so overlapping sequons are both found: NNSS contains N-N-S at 0 and
# N-S-S is not a sequon, but NNTT does contain two starts in general.
SEQUON = re.compile(r"(?=(N[^P][ST]))")
N_BOOT, BOOT_SEED = 4000, 11          # matching 10_analyse_retention_by_class.py

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--sequences", required=True,
                    help="'<out>_sequences.csv' from 08_design.py --save-sequences")
parser.add_argument("--manifest",
                    default="results/manifests/candidate_manifest_dataset.csv")
parser.add_argument("--out-dir", default="glyco_context/results/analysis")
parser.add_argument("--label", default=None,
                    help="output name stem (default: derived from --sequences)")
args = parser.parse_args()

LABEL = args.label or Path(args.sequences).stem.replace("_sequences", "")
OUT = Path(args.out_dir); OUT.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(args.manifest, low_memory=False)
if "scoreable" in manifest.columns:
    manifest = manifest[manifest.scoreable.astype(bool)]
sequences = pd.read_csv(args.sequences, low_memory=False)
sequences["structure_pdb_id"] = sequences.structure_pdb_id.astype(str)
sequences["structure_chain_id"] = sequences.structure_chain_id.astype(str)
MODEL = sorted(str(m) for m in sequences.model.dropna().unique())
SEED = sorted(int(s) for s in sequences.seed.dropna().unique())
print(f"{len(sequences)} designs over "
      f"{sequences.groupby(['structure_pdb_id','structure_chain_id']).ngroups} chains; "
      f"model={MODEL} seed={SEED}")

paths = structure_paths(())
site_rows, chain_rows, dropped = [], [], []

for (pdb_id, chain_id), designs_group in sequences.groupby(
        ["structure_pdb_id", "structure_chain_id"]):
    path = paths.get(pdb_id.upper())
    if path is None:
        dropped.append({"pdb": pdb_id, "chain": chain_id,
                        "reason": "structure_not_cached"})
        continue
    native = next((c for c in _parse_chains(path, pdb_id)
                   if c.chain_id == chain_id), None)
    if native is None:
        dropped.append({"pdb": pdb_id, "chain": chain_id,
                        "reason": "chain_absent_from_parse"})
        continue
    wild = native.sequence

    # Designs must be the same length as the native chain, or every index below
    # points somewhere else. A mismatch is a parse disagreement, not noise.
    designs, wrong_length = [], 0
    for seq in designs_group.sequence:
        if isinstance(seq, str) and len(seq) == len(wild):
            designs.append(seq)
        else:
            wrong_length += 1
    if wrong_length:
        dropped.append({"pdb": pdb_id, "chain": chain_id,
                        "reason": f"design_length_mismatch_x{wrong_length}"})
    if not designs:
        continue

    chain_sites = manifest[(manifest.structure_pdb_id.astype(str) == pdb_id) &
                           (manifest.structure_chain_id.astype(str) == chain_id)]
    accession = (str(chain_sites.accession.iloc[0]) if len(chain_sites)
                 else f"{pdb_id}_{chain_id}")

    # Occupied sites come from the manifest and must survive the index check --
    # this is the guard the 2026-08-25 correction added, and a site that fails
    # it is refused rather than scored at whatever residue the index landed on.
    occupied, refused = {}, set()
    for site in chain_sites.itertuples():
        if site.occupancy_status != "occupied_supported":
            continue
        k = int(site.n_model_index)
        if verify_sequon_index(wild, k, str(site.triplet)):
            occupied[k] = (str(site.accession), int(site.position))
        else:
            # Refused, not reclassified. A site whose index cannot be trusted is
            # not thereby a motif-only site -- it is a site we cannot place, and
            # letting it fall through would move a known-occupied sequon into the
            # comparison arm it is supposed to be contrasted against.
            refused.add(k)
            dropped.append({"pdb": pdb_id, "chain": chain_id,
                            "accession": site.accession, "position": site.position,
                            "n_model_index": k, "triplet": site.triplet,
                            "reason": "verify_sequon_index_failed"})

    all_sequons = {m.start() for m in SEQUON.finditer(wild)}
    unsupported = {int(s.n_model_index) for s in chain_sites.itertuples()
                   if s.occupancy_status != "occupied_supported"}
    motif_only = sorted(all_sequons - set(occupied) - refused)

    covered = {i for k in all_sequons for i in (k, k + 1, k + 2)}
    nonsequon_n = [i for i, a in enumerate(wild) if a == "N" and i not in covered]
    windows = [i for i in range(len(wild) - 2)
               if not covered & {i, i + 1, i + 2}]

    def readings(k):
        """The three readings at one sequon, averaged over this chain's designs."""
        cls = classify_retention(designs, k, k + 1, k + 2)   # pattern, reused
        n_ret = np.mean([d[k] == "N" for d in designs])
        exact = np.mean([d[k:k + 3] == wild[k:k + 3] for d in designs])
        return {"frac_n_retained": float(n_ret),
                "frac_pattern_retained": cls["frac_full_sequon_retained"],
                "frac_exact_retained": float(exact),
                "n_designs_scored": cls["n_designs_scored"]}

    for k, (acc, pos) in sorted(occupied.items()):
        site_rows.append({"accession": acc, "position": pos, "site_class": "occupied",
                          "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                          "n_index": k, "wild_triplet": wild[k:k + 3], **readings(k)})
    for k in motif_only:
        site_rows.append({"accession": accession, "position": pd.NA,
                          "site_class": "motif_only",
                          "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                          "n_index": k, "wild_triplet": wild[k:k + 3],
                          # An `observed_unmodified` manifest site has evidence of
                          # ABSENCE, not absence of evidence. The brief defines
                          # motif-only as "not occupied_supported", so it is
                          # included -- flagged so the stricter reading is still
                          # available without a re-run.
                          "manifest_unsupported": k in unsupported, **readings(k)})

    ident = [[a == b for a, b in zip(wild, d)] for d in designs]
    chain_rows.append({
        "accession": accession, "structure_pdb_id": pdb_id,
        "structure_chain_id": chain_id, "chain_length": len(wild),
        "n_designs": len(designs),
        "n_occupied": len(occupied), "n_motif_only": len(motif_only),
        "frac_nonsequon_n_retained": float(np.mean(
            [[d[i] == "N" for i in nonsequon_n] for d in designs]))
            if nonsequon_n else np.nan,
        "n_nonsequon_n": len(nonsequon_n),
        "frac_control_triplet_exact": float(np.mean(
            [[row[i] and row[i + 1] and row[i + 2] for i in windows] for row in ident]))
            if windows else np.nan,
        "n_control_windows": len(windows),
        "background_mutation_rate": float(1.0 - np.mean(ident)),
    })

sites = pd.DataFrame(site_rows)
chains = pd.DataFrame(chain_rows)
sites.to_csv(OUT / f"{LABEL}_sites.csv", index=False)
chains.to_csv(OUT / f"{LABEL}_chains.csv", index=False)
if dropped:
    pd.DataFrame(dropped).to_csv(OUT / f"{LABEL}_dropped.csv", index=False)


def protein_bootstrap(values_by_protein, n_boot=N_BOOT, seed=BOOT_SEED):
    """Percentile interval, resampling PROTEINS.

    Several sequons on one chain share a single set of designs, so resampling
    sites would treat replicates of one draw as independent observations and
    give an interval far too narrow.
    """
    groups = [np.asarray(v, dtype=float) for v in values_by_protein if len(v)]
    if len(groups) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        draws[i] = np.concatenate([groups[j] for j in pick]).mean()
    return tuple(float(x) for x in np.percentile(draws, [2.5, 97.5]))


def summarise(frame, column, unit="accession"):
    vals = [g[column].dropna().to_numpy() for _, g in frame.groupby(unit)]
    pooled = np.concatenate([v for v in vals if len(v)]) if vals else np.array([])
    low, high = protein_bootstrap(vals)
    return {"mean": float(pooled.mean()) if len(pooled) else float("nan"),
            "ci95": [round(low, 4), round(high, 4)],
            "n_sites": int(len(pooled)),
            "n_proteins": int(frame[unit].nunique())}


occ = sites[sites.site_class == "occupied"]
mot = sites[sites.site_class == "motif_only"]

# ---- 1. N retention by class -------------------------------------------------
by_class = {
    "occupied sequon (Asn)": summarise(occ, "frac_n_retained"),
    "motif-only sequon (Asn)": summarise(mot, "frac_n_retained"),
    "non-sequon asparagine": summarise(chains.dropna(subset=["frac_nonsequon_n_retained"]),
                                       "frac_nonsequon_n_retained"),
}

# ---- 2. within-protein paired contrast (the headline) ------------------------
# Both site classes on one chain share that chain's designs, so the contrast is
# taken per chain and then resampled over proteins.
per_chain = (sites.groupby(["accession", "structure_pdb_id", "structure_chain_id",
                            "site_class"])
             [["frac_n_retained", "frac_pattern_retained", "frac_exact_retained"]]
             .mean().reset_index())
wide = per_chain.pivot_table(
    index=["accession", "structure_pdb_id", "structure_chain_id"],
    columns="site_class",
    values=["frac_n_retained", "frac_pattern_retained", "frac_exact_retained"])
wide.columns = [f"{a}__{b}" for a, b in wide.columns]
paired = wide.dropna().reset_index()
paired.to_csv(OUT / f"{LABEL}_paired.csv", index=False)

contrasts = {}
for reading in ("n_retained", "pattern_retained", "exact_retained"):
    o, m = f"frac_{reading}__occupied", f"frac_{reading}__motif_only"
    diff = paired[o] - paired[m]
    tmp = paired.assign(_d=diff)
    low, high = protein_bootstrap([g["_d"].to_numpy() for _, g in tmp.groupby("accession")])
    contrasts[reading] = {
        "occupied_mean": float(paired[o].mean()),
        "motif_only_mean": float(paired[m].mean()),
        "paired_difference": float(diff.mean()),
        "ci95": [round(low, 4), round(high, 4)],
        "excludes_zero": bool(low > 0 or high < 0),
        "n_chains": int(len(paired)),
        "n_proteins": int(paired.accession.nunique()),
        # Chain-weighted, unlike the other sections, and necessarily so: the
        # chain is the unit of PAIRING -- both classes on it share one set of
        # designs -- so a chain contributes one difference however many sequons
        # it carries. Pooling sites here would break the pairing.
        "weighting": "chain-weighted; the chain is the pairing unit",
    }

# ---- 3. exact triplet against control triplet (exact on BOTH sides) ----------
# Site-weighted point estimate, protein-resampled interval -- the convention in
# 10_analyse_retention_by_class.py and 56_, which pool sites from resampled
# proteins rather than averaging each protein to one value first. Averaging to
# the protein first is a different estimand (it weights a chain carrying one
# sequon the same as a chain carrying six) and shifts this contrast enough to
# change whether its interval covers zero.
tri = occ.merge(chains[["structure_pdb_id", "structure_chain_id",
                        "frac_control_triplet_exact"]],
                on=["structure_pdb_id", "structure_chain_id"], how="inner")
tri = tri.dropna(subset=["frac_exact_retained", "frac_control_triplet_exact"])
tri = tri.rename(columns={"frac_exact_retained": "sequon_exact",
                          "frac_control_triplet_exact": "control_exact"})
tri["_d"] = tri.control_exact - tri.sequon_exact
low, high = protein_bootstrap([g["_d"].to_numpy() for _, g in tri.groupby("accession")])
triplet = {
    "sequon_exact": float(tri.sequon_exact.mean()),
    "control_triplet_exact": float(tri.control_exact.mean()),
    "control_minus_sequon": float(tri["_d"].mean()),
    "ci95": [round(low, 4), round(high, 4)],
    "excludes_zero": bool(low > 0 or high < 0),
    "n_sites": int(len(tri)),
    "n_proteins": int(tri.accession.nunique()),
    "note": "exact against exact; pattern-against-exact is not a fair comparison",
    "weighting": "site-weighted mean, protein-resampled interval, as 10_ and 56_",
    "control_definition": "windows touching NO sequon of any kind. 56_ excluded "
                          "only manifest (occupied) sequon positions, so its "
                          "control pool includes motif-only sequons and this one "
                          "does not -- expect a small difference from its 9.5%.",
}

summary = {
    "model": MODEL, "seed": SEED,
    "designs_per_chain": int(chains.n_designs.median()) if len(chains) else 0,
    "temperature": float(sequences.temperature.dropna().iloc[0])
                   if sequences.temperature.notna().any() else None,
    "bootstrap": {"over": "proteins (accession), not sites",
                  "n_boot": N_BOOT, "seed": BOOT_SEED},
    "attrition": {
        "chains_in_sequence_file": int(sequences.groupby(
            ["structure_pdb_id", "structure_chain_id"]).ngroups),
        "chains_analysed": int(len(chains)),
        "chains_dropped": int(len(dropped)),
        "occupied_sites": int(len(occ)),
        "motif_only_sites": int(len(mot)),
        "motif_only_of_which_manifest_unsupported":
            int(mot.manifest_unsupported.sum()) if "manifest_unsupported" in mot else 0,
        "verify_sequon_index_refusals":
            int(sum(d.get("reason") == "verify_sequon_index_failed" for d in dropped)),
        "chains_with_both_classes": int(len(paired)),
    },
    "n_retention_by_class": by_class,
    "within_protein_paired_contrast": contrasts,
    "exact_triplet_vs_control": triplet,
    "pattern_retention_occupied": summarise(occ, "frac_pattern_retained"),
    "pattern_retention_motif_only": summarise(mot, "frac_pattern_retained"),
    # 56_sequon_retention_rate.py did NOT filter on occupancy_status: it verified
    # every scoreable manifest site, so its 13.6% is over occupied sites AND the
    # 28 internal controls. Reproducing that set is the like-for-like check;
    # comparing it against the strictly-occupied figure above would be comparing
    # two different populations and calling the difference a discrepancy.
    "pattern_retention_all_manifest_sites_56_equivalent": summarise(
        pd.concat([occ, mot[mot.manifest_unsupported == True]]) if len(mot) else occ,
        "frac_pattern_retained"),
    "background_mutation_rate": summarise(chains, "background_mutation_rate"),
    "git": _git_state(),
}
(OUT / f"{LABEL}_summary.json").write_text(json.dumps(summary, indent=2))

print(f"\nchains analysed {len(chains)}; occupied {len(occ)}; motif-only {len(mot)}; "
      f"both classes {len(paired)}; dropped {len(dropped)}")
print(f"\n{'class':32s} {'sites':>6s} {'proteins':>9s} {'N retained':>11s} {'95% CI':>20s}")
print("-" * 84)
for label, s in by_class.items():
    print(f"{label:32s} {s['n_sites']:>6d} {s['n_proteins']:>9d} {s['mean']:>11.4f} "
          f"[{s['ci95'][0]:>+7.4f},{s['ci95'][1]:>+7.4f}]")
print(f"\nwithin-protein paired contrast (occupied - motif-only), "
      f"{contrasts['n_retained']['n_chains']} chains:")
for reading, c in contrasts.items():
    star = "  *" if c["excludes_zero"] else ""
    print(f"  {reading:18s} {c['occupied_mean']:.4f} vs {c['motif_only_mean']:.4f}  "
          f"diff {c['paired_difference']:+.4f} "
          f"[{c['ci95'][0]:+.4f},{c['ci95'][1]:+.4f}]{star}")
print(f"\nexact vs control triplet: sequon {triplet['sequon_exact']:.4f}, "
      f"control {triplet['control_triplet_exact']:.4f}, "
      f"control-sequon {triplet['control_minus_sequon']:+.4f} "
      f"[{triplet['ci95'][0]:+.4f},{triplet['ci95'][1]:+.4f}]")
print(f"pattern retention, occupied only:        "
      f"{summary['pattern_retention_occupied']['mean']:.4f} "
      f"(n={summary['pattern_retention_occupied']['n_sites']})")
_eq = summary["pattern_retention_all_manifest_sites_56_equivalent"]
print(f"pattern retention, 56_-equivalent set:   {_eq['mean']:.4f} "
      f"(n={_eq['n_sites']})  <- compare against 0.136 (56_) / 0.130 (ARC)")
print(f"\nwrote {OUT}/{LABEL}_{{sites,chains,paired,summary}}")
