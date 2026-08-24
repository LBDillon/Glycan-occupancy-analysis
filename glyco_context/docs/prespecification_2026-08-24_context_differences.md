# Pre-specification — occupancy-associated context differences

*Written and committed 2026-08-24, **before** any group difference was computed.
Its purpose is to fix the analysis choices while they can still be made
honestly. If a later analysis departs from this document, the departure is
reported as a departure rather than presented as the plan.*

## Datasets

**Primary:** `context_triplet_core.csv` — every feature in the row describes the
sequon it names (triplet agrees, all three residues located, mapping
continuous). 2,556 sites.

**Sensitivity only:** `context_asn_core.csv` (2,624), used exclusively for
features centred on the asparagine, and only where named in advance. It is never
used for +1 or +2 exposure, secondary structure or geometry.

**Not used for inference:** `context_construct_review.csv` (104). Inspected for
patterns, reported, never tested.

## Reference and comparisons

The reference distribution is **occupied-supported sites** (318 sites, 233
proteins).

Two comparisons, reported separately and never pooled:

1. **Internal controls** (31 sites, 24 proteins) — the more informative label:
   the structure models glycans elsewhere but not here, so sugars demonstrably
   survived preparation and this depositor demonstrably modelled them. Almost no
   power. Pre-specified as **direction-agreement evidence, not an independent
   test.**
2. **Secretory-unannotated** (2,207 sites, 1,098 proteins) — power, but a weak
   and contaminated label. These are never called negatives.

**Bacterial (3,280) and cytosolic (3,543) are diagnostics, not tests.** They are
expected to differ from occupied sites for compartment and composition reasons
that have nothing to do with occupancy. They enter the report to characterise
that confound, and are excluded from the confirmatory family.

Interpretation rule, fixed in advance:

| Internal | Secretory | Reading |
|---|---|---|
| same direction | same direction | more credible occupancy-associated feature |
| null | effect | power, contamination, or population composition |
| effect | null | possibly real, imprecise |
| opposite | opposite | label or population dependence; investigate before claiming |

## The feature family

Tested (continuous, standardised mean difference):

`n_rsa`, `plus1_rsa`, `plus2_rsa`, `loop_run_length`, `n_neighbours_8a`,
`nd2_atoms_8a_same_chain`, `nd2_residues_8a_same_chain`,
`nd2_atoms_8a_other_chain`, `sidechain_neighbour_residues_5a`,
`neighbour_net_charge_8a`, `neighbour_hydrophobic_fraction_8a`,
`neighbour_aromatic_fraction_8a`, `nearest_aromatic_sidechain_nd2`,
`uniprot_residues_after_asn`, `uniprot_residues_after_sequon`,
`distance_to_n_terminus_resolved`, `distance_to_c_terminus_resolved`

Tested (categorical, difference in proportion per level):

`n_ss_coarse`, `plus1_ss_coarse`, `plus2_ss_coarse`, `aromatic_within_8a`

**Excluded, with reasons fixed in advance:**

- `nearest_disulfide_sg_nd2` and `nearest_disulfide_ca` — coverage is 52.5% in
  occupied sites, 15.5% secretory, 9.0% bacterial, 1.6% cytosolic. A gap that
  size means a cross-arm comparison compares missingness, and missingness here
  tracks compartment: disulfides form in the oxidising secretory environment.
  Reported descriptively within the occupied set only, and named as a
  compositional marker rather than a candidate occupancy feature.
- `n_phi`, `n_psi` and the +1/+2 dihedrals — angles are circular, so a
  standardised mean difference is not defined on them. Backbone geometry enters
  as a categorical Ramachandran region instead.
- `neighbour_*_count_8a` where the corresponding fraction is tested — redundant
  with the fraction plus `n_neighbours_8a`.
- `chain_length_resolved` — a property of the deposited construct, not of the
  site.
- `n_residue` (constant), `plus2_residue` (equivalent to the NXS/NXT stratifier).
- All QC and provenance columns. They describe measurement quality and are
  reported beside the estimates, never tested as biology.

## Estimation and uncertainty

Effect size is the standardised mean difference

    SMD = (mean_occ - mean_cmp) / sqrt((var_occ + var_cmp) / 2)

reported as a description of separation, **not** a probability of occupancy.

Sites within a protein are not independent — one protein contributes up to 7
occupied sites and up to 19 secretory ones. All intervals therefore come from a
**cluster bootstrap resampling proteins**, not rows: 2,000 replicates,
percentile 95% interval, with the bootstrap p-value
`2 * min(P(SMD* <= 0), P(SMD* >= 0))`.

## Multiplicity

Benjamini–Hochberg within each comparison, across the whole pre-specified family
in that comparison. Matches the correction used for the occupancy benchmark.

## Missing data

Complete-case per feature, with **coverage reported beside every estimate**.
Missingness is not assumed random: flexible and glycosylated regions are
preferentially unresolved, which is the direction that matters here.
`loop_run_length` is additionally reported split by `loop_run_censored`, because
a run reaching an unresolved boundary is a lower bound.

## Pre-specified stratification and sensitivities

- **NXS versus NXT**, reported separately for every primary feature.
- **Evidence tier** — by `support_count` and by `glycan_modelled_at_site`.
- **Structure choice** — `no_alternative` against `selected_from_alternatives`,
  since selection favoured glycan-bearing structures and is not label-blind.
- **Missingness** — per-feature coverage by population beside each estimate.

## What this analysis cannot conclude

It cannot establish that any feature causes glycosylation, that a
secretory-unannotated site is unoccupied, or that a model represents a feature
because its score correlates with one. Associations here motivate the masking
and perturbation work; they do not substitute for it.
