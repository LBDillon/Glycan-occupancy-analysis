# Analysis options

What can be done with the site tables this module produces, in a recommended
order of work, and what will go wrong if the caveats are ignored.

Option 1 below has since been carried out for ProteinMPNN; see
[`primary_result.md`](primary_result_SUPERSEDED_2026-08-25.md). The rest remain unimplemented, and the
module's job still ends at the evidence tables.

Read [`evidence_sources.md`](../evidence_sources.md) first. Every analysis below
inherits one constraint from it: **no negative class can be derived from
annotation.** A site with no supporting layer is `unknown`, never a negative.

Two distinct things do exist, and conflating them is the main hazard here:

- **Definitive biochemical negatives — none.** Nothing establishes that a sequon
  was examined and found unmodified. Obtaining these is a data-acquisition
  problem — PNGase F / H₂¹⁸O occupancy glycoproteomics — not a coding one.
- **Informative structural internal controls — 32 sites, 25 proteins.** Sequons
  with no modelled glycan in structures that model glycans elsewhere, from hosts
  competent to glycosylate. Absence is informative here without being proven.

The second is what the occupancy comparison actually runs on, and its scarcity is
that analysis's binding constraint: 28 of the 32 are structurally scoreable, and
16 survive matching.

---

## 1. Occupancy-model benchmarking

**Do this first.** It is the use for which auditable provenance is the scarce
ingredient, and it needs the least new machinery.

The positive set is ready to use. `experimental_sites_all.csv` gives 922 sites
in 703 proteins, each with `support_sources` recording which independent layers
support it, and each traceable through `provenance.json` to a hashed input.
Three graded positive sets are available for a sensitivity analysis without any
extra work:

| Set | Sites | Use |
|---|---|---|
| `support_count >= 2` | 430 (109 of them supported by all three layers) | Highest-confidence core |
| `experimental_sites_uniprot_baseline.csv` | 505 | Curated-only, frozen and regression-tested |
| `experimental_sites_all.csv` | 922 | All layers, including 44 sites that fail UniProt policy but carry a structural glycan |

Those 32 sites are especially valuable as a held-out check: they carry a
modelled covalent glycan but fail UniProt's evidence policy, so a model trained
on curated annotation alone has not seen them.

### The negative set is the hard part

The `observed_unmodified` class exists, but it is 32 sites and it is a structural
internal control rather than a proven negative — enough for a matched paired
comparison, nowhere near enough to train against. The tempting move — treat the
3,385 excluded sites as negatives — is wrong, and wrong in a way that will
produce a good-looking and meaningless benchmark.

Those 3,385 sites are `unknown`. They are dominated by proteins nobody has
studied closely. A classifier trained to separate 922 positives from 3,385
"negatives" will learn the features that predict *being well studied* — model
organism, protein family, abundance, medical relevance, structural
tractability — and score highly, because those features really do separate the
two groups. It will be measuring curation effort.

### Use the unsupported sites as an unlabelled set

The defensible framing is **positive-unlabelled (PU) learning**. Sequon-matching
sites with no supporting layer form the unlabelled set `U`: a mixture of true
positives nobody has looked at and true negatives, in unknown proportions.
Standard PU approaches apply — estimating the label frequency, treating `U` as
weighted, or bagging over `U` subsamples.

Reporting rules that follow:

- Never call `U` "negative" in a figure, table, caption or metric name.
- Precision and F1 computed against `U`-as-negative are not interpretable.
  Report ranking metrics — how highly the model ranks known positives within
  `U` — and be explicit that "false positives" are unverified, not wrong.
- If you must approximate a negative set, do it explicitly and defend it: for
  example, sequons in proteins that are heavily annotated overall (many
  experimentally supported sites elsewhere in the same protein) but where this
  particular sequon is unannotated. That at least conditions on someone having
  looked at the protein. It is still an assumption, and it should be stated as
  one and varied in a sensitivity analysis.
- Match the unlabelled comparison set on annotation depth, not just on sequence.
  Sampling `U` uniformly reproduces the annotation-bias confound in section 2.

---

## 2. Sequon-loss context analysis

This is the original evolutionary question the ortholog database was built for:
**of the sequons lost between orthologs, are experimentally occupied ones lost
at a different rate than merely motif-matching ones?**

A sequon is a motif; occupancy is a modification. If occupancy is functionally
constrained, occupied sequons should be lost less readily than unoccupied ones.
The site tables make the comparison expressible for the first time, because
until now "lost a sequon" was the only available statement.

### Shape of the analysis

`site_pair_associations.csv` keys `(accession, position, pair_id)` and carries
`cluster_id`, `neg_accession`, `source` and `homology_qc_bucket`. Join it to
`experimental_sites_all.csv` on `(accession, position)` to label each pair-level
loss event by whether its positive-side site is experimentally supported, then
compare loss rates between supported and unsupported sites.

Restrict to a homology-quality subset before drawing conclusions — start with
`experimental_sites_strict.csv` (strict ortholog-like only), and use
strict+plausible as the sensitivity expansion. Poor-quality "ortholog" pairs
produce apparent losses that are alignment artefacts.

### The annotation-bias confound

This is the analysis most exposed to it, and it must be addressed head-on rather
than mentioned in a limitations paragraph.

Supported and unsupported sites do not differ only in occupancy. Supported sites
sit disproportionately in well-studied proteins, and well-studied proteins are
not a random sample of the proteome: they skew toward human and mouse, toward
secreted and cell-surface proteins, toward medically relevant families, toward
abundant proteins, and toward proteins with structures. Every one of those
properties has its own relationship to evolutionary rate. A difference in loss
rate between supported and unsupported sites is therefore confounded by default,
and the naive comparison will find *something* whether or not occupancy matters.

Ways to address it, in increasing order of effort and persuasiveness:

- **Condition on the protein.** Compare supported and unsupported sequons
  *within the same protein*. This removes every protein-level confound at
  once — species, abundance, family, study intensity — because both sites share
  them. It costs power (only proteins carrying both kinds of site contribute)
  and it is the single most convincing version of the analysis. Do this one
  first.
- **Match on annotation depth.** Use the number of experimentally supported
  sites elsewhere in the same protein as a study-intensity proxy, and either
  match on it or include it as a covariate. It is crude but it is measurable
  directly from these tables.
- **Match on species and function.** `proteins_master.csv` in the ortholog
  database carries lineage and `function_category` fields; the existing control
  generators in `scripts/08_structure.py` already match on `function_category`
  and `lineage_kingdom`, which is a working precedent.
- **Model it.** Regress loss on support status with study-intensity, species and
  family terms, and report how much the support coefficient moves as they are
  added. If it survives, say by how much it shrank.
- **Stratify by evidence layer.** Structural support and curated support have
  different biases: structural evidence requires a crystallisable protein,
  curated evidence requires a motivated curator. If the effect appears under one
  layer and not the other, that is informative about the confound rather than
  about biology. The 32 structure-only sites are the natural probe.

Report the direction and rough size of the residual bias you could not remove.
An honest bounded claim is worth more than a clean unconditional one.

---

## 3. Cluster-aware statistics

**Sites within an ortholog cluster are not independent observations.** They
share ancestry, sequence, structure and, usually, annotation history. Treating
them as independent inflates the effective sample size and makes p-values
arbitrarily small.

`cluster_id` is in `site_pair_associations.csv`, one per pair association.
Collapse it onto sites (a site inherits the cluster of its associations; check
for the multi-cluster case rather than assuming it away) and use it as the unit
of resampling everywhere:

- **Cluster-level bootstrap or permutation.** Resample whole clusters, not
  sites. This is the simplest correct thing and should be the default for every
  confidence interval and every test in sections 1, 2 and 4.
- **Mixed models** with a random intercept on `cluster_id` when a regression is
  wanted rather than a resampling test. Add a nested protein-level random effect
  if a protein contributes many sites, which is common.
- **Cluster-aware cross-validation** for anything in section 1. Split on
  `cluster_id`, never on sites. A random site-level split puts near-identical
  homologues on both sides and produces an optimistic score that will not
  reproduce.

Report both the number of sites and the number of clusters behind every headline
number. The ratio is the honest indication of how much independent information
is present, and readers cannot recover it from a site count alone.

---

## 4. Structure-based analysis

Two structural resources, with different coverage and different meaning.

**The existing structural context table**
(`paths.existing_structural_context`, `pdb_site_structural_context.csv`) carries
per-site relative solvent accessibility (`rsa`, `rsa_bin`), `sasa`,
`secondary_structure`, `mean_bfactor`, neighbour composition within 8 A,
`nearby_glycan_count_6a`, chain and UniProt terminal distances, and window
resolution fractions. It supports the obvious questions: are supported sites
more solvent-exposed than unsupported ones, are they enriched in particular
secondary structure, are they further from termini.

Its ceiling matters. It was generated for a case-control comparison, not for
this module: it spans 1,050 unique positive-side sites, its `status` column
contains only `ok` rows because unresolved sites were filtered out before
writing, and it covers **87 of the 505** baseline sites. Absence of a row
therefore encodes `not_assessed` — never `not_resolved`, and never "no
structure". Do not compute a denominator from it.

**The new structural mapper** (`structure_site_evidence.csv`) assesses every
candidate site against a cached structure and reports the outcome on an explicit
ladder rather than filtering failures out. That is what it adds: a denominator.
From the current build, over all 4,307 candidates:

| `structure_tier` | `structure_detail` | Sites |
|---|---|---|
| `structure_linked_glycan` | (blank) | 172 |
| `structure_residue_resolved` | (blank — genuinely resolved, no glycan) | 332 |
| `structure_residue_resolved` | `low_confidence_chain_match` | 4 |
| `structure_residue_unresolved` | `position_not_in_model` | 252 |
| `structure_not_assessed` | `no_cached_structure` | 3,543 |
| `structure_not_assessed` | `mmcif_linkage_unsupported` | 4 |

It also supplies chain, `resseq`, insertion code and the residue actually
observed, which is what a downstream RSA or burial calculation needs in order to
address the right atom. Combining the two — the mapper for coverage and residue
identity, the existing table for geometry where it overlaps — is the practical
route.

A caution specific to `nearby_glycan_count_6a`: a glycan near a site is not
evidence about *that* site's occupancy. Neighbouring occupied sequons are common
in glycan-dense regions and the correlation is real, but it is a proximity
statistic, not a per-site observation.

---

## 5. Controls

The canonical structural comparison in `scripts/08_structure.py` already
generates three control types, written to
`results/database_current/analyses/pdb_structure_comparison/` (and to the
`homology_qc/strict_ortholog_like/` variant). Each answers a different "compared
to what", and they are not interchangeable.

| Control | Construction | Controls for | Use when |
|---|---|---|---|
| `same_pair_random_position` | Same ortholog pair, a randomly chosen aligned non-sequon position, excluding a window around every sequon in the positive sequence | Everything about the pair — species, divergence, structure quality, alignment method | Asking whether something is special about *the sequon position* rather than about the pair. The tightest control, and the first one to reach for. |
| `same_cluster_preserved_sequon` | A third protein from the same cluster, same function category, glyco-positive, in which the observed positive sequon is **preserved** rather than lost | Cluster, family and function, while varying only whether the sequon survived | Asking whether *losing* the sequon is what matters, as opposed to having one. The direct comparator for section 2. |
| `unlinked_matched` | A protein from a *different* cluster, matched on function category and kingdom (relaxing to same-kingdom, then to any PDB protein, when the pool is under 20), with length ratio at least 0.5 | Broad compositional and technical effects | Establishing a background expectation from unrelated proteins. The loosest, and the right one for asking whether an effect is specific to related proteins at all. |

Notes on using them:

- Each control writes a matching `*_dropped.csv` alongside the kept rows.
  **Read the dropped file.** Controls that fail structural comparison are not
  missing at random — they fail for reasons correlated with structure quality —
  so a control set with heavy dropout is a biased comparator, and the dropout
  rate belongs in the methods.
- The controls are all generated at the **pair** level. Reusing them for a
  site-level analysis reintroduces exactly the duplication this module exists to
  remove: deduplicate onto `(accession, position)` first.
- These controls were built for a structural comparison, not for an occupancy
  comparison. `same_cluster_preserved_sequon` transfers most directly to the
  occupancy question; the other two need their assumptions rechecked against
  whatever is being asked.
- Report all three where they disagree. Disagreement between a tight and a loose
  control is a result about specificity, not a problem to be resolved by picking
  the favourable one.

---

## 6. Limitations

### Annotation bias

The dominant limitation, and the one that most easily becomes an artefact. It is
described in full in section 2. In short: the evidence layers record what has
been looked at, absence of annotation is not evidence of absence, and any
comparison between supported and unsupported sites is confounded by study
intensity until it is shown not to be.

### Taxonomic skew

GlyGen and GlyConnect are heavily skewed toward human and mouse. Their curation
and their upstream mass-spectrometry sources concentrate on those organisms, so
support from those layers is partly a statement about which species a protein
comes from. Any cross-species comparison should either stratify by species or
restrict to the UniProt and structure layers, and species composition should be
reported alongside any per-layer count.

### One cached structure per accession

The structure manifest holds **one** downloaded structure per accession (one
`output_path`, one `pdb_id`), while its own `all_pdb_ids` column — carried over
from `pdb_ids` in `proteins_master.csv` — lists up to **2,422** cross-referenced
entries for a single accession, and 1,474 accessions have more than one. A site
is therefore assessed against one arbitrary structure. If that
entry happens to be an unglycosylated construct, a bacterially expressed
fragment, or a structure that does not span the residue, the site lands on
`structure_residue_resolved`, `structure_residue_unresolved` or
`structure_not_assessed` even though another deposited entry shows the glycan
plainly. Structural support is therefore a **lower bound**, and the 172 linked
sites should be read as "at least 172".

### `observed_unmodified` is unavailable

There is no negative class and no source that can supply one. Every analysis
must be designed around a positive set and an unlabelled set. See section 1.

### Point mutants lose to perfect copies

`SAME_PROTEIN_IDENTITY_TOLERANCE = 0.02` in `structures.py` groups chains within
0.02 identity of the best-matching chain as copies of the same protein, then
prefers a linked chain among that group. When one PDB entry holds **both** a
perfectly matching chain (identity 1.0) and an engineered point mutant at around
0.975 identity, the mutant falls outside the tolerance and is excluded from
chain selection. If the mutant carried the only modelled glycan, that glycan is
missed and the site is reported `structure_residue_resolved`.

This was demonstrated reproducibly during review by construction. It was **not**
observed in the current corpus, so it does not affect the 172 figure as far as
can be determined — but the mechanism is real, and it can only ever remove
evidence.

The failure direction is deliberate: **a missed glycan, never a fabricated
one.** Loosening the tolerance to admit mutants would risk readmitting
cross-protein misattribution, which is the failure this gate exists to
prevent — MHC alpha and beta chains, which co-occur in one entry and are
different proteins, sit at 0.385 identity. Widening the window toward that
number eventually credits one protein's glycan to another. The current setting
accepts a small, one-directional loss to keep every asserted glycan trustworthy.
Treat structural support as a lower bound (see above) rather than tuning this
constant.

### `structure_residue_resolved` carries two meanings

The tier is overloaded, and the two meanings are separated **only** by
`structure_detail`:

| `structure_tier` | `structure_detail` | Actual meaning |
|---|---|---|
| `structure_residue_resolved` | (blank) | The residue was genuinely resolved in a credibly matched chain, with no glycan linkage attached |
| `structure_residue_resolved` | `low_confidence_chain_match` | No chain aligned well enough to identify this protein. This is a best positional guess and asserts nothing about occupancy |

The second case exists so the module never asserts a glycan from a chain it
cannot credibly call this protein — a correct conservative choice — but it lands
on the same tier as the first.

**Any analysis grouping by `structure_tier` alone will conflate them.** Always
group by tier **and** detail. In the current build only 4 sites are affected, so
the error is small today, but it scales with corpus growth and it is silent: a
conflated count looks exactly like a correct one.

---

## 7. Phase-2 work

Ordered by how much each would improve the dataset.

1. **Per-site PDB selection across all cross-referenced entries.** The largest
   single win. Instead of one arbitrary cached structure per accession, consider
   every entry in `all_pdb_ids`, and for each site prefer a structure that spans
   the residue, is eukaryotically expressed, and carries glycan density. This
   directly attacks the "one cached structure" limitation and would raise the
   structural layer from a lower bound toward a real count. It costs storage and
   download time, and needs a defined, recorded selection rule so results stay
   reproducible.

2. **mmCIF `_struct_conn` parsing.** `parse_link_records` handles the
   fixed-column PDB format only, so mmCIF-only entries return
   `structure_not_assessed` with `mmcif_linkage_unsupported` — 4 sites in the
   current build, but a growing share of the PDB overall, and large structures
   are mmCIF-only by necessity. `_struct_conn` carries the same covalent-bond
   information in a parseable form. This becomes urgent if step 1 is done, since
   selecting across all entries will pull in many mmCIF-only ones.

3. **SIFTS-based residue mapping instead of alignment.** SIFTS provides
   authoritative, curated UniProt-to-PDB residue correspondences. Replacing
   pairwise local alignment with SIFTS would eliminate the whole class of
   chain-selection problems above — the identity tolerance, the coverage gates,
   the `low_confidence_chain_match` state and the overloaded
   `structure_residue_resolved` tier all exist only because the mapping is
   inferred rather than looked up. It also removes the alignment cost, which
   dominates runtime. The trade-off is a new external dependency that must be
   versioned and cached like any other input, and a fallback path for entries
   SIFTS does not cover.

Also worth doing, smaller:

- Evaluate the GlyConnect layer on a populated cache. It is fetched and parsed
  but does not confer positivity (`glyconnect_qualifies = false`), on the
  argument that GlyGen already ingests it. Measuring how many of its sites GlyGen
  missed would let that stay a decision rather than an assumption.

Two earlier items on this list are now done: the GlyGen cache has been populated
(1,714 cross-referenced accessions) and rerun, and the resulting enriched counts
are frozen in `tests/snapshots/enriched_2026-08-06.json`, so the full pipeline
now has the same regression protection as the UniProt baseline.
