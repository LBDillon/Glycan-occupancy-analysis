# Evidence sources

What each source can establish, what it cannot, and how this module encodes the
difference.

The module treats evidence as **layers**, not as a single classification. Each
layer is scored independently, and `support_sources` records which ones fired.
No layer can silently overwrite another's call.

| Source | Can establish | Cannot establish | Role here |
|---|---|---|---|
| UniProt CARBOHYD | That a curator asserted glycosylation at an exact residue, with an auditable evidence code | That an unannotated site is unglycosylated | Primary layer |
| GlyGen | That an independent aggregation reports the site, often with an identified glycan | Anything, when the category is `predicted` or the only citation is UniProtKB | Independent layer |
| GlyConnect | That a glycan structure has been attributed to that asparagine | Broad coverage; it is small and largely already inside GlyGen | Corroboration only, off by default |
| PDB `LINK` records | That a glycan was covalently modelled onto that asparagine in a deposited structure | That a site without a linkage was unoccupied | Independent layer |
| GlyCosmos | Nothing per-site here | Per-site occupancy: no usable REST endpoint | Documented cross-reference only |
| OrthoDB / homology QC | Which ortholog subset a site belongs to | Any occupancy claim whatsoever | Subset selection only |

## The rule that governs all of them

**Absence of annotation is not evidence of absence.** Every source above is a
record of what somebody looked at. Well-studied proteins accumulate annotations
across all of them; obscure proteins accumulate none. A site with no supporting
layer is `unknown`, never a negative. Treating unannotated sites as negatives
would measure curation effort and report it as glycobiology.

This is why `occupancy_status` has a third value, `observed_unmodified`, which
**none of the annotation sources below can populate**. Establishing it requires
evidence that a site was *examined and found bare*, and an annotation record
cannot supply that.

It is populated from structures instead, under two conditions that together make
absence informative: the entry models a glycan elsewhere in the same structure,
and the protein was expressed in a host competent to glycosylate. 32 sites across
25 proteins qualify — sequons with no modelled glycan under internal-control
conditions, which is not the same as a proven biochemical negative. See
`structures.py` and [`negative_controls.md`](control_sets.md).

---

## UniProt

`uniprot.py`, `evidence.py`. Read from the dated gzipped TSV snapshot at
`paths.uniprot_tsv`, never from the live API, so a run is reproducible against a
fixed release.

UniProt is the only source consulted for exact-position curated evidence,
because it is the only one whose **site-level features carry auditable
per-feature evidence codes**. Every other source reports a site; UniProt reports
a site plus the reason a curator believes it.

### What counts

Only a CARBOHYD feature that is all of the following:

- at **exactly** the candidate position — nearby features are never accepted;
- a single-residue position. Ranges (`10..12`) and uncertain positions (`?`) are
  parsed to `position = None` with `parse_status = "uncertain_or_range_position"`
  and excluded with that reason, rather than being silently truncated to their
  start coordinate;
- typed `N-linked` **on asparagine**. `_glyco_type` in `uniprot.py` separates
  N-linked asparagine features from `N-linked-other`, from O-, C- and S-linked
  features, and from **glycation** — a non-enzymatic sugar adduct that reads
  superficially like glycosylation in the note text and is not the same thing.

A feature carrying several ECO codes retains all of them; the site is classified
under its strongest tier.

### ECO tiers

Tier names follow the Evidence and Conclusion Ontology. Labels below are the
ontology's; follow the term URL to check any of them.

| Tier in this module | ECO code | Ontology label | Term URL |
|---|---|---|---|
| `manual_experimental` | ECO:0000269 | experimental evidence used in manual assertion | <https://evidenceontology.org/term/ECO%3A0000269/> |
| `manual_combinatorial` | ECO:0007744 | combinatorial computational and experimental evidence used in manual assertion | <https://evidenceontology.org/term/ECO%3A0007744/> |
| `manual_curator_inference` | ECO:0000305 | curator inference used in manual assertion | <https://evidenceontology.org/term/ECO%3A0000305/> |
| `sequence_similarity` | ECO:0000250 | sequence similarity evidence used in manual assertion | <https://evidenceontology.org/term/ECO%3A0000250/> |
| `manual_sequence_model` | ECO:0000255 | match to sequence model evidence used in manual assertion | <https://evidenceontology.org/term/ECO%3A0000255/> |
| `automatic_sequence_model` | ECO:0000256 | match to sequence model evidence used in automatic assertion | <https://evidenceontology.org/term/ECO%3A0000256/> |
| `automatic_sequence_model` | ECO:0000259 | match to InterPro signature evidence used in automatic assertion | <https://evidenceontology.org/term/ECO%3A0000259/> |
| `annotation_without_qualifying_evidence` | — | exact feature present, no recognised code | — |
| `exact_feature_absent` | — | no exact N-linked feature at the candidate position | — |

Precedence is the order in `UNIPROT_TIER_ORDER`; the strongest tier present
wins.

### ECO:0007744 does not mean "seen in a structure"

This matters enough to state on its own.

`ECO:0007744` is **"combinatorial computational and experimental evidence used
in manual assertion"**. It records that a curator combined computational and
experimental lines of evidence — in UniProt practice most often large-scale
glycoproteomics — to make the assertion. It says nothing about a structure. It
is not a PDB cross-reference, it does not imply coordinates exist, and it must
never be read as a structural observation.

**The earlier label in this repository was wrong.** An older script,
`pipeline/identify_sequons.py`, maps `ECO:0007744` to a tier it calls
`pdb_evidence` and documents as "UniProt ECO:0007744 PDB structural evidence";
`analysis/ortholog_sequon_conservation/docs/TERMS_IN_CONTEXT.md` describes the
code as "combinatorial evidence backed by a PDB structure". The design spec
records the same convention under the name "PDB-backed". All of these are
incorrect readings of the ontology term. This module does not use that label
anywhere, and the string appears nowhere in the package except in this paragraph
documenting that the earlier label was wrong. Structural evidence in this module
comes only from the structure layer below, which reads actual coordinate files.

Anyone reconciling this module's counts against older repository outputs should
expect the discrepancy: what the old tier called structural evidence is here
`manual_combinatorial`, and its structural claim is unsupported.

### Default policy and its composition

`policy.qualifying_uniprot_tiers = ["manual_experimental", "manual_combinatorial"]`.
That yields the frozen baseline of **505 sites in 401 proteins**.

Of those 505:

- **498** carry `ECO:0000269`;
- **exactly 7** rest solely on `ECO:0007744`. Because that code is the weaker and
  more often misread of the two, those 7 can be dropped in a sensitivity pass.
  Isolate them by filtering `experimental_sites_uniprot_baseline.csv` on
  `uniprot_tier == "manual_combinatorial"`.

Neither qualifying tier implies structural observation.

`manual_curator_inference` (ECO:0000305) does **not** qualify by default, but is
written to `curator_inferred_sensitivity_sites.csv` so the effect of admitting
it can be measured rather than argued about.

### Exclusion reasons

Every candidate that does not qualify carries one machine-readable reason, so
exclusions are auditable rather than assumed:

| Reason | Meaning |
|---|---|
| `exact_feature_absent` | Accession is in the snapshot; no exact N-linked feature at this position |
| `sequence_model_only` | Strongest tier is a sequence-model match (ECO:0000255/0000256/0000259) |
| `curator_inference_only` | Strongest tier is ECO:0000305 |
| `sequence_similarity_only` | Strongest tier is ECO:0000250 |
| `annotation_without_qualifying_evidence` | Exact feature present, no recognised evidence code |
| `accession_absent_from_snapshot` | Accession is not in the dated UniProt snapshot at all |

`uniprot_exact_n_linked_sites.csv` holds every parsed feature, including ones no
candidate uses, so any exclusion can be checked against the source text rather
than taken on trust.

---

## GlyGen

`glygen.py`. `GET https://api.glygen.org/protein/detail/{accession}/`, cached to
`data/cache/glygen_protein_detail.jsonl`.

GlyGen aggregates glycosylation reports from many upstream sources. It is a
genuinely independent layer — but only for some of its categories, and the
distinction is the whole point of this section.

### The `site_category` taxonomy

Only `n-linked` entries with `start_pos == end_pos` are used; ranges and
non-N-linked entries are dropped. Each remaining entry is classified by
`classify_glygen_entry`:

| GlyGen category | Tier assigned | Qualifies? | Why |
|---|---|---|---|
| `reported_with_glycan` | `glygen_reported_with_glycan` | yes | An identified glycan structure is attached to the site by an independent report |
| `reported`, with at least one non-UniProtKB evidence database | `glygen_reported_independent` | yes | Genuinely independent corroboration |
| `reported`, evidence databases are UniProtKB only | `glygen_reported_uniprot_derived` | **no** | Circular — this is UniProt's own call round-tripped through GlyGen |
| `predicted` / `predicted_with_glycan` | `glygen_predicted` | **no** | This is UniProt's sequence-model rule re-exported |
| anything else | `glygen_unclassified` | no | Unrecognised |

`_UNIPROT_DATABASES` = `{UniProtKB, UniProtKB/Swiss-Prot, UniProtKB/TrEMBL}`.
When a site has several entries, the strongest tier wins and the union of all
categories, evidence databases, PubMed ids and GlyTouCan accessions is retained
in `glygen_site_evidence.csv`.

### Why the two exclusions are not fussiness

Both were measured, over a 60-accession sample of GlyGen protein-detail
payloads:

- **`predicted` is the UniProt rule.** Its evidence records cite UniProtKB in
  **136 of 137** cases (the remaining one cites PDB). Counting `predicted` would
  readmit sequence-model predictions into a dataset whose entire value is that
  it contains none.
- **`reported_with_glycan` is genuinely independent.** Its evidence records cite
  PubMed (46), PDC (40), DOI (25), GlyConnect (6) and Data Submission (3), with
  **zero UniProtKB citations**.
- **`reported` is mixed**, which is why it is split on its evidence databases
  rather than accepted wholesale. A `reported` entry whose only citation is
  UniProtKB tells you nothing UniProt did not already say; treating it as a
  second independent layer would double-count one curator's judgement and
  inflate `support_count`.

### Validation against the UniProt gold standard

50 sites already qualifying under UniProt's experimental policy were queried
against GlyGen. **All 50 came back `reported` or `reported_with_glycan`. None
came back `predicted`, and none were absent.** Where both sources have an
opinion about a site known to be experimentally supported, they agree — which is
what makes the `predicted` exclusion safe rather than merely conservative.

### What GlyGen cannot establish

That a site is unglycosylated. GlyGen's absence of a site means no aggregated
source reported it, which is dominated by how much attention the protein has
received. It is also taxonomically skewed toward human and mouse; see
`analysis_options.md`.

---

## GlyConnect

`glyconnect.py`. Off by default (`layers.glyconnect = false`,
`policy.glyconnect_qualifies = false`).

Retrieval is two steps:

1. `GET https://glyconnect.expasy.org/api/proteins?uniprot={accession}` returns
   a list of hits; the first hit's `id` is taken.
2. `GET https://glyconnect.expasy.org/api/proteins/{protein_id}` returns
   `structures`, each with a `type` and `composition_string` and a `sites` list.

Sites are strings, not integers. The format is `"Asn-80"` — a residue name,
a hyphen, then the 1-indexed UniProt position. `parse_site_label` accepts only
that pattern (case-insensitively) and returns `None` for anything else, so
O-linked serine/threonine site strings and free-text labels are rejected rather
than mis-parsed. Only structures typed `n-linked` are read. Per position the
module records how many glycan structures were attributed to it and the set of
composition strings.

**Why it is corroboration only, and off:**

- **Coverage is low.** Only **183 of 2,591** non-qualifying accessions carry a
  GlyConnect cross-reference at all.
- **It is not independent of GlyGen.** GlyGen already ingests GlyConnect —
  GlyConnect appears as an evidence database behind `reported_with_glycan`
  entries. Enabling both as qualifying layers would count the same underlying
  report twice and inflate `support_count`.

Enable it to check whether a specific site has an attributed glycan structure,
or as a sensitivity analysis. Set `policy.glyconnect_qualifies = true` only
knowing that it partially double-counts GlyGen.

---

## GlyCosmos

Cross-reference only. **Not implemented as a layer, deliberately.**

GlyCosmos has no clean per-site REST endpoint that returns glycosylation
positions for a UniProt accession in the way the GlyGen protein-detail endpoint
does. The browser-facing protein URL redirects into the web application rather
than serving a stable machine-readable payload, so there is nothing to parse
reproducibly and nothing to cache. Scraping the rendered page would produce
site-level claims with no auditable provenance, which is exactly what this
module exists to avoid.

It remains a useful place for a human to check a specific protein by hand. It
contributes nothing to `support_count`.

---

## PDB linkage records

`structures.py`. The only layer that reads direct physical evidence.

A `LINK` record in a PDB coordinate file names a covalent bond explicitly:

```
LINK         ND2 ASN A  86                 C1  NAG A 477A     1555   1555  1.54
```

That record says: atom `ND2` of `ASN 86` in chain `A` is bonded to atom `C1` of
`NAG 477A` in chain `A`, at 1.54 Å. The format is fixed-column, and
`parse_link_records` reads it by column, not by splitting on whitespace:

| Columns (1-indexed) | Field | First partner | Second partner |
|---|---|---|---|
| 13-16 / 43-46 | atom name | `ND2` | `C1` |
| 18-20 / 48-50 | residue name | `ASN` | `NAG` |
| 22 / 52 | chain id | `A` | `A` |
| 23-26 / 53-56 | residue sequence number | `86` | `477` |
| 27 / 57 | insertion code | (blank) | `A` |

`parse_link_records` takes both partners in
either order, and keeps the record when one side is `ASN` and the other is a
glycan residue in `GLYCAN_RESNAMES` = `{NAG, NDG, BGC, GLC, MAN, BMA, FUC, GAL,
XYS}`. Insertion codes are preserved throughout — `477A` in the example above is
real, and dropping the `A` would map the linkage to the wrong residue.

A LINK between an asparagine and a glycan residue is **direct physical evidence
that this asparagine carried a glycan in that structure**. It is the strongest
single statement any of these sources makes, because it is a modelled covalent
bond rather than a curator's or an aggregator's assertion.

Measured over the first 400 cached structures: **88 carry an ASN-glycan LINK
record.**

### Mapping UniProt positions onto chains

The UniProt sequence is locally aligned to each chain's observed residues to get
`uniprot_position → (chain, resseq, icode)`, then the mapped residue is checked
against the linkage set. Chain selection is gated in a specific order, and the
order matters:

1. **Credibility first.** A chain is admitted only if the alignment covers at
   least `MIN_CHAIN_COVERAGE = 0.5` of the chain and spans at least
   `MIN_ALIGNED_RESIDUES = 30` residues. Identity alone is not enough: it is
   measured over the aligned block, so an unrelated chain sharing a short motif
   scores identity 1.0 on a two-residue block and would otherwise pass.
2. **Identity, not coverage, decides "same protein."** Among credible chains,
   those within `SAME_PROTEIN_IDENTITY_TOLERANCE = 0.02` of the best identity are
   treated as copies of the same protein. Homodimer copies align at essentially
   identical identity while differing in coverage.
3. **Prefer a linked chain** among those copies, then rank the rest by coverage.

If nothing aligns credibly, the site is reported as
`structure_residue_resolved` with `structure_detail = "low_confidence_chain_match"`
and **no glycan is ever asserted from it**. The failure direction throughout is
deliberate: a missed glycan, never a fabricated one. See the limitations section
of `analysis_options.md` for two consequences of this design.

### The tier ladder

| `structure_tier` | Meaning | Qualifies? |
|---|---|---|
| `structure_linked_glycan` | ASN-glycan LINK at the mapped residue | yes |
| `structure_residue_resolved` | Residue is present in the model with no linkage attached | no |
| `structure_residue_unresolved` | Position is not in the model (`position_not_in_model`) | no |
| `structure_not_assessed` | No cached structure, mmCIF-only entry, unreadable file, or no UniProt sequence | no |

`structure_detail` carries the specific reason and must be read alongside the
tier.

`structure_observed_residue` records the amino acid actually found at the mapped
position. An N-linked site must map to an asparagine; anything else is a mapping
failure that would otherwise be invisible. Check it before trusting a row.

### Absence of a linkage is not evidence of absence

`structure_residue_resolved` means only that the residue was modelled and no
linkage record pointed at it. It must **never** be read as "observed
unmodified". At least three routine causes produce a bare asparagine in a
deposited structure that is glycosylated in life:

- **Deglycosylation before crystallisation.** Glycans are heterogeneous and
  flexible and frequently obstruct crystallisation, so they are commonly removed
  enzymatically before the experiment.
- **Bacterial expression.** Structures of eukaryotic proteins are often solved
  from *E. coli*-expressed material, which does not perform eukaryotic N-linked
  glycosylation at all.
- **Unmodelled density.** Even when a glycan is present it is often disordered
  beyond the first sugar or entirely, and what cannot be modelled is not
  deposited — and if the depositor did not build the linkage, there is no LINK
  record to read.

The mmCIF case is a fourth, purely technical, gap: `parse_link_records` handles
the fixed-column PDB format only, so mmCIF-only manifest entries return
`structure_not_assessed` with `structure_detail = "mmcif_linkage_unsupported"`
rather than a false negative. Parsing `_struct_conn` is phase-2 work.

---

## OrthoDB and homology QC

**Subset selection only. Never occupancy evidence.**

`orthologs.py` uses the canonical pair table and the homology QC table to decide
which ortholog subset each site belongs to. A site joins the strict subset if at
least one of its analysis-ready associations is `strict_ortholog_like`, and the
strict-plus-plausible subset if at least one is `strict_ortholog_like` or
`plausible_ortholog_like`. Membership is by `any`, and deduplication onto
`(accession, position)` happens after the join.

This is orthology, not glycobiology. That a site sits in a high-quality ortholog
pair says something about how comparable that pair is; it says nothing about
whether a sugar is attached. `assign_subsets` never touches `support_sources`,
`support_count`, `experimental_positive` or `occupancy_status`, and no homology
signal appears in `combine_layers`. Cluster and pair identifiers are kept in
`site_pair_associations.csv` — outside the site tables — precisely so they cannot
leak into an occupancy claim.

Their real analytical use is downstream: `cluster_id` is what makes
cluster-aware resampling possible, because sites within a cluster are not
independent. See `analysis_options.md`.
