# Control sets

Why the comparison sets exist, what each uses as evidence for "no glycan here",
how the diagnostics behave, and what the data audit found. Merges the former
`negative_controls.md`, `diagnostic_controls.md` and `phase1_control_audit.md`.

Terms are defined in [`glossary.md`](glossary.md); current results live in
[`OVERVIEW.md`](OVERVIEW.md) and are deliberately not repeated here.

---

## Why these exist

The resource holds 922 sites with experimental evidence of a glycan, 32 observed
unmodified under an informative experimental context, and 3,353 whose occupancy is
unknown. That is enough to ask whether a protein model has learned the sequon motif or
its biological use — but not enough to answer it convincingly, because 32 negatives
cannot support a null result, and because any comparison between occupied and
unoccupied sites is open to the objection that the two groups differ in something other
than occupancy.

These control sets exist to close that gap. They are **not** part of the dataset. They
are never pooled with the 4,307 candidate sites, never counted in any site total, and
always carry a `control_set` label recording where they came from.

## The problem they solve

Occupied sites are not a random sample of sequons. They sit in secreted and membrane
proteins, in solvent-exposed loops, in well-studied organisms. Any of those properties
can be read by a structure-based model, so a model that separates occupied from
unoccupied sites may be detecting one of them rather than occupancy itself.

The defence is not a better single control. It is several controls whose confounds do
not overlap, so that no one alternative explanation can account for a result that holds
across all of them.

| Control set | Matches the positives on | Confounded by | Approx. sequons |
|---|---|---|---|
| Structural internal controls (the 32) | organism, compartment, experiment | very small; structurally biased | 32 |
| Cytosolic eukaryotic | taxonomy, curation depth | subcellular compartment | 19,337 |
| Bacterial extracytoplasmic | compartment, membrane translocation | taxonomy, fold repertoire | 5,865 |
| Eukaryotic secretory, unannotated | taxonomy **and** compartment | the negative label itself | 4,418 |

The first three were designed to be read together: their confounds do not overlap, so
if a model separates occupied sites from the cytosolic set but not the 32 it has learned
where proteins live, and if it separates them from the bacterial set but not the
cytosolic set it has learned taxonomy.

That triangulation was undermined in practice — corrected scores put the three
comparisons in three different directions with no interpretable ordering, so the gap
between them stopped being a measurement. Set 3 takes a different route: rather than
combining two confounded sets, it removes both confounds at once and accepts a weaker
negative label in exchange. It is the comparison to read first, with its contamination
caveat attached.

## What each set uses as evidence for "no glycan here"

This is the part that matters most, because the four sets do **not** make the
same kind of claim. They are ordered here from strongest negative label to
weakest.

### First, what counts as a positive

A site is `occupied_supported` if **any** of three independent layers supports
it. The layers are not required to agree.

| Layer | What it requires |
|---|---|
| UniProt | a `CARBOHYD` feature at that exact residue whose ECO codes reach `ECO:0000269` (manual experimental) or `ECO:0007744` (manual combinatorial) |
| GlyGen | a reported glycan at that accession and position |
| Structure | a covalent ASN–glycan linkage on that residue in a deposited entry (`LINK` records in PDB, `_struct_conn` in mmCIF) |

`ECO:0000305` (curator inference) is deliberately **not** sufficient. It is
tracked separately and reported as a sensitivity set, because a curator's
inference is not an experiment.

Everything below is about the mirror-image question, which is much harder.

### Tier 1 — observed absence: the 32 internal controls

The only set where absence was actually *looked at*. All four conditions must
hold:

1. the residue is resolved in a deposited structure;
2. it carries no glycan linkage;
3. the **same structure** models at least one glycan at some other residue;
4. the protein was expressed in a host competent to glycosylate.

Conditions 3 and 4 are what make the absence informative. They establish that
sugars survived sample preparation and that this depositor was both able and
willing to model them. A bare asparagine under those conditions is a decision
rather than a silence.

Still not proof: the glycan may have been present, disordered, and left
unmodelled. This is a statement about the deposited model, not the molecule.

### Tier 2 — absence by biology: cytosolic and bacterial

No annotation is consulted for the negative claim at all. The argument is
mechanistic: the enzyme and the substrate never meet.

| | Cytosolic eukaryotic | Bacterial extracytoplasmic |
|---|---|---|
| Inclusion | `GO:0005829` (cytosol), `taxonomy_id:2759` | `taxonomy_id:2`, periplasm / outer membrane / secreted (`KW-0574`, `KW-0998`, `KW-0964`) |
| Exclusion | signal peptide `KW-0732`, transmembrane `KW-0812`, glycoprotein `KW-0325` | glycoprotein `KW-0325`, plus every clade with known machinery |
| Why no glycan | never enters the secretory pathway, so never meets OST | the clade has no OST or N-glycosyltransferase |

The clade exclusions matter and are done by **known machinery**, not by
annotation: Archaea wholesale (AglB is a genuine OST), plus *Campylobacter*
(PglB), *Helicobacter*, *Haemophilus*, *Actinobacillus*, *Yersinia* and
*Kingella* (HMW1C-type N-glycosyltransferases, which act on N-X-S/T). Excluding
these by "not annotated as glycosylated" would repeat exactly the error this
project exists to avoid.

The glycoprotein-keyword exclusion is kept as a second line of defence, never as
the primary argument.

### Tier 3 — absence of annotation: the eukaryotic secretory set

The weakest claim, and the reason it is labelled honestly.

| | |
|---|---|
| Inclusion | `taxonomy_id:2759`, `database:pdb`, and secreted **or** transmembrane **or** signal peptide (`KW-0964`, `KW-0812`, `KW-0732`) |
| Exclusion | glycoprotein `KW-0325` |
| Why no glycan | **nobody has recorded one** |

That is the whole argument. There is no mechanistic reason these sequons cannot
be glycosylated — they are in the right compartment, in the right kind of
organism, on the right kind of protein. The only thing standing behind the
negative label is that UniProt has no glycoprotein keyword for the entry.

Two things partially defend it:

- **A cross-check against our own layers.** Every one of the 4,418 sequons was
  tested against the 922 occupied sites, the GlyGen evidence table and the
  structural glycan linkages. Zero hits. So the negative label does not rest on
  a single source agreeing with itself.
- **The direction of the error.** Contamination makes the groups more alike, so
  it pulls any difference towards zero.

Neither makes it a clean negative, and it should never be described as one.

### Summary

| Set | Evidence for "no glycan" | Strength | Sequons |
|---|---|---|---|
| Internal controls | observed bare, in a structure that shows glycans elsewhere | strongest | 32 |
| Cytosolic eukaryotic | compartment makes contact with OST impossible | strong, confounded | 19,337 |
| Bacterial extracytoplasmic | clade has no glycosylation machinery | strong, confounded | 5,865 |
| Eukaryotic secretory | no annotation exists | weakest, unconfounded | 4,418 |

The trade runs diagonally: the sets with the firmest negative label are the ones
whose populations differ most from the occupied sites, and the set that matches
the occupied sites best has the least defensible label. No set is best on both
axes, which is why more than one is reported.

## Set 1 — cytosolic eukaryotic

Proteins localised to the cytosol, with no signal peptide and no transmembrane region.

N-linked glycosylation happens in the ER lumen, catalysed by oligosaccharyltransferase.
A protein that never enters the secretory pathway never meets OST, so its sequons are
not merely unobserved — they are unreachable. This is positive knowledge from the
biology rather than an inference from missing annotation, which makes these stronger
negatives than anything the annotation databases can supply.

The cost is that the comparison now differs in compartment. Secreted and cytosolic
proteins differ in composition, disulfide content, hydrophobic patterning and fold
family, all visible to a structural model.

**A precision point:** cytosolic and nuclear proteins *are* glycosylated, by O-GlcNAc
transferase, on serine and threonine. That is different chemistry and does not
contaminate an N-linked negative set, but the claim must be written as "never
N-glycosylated", not "not glycosylated".

## Set 2 — bacterial extracytoplasmic

Periplasmic, outer-membrane and secreted bacterial proteins, from clades with no known
N-glycosylation machinery.

This fixes what is wrong with set 1. These proteins cross a membrane, carry signal
peptides, and fold in an oxidising, disulfide-forming compartment — the periplasm is the
closest bacterial analogue to the ER lumen. What differs from the positives is the
presence of an OST, which is much closer to isolating the variable of interest.

The cost is taxonomic: bacterial and eukaryotic secreted proteins differ in fold
repertoire, composition, and evolutionary history with respect to glycans.

### Exclusions, and why each one

Bacterial N-glycosylation is rare but real, and excluding it by annotation alone would
repeat the absence-of-evidence error this project exists to avoid. Clades are therefore
excluded by **known machinery**:

| Excluded | Reason |
|---|---|
| Archaea (all) | Glycosylate via AglB; a genuine OST, not an edge case |
| *Campylobacter* | PglB, the best-characterised bacterial OST |
| *Helicobacter* | PglB-family machinery |
| *Haemophilus* | HMW1C-type cytoplasmic N-glycosyltransferase, acts on N-X-S/T |
| *Actinobacillus* | HMW1C homologue |
| *Yersinia* | HMW1C homologue |
| *Kingella* | HMW1C homologue |

The annotation filter is kept as a second line of defence, not as the primary one.

## Set 3 — eukaryotic secretory, unannotated

Added 2026-08-19, after the internal-control class was ruled out as a growth
route. This set makes the opposite trade from the other two, and it is the
reason it exists.

Sets 1 and 2 each hold the negative label firmly — these proteins genuinely
cannot be N-glycosylated — and pay for it with a confound: compartment for the
cytosolic set, taxonomy for the bacterial one. Set 3 removes **both** confounds
by taking eukaryotic secreted and membrane proteins, the same kind of protein in
the same kind of cell as the occupied sites, and pays for it in the label
instead. These sequons are *not annotated* as glycosylated, which is not the same
as being *annotated unglycosylated*.

**Query:** reviewed, `taxonomy_id:2759`, `database:pdb`, secreted or
transmembrane or signal-peptide (`KW-0964`, `KW-0812`, `KW-0732`), and NOT
`KW-0325` (Glycoprotein).

**Scale:** 3,619 proteins, 4,418 sequons in 1,543 of them, one deposited
structure fetched per protein. Human-dominated (912 human, 283 mouse), which
matches the occupied set's taxonomy closely.

**Purity:** zero overlap with any of this project's own positive evidence — the
922 occupied sites, GlyGen, and structural glycan linkages all return nothing.
The UniProt keyword exclusion is doing real work rather than deferring to a
source we already used.

### The contamination, stated up front

About **half** of all eukaryotic secretory proteins with a solved structure carry
a glycoprotein keyword (4,067 against 3,619). So the unannotated half certainly
contains genuine glycosylation sites that nobody has recorded. This set has false
negatives by construction and there is no way to remove them without the
annotation that is missing by definition.

What saves it is the *direction* of the error. False negatives make the two
groups more alike, so they pull any measured difference **towards zero**. That
cuts two ways and the asymmetry matters:

- The set can support **"no difference detected"** well, because dilution is
  working against finding a difference and one still was not found.
- The set could **not** be used to argue a positive finding away, and equally, a
  positive difference here would be *more* striking than it looks, since it
  survived dilution rather than being created by it.

This should be decided before the result is seen, and it is recorded here for
that reason.

## A second use, independent of the modelling

In an organism with no OST, sequons are under no glycosylation-related selection. The
bacterial set therefore describes the **neutral sequon distribution** — what sequon
density and local context look like with the selective pressure removed.

That makes it a null model for the evolutionary half of the project: whether sequons in
eukaryotic secreted proteins are enriched or depleted relative to neutral expectation, at
matched amino acid composition. This question needs no model scoring at all.

## What the built sets look like

| Set | Proteins | Sequons | Per 1000 residues |
|---|---|---|---|
| Cytosolic eukaryotic | 9,561 | 19,337 | 4.54 |
| Bacterial extracytoplasmic | 1,410 | 5,865 | 8.75 |

Bacterial extracytoplasmic proteins carry nearly twice the sequon density, and that is
composition rather than biology: they are richer in asparagine (0.060 against 0.039) and
threonine (0.069 against 0.052), and poorer in the proline that blocks a sequon.

The useful check is against the density expected from composition alone, treating each
position independently:

| Set | Expected | Observed | Ratio |
|---|---|---|---|
| Cytosolic eukaryotic | 4.55 | 4.54 | 1.00 |
| Bacterial extracytoplasmic | 8.26 | 8.75 | 1.06 |

Both sit essentially at neutral expectation, which is what proteins under no
glycosylation-related selection should do. That is a direct check on the premise these
sets rest on, and it passed. It also sets up the comparison that matters for the
evolutionary question: whether eukaryotic *secreted* proteins deviate from the same
expectation, which would be selection on sequons rather than chance.

## Rules these sets must obey

1. **Never pooled.** Each carries a `control_set` label and lives in its own table.
   They are not sites in the resource and are excluded from every site count.
2. **Never counted as evidence.** They say nothing about the 4,307 candidates.
3. **Experimental structures only.** If positives used deposited structures and controls
   used predicted ones, structure source would confound compartment and taxonomy on top
   of everything else.
4. **Composition reported.** Amino acid usage differs between the sets, so sequon density
   differs for reasons unrelated to glycosylation. Any comparison that depends on it must
   match on it.
5. **Outside the ortholog framework.** These proteins have no ortholog-pair context and
   cannot participate in the sequon-loss analysis.

## Considered and rejected

*Plasmodium falciparum* would control taxonomy and compartment at once, being a eukaryote
with almost no N-glycosylation. Its proteome is extremely AT-rich and low-complexity and
few structures exist, so it trades a clean confound for a messy one.

## How to build them

Stage 1 — sequence inventory and sequon scan. Fast, no structure downloads:

```bash
python -m experimental_glycosylation_sites fetch-controls
```

Writes `results/datasets/negative_control_sites.csv` (one row per sequon, labelled by
`control_set`) and `results/datasets/negative_control_summary.json` (counts, sequon density,
amino acid composition, and the excluded taxa with their reasons). Protein sequences go
to `data/cache/negative_control_proteins.csv.gz`, since they are bulky and only the
feature stage needs them.

**Stage 2 — structural features — is built and has run.** It was separated from
stage 1 because it needs deposited structures for the control proteins, a download on
the order of thousands of entries. The feature code is shared with the main dataset
(`features.residue_features`).

Of 15,107 control sequons, **6,823 have structural features** (3,543 cytosolic, 3,280
bacterial). Of those, **6,092 can be scored by ProteinMPNN** (3,024 and 3,068):
scoreability is settled from the coordinates alone, before matching, because a residue
with an incomplete backbone is silently returned as a non-distribution. See
`pipeline/05_scoreability.py` and `[amendment_1]` in `config/scoring_frozen.toml`.

```bash
python pipeline/04_build_candidate_manifest.py controls
python pipeline/05_scoreability.py results/manifests/candidate_manifest_controls.csv \
                               results/manifests/scoreability_controls.csv
python pipeline/06b_match_diagnostics.py
```

Results are in [`diagnostic_controls.md`](control_sets.md). Both sets remain
**diagnostics**: each is confounded by construction, and neither substitutes for the 32
internal controls, which are the only ones holding organism, compartment and
experimental context roughly constant.

---

## The diagnostics in detail — bacterial and cytosolic

These two sets are confounded on purpose and in opposite directions. They are
informative about how the measurement behaves and are not a valid answer to the
occupancy question.

**These are not the primary result and cannot substitute for it.** The primary

> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](OVERVIEW.md) for current numbers.

comparison is [`primary_result.md`](archive/primary_result_SUPERSEDED_2026-08-25.md).

**A fourth control set now exists** — eukaryotic secretory, unannotated — which is
*not* a diagnostic. It matches the occupied sites on both taxonomy and
compartment, so it carries neither confound described below, and it supplies 262
matched pairs against the primary's 16. It gives +0.073 SD, CI [−0.056, +0.346].
Its trade is a weaker negative label rather than a confounded population; see
[`negative_controls.md`](control_sets.md) and
[`primary_result.md`](archive/primary_result_SUPERSEDED_2026-08-25.md). The two sets below remain diagnostics.

Both control sets are confounded by construction, and deliberately so. The
cytosolic set matches taxonomy and differs in subcellular compartment; the
bacterial set matches compartment and differs in taxonomy. They were built to be
read against one another, on the reasoning that a signal surviving both
orthogonal confounds is less likely to be produced by either.

They were also built because there are only 28 scoreable internal controls. That
scarcity is the real constraint, and a large well-matched comparison against the
wrong population does not relieve it.

### Matching

Same discipline as the primary: scoreable pool fixed before matching, 1:1
without replacement, exact NXS/NXT, caliper 0.25 on relative accessibility,
neighbour count and hydrophobic fraction.

| | Bacterial | Cytosolic |
|---|---|---|
| Scoreable controls available | 3,068 | 3,024 |
| Matched pairs | 278 | 270 |
| Occupied proteins / ortholog clusters | 210 / 95 | 205 / 96 |
| Resampling units | 58 | 75 |
| Most-reused control protein | 4 contrasts | 3 contrasts |
| Worst \|SMD\| after matching | 0.003 | 0.009 |

Control-protein reuse is why the resampling unit is the connected component of
the graph joining occupied ortholog clusters to shared control proteins. With 278
contrasts across only 58 such units, a site-level interval would be far too
narrow — as the comparison below shows.

### Results

All on the common reference scale (SD = 1.3316 log-odds, from dataset sites only).

Primary figures are from the deterministic optimal matching (Amendment 2).

| Comparison | n | Mean | 95% CI (cluster-aware) | Site-level CI (not used) | Verdict |
|---|---|---|---|---|---|
| **Internal control** (primary) | 16 | **+0.458 SD** | [−0.227, +1.098] | [−0.200, +1.121] | inconclusive |
| Bacterial extracytoplasmic | 278 | −0.174 SD | [−0.459, −0.027] | [−0.334, −0.015] | directional, magnitude undetermined |
| Cytosolic eukaryotic | 270 | +0.062 SD | [−0.145, +0.270] | [−0.112, +0.233] | inconclusive |

Robustness:

| | Bacterial | Cytosolic |
|---|---|---|
| Occupied higher | 125 of 278 | 151 of 270 |
| Median | −0.191 | +0.214 |
| Sign test | p = 0.105 | p = 0.059 |
| Wilcoxon | p = 0.044 | p = 0.228 |

### What this does and does not support

**The earlier "gradient" is withdrawn.** The previous version of this analysis
reported all three comparisons as negative, shrinking as matching improved
(−0.237, −0.145, −0.057 SD), and read that ordering as evidence that the
apparent effects came from compartment and taxonomy rather than occupancy. That
pattern was an artefact of the corrupted scores. It does not survive correction.

What replaces it is less tidy. The three comparisons now point in **different
directions**: positive but not significant against internal controls, negative
against bacterial controls, and indistinguishable from zero against cytosolic
controls. There is no monotone ordering to interpret.

That is a genuine inconsistency and it should be treated as one. The most
economical reading is that the two diagnostic contrasts are dominated by their
respective confounds — bacterial folds and cytosolic environments differ from
secreted eukaryotic ones in ways ProteinMPNN can see directly — and that neither
tells us much about occupancy. But this is the reading that was *not* available
before correction, when the confounds appeared to be pushing the same way, and
it should not be presented as a finding.

The primary comparison is the only one that holds organism, compartment and
experimental context roughly constant. It has 16 pairs. Nothing in this appendix
changes that, and nothing here should be quoted as support for the primary
estimate.

### Subtype

| Comparison | NXS | NXT |
|---|---|---|
| Bacterial | −0.167 SD (n=132) | −0.181 SD (n=146) |
| Cytosolic | −0.094 SD (n=135) | +0.218 SD (n=135) |

The cytosolic subtypes disagree in sign at n=135 each, which is not a sample-size
problem. It is unexplained and recorded as such.

### Artefacts

`results/matching/matched_pairs_{bacterial,cytosolic}.csv`,
`results/matching/matching_balance_{bacterial,cytosolic}.json`,
`results/analysis/contrasts_{bacterial,cytosolic}.csv`,
`results/analysis/analysis_{bacterial,cytosolic}.json`,
`results/scores/scores_controls.csv`, `results/manifests/manifest_matched_controls.csv`,
`results/manifests/scoreability_controls.csv`.

Reproduce with `pipeline/06b_match_diagnostics.py`, then `pipeline/07_score.py` on
`results/manifests/manifest_matched_controls.csv`, then
`pipeline/09_analyse_scores.py bacterial` and `... cytosolic`.

---

## Data audit — what verification found



> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](OVERVIEW.md) for current numbers.

> **Historical: a pre-scoring snapshot.** This records the state of the control
> data at the gate, before any ProteinMPNN scoring ran, and is kept as the record
> that the gate was passed. Its counts predate the corrections of 18 August 2026:
> sites ProteinMPNN cannot decode had not yet been identified, so the totals here
> are larger than the scoreable sets the analysis actually used.
>
> Current figures: [`primary_result.md`](archive/primary_result_SUPERSEDED_2026-08-25.md) and
> [`diagnostic_controls.md`](control_sets.md). What changed and why:
> [`correction_2026-08-18_SUPERSEDED_2026-08-31.md`](archive/correction_2026-08-18_SUPERSEDED_2026-08-31.md).

Gate for model scoring. Nothing in Phase 2 onward may proceed on data that fails
these checks.

### Verdict

**Pass, after two defects were found and repaired.** Both were real and both
would have invalidated the comparisons had scoring started on the previous
state.

### Defects found

#### 1. Provenance scrambled — 43.5% of rows mislabelled

`build_features` re-sorts its output by `(accession, position)`. The caller
attached `control_set` afterwards positionally, pairing labels from input order
against rows in sorted order.

| | bacterial (source) | cytosolic (source) |
|---|---|---|
| **labelled bacterial** | 3,258 | 2,608 |
| **labelled cytosolic** | 2,607 | 3,505 |

5,215 of 11,978 rows carried the wrong set, so the bacterial and cytosolic
matched comparisons were thoroughly mixed. `occupancy_status`, which travelled
with the row rather than being reattached, was correct — which is what localised
the fault.

**Repair.** `build_features` now takes `carry_columns` and copies provenance per
row, so no caller can reattach positionally. The joining script uses a key-based
merge with `validate="one_to_one"`. A regression test supplies input in
reverse-sorted order, which fails against the old behaviour.

#### 2. The cytosolic set was not eukaryotic

`go:0005829` carries no taxonomic restriction, and none was applied. The set
contained at least **1,929 bacterial** and **302 archaeal** proteins; *E. coli*
was its third most common organism.

The archaeal contamination is the serious part. Archaea N-glycosylate via AglB,
so those were potentially genuine glycoproteins inside a set whose definition is
that its members cannot be glycosylated. It also explains the 13 accessions that
appeared in *both* control sets: bacterial proteins annotated cytosolic and
periplasmic satisfied both queries.

**Repair.** `taxonomy_id:2759` added to the cytosolic query; inventory and all
downstream files rebuilt.

### Assertions

| Check | Result |
|---|---|
| `control_set` agrees with source inventory | **0** disagreeing rows |
| `occupancy_status` equals `control_<control_set>` | **0** disagreeing rows |
| Each matched comparison contains only its intended controls | **0** violations, across 3 comparisons × 5 seeds |
| No control protein in both biological control sets | **0** accessions |
| Cytosolic set is eukaryotic | **0** bacterial or archaeal organisms |

### Final counts

**Control inventories**

| Set | Proteins | Sequons | With structural features |
|---|---|---|---|
| Cytosolic eukaryotic | 4,764 | 16,386 | 3,543 sequons / 1,719 proteins |
| Bacterial extracytoplasmic | 1,141 | 5,865 | 3,280 sequons / 1,033 proteins |

Cytosolic fell from 19,337 to 16,386 sequons when the non-eukaryotic proteins
were removed. Bacterial was unaffected, as expected.

**Dataset sites with structural features**

| Class | Sites | Proteins |
|---|---|---|
| `occupied_supported` | 332 | 245 |
| `observed_unmodified` | 32 | — |

**Matched sets** (seed 0; 1:1 for the scarce comparison, 1:5 otherwise)

| Comparison | Controls available | Occupied sites matched | Pairs |
|---|---|---|---|
| vs observed-unmodified | 32 | 24 | 24 |
| vs bacterial control | 3,280 | 280 | 1,180 |
| vs cytosolic control | 3,543 | 275 | 1,196 |

#### Reading the observed-unmodified comparison correctly

It contains **24 occupied sites matched to 24 observed-unmodified controls**.
Those 24 are drawn from the **332 structurally scoreable occupied sites** — they
are not "24 of 32". The 32 is the size of the available control pool, of which 24
found a partner inside the caliper. The statistical unit is the occupied site,
and there are 24 of them in this comparison.

### Weighted balance

Each occupied site carries total weight one; the controls matched to it share
weight one between them. Without this, a case matched to five controls would
contribute five times as much to the control mean as a case matched to one, so
the statistic would describe the matching's bookkeeping rather than the
comparison. Verified: total control weight equals the number of matched cases.

Standardised mean differences, seed 0, with the range across seeds 0–4:

| Comparison | RSA | Neighbours | Hydrophobic fraction |
|---|---|---|---|
| vs observed-unmodified | +0.017 (−0.049..+0.017) | 0.000 | 0.000 |
| vs bacterial | +0.011 (−0.005..+0.021) | 0.000 | −0.000 |
| vs cytosolic | +0.020 (+0.020..+0.030) | 0.000 | −0.000 |

Worst absolute value across all comparisons, features and seeds: **0.030**,
comfortably inside the conventional 0.1 threshold. Seed-to-seed spread is at most
0.066 and mostly far smaller, so the matching is not seed-sensitive.

**The exact zeros are real, not a bug.** `n_neighbours_8a` is integer-valued and
`hydrophobic_fraction_8a` is a ratio over 5–6 neighbours, so both are effectively
discrete and the matcher finds exact partners for 100% of pairs. Confirmed at
pair level; RSA, being continuous, matches to within about 0.01. It was checked
precisely because a reported SMD of exactly zero usually indicates an error.

### Known limitations carried forward

- Structures above 20 MB are skipped at parsing, mapping and linkage detection,
  so ribosomal and proteasomal control proteins largely lack features. This is
  reported as `structure_too_large`, not silently dropped, but it biases the
  control subset against large-complex components.
- Fewer than half of control sequons have features (6,823 of 15,107), mostly
  through unmapped residues and the size guard.
- The 32 observed-unmodified sites exist only where someone solved a glycoprotein
  structure, so they represent well-studied secreted and membrane proteins rather
  than a random draw from the unlabelled pool.
- 24 pairs cannot support a null result on their own. Any equivalence claim needs
  a smallest effect size of interest fixed before the labelled contrasts are seen.
