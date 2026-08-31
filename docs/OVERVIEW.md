# Overview — where this stands

> **⚠ Correction, 2026-08-25 — ProteinMPNN's sequon indexing.** The manifest's
> `model_index` counts observed residues; ProteinMPNN's parser walks the author
> numbering and inserts a placeholder for every absent number, so ProteinMPNN was
> read at the wrong residue for **25.3%** of sites. Corrected, its conditional
> score on the secretory comparison moves from **+0.090 SD (BH 0.30)** to
> **+0.282 SD (BH 0.031)**, and the claim that *ProteinMPNN does not distinguish
> occupied sequons* does not survive. Retention moves the same way (+0.0423 ->
> +0.0700) but does not clear correction. **ESM-IF and ESMC are unaffected** —
> ESM-IF reproduces +0.431 SD exactly. Any ProteinMPNN number below this line
> predates the fix. Full account:
> [`correction_2026-08-25_sequon_indexing.md`](correction_2026-08-25_sequon_indexing.md).

*Current as of 2026-08-31. This is the one document kept in step with the
results, and it now carries the short account as well. Where it disagrees with
another doc, this one is right; superseded prose lives in
[`archive/`](archive/README.md).*

---

## The question

> Does a protein design model treat an experimentally **occupied** N-X-S/T
> sequon differently from a structurally matched sequon that carries **no
> glycan**?

Two outcomes answer it, and they answer different things:

- **Conditional sequon score** — what probability the model holds at the three
  sequon residues, read off the native structure without generating anything.
- **Retention** — whether the motif survives when the model actually redesigns
  the backbone.

The unit of analysis is one UniProt accession plus one 1-indexed residue. Pair
rows from the ortholog database are kept as context but never as observations,
because one asparagine appearing against a dozen orthologs is still one
biochemical fact — 13,816 pair rows collapse to 4,307 sites.

## Why it matters

A null here is not a disappointment. It establishes the baseline that
glycan-aware modelling has to beat: if current models are indifferent to
occupancy, then any future model that is *not* indifferent has demonstrably
learned something. That is only a meaningful baseline if it holds across models,
which is why there are now seven scored configurations rather than one.

## What is being compared

Matching is built from RSA, neighbour counts and hydrophobic fraction —
**never** from model output — so every model is scored on exactly the same
pairs, and a disagreement between models is a disagreement about the same
comparisons.

| Comparison | Pairs | Role |
|---|---|---|
| `optimal` | 16 | **PRIMARY** — internal controls: no modelled glycan under conditions where glycans *were* modelled elsewhere in the same structure |
| `secretory` | 262 | **Best powered** — eukaryotic, secreted/membrane, no glycoprotein keyword |
| `bacterial` | 280 | Diagnostic |
| `cytosolic` | 273 | Diagnostic |

The diagnostics exist to detect confounding, not to support a claim. See
[`control_sets.md`](control_sets.md).

## The models, and what each conditions on

Six models, seven scored configurations — ESM3 appears twice because it is the
only one that can switch its own structure track off.

| Model | Sees | Conditional | Outcomes |
|---|---|---|---|
| **ProteinMPNN** v_48_020 | backbone + all other native residues | bidirectional, 8 decoding orders averaged | score + retention |
| **ESM-IF1** | backbone + native prefix | autoregressive, single pass | score + retention |
| **CARBonAra** s_v6_4 | backbone + all other native residues | one-shot, single pass per position | score + retention |
| **ESM3-open** (struct) | backbone + sequence, structure track **on** | masked position, single pass | score + retention |
| **ESM3-open** (seq) | sequence only, structure track **off** | masked position, single pass | score |
| **ESMC 300M** | **sequence only** | masked position, single pass | score |
| **ProGen2-base** | **sequence only** | autoregressive, single pass | score |

Laid out by what each conditions on, the grid is now complete in both
directions:

|  | masked / bidirectional | causal | one-shot |
|---|---|---|---|
| **structure + sequence** | ProteinMPNN, ESM3 (struct) | ESM-IF | CARBonAra |
| **sequence only** | ESMC, ESM3 (seq) | ProGen2 | — |

This is a **conditioning spectrum**, not seven attempts at one number. Raw score
magnitudes are not comparable across models; the SD-standardised matched-pair
contrast is, and that is what the analysis rests on. Details:
[`models.md`](models.md).

---

## Current results — scores

Effect is the occupied-minus-control paired difference in SD units; positive
means the model scores occupied sequons higher.

### The three models with the full comparison set

All four comparisons, after the alphabet correction and — for ProteinMPNN — the
sequon-indexing correction. ProteinMPNN's n is lower than the others' because the
indexing guard refuses sites whose two parses cannot be reconciled.

| Comparison | ProteinMPNN | ESM-IF1 | ESMC single | ESMC joint |
|---|---|---|---|---|
| **optimal** | −0.386 (n=15) | −0.337 (n=16) | **+0.792** * | +0.711 |
| **secretory** | **+0.282** *** (n=232) | **+0.431** *** (n=262) | **+0.261** *** | −0.113 |
| bacterial | **−0.528** *** (n=251) | **−0.265** ** | **+0.416** *** | +0.050 |
| cytosolic | −0.002 (n=237) | +0.121 * | −0.132 | **−0.529** *** |

`*` p<0.05, `**` p<0.01, `***` p<0.001 (Wilcoxon, uncorrected). Reference SDs
differ per model (0.95 / 1.51 / 1.55 / 1.37).

**Multiplicity correction spans the confirmatory comparisons only** — internal
control and eukaryotic secretory, across both outcomes. Bacterial and cytosolic
sequons cannot be occupied, so they are reported as diagnostics without
correction. Including them had put four guaranteed-significant tests in the
family, which under Benjamini–Hochberg loosened the threshold for everything
else: ProteinMPNN's secretory score reads BH 0.031 corrected across the four
confirmatory tests, against 0.021 when the diagnostics were included.

Only ProteinMPNN and ESM-IF have a correction family. ESMC's column is shown
here because it has the full set of comparisons, but its stars are uncorrected
Wilcoxon like every other cell, and it has no BH-adjusted p-value anywhere.

### The four later configurations, secretory only

These were run against the best-powered comparison only. They have bootstrap
intervals and Wilcoxon p-values but **no BH-corrected family and no diagnostic
arms**, so they do not sit on the same inferential footing as the table above and
should not be read as though they did.

| Model | Effect (SD) | 95% CI | n | Wilcoxon p | Verdict |
|---|---|---|---|---|---|
| **CARBonAra** | **+0.288** | [+0.130, +0.543] | 262 | 0.0002 | directional, magnitude undetermined |
| **ESM3 structure** | +0.071 | [−0.178, +0.273] | 218 | 0.14 | inconclusive |
| **ESM3 sequence** | −0.060 | [−0.343, +0.150] | 218 | 0.61 | inconclusive |
| **ProGen2** | +0.028 | [−0.217, +0.324] | 261 | 0.49 | inconclusive |

### What this says

**On the best-powered comparison, four of seven configurations find an effect.**
ESM-IF gives +0.431 SD (BH 8e-04) and ProteinMPNN +0.282 SD (BH 0.031) — the two
that sit inside a correction family. CARBonAra (+0.288) and ESMC-single (+0.261)
are significant on the uncorrected Wilcoxon (p = 2e-04 and 4e-04) but have never
entered a BH family, so no corrected p-value exists for either. ESM-IF's,
ProteinMPNN's and CARBonAra's intervals sit entirely above the ±0.2 SD
equivalence margin.

**The two causal sequence-only readings find nothing.** ProGen2 lands on +0.028
with an interval spanning zero, and ESM3 with its structure track off gives
−0.060. That ESMC finds +0.261 from sequence alone while ProGen2 finds nothing
is a difference in *conditional*, not in information available: ESMC is
bidirectional and sees the downstream S/T when scoring the asparagine, where a
causal decoder sees only the prefix.

**The diagnostics have gone incoherent, and that is informative.** Bacterial:
ProteinMPNN −0.528, ESM-IF −0.265, ESMC +0.416 — same sites, opposite signs, all
highly significant. A control set producing large significant effects in
*contradictory directions* across architectures is not measuring occupancy; it is
measuring how each model responds to compositional differences between bacterial
and eukaryotic secretory proteins. That is strong evidence the bacterial and
cytosolic sets cannot support any claim about glycosylation.

### The primary comparison cannot settle anything

`optimal` is 16 pairs over 12 resampling units, and the models disagree in sign
(−0.386, −0.337, +0.792, +0.711). That is what underpowered looks like, not a
contradiction — but it means the primary comparison contributes little, and the
weight sits on `secretory`, which pays with a weaker label (absence of
annotation, not annotated absence) and is contaminated by construction.

## The masking result — and what it settles

Both arms carry the motif, so motif recognition alone cannot separate them.
Hiding the whole sequon asks what the surroundings alone say. Positive change
means the preference **shrinks** when the motif is hidden.

| | motif visible | motif hidden | change | 95% CI | p |
|---|---|---|---|---|---|
| *sequence only* | | | | | |
| ESMC | +0.404 | −0.155 | **+0.558** | [+0.447, +0.724] | 0.0001 |
| ProGen2 | +0.061 | −0.071 | **+0.133** | [+0.089, +0.217] | 0.0001 |
| ESM3 (seq) | −0.064 | −0.052 | −0.012 | [−0.062, +0.050] | 0.65 |
| *structure-conditioned* | | | | | |
| ESM-IF | +0.653 | +0.610 | **+0.043** | [+0.010, +0.073] | 0.016 |
| ProteinMPNN | +0.266 | +0.301 | **−0.035** | [−0.088, −0.008] | 0.008 |
| ESM3 (struct) | +0.085 | +0.118 | −0.033 | [−0.096, +0.010] | 0.12 |

*Log odds. CARBonAra has no row: it is one-shot, so it has no motif-visible arm
to contrast. ESM-IF and ProGen2 are marginalised over the hidden residues rather
than shown an off-distribution `<mask>` token, because a causal decoder cannot
hide an upstream residue from what follows it.*

**This answers the question the 26/08 note left open.** That note recorded ESMC's
collapse and ProteinMPNN's survival and could not say which pattern ESM-IF would
follow — it was listed as the highest-value open experiment. It has now run, and
the split is by conditioning, not by architecture:

- **Every sequence-only reading depends on seeing the motif.** ESMC loses its
  entire preference and ProGen2 loses more than it had. What looked like a
  sequence model recognising occupancy is largely **the motif reinforcing
  itself** — with native N and X visible the model is more confident about S/T at
  occupied sites than at controls, and blind to them the difference is gone.
- **No structure-conditioned reading does.** ESM-IF keeps +0.398 of +0.431.
  ProteinMPNN and ESM3-struct get very slightly *stronger* when blinded. The
  shrinkage in ESM-IF's case is statistically real but roughly a tenth of the
  effect.

So the earlier worry — that ProteinMPNN's and ESM-IF's effects might be the same
artefact as ESMC's — is not supported. Whatever the structure-conditioned models
are responding to, it is not the visible remainder of the motif.

**One caveat bounds the reading, and it is not small.** A structure-conditioned
model **cannot be fully blinded to the sequon**: joint masking hides the other two
residues' identities while the backbone stays, and backbone geometry at a sequon
is itself informative. Occupied sequons really do sit in distinguishable local
geometry — solvent-exposed, often in loops and turns — and matching controls RSA,
neighbour count and hydrophobic fraction only imperfectly. So "survives masking"
means "does not need the neighbouring residue identities". It does **not**
distinguish *the model has learned something about glycosylation* from *occupied
sequons sit on backbone geometry this model can already read*.

## The ESM3 structure contrast

ESM3-open is the only model here that can be run with its structure track on and
off from **one set of weights**, on the same sites. That removes the confound
every cross-model comparison carries, where a difference in conditioning is also
a difference in training data, scale and architecture.

| Arm | structure kept | structure withheld | difference | 95% CI | p |
|---|---|---|---|---|---|
| motif visible | +0.085 | −0.064 | **+0.149** | [+0.037, +0.292] | 0.012 |
| motif hidden | +0.114 | −0.048 | **+0.162** | [+0.056, +0.320] | 0.003 |

**Structure carries occupancy-associated information that sequence alone does
not.** The contrast holds in both masking arms and is slightly *larger* with the
motif hidden — the same direction as the masking table, from within one model.

Two things keep this from being decisive on its own. Neither ESM3 arm reaches
significance as a standalone effect (+0.071 and −0.060, both intervals spanning
zero); it is the *difference between the arms* that is significant, which is a
weaker claim than either arm being individually established. And the caveat above
applies here too: the structure track supplies backbone geometry, so this
identifies structure as the channel without establishing what about the structure
is doing the work.

## Current results — retention

Now run for four models, and it agrees with the direction of the scores.

### Eukaryotic secretory, all four models

| Model | occupied | control | difference | 95% CI | n pairs | Wilcoxon p |
|---|---|---|---|---|---|---|
| **ESM-IF1** | 15.1% | 5.9% | **+0.0925** | [+0.033, +0.152] | 245 | 0.0001 |
| ProteinMPNN | 12.1% | 5.1% | +0.0700 | [−0.004, +0.119] | 216 | 0.0033 |
| CARBonAra | 14.9% | 9.4% | +0.0552 | [−0.006, +0.089] | 245 | 0.028 |
| ESM3 (struct) | 12.3% | 8.8% | +0.0353 | [−0.024, +0.090] | 206 | 0.099 |

**Only ESM-IF's interval excludes zero.** All four point estimates are positive
and they order the same way the scores do, but three of the four are compatible
with no effect. ESMC, ProGen2 and ESM3-seq have no retention row at all: a model
that conditions on sequence rather than a backbone has nothing to redesign.

### The other comparisons, for the two models with a full set

| Comparison | ESM-IF1 | ProteinMPNN |
|---|---|---|
| eukaryotic secretory | **+0.0925** (BH 0.015) | +0.0700 (BH 0.225) |
| internal control | −0.023 (BH 1.00) | −0.081 (BH 0.67) |
| bacterial *(diagnostic)* | −0.054 (BH 0.003) | −0.084 (BH 0.002) |
| cytosolic *(diagnostic)* | +0.022 (BH 0.75) | +0.019 (BH 0.89) |

**ESM-IF clears correction on both outcomes**; ProteinMPNN on the score only. Its
retention moved the same way after the indexing correction (+0.0423 to +0.0700)
without clearing.

The two outcomes are methodologically independent — scoring reads probabilities
off the native backbone without generating anything, retention samples 32
unconstrained designs and counts survivals — so they cannot share an artefact
through the scoring path.

This is also the first test of the earlier lead (occupied sequons retained more
in eukaryotic-secretory context, not at all in the confounded sets). That lead
was generated on ProteinMPNN with the broken alphabet; it now holds clearly for
**ESM-IF**, directionally for ProteinMPNN and CARBonAra, and weakly for ESM3.

One oddity worth recording: CARBonAra's internal-control retention is +0.170 with
an interval excluding zero, the only positive internal-control result anywhere in
the benchmark. On 16 pairs over 12 resampling units, with Wilcoxon p = 0.09, this
is noise until something replicates it.

## What went wrong, and what it cost

Four defects, all found by checking rather than by reasoning:

1. **ProteinMPNN's token alphabet** was a three-letter lookup table from inside
   `parse_PDB_biounits`. `p_asn_at_n` read P(aspartate). Decoding the model's own
   `S` tensor reproduced the native sequence 19.97% of the time with the wrong
   string and 99.53% with the right one. The test asserted the constant against a
   copy of itself, locking the defect in.
   → [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md)
2. **ESM-IF's residue indices** came from biotite; the manifest's came from
   Biopython. They disagreed on ~5% of matched-set sites.
3. **ESMC's token offset** — obvious, and exactly the kind of thing that is
   silently wrong when a tokenizer changes.
4. **A retention run reported success while dropping 21% of its sites** —
   `08_design.py` records failures to a side file and carries on, so a short table
   looked complete.

The first three are one defect: *an assumption about how a model represents its
own input, believed rather than checked.* The fix is the same each time — round
trip the model's representation back to something you already know. All three
checks now run in code. The fourth is a reporting failure, and the shard merge
now refuses to stay quiet about coverage.

## Benchmark state is frozen

`benchmark_frozen/2026-08-23/` holds the analysis outputs and matched pairs
verbatim, plus SHA-256 for every score, design and manifest table.
`pipeline/40_freeze_benchmark.py --check` reports any drift. This exists because
`results/` is gitignored, so roughly twenty hours of ARC and laptop compute lived
only in untracked files — and because stage 10 was found writing
`retention_by_class.json` to a fixed path, so running it for a second model
silently overwrote the first.

→ [`running_on_arc.md`](running_on_arc.md)

## Current results — glyco-site context

The fixed-sequon context-retention test asks whether protecting the motif during
redesign also protects its surroundings. 285 occupied sites, 207 proteins, 32
ProteinMPNN designs each with the sequon held fixed.

| Quantity | Mean | 95% CI | p |
|---|---|---|---|
| design − wild type | **+0.050** | [0.033, 0.067] | 0.0005 |
| random control − wild type | +0.017 | [0.005, 0.031] | 0.006 |
| **design − random** | **+0.033** | [0.020, 0.046] | 0.0005 |

**Protecting the sequon does not protect its environment**, and ProteinMPNN
drifts further from natural occupied context than changing the same number of
residues at random.

A composition control qualifies the mechanism: measuring the same classes at
non-sequon positions of the same designed chains shows the drift is a **global
composition preference, not a local disregard for glycosylation context**. Of the
seven chemical classes, only two show a significant difference between the
flanking window and the rest of the chain — hydrophobics, depleted near the
sequon relative to elsewhere (−0.016, p = 0.026), and cysteine, slightly enriched
(+0.004, p = 0.036). Proline is added across the whole chain and somewhat less
near the sequon; glycine is added slightly *more* near it, but not significantly.
The headline is that ProteinMPNN reshapes composition everywhere, and the sequon
neighbourhood is part of everywhere.

Directly relevant to SugarFix preserve mode: sequon retention as a metric does
not capture what redesign does around a protected site.

### Is a sequon lost faster than anything else?

No. Unconstrained redesign over 318 occupied sites, 220 proteins:

| Quantity | Mean | 95% CI |
|---|---|---|
| Sequon retained, exact triplet | 7.6% | [4.9, 10.4] |
| Control triplet retained, exact | 9.5% | [9.0, 10.0] |
| **control − sequon** | **1.8 pp** | **[−1.0, +4.6]** |

**No detectable excess loss**: ProteinMPNN shows no evidence of selectively
protecting occupied sequons, and their loss rate is broadly consistent with its
overall sequence recovery (background mutation rate 55.7%). Note the interval
reaches +4.6 pp, so this is *not* a demonstration that the rates are equal.

It is also not evidence against glycan blindness — a blind model would be
expected to lose the motif at the ordinary rate, so the two are not alternatives.
What it establishes is that a known biological requirement is treated as ordinary
mutable sequence unless the model is told otherwise.

Validated against the ARC retention run (13.6% on the pattern definition against
13.0% there, correlated 0.993); sequon residues are not unusually vulnerable
(N/S/T average 41.7% retention against 45.5% overall). Sequons are lost in
roughly 92% of designs.

Two caveats carried with this result: the random control was corrected on 26/08
after it was found not to be mutation-count matched (the result moved +0.032 to
+0.033), and disulfide cysteines were never held fixed despite the
pre-specification saying they would be — recorded as a deviation, not correctable
without regenerating the designs. The cysteine enrichment noted above sits in the
same area and is worth remembering when that deviation is eventually closed.

Full accounts:
[`sequon retention rate`](../glyco_context/docs/findings_2026-08-26_sequon_retention_rate.md),
[`context retention`](../glyco_context/docs/findings_2026-08-26_context_retention.md).

## Next, in priority order

1. ~~**Finish retention** for both models.~~ **Done**, 25/08.
2. ~~**Regenerate ProteinMPNN's joint-mask scores.**~~ **Done**, 26/08.
3. ~~**Run ESM-IF joint masking.**~~ **Done**, 31/08. Its preference survives
   (+0.398 of +0.431), which resolves the open question above: motif
   reinforcement explains the sequence-only effects and not the
   structure-conditioned ones.
4. ~~**Regenerate figures and results docs.**~~ **Done**, 31/08.
5. ~~**ESM3-open**, for the sequence-only / structure-conditioned contrast from
   one set of weights.~~ **Done**, 30/08. See the structure contrast above.
6. **Give the four later configurations the full comparison set.** CARBonAra,
   ProGen2 and both ESM3 arms have secretory only. Until they have internal
   control, bacterial and cytosolic arms and enter the BH family, the two results
   tables cannot be merged, and CARBonAra's +0.288 in particular is doing more
   rhetorical work than its inferential status supports.
7. **Separate backbone geometry from learned glycan knowledge.** This is now the
   binding limitation on the headline. Both the masking result and the ESM3
   structure contrast identify structure as the channel, and neither can say
   whether the models learned anything about glycosylation or are reading local
   geometry that happens to correlate with occupancy. Tightening the geometric
   matching, or scoring a geometry-matched set with occupancy shuffled, would
   bound it.
8. **Pre-register before the next analysis.** The secretory result was not
   pre-specified and is the same set that generated the retention lead, so it is
   not independent confirmation of it. Freezing the contrasts and the correction
   family before the next numbers land is the cheapest available guard.

## The comparison this is all building toward

Frozen pretrained model versus the **same** model after glycosylation-specific
training, on identical held-out proteins and matched controls. Comparing
unrelated architectures confounds training data with architecture and scale; the
same model before and after does not. The ESM3 structure contrast is a rehearsal
of exactly this logic — one set of weights, one variable changed — and it is
noticeably cleaner to interpret than any cross-model row in this document.

The benchmark is already shaped for it: everything downstream of scoring is
model-agnostic, and the resampling unit is the ortholog cluster, which is also
the right unit for separating training from evaluation proteins.

The eukaryotic secretory controls must stay **evaluation** controls. They are
unannotated, not annotated-negative, and training on them as negatives would
teach a model curation patterns.

---

## Where to read next

| Document | For |
|---|---|
| [`methodology_explainer.md`](methodology_explainer.md) | **the whole project: data, assumptions, maths, limits** |
| [`models.md`](models.md) | how one model became seven configurations, what each conditions on, how to run them |
| [`glossary.md`](glossary.md) | every term, in plain language |
| [`control_sets.md`](control_sets.md) | what evidence stands behind each control set |
| [`evidence_sources.md`](evidence_sources.md) | what each database can and cannot establish |
| [`figures_and_captions.md`](figures_and_captions.md) | every figure, and what it does not show |
| [`methods_sequon_indexing.md`](methods_sequon_indexing.md) | how a site becomes an index a model can be read at |
| [`running_on_arc.md`](running_on_arc.md) | job arrays, sharding, merging |
| [`why_glyco_site_context_analysis.md`](../glyco_context/docs/why_glyco_site_context_analysis.md) | why the context branch exists, and what was abandoned |
| [`correction_2026-08-25_sequon_indexing.md`](correction_2026-08-25_sequon_indexing.md) | the defect that reversed the ProteinMPNN headline |
| [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) | the alphabet defect in full |
| [`archive/README.md`](archive/README.md) | superseded documents, and why each was retired |
| [`../README.md`](../README.md) | how to run everything |

The two correction notes are dated records that are never updated, but they
explain numbers that are still current, so they are read alongside this document
rather than instead of it. Everything else that went stale is in
[`archive/`](archive/README.md) — including the 23/08 and 26/08 summaries, whose
role this document has taken over.
