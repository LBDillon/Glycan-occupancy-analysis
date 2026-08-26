# Overview — where this stands

> **⚠ Correction, 2026-08-25 — ProteinMPNN's sequon indexing.** The manifest's
> `model_index` counts observed residues; ProteinMPNN's parser walks the author
> numbering and inserts a placeholder for every absent number, so ProteinMPNN was
> read at the wrong residue for **25.3%** of sites. Corrected, its conditional
> score on the secretory comparison moves from **+0.090 SD (BH 0.30)** to
> **+0.282 SD (BH 0.021)**, and the claim that *ProteinMPNN does not distinguish
> occupied sequons* does not survive. Retention moves the same way (+0.0423 ->
> +0.0700) but does not clear correction. **ESM-IF and ESMC are unaffected** —
> ESM-IF reproduces +0.431 SD exactly. Any ProteinMPNN number below this line
> predates the fix. Full account:
> [`correction_2026-08-25_sequon_indexing.md`](correction_2026-08-25_sequon_indexing.md).

*Current as of 2026-08-23. This is the one document kept in step with the
results. Where it disagrees with another doc, this one is right and the other is
flagged with a staleness banner.*

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
which is why there are now three.

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
[`negative_controls.md`](control_sets.md).

## The three models, and what each conditions on

| Model | Sees | Conditional | Outcomes |
|---|---|---|---|
| **ProteinMPNN** v_48_020 | backbone + all other native residues | bidirectional, 8 decoding orders averaged | score + retention |
| **ESM-IF1** | backbone + native prefix | autoregressive, single pass | score + retention |
| **ESMC 300M** | **sequence only** | masked position, single pass | score |

This is a **conditioning spectrum**, not three attempts at one number. Raw score
magnitudes are not comparable across models; the SD-standardised matched-pair
contrast is, and that is what the analysis rests on. Details:
[`second_model_esm_if.md`](models.md),
[`third_model_esmc.md`](models.md),
[`adding_models_explainer.md`](models.md).

---

## Current results — scores

All four comparisons, all models, after the alphabet correction and — for
ProteinMPNN — the sequon-indexing correction. Effect is the occupied-minus-control
paired difference in SD units; positive means the model scores occupied sequons
higher. ProteinMPNN's n is lower than the other models' because the indexing
guard refuses sites whose two parses cannot be reconciled.

| Comparison | ProteinMPNN | ESM-IF1 | ESMC single | ESMC joint |
|---|---|---|---|---|
| **optimal** | −0.386 (n=15) | −0.337 (n=16) | **+0.792** * | +0.711 |
| **secretory** | **+0.282** * (n=232) | **+0.431** *** (n=262) | **+0.261** *** | −0.113 |
| bacterial | **−0.528** *** (n=251) | **−0.265** ** | **+0.416** *** | +0.050 |
| cytosolic | −0.002 (n=237) | +0.121 * | −0.132 | **−0.529** *** |

`*` p<0.05, `**` p<0.01, `***` p<0.001 (Wilcoxon, uncorrected). Reference SDs
differ per model (1.36 / 1.51 / 1.55 / 1.37).

**Multiplicity correction spans the confirmatory comparisons only** — internal
control and eukaryotic secretory, across both outcomes. Bacterial and cytosolic
sequons cannot be occupied, so they are reported as diagnostics without
correction. Including them had put four guaranteed-significant tests in the
family, which under Benjamini–Hochberg loosened the threshold for everything
else: ProteinMPNN's secretory score reads BH 0.031 corrected across the four
confirmatory tests, against 0.021 when the diagnostics were included.

### What this says

**On the best-powered comparison, three of three models find a real effect.**
ESM-IF gives +0.431 SD (BH p = 2e-06), ESMC-single +0.261 SD (BH p = 8e-04), and
ProteinMPNN — once its sequon indexing is corrected — +0.282 SD (BH 0.031),
both with intervals entirely above the ±0.2 SD equivalence margin. ProteinMPNN
does not resolve it (+0.090, inconclusive).

**But sequence alone reproduces most of it.** ESMC sees no structure and still
finds +0.261 SD. So the honest reading of ESM-IF's result is not "a
structure-conditioned model recognises occupancy" but "occupied sequons sit in
distinguishable **sequence** contexts, and a structure-conditioned model sees
that too."

**And the joint-masking sensitivity mostly dissolves it.** With all three sequon
positions masked so the model cannot see the rest of the motif, ESMC's secretory
effect collapses from +0.261 to −0.113 (p = 0.67), and its bacterial effect from
+0.416 to +0.050. The sequence-only signal is therefore **largely the motif
reinforcing itself** — when the model can see native N and X it is more confident
about S/T at occupied sites than at controls, and when it cannot, the difference
is gone.

**This is the most important open question.** ProteinMPNN and ESM-IF both
condition on the rest of the native sequon too, so their effects may be the same
artefact. Neither has a joint-masking variant yet. ProteinMPNN could gain one
cheaply — `conditional_probs` already accepts a position set.

**The diagnostics have gone incoherent, and that is informative.** Bacterial:
ProteinMPNN −0.396, ESM-IF −0.265, ESMC +0.416 — same sites, opposite signs, all
highly significant. A control set producing large significant effects in
*contradictory directions* across architectures is not measuring occupancy; it is
measuring how each model responds to compositional differences between bacterial
and eukaryotic secretory proteins. That is strong evidence the bacterial and
cytosolic sets cannot support any claim about glycosylation.

### The primary comparison cannot settle anything

`optimal` is 16 pairs over 12 resampling units, and the models disagree in sign
(+0.640, −0.337, +0.792, +0.711). That is what underpowered looks like, not a
contradiction — but it means the primary comparison contributes little, and the
weight sits on `secretory`, which pays with a weaker label (absence of
annotation, not annotated absence) and is contaminated by construction.

## Current results — retention

Now complete for both models, and it agrees with the scores.

| Comparison | ESM-IF1 | ProteinMPNN |
|---|---|---|
| eukaryotic secretory | **+0.0925** (BH 0.015) | +0.0700 (BH 0.225) |
| internal control | −0.023 (BH 1.00) | −0.045 (BH 0.87) |
| bacterial *(diagnostic)* | −0.054 (BH 0.003) | −0.070 (BH 0.002) |
| cytosolic *(diagnostic)* | +0.022 (BH 0.75) | −0.003 (BH 0.89) |

**ESM-IF clears correction on both outcomes**; ProteinMPNN on the score only.
Its retention moved the same way after the indexing correction (+0.0423 to
+0.0700) without clearing. Retention
15.1% vs 5.9% on the secretory comparison, 93 informative pairs of 245.

The two outcomes are methodologically independent — scoring reads probabilities
off the native backbone without generating anything, retention samples 32
unconstrained designs and counts survivals — so they cannot share an artefact
through the scoring path.

This is also the first test of the earlier lead (occupied sequons retained more
in eukaryotic-secretory context, not at all in the confounded sets). That lead
was generated on ProteinMPNN with the broken alphabet; it now **holds for ESM-IF
and not for ProteinMPNN**.

## The masking result

Both arms carry the motif, so motif recognition alone cannot separate them.
Hiding the whole sequon asks what the surroundings alone say.

| | gap, motif visible | gap, motif hidden | change |
|---|---|---|---|
| ESMC (sequence only) | +0.121 | +0.007 | **+0.113** [+0.071, +0.147] |
| ProteinMPNN (structure) | +0.022 ⚠ | +0.029 ⚠ | −0.007 ⚠ |

> **⚠ The ProteinMPNN row is stale and must not be interpreted.** Its joint-mask
> scores were computed before the 2026-08-25 sequon-indexing correction, which
> moved that model's unmasked secretory result from +0.090 to +0.282 SD. The
> masked variant has not been regenerated, so the comparison between the two
> columns is not currently meaningful for ProteinMPNN. ESMC is unaffected — its
> adapter validated the triplet and never used the raw index.

The sequence-only model's discrimination depends entirely on seeing the motif in
context. Whether the structure-conditioned model's does is currently unknown,
pending that regeneration. ESM-IF has no joint variant
— `<mask>` is off-distribution in its decoder prefix (it sends 93% of the
probability mass onto aromatics), so the honest version needs marginalisation
rather than substitution. Not built.

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
| random control − wild type | +0.018 | [0.006, 0.031] | 0.002 |
| **design − random** | **+0.033** | [0.020, 0.046] | 0.0005 |

**Protecting the sequon does not protect its environment**, and ProteinMPNN
drifts further from natural occupied context than changing the same number of
residues at random.

A composition control qualifies the mechanism: measuring the same classes at
non-sequon positions of the same designed chains shows ProteinMPNN adds proline
and glycine and removes aromatics **across the whole chain**, slightly less near
the sequon than elsewhere. The drift is a global composition preference, not a
local disregard for glycosylation context.

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
overall sequence recovery. Note the interval reaches +4.6 pp, so this is *not* a
demonstration that the rates are equal.

It is also not evidence against glycan blindness — a blind model would be
expected to lose the motif at the ordinary rate, so the two are not alternatives.
What it establishes is that a known biological requirement is treated as ordinary
mutable sequence unless the model is told otherwise.

Validated against the ARC retention run (13.0% against 13.6% on the pattern
definition, correlated 0.993); sequon residues are not unusually vulnerable
(N/S/T average 41.7% retention against 45.5% overall). Sequons are lost in
roughly 92% of designs.

Full account:
[`../glyco_context/docs/findings_2026-08-26_sequon_retention_rate.md`](../glyco_context/docs/findings_2026-08-26_sequon_retention_rate.md).

Two caveats carried with this result: the random control was corrected on 26/08
after it was found not to be mutation-count matched (the result moved +0.032 to
+0.033), and disulfide cysteines were never held fixed despite the
pre-specification saying they would be — recorded as a deviation, not correctable
without regenerating the designs.

Full account, including what the result does not license:
[`../glyco_context/docs/findings_2026-08-26_context_retention.md`](../glyco_context/docs/findings_2026-08-26_context_retention.md).

## Next, in priority order

1. ~~**Finish retention** for both models.~~ **Done**, 25/08 — merged from ARC
   and rerun through stages 10, 10b and 11 with the corrected indexing.
2. **Regenerate ProteinMPNN's joint-mask scores**, which predate the
   sequon-indexing correction and are the last stale numbers in this document
   (see the masking section, where that row is flagged).
3. **Add joint masking to ProteinMPNN.** This is the cheapest high-value
   experiment available: it tests whether the structure-conditioned effects are
   the same motif-self-reinforcement artefact that ESMC's sensitivity exposed. If
   they survive joint masking and ESMC's does not, that is a real
   structure-specific finding. If they collapse too, the effect is a property of
   sequon context.
4. ~~**Regenerate figures and results docs.**~~ **Done**, 26/08 for everything
   except the joint-mask row above.
5. **ESM3-open**, for sequence-only / structure-only / sequence+structure from one
   set of weights — the cleanest possible conditioning contrast. Gated behind a
   licence and a token; leave the function and annotation tracks empty or they can
   leak the label.
6. **Pre-register before the next analysis.** The secretory result was not
   pre-specified and is the same set that generated the retention lead, so it is
   not independent confirmation of it. Freezing the contrasts and the correction
   family before the ARC numbers land is the cheapest available guard.

## The comparison this is all building toward

Frozen pretrained model versus the **same** model after glycosylation-specific
training, on identical held-out proteins and matched controls. Comparing
unrelated architectures confounds training data with architecture and scale; the
same model before and after does not. The benchmark is already shaped for it —
everything downstream of scoring is model-agnostic, and the resampling unit is
the ortholog cluster, which is also the right unit for separating training from
evaluation proteins.

The eukaryotic secretory controls must stay **evaluation** controls. They are
unannotated, not annotated-negative, and training on them as negatives would
teach a model curation patterns.

---

## Where to read next

| Document | For |
|---|---|
| [`summary_2026-08-26.md`](summary_2026-08-26.md) | **the short account, current** |
| [`methodology_explainer.md`](methodology_explainer.md) | **the whole project: data, assumptions, maths, limits** |
| [`why_glyco_site_context_analysis.md`](../glyco_context/docs/why_glyco_site_context_analysis.md) | why the context branch exists, and what was abandoned |
| [`glossary.md`](glossary.md) | every term, in plain language |
| [`models.md`](models.md) | how one model became three, what each conditions on, how to run them |
| [`control_sets.md`](control_sets.md) | what evidence stands behind each control set |
| [`figures_and_captions.md`](figures_and_captions.md) | every figure, and what it does not show |
| [`methods_sequon_indexing.md`](methods_sequon_indexing.md) | how a site becomes an index a model can be read at |
| [`correction_2026-08-25_sequon_indexing.md`](correction_2026-08-25_sequon_indexing.md) | the defect that reversed the ProteinMPNN headline |
| [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) | the alphabet defect in full |

Dated records — never updated, kept as accounts of a moment:
`summary_2026-08-23.md`, the four `correction_*` notes,
`walkthrough_annotated_2026-08-19.md`, and `archive/`.
| [`running_on_arc.md`](running_on_arc.md) | job arrays, sharding, merging |
| `../README.md` | how to run everything |

Docs carrying a staleness banner have sound arguments and superseded numbers.
