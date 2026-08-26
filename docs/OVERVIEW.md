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
[`negative_controls.md`](negative_controls.md).

## The three models, and what each conditions on

| Model | Sees | Conditional | Outcomes |
|---|---|---|---|
| **ProteinMPNN** v_48_020 | backbone + all other native residues | bidirectional, 8 decoding orders averaged | score + retention |
| **ESM-IF1** | backbone + native prefix | autoregressive, single pass | score + retention |
| **ESMC 300M** | **sequence only** | masked position, single pass | score |

This is a **conditioning spectrum**, not three attempts at one number. Raw score
magnitudes are not comparable across models; the SD-standardised matched-pair
contrast is, and that is what the analysis rests on. Details:
[`second_model_esm_if.md`](second_model_esm_if.md),
[`third_model_esmc.md`](third_model_esmc.md),
[`adding_models_explainer.md`](adding_models_explainer.md).

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

### What this says

**On the best-powered comparison, three of three models find a real effect.**
ESM-IF gives +0.431 SD (BH p = 2e-06), ESMC-single +0.261 SD (BH p = 8e-04), and
ProteinMPNN — once its sequon indexing is corrected — +0.282 SD (BH 0.021),
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
| ProteinMPNN (structure) | +0.022 | +0.029 | −0.007 [−0.016, −0.002] |

The sequence-only model's discrimination depends entirely on seeing the motif in
context; the structure-conditioned model's does not. ESM-IF has no joint variant
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

## Next, in priority order

1. **Finish retention** for both models, merge, and rerun stages 10, 10b, 11.
   Until this lands the picture is half-complete and the live lead is untested.
2. **Add joint masking to ProteinMPNN.** This is the cheapest high-value
   experiment available: it tests whether the structure-conditioned effects are
   the same motif-self-reinforcement artefact that ESMC's sensitivity exposed. If
   they survive joint masking and ESMC's does not, that is a real
   structure-specific finding. If they collapse too, the effect is a property of
   sequon context.
3. **Regenerate figures and the results docs** from corrected numbers, and retire
   the staleness banners.
4. **ESM3-open**, for sequence-only / structure-only / sequence+structure from one
   set of weights — the cleanest possible conditioning contrast. Gated behind a
   licence and a token; leave the function and annotation tracks empty or they can
   leak the label.
5. **Pre-register before the next analysis.** The secretory result was not
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
| [`why_glyco_site_context_analysis.md`](../glyco_context/docs/why_glyco_site_context_analysis.md) | **why this work exists, and what the benchmark could not answer** |
| [`adding_models_explainer.md`](adding_models_explainer.md) | how one model became three, and what broke |
| [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) | the alphabet defect in full |
| [`second_model_esm_if.md`](second_model_esm_if.md) | ESM-IF's conditional, index mapping, batched generation |
| [`third_model_esmc.md`](third_model_esmc.md) | ESMC, the two masking schemes, the environment split |
| [`negative_controls.md`](negative_controls.md) | what evidence stands behind each control set |
| [`glossary.md`](glossary.md) | the terms, in plain language |
| [`rationale_and_progress.md`](rationale_and_progress.md) | why the module exists |
| [`running_on_arc.md`](running_on_arc.md) | job arrays, sharding, merging |
| `../README.md` | how to run everything |

Docs carrying a staleness banner have sound arguments and superseded numbers.
