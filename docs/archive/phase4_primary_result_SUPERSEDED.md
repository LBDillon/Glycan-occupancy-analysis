> ## ⚠️ SUPERSEDED — the numbers in this document are withdrawn
>
> Every ProteinMPNN score quoted below was computed with a scorer that returned
> twenty-one ones for any residue with an incomplete backbone, giving
> P(S)+P(T) = 2 and a score near +13.8. This corrupted 105 of 2,564 sites and
> inflated the reference SD from 1.33 to 2.62, so all standardised effects here
> are wrong — including the primary estimate, whose sign reverses once the
> invalid rows are removed.
>
> **Current result: [`primary_result.md`](primary_result_SUPERSEDED_2026-08-25.md).**
> Correction record: `config/scoring_frozen.toml`, section `[amendment_1]`.

# Phase 4 — occupied versus unoccupied sequons under ProteinMPNN

Zero-shot conditional scoring, ProteinMPNN v_48_020, against the frozen
configuration in `config/scoring_frozen.toml`. The configuration, including the
equivalence margin, was written and committed before any labelled contrast was
computed.

## The question

Does ProteinMPNN assign a higher site-level probability to an intact N-X-S/T
sequon when that site is experimentally occupied?

## Reference scale

Estimated across all 350 structurally scoreable dataset sites, pooled, without
consulting their labels, as the frozen rule requires.

| | Value |
|---|---|
| Reference SD | 2.6169 log-odds |
| Equivalence margin (±0.2 SD) | ±0.5234 log-odds |

The margin is an **exploratory statistical threshold, not a biologically
validated one**. No claim is made that a 0.2 SD difference is or is not
biologically meaningful.

## Primary result — inconclusive

Occupied sites against observed-unmodified sites: the comparison matched on
organism, compartment and experimental context.

| | Value |
|---|---|
| Sites / proteins / ortholog clusters | 22 / 22 / 19 |
| Mean contrast | −0.148 log-odds |
| 95% CI (cluster bootstrap) | [−2.808, +2.395] |
| Standardised effect | −0.057 SD |
| 95% CI, standardised | [−1.073, +0.915] |
| P(occupied scores higher) | 0.545 |

**Verdict: inconclusive.** The interval spans roughly ±1 SD, about five times
the width of the equivalence margin. It contains zero, and it contains effects
far larger than the margin in both directions.

This does **not** show that ProteinMPNN scores occupied and unoccupied sequons
alike. It shows that 22 matched pairs cannot distinguish those possibilities.
The point estimate sits near zero, but a point estimate without a usable
interval is not a result.

## Secondary comparisons

Diagnostic only. Each is confounded by construction — the cytosolic set differs
in subcellular compartment, the bacterial set in taxonomy — so they are read
against one another, not as substitutes for the primary comparison.

| Comparison | Sites | Standardised | 95% CI (SD) | P(occupied higher) |
|---|---|---|---|---|
| vs bacterial extracytoplasmic | 269 | −0.145 | [−0.267, −0.021] | 0.461 |
| vs cytosolic eukaryotic | 262 | −0.237 | [−0.383, −0.084] | 0.450 |

Both intervals exclude zero but straddle the margin, so neither establishes
equivalence nor a difference beyond the margin. They are too imprecise to call
either way, which is a different statement from either finding.

## What the three comparisons say together

The effect shrinks as the comparison gets better matched:

| Comparison | Differs from the occupied set in | Standardised effect |
|---|---|---|
| Cytosolic | subcellular compartment | −0.237 SD |
| Bacterial | taxonomy and fold repertoire | −0.145 SD |
| Observed-unmodified | neither — same organism, compartment, experiment | −0.057 SD |

That gradient is what the orthogonal-confound design was built to detect. It is
consistent with the apparent differences in the control comparisons arising from
compartment and taxonomy rather than from occupancy, since the residual effect is
smallest precisely where those two confounds are removed. It is consistent with,
not proof of: the best-matched comparison is also the least precise, so the
ordering could equally reflect its wide interval.

Every point estimate is negative. Across all three comparisons ProteinMPNN
assigns slightly **lower** conditional probability to the motif-forming residues
at occupied sites than at matched unoccupied ones. Nothing here supports the
model preferring occupied sequons.

## Sequon subtype

| Comparison | NXS | NXT |
|---|---|---|
| vs bacterial | −0.111 SD (n=129) | −0.175 SD (n=140) |
| vs cytosolic | −0.280 SD (n=125) | −0.197 SD (n=137) |
| vs observed-unmodified | −0.339 SD (n=8) | +0.105 SD (n=14) |

The subtypes point in opposite directions in the primary comparison, but at
n=8 and n=14 that is noise, recorded for completeness rather than interpreted.

## How to describe this result

Any effect here is an **occupancy-associated statistical preference**, not
evidence that the model represents glycosylation mechanistically. Equally, the
inconclusive primary result is not evidence that ProteinMPNN contains no
glycosylation-related information — the study lacks the precision to say.

There is a further reason to expect little: ProteinMPNN parses only ATOM
records, so glycans are HETATM and structurally invisible to it. The model has
never seen a glycan in training or at inference. A null is close to the prior
expectation, which is what makes this a negative control for glycan-aware
modelling rather than an open question resolved.

## Limitations

- **22 pairs.** The binding constraint. Growing the observed-unmodified class —
  realistically through PNGase F / H₂¹⁸O occupancy glycoproteomics — is what
  would make the primary comparison decisive.
- Control comparisons carry the confounds they were designed to carry.
- 114 of 2,640 sites (4.3%) could not be scored because ProteinMPNN reads only
  PDB-format files and their structures are mmCIF-only. This did not touch the
  primary comparison: occupied and observed-unmodified sites were 100% scored,
  and only controls lost ~5%.
- The observed-unmodified sites exist only where a glycoprotein structure was
  solved, so they are well-studied secreted and membrane proteins rather than a
  random draw.
- Conditional scores depend on decoding order; 8 seeded orders were used after a
  blinded convergence check (median 8-vs-16 difference 0.0022 SD).

## Artefacts

`results/mpnn_conditional_scores.csv`, `results/mpnn_conditional_scores_unmatched.csv`,
`results/phase4_primary_analysis.json`, `results/contrasts_*.csv`,
`results/convergence_check.json`, `config/scoring_frozen.toml`.
