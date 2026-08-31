> **⚠ SUPERSEDED.** The numbers here predate the 2026-08-25 sequon-indexing
> correction, and the ProteinMPNN figures also predate the 2026-08-20 alphabet
> correction. The reasoning and method stand; the quantities do not. Current
> results live in [`../OVERVIEW.md`](../OVERVIEW.md), which is the single place
> they are maintained. Kept as a dated record of what was concluded at the time.

# The primary comparison


> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](../correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](../OVERVIEW.md) for current numbers.

> Given the native protein backbone and the surrounding native sequence, does
> ProteinMPNN assign a higher conditional sequon score to experimentally
> occupied sites than to structurally matched sites with no modelled glycan?

One model, one score, one comparison. Everything else is secondary and reported
separately. **The ProteinMPNN analysis is frozen as of 18 August 2026.**

## Conclusion

**Occupied sites tend to score higher than matched sites with no modelled
glycan, but 16 pairs do not establish a precise or statistically robust
difference.**

Both halves of that sentence carry weight, and the evidence for each is
different in kind.

*Tend to score higher* is the stable part. Every point estimate computed under
every matching we tried is positive — the deterministic optimum, the earlier
greedy seed, and all 200 seeds of a greedy sweep, spanning +0.286 to +0.699 SD.
Nothing we varied produced a negative estimate.

*Not precise or statistically robust* is equally well supported. The primary
interval includes zero. Occupied sites score higher in only 9 of 16 pairs. And
across the 200-seed sweep the interval excludes zero in just 38% of cases, so
whether the result reaches significance is decided by an arbitrary choice inside
the matching algorithm rather than by the data.

## Result

![primary comparison](../../results/figures/primary_optimal.png)

Matching is deterministic: the assignment maximising the number of admissible
pairs and, among those, minimising total distance. No seed, no ordering.

| | Value |
|---|---|
| Matched pairs | 16 |
| Occupied proteins / ortholog clusters / resampling units | 16 / 13 / 12 |
| Mean paired difference | **+0.610 log-odds (+0.458 SD)** |
| 95% CI, cluster-aware | [−0.302, +1.463] → **[−0.227, +1.098] SD** |
| Equivalence margin | ±0.266 log-odds (±0.2 SD) |
| Verdict | **inconclusive** |

The interval spans more than five times the equivalence margin, and it is not
symmetric: it reaches to +1.098 SD but only to −0.227 SD. So it is consistent
with no difference, with a difference well beyond the margin in the positive
direction, and only with small negative differences — a large effect favouring
the controls is the one thing these 16 pairs do exclude.

### Robustness

| | Value |
|---|---|
| Occupied scores higher in | 9 of 16 pairs |
| Sign test | p = 0.80 |
| Wilcoxon signed-rank | p = 0.19 |
| Median | +0.777 log-odds |
| 10% trimmed mean | +0.618 |
| Leave-one-out mean | stays within [+0.377, +0.851] |

The mean is positive because the negative contrasts are small and the positive
ones large, not because most pairs point the same way. The effect lives in
magnitude, not consistency. It is not an outlier artefact — no single pair can
be dropped to make the mean negative — but nor is it a difference that a
direction-based test can detect.

### Matching sensitivity

The reason matching is now deterministic.

| Matching | Contrasts | Mean | 95% CI (SD) | Excludes zero |
|---|---|---|---|---|
| **Optimal (primary)** | 16 | +0.458 SD | [−0.227, +1.098] | no |
| Greedy, seed 0 | 16 | +0.649 SD | [+0.075, +1.243] | yes |
| Greedy, k=5 | 14 cases (16 control links) | +0.570 SD | [+0.051, +1.165] | yes |
| Greedy, 200 seeds | 16 | +0.286 to +0.699 SD | — | 75/200 (38%) |

At k=5 the 16 control links collapse to 14 contrasts, because the unit of
analysis is the occupied site and two cases absorb more than one control each.
Under 1:1 matching the pair count is 16 for every seed, so the **caliper**
limits it, not the algorithm: 12 of the 28 controls have no admissible partner within 0.25. The
optimal matching finds tighter pairs than greedy (mean distance 0.099 vs 0.113)
and better balance, but not more of them.

The seed-0 result was one draw from this distribution reported as a value. That
is why it is now a sensitivity analysis rather than the headline.

### Sequon subtype

Matched exactly, so like-for-like.

| Subtype | n | Mean |
|---|---|---|
| NXS | 7 | +0.721 log-odds (+0.542 SD) |
| NXT | 9 | +0.524 log-odds (+0.394 SD) |

Both positive; at these counts the difference between them is not interpretable.

## What the controls are

The internal-control class is **sequons with no modelled glycan under
internal-control conditions**: sites in solved glycoprotein structures where a
glycan would have been visible had one been present. They are the only controls
holding organism, compartment and experimental context roughly constant, and
there are 28 that can be scored.

They are **not** demonstrated to be chemically unmodified. Absence of a modelled
glycan is weaker evidence than absence of glycosylation.

## Reference scale

| | Value |
|---|---|
| Scoreable dataset sites | 342 (314 occupied, 28 internal control) |
| Reference SD | **1.3316** log-odds |
| Equivalence margin (±0.2 SD) | ±0.2663 log-odds |

Pooled across dataset sites only, without consulting labels. Control sites are
excluded by construction, so the scale cannot move when a control pool is rebuilt.

## Corrections behind this result

Two amendments, both recorded in `config/scoring_frozen.toml`.

| | Before | Now |
|---|---|---|
| Invalid probability vectors | 105 sites scored at ~+13.8 | rejected at two independent guards |
| Reference SD | 2.6169 (inflated by those rows) | 1.3316 |
| Scoreability | established after matching | before, from coordinates alone |
| Sequon subtype | ~45% of pairs matched NXS to NXT | matched exactly |
| Control-protein reuse | ignored in the bootstrap | resampling unit is the connected component |
| Matching | greedy, seed-dependent | deterministic optimum |
| `matched_pairs.csv` | no runner produced it | `pipeline/06_match_primary.py` |
| Unmatched-site scoring | depended on a `/tmp` file | `pipeline/04_build_candidate_manifest.py` |

The originally reported −0.057 SD is withdrawn: it was computed against the
inflated reference SD, from scores including invalid rows.

## Limitations

- **16 pairs.** The binding constraint, and correcting the analysis tightened it.
  Growing the internal-control class — realistically through PNGase F / H₂¹⁸O
  occupancy glycoproteomics — is what would make this decisive.
- The caliper discards 12 of 28 controls. Relaxing it would trade matching
  quality for pairs; not attempted here.
- Matching holds local structure and subtype constant, **not** organism.
- Internal controls are sites with no modelled glycan, not proven negatives.
- The margin is an exploratory statistical threshold, not a biological one.
- ProteinMPNN parses only ATOM records, so glycans are invisible to it. Any
  effect is an occupancy-associated statistical preference, not evidence that
  the model represents glycosylation mechanistically.

## The parallel comparison — better powered, same answer

Because the internal-control class is not being grown, a fourth control set was
built: eukaryotic, secreted or membrane, with a solved structure and **no**
glycoprotein annotation. It matches the occupied sites on taxonomy *and*
compartment — neither confound the diagnostics carry — at the cost of a weaker
negative label, since absence of annotation is not annotated absence.

| | Value |
|---|---|
| Matched pairs | **262** |
| Mean paired difference | **+0.097 log-odds (+0.073 SD)** |
| 95% CI, cluster-aware | **[−0.056, +0.346] SD** |
| Verdict | inconclusive |

The point estimate is *inside* the equivalence margin; only the upper bound
escapes it. This is the narrowest interval in the study — four times tighter than
the primary — and it sits essentially on zero. Roughly twice the pairs would
likely settle the question.

Contamination is real and quantified: about half of eukaryotic secretory proteins
with a structure carry a glycoprotein keyword, so the unannotated half certainly
holds unrecorded sites. A name audit found 32 of 1,543 control proteins suspect,
8 of which reach the matched pairs; removing them changes nothing (+0.073 →
+0.074 SD). See `pipeline/13_name_audit.py` and
[`negative_controls.md`](../control_sets.md).

## Significance

Eight tests — four control sets × two outcomes — with a cluster-level permutation
test rather than Wilcoxon, because the pairs are not independent.
**No test survives correction**: smallest raw p 0.030, smallest Holm 0.237,
smallest Benjamini–Hochberg 0.120. Full detail in
[`significance.md`](significance_SUPERSEDED_2026-08-25.md).

The pre-specified inference was never a p-value but the equivalence assessment
above, which is unaffected.

## Secondary work

Rerun on validated scores; neither is part of the conclusion above.

- **[Diagnostic controls](../control_sets.md)** — bacterial (278 pairs,
  −0.174 SD) and cytosolic (270 pairs, +0.062 SD). The three comparisons point
  in different directions, so the earlier "gradient" reading is withdrawn.
- **Design retention** — sweep complete, scoreable sites only (2,423 of 2,526):
  mean full retention 0.072, **81.6%** of sites lose the sequon in all 32 designs,
  and the conditional score predicts retention at Spearman **+0.547** with
  monotonic quintiles (0.000 in the lowest, 0.301 in the highest). 114 of 2,640
  manifest sites are uncovered, all from logged chain parse failures, skewed
  towards large recent depositions. `results/retention_analysis.json`.

The earlier retention contrast (+0.084, CI [+0.007, +0.176]) is **withdrawn**: it
used a site-level bootstrap and its pairs no longer exist after rematching.

## Reproducing

```
pipeline/04_build_candidate_manifest.py dataset
pipeline/05_scoreability.py results/manifests/candidate_manifest_dataset.csv results/manifests/scoreability_dataset.csv
pipeline/06_match_primary.py                     # optimal + both greedy sensitivities
pipeline/07_score.py results/manifests/candidate_manifest_dataset.csv results/scores/scores_dataset.csv
pipeline/09_analyse_scores.py optimal          # PRIMARY
pipeline/12_matching_sensitivity.py              # 200-seed sweep
pipeline/23_figures_primary.py optimal
```

Artefacts: `results/analysis/analysis_optimal.json`, `results/analysis/contrasts_optimal.csv`,
`results/matching/matched_pairs_optimal.csv`, `results/matching/matching_balance_optimal.json`,
`results/analysis/matching_sensitivity.json`, `results/analysis/matching_seed_sweep.csv`,
`results/scores/scores_dataset.csv`, `results/figures/primary_optimal.png`,
`config/scoring_frozen.toml` (`[amendment_1]`, `[amendment_2]`).

## Still open

- The retention contrast has not been recomputed on the current 16 pairs. The
  earlier one is withdrawn, not replaced.
- Caliper sensitivity not run. It is the one analysis that could change the pair
  count, since 12 of the 28 scoreable controls fail the caliper rather than any
  other criterion.
- No second model. The ProteinMPNN analysis is frozen; nothing further runs
  until the pair count improves.

`docs/phase4_primary_result.md` was superseded by this document and has been
moved to `docs/archive/phase4_primary_result_SUPERSEDED.md`.
`docs/methodology_explainer.md` has been rewritten against the current result.
