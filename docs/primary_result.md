# The primary comparison

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

![primary comparison](../results/figures/primary_optimal.png)

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
| `matched_pairs.csv` | no runner produced it | `runners/match_primary.py` |
| Unmatched-site scoring | depended on a `/tmp` file | `runners/build_candidate_manifest.py` |

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

## Secondary work

Rerun on validated scores; neither is part of the conclusion above.

- **[Diagnostic controls](diagnostic_controls.md)** — bacterial (278 pairs,
  −0.174 SD) and cytosolic (270 pairs, +0.062 SD). The three comparisons point
  in different directions, so the earlier "gradient" reading is withdrawn.
- **Design retention** — frozen snapshot, scoreable sites only: mean full
  retention 0.078, 80.4% of sites lose the sequon in all 32 designs, conditional
  score predicts retention at Spearman +0.559 with monotonic quintiles. Coverage
  is alphabetical, not a random sample. `results/retention_analysis.json`.

The earlier retention contrast (+0.084, CI [+0.007, +0.176]) is **withdrawn**: it
used a site-level bootstrap and its pairs no longer exist after rematching.

## Reproducing

```
runners/build_candidate_manifest.py dataset
runners/scoreability.py results/candidate_manifest_dataset.csv results/scoreability_dataset.csv
runners/match_primary.py                     # optimal + both greedy sensitivities
runners/score_all.py results/candidate_manifest_dataset.csv results/scores_dataset.csv
runners/primary_analysis.py optimal          # PRIMARY
runners/matching_sensitivity.py              # 200-seed sweep
runners/primary_plot.py optimal
```

Artefacts: `results/analysis_optimal.json`, `results/contrasts_optimal.csv`,
`results/matched_pairs_optimal.csv`, `results/matching_balance_optimal.json`,
`results/matching_sensitivity.json`, `results/matching_seed_sweep.csv`,
`results/scores_dataset.csv`, `results/figures/primary_optimal.png`,
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
