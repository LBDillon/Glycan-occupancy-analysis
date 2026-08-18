# The primary comparison

> Given the native protein backbone and the surrounding native sequence, does
> ProteinMPNN assign a higher conditional sequon score to experimentally
> occupied sites than to structurally matched sites with no modelled glycan?

One model, one score, one comparison. Everything else in this module is
secondary to this question and is reported separately.

**Status: provisional.** The estimate below rests on 16 matched pairs and
replaces an earlier result computed with a defective scorer. It should not be
quoted as settled.

## What changed, and why the earlier number is void

An external review found that ProteinMPNN's `conditional_probs` returns a row
of zeros for any residue whose backbone is incomplete — a missing N, CA, C or
O. Exponentiating that row gives twenty-one ones: P(N) = 1, P(S)+P(T) = 2, and
a score near +13.8. The manifest had confirmed that all three sequon residues
had coordinates, which is a different claim from the model having accepted
them, and nothing checked that the returned vectors were probability
distributions.

This corrupted 105 of 2,564 scored sites, 8 of them dataset sites, and inflated
the reference standard deviation from 1.33 to 2.62. Every standardised effect
computed against it was wrong, including the reported primary estimate of
−0.057 SD. That number is withdrawn.

The scorer now refuses undecoded positions and refuses any row that is not a
distribution. The two checks are independent. A regression test drives real
backbone geometry with one sequon oxygen deleted through the model and asserts
the all-ones row appears, so the defect cannot return silently.

## What the controls are

The 32-site internal-control class is **sequons with no modelled glycan under
internal-control conditions**: sites in solved glycoprotein structures where a
glycan would have been visible had one been present. They are the only controls
that hold organism, compartment and experimental context roughly constant, and
they are the binding constraint on this study.

They are **not** demonstrated to be chemically unmodified. Absence of a modelled
glycan is weaker evidence than absence of glycosylation, and the earlier name
"observed-unmodified" claimed more than the data support.

## Reference scale

| | Value |
|---|---|
| Scoreable dataset sites | 342 (314 occupied, 28 internal control) |
| Reference SD | **1.3316** log-odds (was 2.6169) |
| Equivalence margin (±0.2 SD) | ±0.2663 log-odds |

Pooled across all scoreable dataset sites without consulting their labels, as
the frozen rule requires. The rule is unchanged; only the set it applies to is
corrected.

## Result

![primary comparison](../results/figures/primary_primary.png)

| | Value |
|---|---|
| Matched pairs | 16 |
| Occupied proteins / ortholog clusters / resampling units | 16 / 13 / 12 |
| Mean paired difference | **+0.864 log-odds (+0.649 SD)** |
| 95% CI, cluster-aware | [+0.100, +1.656] → **[+0.075, +1.243] SD** |
| Verdict | **directional, magnitude undetermined** |

The interval excludes zero, so a positive difference is indicated. It also
reaches from inside the equivalence margin to five times beyond it, so the size
of that difference is not established. Neither "the model prefers occupied
sequons by a meaningful amount" nor "the two are equivalent" is supported.

### The direction is not consistent

This is the most important qualification and it is not visible in the mean.

| | Value |
|---|---|
| Occupied scores higher in | 9 of 16 pairs |
| Sign test | p = 0.80 |
| Wilcoxon signed-rank | p = 0.058 |
| Median | +0.978 log-odds |
| 10% trimmed mean | +0.829 |
| Leave-one-out mean | stays within [+0.694, +1.001] |

Nearly half the pairs point the other way. The mean is positive because the
seven negative contrasts are all small — none below −1.19 — while the nine
positive ones reach +3.41. The effect is carried by **magnitude asymmetry, not
by direction**, and a test that reads only direction finds nothing.

It is not an outlier artefact: no single pair can be dropped to move the mean
outside [+0.69, +1.00]. But a bootstrap of the mean excluding zero, a Wilcoxon
at p = 0.058 and a sign test at p = 0.80 are genuinely mixed evidence at n=16.

### Sequon subtype

Matched exactly, so these are like-for-like.

| Subtype | n | Mean |
|---|---|---|
| NXS | 7 | +1.141 log-odds (+0.857 SD) |
| NXT | 9 | +0.649 log-odds (+0.487 SD) |

Both positive. At n=7 and n=9 the difference between them is not interpretable.

### Sensitivity: k = 5

The pre-specified matching is 1:1. At k=5 — the earlier setting — the result is
qualitatively unchanged: 14 contrasts, mean +0.703 log-odds, median +0.768,
Wilcoxon p = 0.119, same verdict. The conclusion does not depend on that choice.

## What was fixed to get here

| Issue | Before | Now |
|---|---|---|
| Invalid probability vectors | 105 sites scored at ~+13.8 | rejected at two independent guards |
| Scoreability | established after matching | established before, from coordinates alone |
| Sequon subtype | ~45% of pairs matched NXS to NXT | matched exactly |
| Control-protein reuse | ignored in the bootstrap | resampling unit is the connected component |
| `matched_pairs.csv` | no runner produced it | `runners/match_primary.py` |
| Unmatched-site scoring | depended on a `/tmp` file | `runners/build_candidate_manifest.py` |

## Limitations

- **16 pairs.** The binding constraint, and it has tightened: requiring exact
  subtype and dropping unscoreable sites cost pairs that the earlier count
  should not have had. Growing this class — realistically through PNGase F /
  H₂¹⁸O occupancy glycoproteomics — is what would make the comparison decisive.
- Matching holds local structure and subtype constant. It does **not** match
  organism, and the matched pairs span species.
- Internal controls are sites with no modelled glycan, not proven negatives.
- The margin is an exploratory statistical threshold, not a biologically
  validated one.
- ProteinMPNN parses only ATOM records, so glycans are invisible to it. Any
  effect is an occupancy-associated statistical preference, not evidence that
  the model represents glycosylation mechanistically.
- The k=1 choice was pre-specified before any corrected score existed, but not
  blind: the defective run's estimate had been seen.

## Reproducing

```
runners/build_candidate_manifest.py dataset     # map candidates to structures
runners/scoreability.py  ... candidate_manifest_dataset.csv scoreability_dataset.csv
runners/match_primary.py                        # matched sets + balance, k=1 and k=5
runners/score_all.py     ... candidate_manifest_dataset.csv scores_dataset.csv
runners/primary_analysis.py primary             # contrast, CI, equivalence
runners/primary_plot.py primary                 # the figure above
```

Artefacts: `results/candidate_manifest_dataset.csv`, `results/scoreability_dataset.csv`,
`results/matched_pairs_primary.csv`, `results/matching_balance_primary.json`,
`results/scores_dataset.csv`, `results/contrasts_primary.csv`,
`results/analysis_primary.json`, `results/figures/primary_primary.png`,
`config/scoring_frozen.toml` (see `[amendment_1]`).

## Secondary work, reported separately

Both have been rerun on validated scores. Neither is part of the conclusion above.

- **[Diagnostic controls](diagnostic_controls.md)** — bacterial (278 pairs,
  −0.174 SD) and cytosolic (270 pairs, +0.062 SD). The three comparisons now
  point in different directions, so the earlier "gradient" reading is withdrawn.
  These are confounded by construction and do not corroborate the primary result.
- **Design retention** — on a frozen snapshot, scoreable sites only: mean full
  retention 0.078, 80.4% of sites lose the sequon in all 32 designs, and the
  conditional score predicts retention at Spearman +0.559 with monotonic
  quintiles. The sweep was interrupted and covers chains alphabetically, so its
  coverage is not a random sample. `results/retention_analysis.json`.

The earlier claim that retention differed between occupied and control sites
(+0.084, CI [+0.007, +0.176]) is **withdrawn**: it used a site-level bootstrap,
and the matched pairs it rested on no longer exist after rematching.

## Still open

- `docs/phase4_primary_result.md` and `docs/methodology_explainer.md` are
  bannered as superseded but not rewritten.
- The retention contrast has not been recomputed on the new 16 matched pairs.
- No second model. Per the agreed scope, nothing else runs until this result is
  stable.
