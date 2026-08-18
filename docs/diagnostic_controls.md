# Appendix: bacterial and cytosolic diagnostics

**These are not the primary result and cannot substitute for it.** The primary
comparison is [`primary_result.md`](primary_result.md).

Both control sets are confounded by construction, and deliberately so. The
cytosolic set matches taxonomy and differs in subcellular compartment; the
bacterial set matches compartment and differs in taxonomy. They were built to be
read against one another, on the reasoning that a signal surviving both
orthogonal confounds is less likely to be produced by either.

They were also built because there are only 28 scoreable internal controls. That
scarcity is the real constraint, and a large well-matched comparison against the
wrong population does not relieve it.

## Matching

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

## Results

All on the common reference scale (SD = 1.3316 log-odds, from dataset sites only).

| Comparison | n | Mean | 95% CI (cluster-aware) | Site-level CI (not used) | Verdict |
|---|---|---|---|---|---|
| **Internal control** (primary) | 16 | **+0.649 SD** | [+0.075, +1.243] | [+0.123, +1.202] | directional, magnitude undetermined |
| Bacterial extracytoplasmic | 278 | −0.174 SD | [−0.459, −0.027] | [−0.334, −0.015] | directional, magnitude undetermined |
| Cytosolic eukaryotic | 270 | +0.062 SD | [−0.145, +0.270] | [−0.112, +0.233] | inconclusive |

Robustness:

| | Bacterial | Cytosolic |
|---|---|---|
| Occupied higher | 125 of 278 | 151 of 270 |
| Median | −0.191 | +0.214 |
| Sign test | p = 0.105 | p = 0.059 |
| Wilcoxon | p = 0.044 | p = 0.228 |

## What this does and does not support

**The earlier "gradient" is withdrawn.** The previous version of this analysis
reported all three comparisons as negative, shrinking as matching improved
(−0.237, −0.145, −0.057 SD), and read that ordering as evidence that the
apparent effects came from compartment and taxonomy rather than occupancy. That
pattern was an artefact of the corrupted scores. It does not survive correction.

What replaces it is less tidy. The three comparisons now point in **different
directions**: positive against internal controls, negative against bacterial
controls, and indistinguishable from zero against cytosolic controls. There is
no monotone ordering to interpret.

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

## Subtype

| Comparison | NXS | NXT |
|---|---|---|
| Bacterial | −0.167 SD (n=132) | −0.181 SD (n=146) |
| Cytosolic | −0.094 SD (n=135) | +0.218 SD (n=135) |

The cytosolic subtypes disagree in sign at n=135 each, which is not a sample-size
problem. It is unexplained and recorded as such.

## Artefacts

`results/matched_pairs_{bacterial,cytosolic}.csv`,
`results/matching_balance_{bacterial,cytosolic}.json`,
`results/contrasts_{bacterial,cytosolic}.csv`,
`results/analysis_{bacterial,cytosolic}.json`,
`results/scores_controls.csv`, `results/manifest_matched_controls.csv`,
`results/scoreability_controls.csv`.

Reproduce with `runners/match_diagnostics.py`, then `runners/score_all.py` on
`results/manifest_matched_controls.csv`, then
`runners/primary_analysis.py bacterial` and `... cytosolic`.
