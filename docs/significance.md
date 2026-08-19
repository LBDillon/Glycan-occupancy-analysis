# Significance testing

Eight tests: four control sets × two outcomes. Run with
`python runners/significance.py`. Raw numbers in `results/significance.csv`.

## Headline

**No test survives correction for the number of comparisons made.** The smallest
raw p is 0.030; the smallest after Benjamini–Hochberg is 0.120, and after Holm
0.237. Zero of eight tests reach 0.05 corrected.

## Results

| Outcome | Comparison | Pairs | Units | Effect | p (perm) | p (Wilcoxon) | p Holm | p BH |
|---|---|---|---|---|---|---|---|---|
| Score | internal control | 16 | 12 | +0.458 SD | 0.251 | 0.193 | 1.000 | 0.436 |
| Score | eukaryotic secretory | 262 | 72 | +0.073 SD | 0.327 | 0.157 | 1.000 | 0.436 |
| Score | bacterial | 280 | 54 | −0.157 SD | **0.030** | 0.063 | 0.237 | 0.120 |
| Score | cytosolic | 273 | 77 | +0.067 SD | 0.285 | 0.208 | 1.000 | 0.436 |
| Retention | internal control | 16 | 12 | +0.152 | 0.125 | 0.063 | 0.747 | 0.332 |
| Retention | eukaryotic secretory | 245 | 73 | +0.043 | **0.030** | 0.008 | 0.237 | 0.120 |
| Retention | bacterial | 254 | 55 | +0.002 | 0.816 | 0.990 | 1.000 | 0.821 |
| Retention | cytosolic | 251 | 73 | +0.003 | 0.821 | 0.819 | 1.000 | 0.821 |

## Why a permutation test rather than Wilcoxon

Wilcoxon and the sign test assume the pairs are independent. They are not.
Occupied sites in the same ortholog cluster are near-copies of one another, and
a single control protein can be matched to several occupied cases. Those
dependencies were already handled in the confidence intervals — the bootstrap
resamples connected components — but the rank tests quoted alongside them were
not.

The test used here flips the sign of every contrast within a whole component at
once. Under the null, a component's contrasts are exchangeable in sign;
individual contrasts are not.

The difference is not cosmetic. For the eukaryotic secretory retention contrast:

| | p |
|---|---|
| Wilcoxon, treating 245 pairs as independent | 0.008 |
| Cluster permutation, 73 resampling units | 0.030 |

The effective sample size is **73**, not 245. Roughly a quarter of what a naive
test assumes.

## Why correct across all eight

These were not eight pre-planned tests. The eukaryotic secretory set was added
after seeing that the internal-control comparison was too small, and the paired
retention analysis was run after seeing the class averages. That is a reasonable
way to work, but it means the family of tests grew in response to results, which
is exactly the situation multiple-comparison correction exists for.

Holm controls the family-wise error rate; Benjamini–Hochberg the false discovery
rate. Both are reported. Neither leaves anything below 0.05.

## What this does and does not change

**Unchanged.** The pre-specified primary inference was never a p-value. It is an
equivalence assessment of the conditional score against a ±0.2 SD margin, and it
remains inconclusive in every comparison — the intervals are too wide to
establish equivalence and too close to zero to establish a difference.

| Comparison | 95% CI (SD) | Verdict |
|---|---|---|
| internal control | [−0.227, +1.098] | inconclusive |
| eukaryotic secretory | [−0.056, +0.346] | inconclusive |
| bacterial | [−0.494, −0.058] | directional, magnitude undetermined |
| cytosolic | [−0.056, +0.298] | inconclusive |

**Changed.** The paired retention difference against the eukaryotic secretory
set was previously described as excluding zero with the confidence interval,
sign test and Wilcoxon all agreeing. Two of those three were computed under an
independence assumption that does not hold, and none of them accounted for the
number of comparisons run. The corrected reading is that it is a **suggestive
lead that does not reach significance**, not a finding.

The pattern it sits in remains interesting and is worth testing properly: the
effect appears in both comparisons matched on eukaryotic secretory context and
in neither confounded set, which is the opposite of what confounding by
compartment or taxonomy would produce. That is a hypothesis for a
pre-registered test on data not used to generate it — not a result.
