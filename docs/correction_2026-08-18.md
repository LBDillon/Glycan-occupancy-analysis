# What changed on 18 August 2026, and why


> **Note (2026-08-20).** This documents the *August 18th* correction. Its
> numbers are the historical record of that moment and are deliberately left as
> they were. A later, larger correction — the ProteinMPNN token alphabet —
> supersedes them for any current purpose. See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md).

A short account of two rounds of correction on the same day. The first
(Amendment 1) fixed a defect in the scorer; the second (Amendment 2, at the end
of this document) removed an arbitrary dependence in the matching. The detailed
result is in [`primary_result.md`](archive/primary_result_SUPERSEDED_2026-08-25.md); the machine-readable
record is `config/scoring_frozen.toml`, sections `[amendment_1]` and
`[amendment_2]`.

## What prompted it

A review of the scoring code and its output rows found that some scored sites
were not scores at all. The check that caught it was simple: a probability
vector has to sum to one, and some of these summed to twenty-one.

## What was wrong

ProteinMPNN's `conditional_probs` creates an empty table of log-probabilities
and fills in only the residues it actually decodes. It refuses to decode any
residue whose backbone is incomplete — a missing N, CA, C or O — and leaves that
row as zeros. Exponentiating zeros gives twenty-one ones. Read as probabilities
that says P(asparagine) = 1 and P(serine or threonine) = 2, which is impossible,
and it produced a site score of about +13.8 where real scores run from −5 to +1.

The manifest had confirmed that all three sequon residues had coordinates and a
position in the chain. That is a genuinely different claim from the model having
accepted them, and nothing downstream checked the difference.

**Why it mattered more than the count suggests.** Only 105 of 2,564 sites were
affected, and only 8 of those were dataset sites. But those few enormous values
inflated the reference standard deviation from 1.33 to 2.62. Every standardised
effect in the study was divided by that number, so a 4% data problem became a
100% error in the scale everything was measured against — and the primary
estimate reversed sign once the invalid rows were removed.

## What I changed, and why

**Two independent guards in the scorer.** It now refuses any position the model
declined to decode, and separately refuses any row that is not a probability
distribution. Two checks rather than one because either alone can be defeated:
the first depends on correctly reading the model's mask, the second on the
arithmetic. Neither is a single point of failure.

**A regression test built from real geometry.** The existing small fixture has
collinear coordinates, which makes the model behave strangely for unrelated
reasons. So the new fixture is a real 21-residue backbone with one sequon oxygen
deleted. It reproduces the all-ones row through the actual model, which means
the test would fail if the bug ever returned.

**Scoreability is now settled before matching, not after.** This was the
structural fix rather than the surface one. Whether a site can be scored depends
only on its coordinates, so it can be determined without running the model at
all. Deciding it afterwards meant sites were matched into balanced sets and then
dropped, quietly unbalancing the sets that matching had just balanced. Checked
against the existing scores, the advance screen predicted every corrupted row
with no false alarms.

**Sequon subtype is matched exactly.** Roughly 45% of pairs had been matching an
occupied NXS site against an unoccupied NXT one. NXT is occupied more often and
the two are chemically distinct, so subtype was confounded with occupancy — the
comparison was partly measuring the wrong thing.

**One control per occupied site.** With 28 controls against 314 occupied sites,
allowing up to five controls each let the first few cases absorb the entire pool,
leaving 14 usable cases out of 314. Pairing one-to-one gives every control its
own case. I fixed this before any corrected score existed, but the earlier
defective estimate had already been seen, so it is recorded as a deviation rather
than a blind choice. Running it the old way changes nothing material.

**Uncertainty now accounts for reused control proteins.** The frozen
configuration had promised this and the code did not do it. Two things tie
contrasts together: occupied sites in the same ortholog cluster are near copies,
and one control protein can serve several occupied cases. Resampling on either
alone leaves the other unhandled, so the unit is now the connected group of
clusters and controls that share members.

**Missing runners written.** Some reported stages could not be regenerated from
committed code — the matched pairs had no runner, and one script depended on a
file in `/tmp`. Both are now reproducible from the repository.

**The controls renamed.** The 32-site class was called "observed-unmodified",
which asserts the sites are chemically unmodified. They are not shown to be. They
are sites with **no modelled glycan under internal-control conditions**: strong
evidence of absence by structural standards, but still a statement about the
deposited model rather than the molecule. They remain the best controls available.

## What the results were after Amendment 1

Correcting the scorer reversed the sign. Against the withdrawn −0.057 SD, the
comparison now gave **+0.649 SD, 95% CI [+0.075, +1.243]** on 16 matched pairs,
an interval that excluded zero.

That reading did not survive Amendment 2. See the end of this document.

**The diagnostic controls no longer tell a tidy story.** They previously appeared
to shrink toward zero as matching improved, which was read as evidence that the
apparent effects came from compartment and taxonomy. Corrected, the three
comparisons point in different directions: +0.649 internal, −0.174 bacterial,
+0.062 cytosolic. That ordering is gone and I have not replaced it with another
interpretation.

**Design retention** was affected by the same underlying problem through a
different route: an undecodable residue keeps its native identity in every
generated design, so it appears preserved without the model ever having a choice.
Excluding those sites, retention still tracks the conditional score closely
(Spearman +0.547 over 2,423 sites, rising monotonically across score quintiles), which is the
bridge to the earlier scoping analysis. It is reported separately and forms no part of the
primary conclusion.

## What this round does not settle

Sixteen pairs is the binding constraint, and correcting the analysis tightened it
rather than loosening it — the earlier count included pairs it should not have.
Growing the internal-control class, realistically through occupancy
glycoproteomics, is what would make this comparison decisive. No second model
should run until this one is stable.

The honest summary is that the study now measures what it claims to measure, and
that what it measures is not yet precise enough to answer the question.


---

# Amendment 2: the matching was making the decision

## What was wrong

Nothing in the data. The problem was the algorithm that pairs occupied sites
with controls.

Greedy nearest-neighbour matching walks the cases in a random order and gives
each one its nearest unused control. With 28 controls for 314 occupied sites,
that order decides things: an early case can take the only control that was
admissible for a later one. The order comes from a seed, and the reported result
came from seed 0. It was one draw from a distribution, written down as a value.

## What I did

Replaced it with the assignment that maximises the number of admissible pairs
and, among those, minimises total distance. It is solved directly rather than
approached greedily, so there is no seed and no ordering. Rows are sorted by
accession and position first, so even exact ties break the same way regardless
of how the data arrived. Features, caliper, the exact NXS/NXT requirement and
one-control-per-case are all unchanged — only the pairing algorithm differs.

Then I ran the greedy version 200 times to see how much the old answer had
depended on its seed.

## What that showed

| | |
|---|---|
| Pairs, every seed | 16 |
| Mean contrast across seeds | +0.286 to +0.699 SD |
| Point estimates that were positive | **200 of 200** |
| Intervals that excluded zero | **75 of 200 (38%)** |
| Deterministic optimal matching | +0.458 SD, CI [−0.227, +1.098], inconclusive |

The direction is robust and the significance is not. Every matching we tried put
occupied sites higher; fewer than half produced an interval that excluded zero.
Whether the earlier result "reached significance" was settled by an arbitrary
choice inside the matching rather than by the data, which is exactly why that
choice should not be left to a seed.

The pair count is 16 under every seed and under the optimum, so the caliper, not
the algorithm, is what limits it: 12 of the 28 controls have no partner within
0.25. The optimum finds *better* pairs, not more — mean distance 0.099 against
greedy's 0.113, with tighter balance.

## The headline now

> Occupied sites tend to score higher than matched sites with no modelled
> glycan, but 16 pairs do not establish a precise or statistically robust
> difference.

This supersedes the Amendment 1 reading that the interval excluded zero.

The interval is worth reading asymmetrically. It runs from −0.227 to +1.098 SD,
so it is consistent with no difference and with a substantial positive one, and
it excludes only large differences favouring the controls. That is a weak
conclusion, but it is the one the data support.

## Housekeeping in the same pass

- `primary_analysis.py` now loads only what the requested comparison needs, so an
  analysis of 16 dataset pairs no longer depends on control files it never reads.
- The reference population is reported from dataset sites alone, so the scale
  cannot shift when a control pool is rebuilt.
- Contrast construction, the resampling unit and the bootstrap moved to a shared
  module used by both the primary analysis and the sensitivity sweep — a
  sensitivity check that computed its interval differently from the result it
  tests would look like corroboration while being nothing of the kind.
- `score_unmatched.py` removed: superseded, and it read a file from `/tmp` that
  no one else could reproduce.

**The ProteinMPNN analysis is frozen here.** The constraint is 16 pairs, and no
further modelling changes that.
