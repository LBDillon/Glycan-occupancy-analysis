# What changed on 18 August 2026, and why

A short account of one round of corrections. The detailed result is in
[`primary_result.md`](primary_result.md); the machine-readable record is
`config/scoring_frozen.toml`, section `[amendment_1]`.

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

## What the results are now

The primary comparison gives **+0.649 SD, 95% CI [+0.075, +1.243]**, on 16
matched pairs — the opposite sign to the withdrawn −0.057 SD.

It should be read carefully. The interval excludes zero, so a positive difference
is indicated, but it stretches from inside the equivalence margin to five times
past it, so the size is not established. More importantly, occupied sites score
higher in only 9 of 16 pairs. The mean is positive because the seven negative
contrasts are all small while the nine positive ones are large. The effect lives
in magnitude, not in consistency, and a test that reads only direction finds
nothing at all.

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
(Spearman +0.559, rising monotonically across score quintiles), which is the
bridge to the preprint. It is reported separately and forms no part of the
primary conclusion.

## What this round does not settle

Sixteen pairs is the binding constraint, and correcting the analysis tightened it
rather than loosening it — the earlier count included pairs it should not have.
Growing the internal-control class, realistically through occupancy
glycoproteomics, is what would make this comparison decisive. No second model
should run until this one is stable.

The honest summary is that the study now measures what it claims to measure, and
that what it measures is not yet precise enough to answer the question.
