# What the terms mean


> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](OVERVIEW.md) for current numbers.

Answers to the questions raised on 2026-08-19, in plain language. Grouped by
theme rather than by the order asked.

---

## The score itself

### "The score averages the log odds of N at the first position and of S/T at the third" — what does that mean?

Take one sequon, say **N-I-T** at positions 100–102. Show ProteinMPNN the
backbone and the entire native sequence *except* it re-derives what each residue
should be. At position 100 it outputs a probability for each of the 20 amino
acids: maybe 12% asparagine, 20% aspartate, 8% serine, and so on. We take the
12%.

**Log odds** converts a probability into a symmetric scale. p/(1−p) is the odds;
taking the log makes it symmetric around zero. So p = 0.5 → 0; p = 0.12 → −2.0;
p = 0.88 → +2.0. Without this, probabilities bunch up near zero and averages get
distorted.

We do the same at position 102, but there we add P(serine) + P(threonine),
because either satisfies the motif.

Then average those two log-odds values. That is the score for that site. It runs
from about −5 (the model thinks this sequon is very unlikely) to +1 (quite
likely). All the scores in this study are negative on average, which just means
ProteinMPNN generally does not favour these residues.

### Is this at the sequon level or the protein level?

Sequon level — three residues. We deliberately did not use a whole-protein score
because three residues out of several hundred would be swamped.

### Why exclude the middle residue?

Any residue except proline works there, so a preference at the middle position is
not a preference for the motif — it would be measuring something else. Proline is
recorded separately because a proline there abolishes the sequon.

---

## The control sets

### What are the 32 "internal" controls, and why "internal"?

They are sequons with **no modelled glycan** in a structure that *does* model
glycans somewhere else on the same protein, from an organism that can
glycosylate.

The logic: normally a bare asparagine in a structure proves nothing — the glycan
may have been trimmed off before crystallisation, or the protein expressed in
bacteria, or the sugar too floppy to see. But if the crystallographer clearly
*could* see glycans, because they modelled some, and this particular asparagine
is bare, that absence is a decision rather than a silence.

"Internal" means they come from **inside the same set of structures** as the
occupied sites — same organisms, same secretory compartment, same kind of
experiment. That is what makes them the cleanest comparison. The trade-off is
there are only 32.

They are **not** proven negatives. That is why the wording changed to "no
modelled glycan under internal-control conditions".

### What does "diagnostic" mean for the bacterial and cytosolic sets?

It means: informative about how the measurement behaves, but not a valid answer
to the question.

Both are confounded on purpose, and in opposite ways:

- **Cytosolic** proteins are eukaryotic like the occupied ones, but they live in
  the cytosol and never meet the glycosylation machinery. Confound: compartment.
- **Bacterial** secreted and membrane proteins are the right *kind* of protein,
  but from organisms with no N-glycosylation. Confound: kingdom.

The original hope was that using both would triangulate: if a signal survives
both confounds, it is more likely real. That reasoning is described below.

### What data underlie the comparisons in Figure 2?

The same ProteinMPNN conditional sequon score, computed the same way on all
sites. The only thing that differs is which control group the occupied sites are
being compared against.

### "The logic behind the controls did not survive the corrections" — what does that mean?

Before the bug was fixed, the three comparisons lined up neatly:

| | Old (buggy) |
|---|---|
| Cytosolic | −0.237 |
| Bacterial | −0.145 |
| Internal | −0.057 |

All negative, and shrinking steadily as the control set got better matched. That
looked like meaningful evidence — it suggested the apparent differences were
caused by compartment and kingdom, and that once you removed those confounds
almost nothing was left.

After the fix:

| | Corrected |
|---|---|
| Cytosolic | +0.062 |
| Bacterial | −0.174 |
| Internal | +0.458 |

No ordering, no pattern, three different directions. So the triangulation
argument no longer has anything to stand on. That is all "did not survive"
means — the reasoning was sound, but the numbers it rested on were wrong, and
the corrected numbers do not support it.

---

## Matching

### What is a "match"?

Pairing one occupied site with one control site that sits in a similar local
structural environment, so that a score difference cannot be blamed on the
environment.

### Which metrics?

Three, all computed from the structure around the asparagine:

- **Relative solvent accessibility (RSA)** — how exposed the residue is, from 0
  (fully buried) to 1 (fully exposed). Computed with a rolling-ball algorithm
  and divided by the maximum possible for that amino acid type.
- **Neighbour count within 8 Å** — how many other residues have an atom within
  8 ångströms. A crude density measure: high means tightly packed.
- **Hydrophobic fraction within 8 Å** — of those neighbours, what proportion are
  hydrophobic. Describes the chemical character of the pocket.

A pair is allowed only if it is within a **caliper** — a maximum distance in
these three dimensions combined, set at 0.25 pooled standard deviations. Beyond
that, no match is made rather than a bad match.

### Why must NXS pair with NXS and NXT with NXT?

NXT sequons are glycosylated more often than NXS in nature, and threonine and
serine are chemically different, so the model may score them differently for
reasons having nothing to do with occupancy.

Before this rule, about 45% of pairs had an occupied NXS against an unoccupied
NXT. Any difference measured could then have been a subtype difference wearing
occupancy's clothes. Requiring them identical removes that, at the cost of pairs.

### Greedy matching, seeds, deterministic — what are these?

**Greedy matching** walks through the occupied sites one at a time and gives each
one its nearest unused control. It's called greedy because it takes the best
option available at each step without considering the consequences for later
steps.

**The seed** sets the order it walks through them. It is a number fed to a random
number generator; the same seed always gives the same order.

With only 28 controls this matters enormously. Suppose control C is the only
admissible partner for occupied site B, but site A is also near C and gets
processed first. A takes C, and B is left with nothing. Change the order and both
get matched. So the seed changes both *how many* pairs form and *which* ones.

**Deterministic** matching does not walk in any order. It solves the whole
assignment at once, finding the arrangement that produces the most pairs and,
among those, the tightest total distance. No seed, no ordering, one answer.

### What is the 200-seed sweep for?

To find out how much the old answer depended on that arbitrary seed. Run greedy
matching 200 times with 200 different seeds, and look at the spread.

Result: all 200 gave a positive estimate, but only 75 produced an interval
excluding zero. So the *direction* was robust and the *significance* was being
decided by a random number. That is the argument for using the deterministic
version.

### Does deterministic matching allow multiple matches?

No — still one control per occupied site, each control used once. "Deterministic"
refers only to removing the randomness in how pairs are chosen, not to how many
are allowed. (A variant allowing up to five controls per site was tested as a
sensitivity check; it gives a similar answer with fewer independent cases.)

### What happens without matching?

Worth knowing, because it is the natural question:

| | Occupied | Internal control | Difference |
|---|---|---|---|
| **Unmatched** (all sites) | n=314 | n=28 | +0.213 SD, p = 0.25 |
| **Matched** (pairs) | n=16 | n=16 | +0.458 SD, CI includes 0 |

Both positive, neither significant. Notably the two groups were already fairly
similar structurally before matching (SMD 0.11–0.15, well under the 0.1–0.2
threshold usually considered concerning), so matching bought less here than it
did for the bacterial and cytosolic sets, where the imbalance was 0.5+.

---

## Clusters and uncertainty

### What is a cluster, and what does "cluster-aware" mean?

An **ortholog cluster** is a group of the same protein across different species —
human, mouse, zebrafish versions of the same gene. They come from the ortholog
database this module draws its sites from.

Sites in the same cluster are not independent observations. If human and mouse
serum albumin both have an occupied sequon at the equivalent position, that is
close to one fact observed twice, not two facts.

**Cluster-aware** statistics take this into account. If you treat 16 correlated
observations as 16 independent ones, your confidence interval comes out too
narrow and you will claim more certainty than you have.

### What is the cluster bootstrap?

A bootstrap estimates uncertainty by resampling your data many times with
replacement and seeing how much the answer moves.

A naive bootstrap resamples individual sites. A **cluster** bootstrap resamples
whole groups — either all of a cluster is in a given resample or none of it is.
That respects the fact that members of a group move together.

Here the groups are slightly more complicated than ortholog clusters alone,
because a single control protein can be matched to several different occupied
sites. So the resampling unit is the connected group formed by linking occupied
clusters to the control proteins they share. The 16 pairs collapse into 12 such
groups, which is why the interval is wider than 16 independent pairs would give.

### Why can't an ordinary significance test answer this?

Because we expect, and largely want to be able to state, a **null** — that the
model treats occupied and unoccupied sequons the same.

An ordinary test can only ever reject the null or fail to reject it. Failing to
reject means "we did not find a difference", which is compatible with both "there
is no difference" and "there might be a big difference but our sample is too
small to see it". With 16 pairs, the second is very much live. So a
non-significant p-value here would tell you almost nothing.

**Equivalence testing** turns it around: state in advance how small a difference
would count as "effectively no difference", then ask whether the confidence
interval fits entirely inside that band. If it does, you have positive evidence
of no meaningful effect. If it doesn't, you are honestly inconclusive.

### How was the ±0.2 SD margin chosen?

Somewhat arbitrarily, and it is labelled as such throughout — an *exploratory
statistical threshold, not a biologically validated one*. 0.2 standard deviations
is a conventional "small effect" in the social-science effect-size literature
(Cohen's d), and nothing in glycobiology tells us what a meaningful shift in
ProteinMPNN log-odds would be.

What matters more than the exact number is that it was written down **before**
any comparison was computed. Choosing a margin after seeing the differences is
not a test.

---

## Retention

### What are the quintiles in Figure 4?

Take all 2,423 sites, sort them by their conditional score, and cut into five
equal-sized groups: the lowest-scoring fifth, the next fifth, and so on. Then ask
for each group: on average, what fraction of designs kept the sequon?

### Where does that score come from — is it the same score?

Yes, the same site-level conditional sequon score described at the top. So the
figure is asking whether the number we measure predicts the behaviour we care
about.

### Why does that graph matter?

Because the two analyses measure different things and could easily have been
unrelated.

The conditional score asks *what probability the model holds at this site while
looking at the native sequence*. Retention asks *what the model actually writes
when generating a fresh sequence*. Generation involves hundreds of interacting
decisions, so a sequon could survive or vanish for reasons unrelated to the
model's opinion about that one position.

The figure shows they track each other closely and monotonically. That means the
abstract number is a valid proxy for the behaviour, and the two branches of the
project are describing one underlying quantity rather than two unrelated ones.

### 2,423 sites × 32 designs — is that 77,000?

Yes, about 77,500 designed sequences. Designs are generated per chain and every
sequon on that chain is read from the same 32 designs, so the cost is per chain
rather than per site.

### What determined "scoreable"?

Whether ProteinMPNN will process all three sequon residues — which comes down to
whether the backbone atoms (N, CA, C, O) are all present in the structure. See
the Figure 6 section below.

### Does the model remove non-glycosylated motifs at the same rate?

This changed once a properly matched control set existed, and the change is
instructive.

**Class averages say yes.** Comparing whole populations:

| Class | n | Mean retention | Lost in every design |
|---|---|---|---|
| Occupied | 285 | 0.080 | 78.6% |
| Bacterial control | 1,096 | 0.077 | 79.5% |
| Cytosolic control | 1,021 | 0.066 | 84.4% |
| Eukaryotic secretory control | 255 | 0.036 | 85.5% |
| Internal control | 21 | 0.006 | 90.5% |

**But a class average is the wrong test.** The control sets differ from one
another — the eukaryotic secretory set sits at half the retention of the
bacterial and cytosolic sets — so averaging across populations mixes in whatever
else distinguishes them. Every control was matched site-by-site to an occupied
site, so the paired contrast is available and is the right comparison.

**Paired, the answer is less clean:**

| Comparison | Pairs | Occupied | Control | Difference | p (cluster permutation) |
|---|---|---|---|---|---|
| Internal control | 16 | 0.160 | 0.008 | +0.152 | 0.125 |
| Eukaryotic secretory | 245 | 0.081 | 0.037 | +0.043 | 0.030 |
| Bacterial | 254 | 0.079 | 0.076 | +0.002 | 0.816 |
| Cytosolic | 251 | 0.078 | 0.075 | +0.003 | 0.821 |

Occupied sequons are retained *more* than their matched partners in both
comparisons that hold eukaryotic secretory context constant, and not at all in
the two confounded sets. That ordering is the opposite of what confounding by
compartment or taxonomy would produce, which is what makes it interesting.

**It does not reach significance.** Across the eight tests run (four control
sets × two outcomes), nothing survives multiple-comparison correction — the
smallest corrected p is 0.120. See [`significance.md`](significance.md), which
also explains why the Wilcoxon p of 0.008 originally quoted for the secretory
contrast was too optimistic: it treated 245 pairs as independent when the
effective sample size is 73 clusters.

So the honest answer is: **no detectable difference that survives correction**,
with a suggestive pattern worth testing on data that did not generate it.

### Is this different from the background mutation rate at non-motif positions?

**We cannot say — that number was not recorded.** The retention runner stored only
the per-site classification, not the full designed sequences, so overall sequence
recovery is not recoverable without regenerating the designs. It is a fair
question and a real gap. ProteinMPNN's published sequence recovery is roughly
40–50% on native backbones, which would suggest sequons are lost far more often
than residues in general, but that is a literature comparison rather than a
matched internal one.

---

## The bug (Figure 6)

### Why won't ProteinMPNN decode a residue with an incomplete backbone?

Because it builds its representation of a residue from the geometry of its four
backbone atoms — N, CA, C, O. If one is missing from the deposited structure
(common at flexible loops and chain ends), the geometry is undefined and the
model marks that position as unusable.

### What is the "21 ones" business?

Internally the model fills in a table with one row per residue and 21 columns
(20 amino acids plus an unknown token). It creates the table full of **zeros**,
then writes real values only into the rows it decoded. Rows it skipped stay zero.

Those stored values are *log* probabilities, so the final step exponentiates
them: e^value. And e⁰ = 1. So a skipped row comes out as twenty-one 1s.

Read as probabilities that says P(asparagine) = 1, and P(serine) + P(threonine) =
1 + 1 = 2. A probability of 2 is meaningless — that is the tell. Pushed through
the log-odds formula it produced a score of about +13.8, where real scores span
−5 to +1.

### How was it fixed — just deciding scoreability before matching?

Two separate things, and both were needed:

1. **The guards.** The scorer now refuses to score a site if the model reports
   it did not decode that position, *and* separately refuses any row whose
   probabilities do not sum to 1. Two independent checks, so neither is a single
   point of failure. This stops bad values existing at all.

2. **The ordering.** Scoreability is now determined *before* matching rather than
   after. This is not about correctness of individual scores — it is about study
   design. Previously sites were matched into balanced pairs and then some were
   dropped for being unscoreable, which unbalanced the sets that matching had
   just carefully balanced. Doing it first means the matching only ever sees
   sites that will survive.

Conveniently, scoreability depends only on the coordinates, so it can be checked
without running the model at all.
