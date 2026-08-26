# What we are doing, and why

*A plain-language account of the reasoning behind the method, written so the

> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](OVERVIEW.md) for current numbers.

logic can be checked independently of the code. Current as of 2026-08-18. The
ProteinMPNN analysis is frozen; see [`primary_result.md`](archive/primary_result_SUPERSEDED_2026-08-25.md) for
the result table and [`correction_2026-08-18.md`](correction_2026-08-18.md) for
what was corrected along the way.*

## The question underneath everything

Protein design models such as ProteinMPNN are used to redesign real proteins,
many of which are glycosylated. The models are trained on structures stripped of
glycans, so a natural worry is that they treat a glycosylation site as an
ordinary patch of surface and quietly remove it. An earlier scoping analysis on
a handful of proteins suggested they frequently do remove natural sequons; this
module is what tests that properly.

That leaves a sharper question. When ProteinMPNN keeps or destroys an N-X-S/T
motif, is it responding to anything biological? Specifically: does it treat a
sequon that is *actually glycosylated* differently from one that merely matches
the motif?

Answering that needs something the field mostly lacks — a set of sequons whose
occupancy is known, one way or the other, on evidence rather than assumption.

## Why building the dataset took most of the effort

The ortholog database already knew where sequons are. It did not know which are
used. Two problems had to be solved before any model could be asked anything.

**A sequon is a motif, not a modification.** N-X-S/T is necessary for N-linked
glycosylation but nowhere near sufficient. So "this protein has a sequon" and
"this protein is glycosylated here" are different claims, and only the first was
in the database.

**The database counts comparisons, not sites.** It stores one row per orthologous
comparison, so a single asparagine appears once for every ortholog it was
compared against — 13,816 rows collapse to 4,307 distinct sites. Counting rows
would have multiplied one biochemical fact into dozens of apparent observations.

Re-indexing around one protein and one residue, then gathering evidence from
UniProt, GlyGen and glycan linkages read out of deposited structures, produced
**922 sites with experimental evidence of a glycan**.

## The part that is easy to get wrong

The tempting move is to treat every unannotated sequon as unoccupied and compare
the two groups. That is wrong, and wrong in a way that would look like a result.
Absence of annotation overwhelmingly means nobody looked. Well-studied proteins
accumulate annotations; obscure ones do not. A study built that way measures
curation effort and reports it as glycobiology.

So the dataset carries three states rather than two: **occupied** (922),
**no modelled glycan under internal-control conditions** (32), and **unknown**
(3,353). Most sites are honestly unknown, and the design refuses to convert that
into a negative.

The 32 exist only because of a narrow argument. A bare asparagine in a crystal
structure normally proves nothing — glycans are routinely trimmed before
crystallisation, expressed away in bacteria, or too mobile to model. But if the
*same structure* models glycans at other residues, and the protein came from a
host that can glycosylate, then sugars demonstrably survived and this
crystallographer demonstrably modelled them. A bare asparagine there is a
decision, not a silence.

**They are still not proven negatives.** They are the most informative internal
controls available, and the earlier name for them — "observed unmodified" —
claimed more than the evidence supports. Absence of a modelled glycan is a
statement about the deposited model, not about the molecule.

## Why there are diagnostic control sets as well

Thirty-two internal controls cannot carry a null result. So two further sets of
sequons that *cannot* be N-glycosylated were assembled: cytosolic eukaryotic
proteins, which never enter the secretory pathway and so never meet the enzyme,
and bacterial periplasmic and outer-membrane proteins, whose clades have no
equivalent machinery. After structural feature extraction and the scoreability
screen, 3,024 and 3,068 of these can be scored.

Neither is a clean substitute. The cytosolic set differs from the occupied sites
in subcellular compartment; the bacterial set differs in taxonomy. Their
confounds do not overlap, which was the design intent: a model that separates
occupied sites from the cytosolic set but not from the internal controls has
learned where proteins live; one that separates them from the bacterial set but
not the cytosolic set has learned taxonomy.

**That reading did not survive the corrections.** The three comparisons now point
in different directions rather than forming an interpretable ordering, and the
gradient once read off them is withdrawn. They are reported as diagnostics, in an
appendix, and they do not corroborate the primary result.

## Why matching was necessary before any scoring

Occupied sites are not a random sample of sequons: they sit disproportionately in
exposed loops, because the enzyme has to reach them. Structure-based models can
see exposure. An unmatched comparison of model scores would therefore mostly
restate that occupied sites are exposed, which is already known and says nothing
about glycosylation.

Each occupied site is matched to an internal control of comparable local
environment — accessibility, packing, neighbourhood composition — and, since the
corrections, with **NXS/NXT required to be identical**. Around 45% of pairs had
previously matched an occupied NXS against an unoccupied NXT, which confounded
subtype with occupancy.

Matching is now **deterministic**: the assignment maximising the number of
admissible pairs and, among those, minimising total distance. The earlier greedy
matcher walked the cases in a seeded random order, which mattered enormously with
only 28 controls — an early case could take the only admissible partner for a
later one.

## What is being measured, and what was deliberately not

The primary measurement is the **conditional probability** the model assigns to
the motif-forming residues at their own positions, given the backbone and the
rest of the native sequence. Nothing is generated and the sequon is never
altered. The score averages the log odds of asparagine at the first position and
of serine-or-threonine at the third.

The middle residue is excluded: any residue except proline permits a sequon, so a
preference there is not a preference for the motif. Proline is recorded
separately as the residue whose presence would abolish it.

The whole-protein score was rejected as the primary measure. Three residues
contribute almost nothing to an average over several hundred, so a protein-level
number cannot answer a site-level question.

## One defect worth understanding, because it inverted the answer

ProteinMPNN decodes only residues whose backbone is complete. For any residue
missing an N, CA, C or O it returns a row of zeros, which exponentiates to
twenty-one ones — P(asparagine) = 1, P(serine or threonine) = 2, a score near
+13.8 where real scores run from −5 to +1.

This affected 105 of 2,564 sites, only 8 of them dataset sites. But those few
enormous values inflated the reference scale from 1.33 to 2.62, and every
standardised effect was divided by it. A 4% data problem became a 100% error in
the units, and the first reported result had the wrong sign.

The fix has two independent guards, and scoreability is now settled **before**
matching rather than after — it depends only on the coordinates, so no model pass
is needed. Establishing it afterwards had let unscoreable sites into matched sets
and removed them later, unbalancing the sets matching had just balanced.

---

# The reasoning, in order

## Hypothesis

If ProteinMPNN has learned anything about the biological use of the N-X-S/T
motif, rather than the motif's mere appearance, it should assign a higher
site-level probability to sequons that are actually glycosylated than to matched
sequons that are not.

## The cleanest experiment

Compare **occupied** sites against **internal controls** — sequons with no
modelled glycan under conditions where a glycan would have been visible — matched
on local structural context and sequon subtype.

This is the cleanest because the two groups share nearly everything that is not
occupancy: same kind of organism, same subcellular compartment, same kind of
experiment, and — after matching — the same solvent accessibility, packing,
neighbourhood composition and NXS/NXT identity.

Its weakness is size, and the attrition is steep: 32 internal controls exist, 28
can be scored, and **16** find a partner inside the matching caliper.

## Why the additional pieces were added

**Two diagnostic control sets**, because 16 pairs cannot support a null. They buy
statistical power at the cost of a confound each, deliberately chosen so the
confounds are orthogonal. In the event they disagreed with each other and with
the primary comparison, so they now serve mainly as a caution.

**A frozen configuration**, written before any labelled contrast was computed,
fixing the score definition, the model checkpoint, the seeds, the rule for
estimating the reference scale, and the equivalence margin. Choosing a margin
after seeing the differences is not the same test. Two amendments are recorded
against it, both because something had to be corrected after results had been
seen — which is exactly when a written record is worth having.

**An equivalence margin of ±0.2 standard deviations**, because the expected
answer is "no difference", and an ordinary significance test cannot deliver that.
It is an exploratory statistical threshold, not a biologically validated one.

**A cluster bootstrap over connected components**, because two separate
dependencies run through these contrasts: occupied sites in the same ortholog
cluster are near copies, and one control protein can serve several occupied
cases. Resampling on either alone leaves the other unhandled. One contrast per
occupied site, never one per matched row.

**A blinded convergence check**, because the model's conditional probabilities
depend on a sampled decoding order. Eight orders were adopted only after showing
the 8-versus-16 difference was negligible on 50 sites chosen without reference to
their labels.

**A 200-seed matching sweep**, added after the deterministic matcher revealed how
much the old answer had depended on its seed. This is the piece that changed the
conclusion most.

**Design retention**, because probability and behaviour are different things, and
the earlier scoping analysis was about behaviour. It also provides the bridge: if the
conditional score predicts retention, the two analyses describe one underlying
quantity.

## The test, and where the results sit

**Primary conditional-score test — inconclusive.** Occupied versus internal
controls, 16 pairs: **+0.458 SD, 95% CI [−0.227, +1.098]**. The interval includes
zero and spans more than five times the equivalence margin.

**But the direction is stable.** Every point estimate under every matching tried
is positive — the deterministic optimum, the earlier greedy seed, and all 200
seeds of the sweep, spanning +0.286 to +0.699 SD. Nothing produced a negative
estimate.

**And the significance is not.** Across the 200-seed sweep the interval excludes
zero in only 38% of cases. Whether this result "reaches significance" was being
decided by an arbitrary choice inside the matching algorithm rather than by the
data — which is why that choice is no longer left to a seed.

**The effect is in magnitude, not consistency.** Occupied sites score higher in
just 9 of 16 pairs (sign test p = 0.80). The mean is positive because the
negative contrasts are all small while the positive ones are large. A test that
reads only direction finds nothing.

**Reading the interval properly.** It runs from −0.227 to +1.098, so it is
consistent with no difference and with a substantial positive one, and it
excludes only large differences *favouring the controls*. That asymmetry is the
most that 16 pairs support.

**A fourth, better-powered comparison was then added.** The internal-control
class is not being grown, so a eukaryotic secretory set was built that matches
the occupied sites on taxonomy *and* compartment — removing both confounds at
once — and accepts a weaker negative label instead. 262 matched pairs:
**+0.073 SD, CI [−0.056, +0.346]**, the narrowest interval in the study, sitting
essentially on zero with its point estimate inside the equivalence margin.

**Diagnostics disagree with each other.** Bacterial −0.157 SD; cytosolic +0.067
SD. With the primary +0.458 and the parallel +0.073, the four do not form an
interpretable ordering. No interpretation is offered.

**Nothing reaches significance.** Eight tests across four control sets and two
outcomes, using a cluster-level permutation test rather than Wilcoxon because the
pairs are not independent. Smallest raw p 0.030; nothing survives Holm or
Benjamini–Hochberg correction. See [`significance.md`](archive/significance_SUPERSEDED_2026-08-25.md).

**The bridge holds.** Across 2,423 scoreable sites the conditional score predicts
retention at Spearman **+0.547**, monotonically: sites in the lowest score
quintile retain their sequon in 0% of designs, the highest in 30%. The score
measures something real about what the model will do.

**Retention itself is stark.** **81.6%** of sites lose the sequon in every one of
32 designs; overall retention 0.072. The scoping analysis's impression holds up at
site level and at scale, and its 8-design setting proves unbiased — mean 0.0722 against 0.0721 at
32 designs, correlating 0.98. It was noisier per site, not wrong.

These figures are from the completed sweep. An earlier interim pass at 57%
coverage gave 0.078, 80.4% and +0.559 — close enough that nothing turned on it,
but mildly optimistic in every direction, which is the usual reason not to quote
a partial run.

**One earlier claim withdrawn.** A reported retention difference between occupied
and control sites (+0.084, CI [+0.007, +0.176]) used a site-level bootstrap that
ignored clustering, and rested on matched pairs that no longer exist after
rematching. It is not replaced.

## What would settle it

The binding constraint throughout is 16 pairs, and the corrections tightened it
rather than loosening it — the earlier count included pairs it should not have
had. Notably, the pair count is 16 under *every* matching, so the caliper rather
than the algorithm is what limits it: 12 of the 28 scoreable controls have no
admissible partner.

Growing the internal-control class is the single change that would make this
decisive, and the realistic route is occupancy glycoproteomics — a PNGase F
digest in heavy-oxygen water converts occupied asparagines to labelled aspartate,
so a sequon peptide detected with the asparagine intact is a genuine, quantified
negative.

That is data acquisition rather than analysis. It is also why no second model has
been run: ESM-IF or TriFlow on 16 pairs would produce several imprecise answers
instead of one, and none of them would address the constraint.
