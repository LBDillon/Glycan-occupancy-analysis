# What we are doing, and why

*A plain-language account of the reasoning behind the method, written so the
logic can be checked independently of the code. Current as of 2026-08-18, with
the design-retention sweep still running.*

## The question underneath everything

Protein design models such as ProteinMPNN are used to redesign real proteins,
many of which are glycosylated. The models are trained on structures stripped of
glycans, so a natural worry is that they treat a glycosylation site as an
ordinary patch of surface and quietly remove it. The preprint established that
they frequently do remove natural sequons.

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
**observed unmodified** (32), and **unknown** (3,353). Most sites are honestly
unknown, and the design refuses to convert that into a negative.

The 32 exist only because of a narrow argument. A bare asparagine in a crystal
structure normally proves nothing — glycans are routinely trimmed before
crystallisation, expressed away in bacteria, or too mobile to model. But if the
*same structure* models glycans at other residues, and the protein came from a
host that can glycosylate, then sugars demonstrably survived and this
crystallographer demonstrably modelled them. A bare asparagine there is a
decision, not a silence.

## Why there are negative control sets as well

Thirty-two negatives cannot carry a null result. So two further sets of sequons
that *cannot* be N-glycosylated were assembled: **19,337** in cytosolic
eukaryotic proteins, which never enter the secretory pathway and so never meet
the enzyme, and **5,865** in bacterial periplasmic and outer-membrane proteins,
whose clades have no equivalent machinery.

Neither is a clean substitute. The cytosolic set differs from the occupied sites
in subcellular compartment; the bacterial set differs in taxonomy. But their
confounds do not overlap, and that is the point. A model that separates occupied
sites from the cytosolic set but not from the 32 has learned where proteins live.
One that separates them from the bacterial set but not the cytosolic set has
learned taxonomy. **The gap between the comparisons is the measurement**, not any
one of them alone.

## Why matching was necessary before any scoring

Occupied sites are not a random sample of sequons: they sit disproportionately in
exposed loops, because the enzyme has to reach them. Structure-based models can
see exposure. An unmatched comparison of model scores would therefore mostly
restate that occupied sites are exposed, which is already known and says nothing
about glycosylation.

Each occupied site was matched to unoccupied sites of comparable local
environment — accessibility, packing, neighbourhood composition. After matching,
the largest imbalance across every comparison and feature is 0.03 standard
deviations, comfortably inside the conventional 0.1 threshold.

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

## What is running now

The design-retention sweep: unconstrained ProteinMPNN designs on the same
structures, recording at every original sequon whether the motif survives. This
asks what the model *does* when generating, as opposed to what probability it
holds internally — a different question, because sampling combines many residue
decisions at once.

It uses the preprint's own conditions (temperature 0.1, checkpoint v_48_020) with
the sample count raised from 8 to 32, since 8 designs give a per-site error near
0.18 — far too coarse. The first 8 of each run are still a valid
preprint-condition sample, so both are reported from one pass.

---

# The reasoning, in order

## Hypothesis

If ProteinMPNN has learned anything about the biological use of the N-X-S/T
motif, rather than the motif's mere appearance, it should assign a higher
site-level probability to sequons that are actually glycosylated than to matched
sequons that are not.

## The cleanest experiment

Compare **occupied** sites against **observed-unmodified** sites, matched on
local structural context.

This is the cleanest because the two groups share nearly everything that is not
occupancy: same organisms, same subcellular compartment, same kind of experiment,
and — after matching — the same solvent accessibility, packing and neighbourhood
composition. Very little is left that could explain a difference except the
presence of the glycan itself.

Its weakness is size. Only 32 observed-unmodified sites exist, and 22 survive
matching and exact three-residue mapping.

## Why the additional pieces were added

**Two mechanistic control sets**, because 22 pairs cannot support a null. They
buy statistical power at the cost of a confound each, deliberately chosen so the
confounds are orthogonal and can be read against one another.

**A frozen configuration**, written before any labelled contrast was computed,
fixing the score definition, the model checkpoint, the seeds, the rule for
estimating the reference scale, and the equivalence margin. Choosing a margin
after seeing the differences is not the same test.

**An equivalence margin of ±0.2 standard deviations**, because the expected
answer is "no difference", and an ordinary significance test cannot deliver that.
It is an exploratory statistical threshold, not a biologically validated one.

**A cluster bootstrap**, because sites within an ortholog cluster are not
independent and treating them as such would make every interval too narrow. One
contrast per occupied site, never one per matched row.

**A blinded convergence check**, because the model's conditional probabilities
depend on a sampled decoding order. Eight orders were adopted only after showing
that the 8-versus-16 difference was negligible on 50 sites chosen without
reference to their labels.

**Design retention**, because probability and behaviour are different things, and
the preprint's finding is about behaviour. It also provides the bridge: if the
conditional score predicts retention, the two analyses are describing one
underlying quantity.

## The test, and where the results sit

**Primary conditional-score test — inconclusive.** Occupied versus
observed-unmodified: −0.057 SD, 95% CI [−1.073, +0.915]. The interval is about
five times the equivalence margin and contains zero. This does not show the model
treats the two alike; it shows 22 pairs cannot tell. Reported as inconclusive,
not as evidence of no effect.

**Control comparisons — small, negative, imprecise.** Against bacterial controls
−0.145 SD [−0.267, −0.021]; against cytosolic −0.237 SD [−0.383, −0.084]. Both
exclude zero but straddle the margin, so neither establishes equivalence nor a
difference beyond it.

**The pattern across the three is the informative part.** The effect shrinks as
matching improves: cytosolic −0.237, bacterial −0.145, observed-unmodified
−0.057. That is what the orthogonal-confound design was built to detect, and it
is consistent with the control differences arising from compartment and taxonomy
rather than occupancy. Consistent with, not proof of — the best-matched
comparison is also the least precise.

**Every conditional-score estimate is negative.** Nothing supports the model
preferring occupied sequons.

**The bridge holds.** Conditional score predicts retention with Spearman
ρ ≈ +0.50 over 1,186 sites, monotonically: sites in the lowest score quintile
retain their sequon in 0% of designs, the highest in 32%. The score is measuring
something real about what the model will do.

**Retention itself is stark.** 81% of sites lose the sequon in every one of 32
designs; the overall retention rate is 0.076. The preprint's finding replicates
at site level.

**One tension worth watching.** In the 22 matched pairs, retention is *higher*
for occupied sites: +0.084, 95% CI [+0.007, +0.176], excluding zero. That points
the opposite way to the conditional score. It rests on 8 informative pairs — 14
of 22 are tied at zero retention on both sides — so it is fragile, and it is a
secondary outcome with no pre-specified margin. It should be treated as something
to test properly, not as a finding.

## What would settle it

The binding constraint throughout is 22 pairs. Growing the observed-unmodified
class is the single change that would make the primary comparison decisive, and
the realistic route is occupancy glycoproteomics — a PNGase F digest in
heavy-oxygen water converts occupied asparagines to labelled aspartate, so a
sequon peptide detected with the asparagine intact is a genuine, quantified
negative. That is data acquisition rather than analysis, and it is what the next
stage of this work needs.
