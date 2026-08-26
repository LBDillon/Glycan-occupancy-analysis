# What we are doing, and why

*A plain-language account of the whole project — the data, the assumptions, the
maths, and the limits of each conclusion. Written for someone picking the work
back up, including its author.*

*It states almost no figures. Everything measured is maintained in
[`OVERVIEW.md`](OVERVIEW.md), and duplicating numbers here is how six documents
came to assert the same result and go stale together. Terms are defined in
[`glossary.md`](glossary.md).*

*Written so the logic can be checked independently of the code. What was
corrected along the way, and why it kept happening, is at the end.*

---

# Part 1 — the occupancy benchmark

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

## The test, and how to read it

*Figures are in [`OVERVIEW.md`](OVERVIEW.md). What follows is how to read them,
which is the part that does not change when they do.*

**The primary comparison is underpowered by design.** Occupied against internal
controls is the cleanest label available and yields sixteen pairs. Its interval
spans several times the equivalence margin, so it cannot settle anything on its
own — and that was known before it was run, not discovered afterwards.

**Direction was stable while significance was not.** Every point estimate under
every matching tried came out positive — the deterministic optimum, the earlier
greedy seed, and all 200 seeds of a sweep. But across that sweep the interval
excluded zero in only 38% of cases. Whether the result "reached significance" was
being decided by an arbitrary number inside the matching algorithm rather than by
the data, which is why matching is now deterministic and no seed is involved.

**Magnitude and consistency are different claims.** In the primary comparison
occupied sites scored higher in only nine of sixteen pairs. The mean was positive
because the negative contrasts were small and the positive ones large. A test
reading only direction finds nothing there; a test reading magnitude finds
something. Reporting one without the other would mislead in either direction.

**Read the interval, not the verdict.** An interval running from a small negative
to a large positive is consistent with no difference *and* with a substantial
one, and excludes only large differences favouring the controls. That asymmetry
is often the most a small comparison supports, and it is more informative than
"not significant".

**Why a fourth comparison was added.** The internal-control class cannot be
grown — it is limited by how many depositors modelled some glycans and not
others. So a eukaryotic secretory set was built that matches on taxonomy *and*
compartment, removing both confounds at once, and accepts a weaker negative label
in exchange for two orders of magnitude more pairs. It is the best-powered
comparison and the one the headline rests on.

**The diagnostics do not form an ordering.** Before the corrections, the four
comparisons lined up neatly enough to suggest a triangulation argument: if a
signal survives both the compartment confound and the kingdom confound, it is
more likely real. After the corrections they point in different directions with
no pattern. The reasoning was sound; the numbers it rested on were wrong, and the
corrected ones do not support it. No interpretation is offered in its place.

**Significance is tested by permutation, not Wilcoxon.** Wilcoxon assumes the
pairs are independent and they are not — sites in one ortholog cluster are near
copies, and one control protein can serve several occupied cases. The test flips
the sign of every contrast within a whole resample unit at once, so the null
respects the dependency. Wilcoxon p-values are reported alongside and are
systematically smaller; that gap is the dependency, not extra evidence.
Correction spans the **confirmatory** comparisons only: internal control and
eukaryotic secretory, across both outcomes. Bacterial and cytosolic sequons
cannot be occupied in any compartment sense, so a test against them can never
answer the question and does not belong in a family whose purpose is to price
the chances a real answer had to appear.

That change is not free, and not in the direction one might expect. Those two
diagnostics were reliably and strongly significant, so under Benjamini–Hochberg
they occupied the lowest ranks and loosened the threshold for everything above
them — the secretory result read BH 0.021 in a family of eight and reads 0.031
in a family of four. The help was illegitimate: BH controls the false-discovery
proportion, and stuffing a family with guaranteed discoveries inflates its
estimate of how many are real. Under Holm, which has no such mechanism, the
narrower family is straightforwardly better.

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


---

# Part 2 — the glyco-site context analysis

The benchmark measures *whether* a model discriminates. It cannot say *what it
discriminates on*, because "the surroundings" is not a measurement. The second
branch exists to turn that phrase into numbers.

## What was tried first, and why it was abandoned

The obvious version: describe the environments of occupied sites, describe the
environments of unoccupied ones, and report where the distributions differ.

That ran, and the answer was a confound rather than a finding. Occupied sites and
their secretory-unannotated comparison set **share no proteins and no chains at
all**. So every difference between them is also a difference between two
different sets of proteins — fold class, size, expression system, depositor,
resolution. Controlling composition by matching removes nearly all of it.

This is not the contamination problem. Contamination of the unannotated set
biases toward the null and cannot manufacture an effect. Between-protein
composition can, and did.

The within-protein alternative — occupied and unoccupied sequons in the *same*
protein — has no such confound and 31 sites. Neither route justified more
machinery, so the analysis is archived as a negative result with its reasoning
intact rather than built upon.

## What replaced it

A question that needs no negative label at all:

> When a model redesigns a protein with a naturally occupied sequon held fixed,
> does it rebuild an environment like the ones glycosylation actually occurs in?

The comparison is **paired by construction**: same protein, same backbone, same
position, only the sequence changes. The confound that stopped the first version
cannot arise, because there is no comparison group — the natural distribution is
a yardstick, not a control arm.

## What can and cannot move

ProteinMPNN is fixed-backbone. Solvent accessibility, secondary structure,
distance to the termini and disulfide geometry are inherited from the input
structure and **cannot change**. A null on them would be arithmetic, not
evidence, so they are excluded from the outcome and used to describe the
reference instead.

What fixed-backbone design determines is the local chemistry: which residues sit
in the flanking window and in the shell around ND2. That is the whole outcome
panel — fifteen numbers, deliberately small and interpretable.

## The maths, in order

**The panel.** For one site, the fraction of residues in each of seven chemical
classes in two regions: the ±5 sequence window either side of the sequon
(excluding the sequon), and the residues with any heavy atom within 8 Å of ND2.
Plus formal charge in the shell. Fractions rather than counts, because window and
shell sizes vary.

**The distance.** A site's environment is summarised as one number:

    D = median over the panel of | (x - mu_ref) / sigma_ref |

The median rather than the mean, so a single feature with a small reference
spread cannot dominate. The absolute value because either direction is a
departure.

**The reference.** `mu_ref` and `sigma_ref` come from natural occupied sites in
*other proteins* — leave-one-protein-out, not leave-one-row-out. A protein with
several sequons would otherwise contribute to the reference it is scored
against, and its wild type would look more natural than it is.

**The paired quantity.**

    dD = D(design) - D(wild type)

Positive means the design has moved away from natural occupied context. Designs
of one chain are replicates of one draw, so they are averaged within site before
anything is resampled.

**Uncertainty.** Percentile intervals from a bootstrap that resamples
**proteins**, not sites. The tested property is that duplicating a protein's
sites does not narrow the interval, which row resampling would.

**Multiplicity.** Benjamini–Hochberg across the fifteen per-feature tests.

## The two controls, and why both were needed

**The random control** changes the same number of positions as the design, drawing
replacements from the wild-type chain's own amino-acid frequencies. It exists
because a wild type *is* a natural occupied site and therefore sits inside the
reference by construction — so any perturbation moves it outward, and without a
control that arithmetic would be read as a finding. Drawing uniformly over the
twenty residues instead would have destroyed composition trivially and flattered
the model for no reason.

**The composition control** compares the same classes in the flanking window
against every other designable position of the same designed chains. It answers
what the random control structurally cannot: the random control preserves the
overall amino-acid mix and ProteinMPNN does not, so a difference between them
could be a global preference rather than anything local.

Running only the first would have supported a conclusion the second does not.
That is the argument for both.

## A second experiment that needs no labels

The fixed-sequon test asks what happens to a *protected* motif's surroundings.
The complementary question is what happens to an *unprotected* one, and it is
the cleaner question because it needs no occupancy labels at all:

> When unconstrained design is applied to proteins with known occupied sequons,
> how often are those sequons lost relative to the overall mutation rate and to
> comparable non-glycan three-residue motifs?

"The model destroys sequons" is true and uninformative on its own — the model
changes most residues. The control carries the experiment: every three-residue
window in the same chains that does not touch a sequon, measured on the same
designs. With sequence recovery around 44%, three consecutive residues surviving
is unlikely whatever they spell, and comparing sequon loss against the
*per-residue* rate would make an ordinary outcome look alarming.

Two readings are reported because they differ. Retention as the benchmark defines
it is a **pattern** — asparagine, not proline, then serine or threonine — so a
sequon can survive while its residues change. An arbitrary triplet has no such
latitude, so only exact-against-exact is like for like.

The answer is no detectable excess loss. It is worth being careful about what
that licenses. It does **not** show the loss reflects "general redesign rather
than glycan blindness", because those are not alternatives: a glycan-blind model
would be expected to treat the motif as ordinary sequence and lose it at the
ordinary rate, so the ordinary rate is consistent with blindness rather than
evidence against it. What it shows is that the model does not *selectively
protect* the motif — it treats a known biological requirement as ordinary
mutable sequence. And an interval spanning zero is not equivalence: an excess of
several percentage points remains compatible with the data.

## What the result is, and what it is not

Designs drift away from natural occupied context with the motif fully protected,
and further than composition-preserving random change. **But the composition
shift driving it is global to the chain** — ProteinMPNN adds proline and glycine
and removes aromatics everywhere, slightly *less* near the sequon than elsewhere.

So the supported claim is narrow: *fixing three residues does not hold their
surroundings in place, because redesign changes composition everywhere and the
sequon neighbourhood is part of everywhere.* The claim that a model treats
glycosylation sites particularly badly is **not** supported.

For a design tool the practical consequence survives either way — retention as a
metric does not capture what redesign does around a protected site.

Figures: [`../glyco_context/docs/findings_2026-08-26_context_retention.md`](../glyco_context/docs/findings_2026-08-26_context_retention.md).

---

# The assumptions, stated plainly

Each of these could be wrong, and each would change something specific.

**That absence of annotation is not evidence of absence.** The whole design rests
on it. If unannotated sequons were mostly genuinely unoccupied, a far simpler and
better-powered study would be available, and this one is needlessly conservative.

**That a structure's modelled glycans reflect what was there.** The 32 internal
controls assume that a depositor who modelled some glycans would have modelled
this one. Reasonable, not certain.

**That matched pairs are comparable.** Matching uses three variables. Anything
that differs between occupied and control sites and is not correlated with those
three is uncontrolled.

**That the model's conditional probability means something about biology.** It is
a statement about what a network trained on structures expects, not about what a
cell does. Every effect here is a fact about a model.

**That the fifteen-feature panel captures relevant context.** It captures local
chemistry, because that is what fixed-backbone design can change. Glycosylation
also depends on translation rate, trafficking, competition with folding, and
enzyme abundance — none of which appear anywhere in this project.

**That ortholog clusters and shared control proteins capture the dependence.**
The resample unit is built from those two links. Dependence from any other source
— homology below the clustering threshold, shared structural genomics pipelines —
is not accounted for.

# What this project cannot tell you

- **Whether any site is occupied.** Nothing here is an occupancy predictor and no
  score should be reported as one.
- **Whether a designed protein folds, expresses, or is glycosylated.** No
  experimental validation has been done.
- **Whether any feature causes occupancy.** Every association here is
  observational; the causal test is perturbation, which has not been run.
- **Whether models improve with glycan-aware training.** The comparison that
  would answer it is the same model before and after such training, on the same
  frozen pairs. Comparing different architectures confounds training data with
  architecture and scale.

# What the corrections taught, since it keeps recurring

Four defects in nine days, and all four shared a shape: **an assumption about how
a component represents its own input, believed rather than checked.** A token
alphabet taken from the wrong table. Residue indices from one parser used against
another's output. A guard parameter accepted and ignored. None of them crashed;
each produced plausible numbers for the wrong thing.

The countermeasures that worked were not review but assertion — checking that a
decoded sequence really reads N-X-S/T at the index claimed, that a merge has
every shard, that a mapped position matches the other parser residue for residue.
Each of those now fails loudly instead of proceeding.

The one that would have caught all four earliest is the cheapest: **read the
value back out and check it is what you asked for.**
