# Annotated walkthrough — 19 August 2026

Laura's own account of the week, with the approach and evidence attached to each
point, errors corrected, and every open question answered. Her text is
quoted; annotations follow each block.

---

## Part 1 — the account

> As I understand it: the ProteinMPNN evaluation glycosylation occupancy
> experiment was intended to make the evidence baseline clearer around how models
> are currently sitting on the spectrum of glycosylation relevance awareness.

**Correct.** That framing is the one the module is built around. The design
question is deliberately narrow: not "does ProteinMPNN understand glycosylation"
but "does it treat a sequon that is *actually occupied* differently from one that
merely matches the motif". A null answer is informative because it establishes a
baseline that glycan-aware modelling would have to beat.

> The idea is to take a bunch of proteins that we have experimental evidence for
> having glycans attached at those sequons, and then look at proteins that have
> evidence of the sequence motif being there but not a glycosylation site. Then
> do an analysis with ProteinMPNN of the way the sequons are recognised/treated
> distinctly between them, or not. Both in the log probability for the sequon
> region and in the analysis of how the model designs the proteins given the
> backbones.

**Correct, and this is exactly the two-branch structure.**

| Branch | What it measures | Where |
|---|---|---|
| Conditional score | log-odds the model assigns to N at position 1 and S/T at position 3, given the backbone and the rest of the native sequence | `pipeline/07_score.py` |
| Design retention | whether the sequon survives in 32 unconstrained redesigns of the same backbone | `pipeline/08_design.py` |

They answer different questions — what the model *believes* versus what it
*writes* — and are analysed separately. Figure 4 shows they track each other
(Spearman +0.547), which is what licenses treating them as two views of one
quantity rather than two unrelated measurements.

> No clean one-to-one comparison: the controls matter here because it's tricky to
> get a clean set of negative controls. The proteins that we know for certain are
> not glycosylated also differ from the experimental set in the types of proteins
> they are, cytosolic vs membrane or secreted, and this means they will be
> evaluated differently at a baseline.

**Correct, and quantified.** Before matching, the cytosolic controls differed
from occupied sites by **+0.510 SD** in relative solvent accessibility and the
bacterial by **+0.571 SD**. Both are large — anything above ~0.1 is normally
considered a problem. Structure-based models see exposure directly, so an
unmatched comparison would largely restate that occupied sites sit in exposed
loops. Balance reports: `results/matching/matching_balance_*.json`.

> One approach is to take proteins that are secreted and membrane proteins but
> from species that we know do not have glycosylation mechanisms, so largely from
> bacteria… However, this introduces the second confounding factor… The idea is
> to take both… and use them in combination as the negative controls. Using both
> might help with the confounding factors.

**Correct as a description of the original design, and it is worth knowing that
this reasoning did not survive contact with the corrected data.**

The triangulation logic was: if the confounds are orthogonal, a signal surviving
both is unlikely to be caused by either. That required the three comparisons to
form an interpretable pattern. Before the scoring bug was fixed they did —
cytosolic −0.237, bacterial −0.145, internal −0.057, shrinking neatly as matching
improved. After the fix they point in different directions with no ordering. The
argument was sound; the numbers under it were wrong. See
`docs/concepts.md` § "the logic behind the controls did not survive".

The clade exclusions are worth noting as a point of rigour: bacteria were not
excluded by *annotation* but by **known machinery** — Archaea wholesale (AglB),
plus *Campylobacter* (PglB), *Helicobacter*, *Haemophilus*, *Actinobacillus*,
*Yersinia*, *Kingella* (HMW1C-type, which acts on N-X-S/T). Excluding them for
"not being annotated as glycosylated" would have repeated the exact error the
project exists to avoid.

> The other avenue is to have a higher threshold for uncertainty about the
> annotation reliability… In order to have a true negative control with a similar
> evidence check the sequence motifs would have to have clear "no glycan here"
> annotations, rather than just an absence of a "glycan-here" annotation…

**Correct, and this is precisely the asymmetry the whole design turns on.** A
positive needs one experiment. A negative needs an experiment that *looked and
found nothing*, and almost nothing in the databases records that.

> One thing I thought of for this was to take the proteins that we know have
> really good records for glycosylation, such that they have sequons that are
> properly annotated, and look in the sequences for motifs that do not have the
> annotations and infer that there are no glycans there, but this is still a bit
> messy…

**This is what the 32 internal controls are** — you had already built it. The
implementation is stricter than "well-annotated protein": all four of these must
hold.

1. the residue is resolved in a deposited structure;
2. it carries no glycan linkage;
3. the **same structure** models a glycan at some other residue;
4. the protein was expressed in a host competent to glycosylate.

Conditions 3 and 4 do the work — they establish that sugars survived preparation
and that this depositor was willing and able to model them, so a bare asparagine
is a decision rather than a silence. Code: `structures.py`, `assess_site`.

Your instinct that it would be "small sample size" was right: 32 sites, 28
scoreable, 16 matched.

> …because it is the comparison that is the critical unit of analysis and the
> sample size is rather small, having a lower chance of false negatives is
> something worth considering a lot… Mainly because the result I expect is a lack
> of distinguishing, and I don't want that made messy by the data that was meant
> to be distinguished being more mixed up than required for signal.

**Right instinct, and worth making the direction explicit, because it cuts one
way only.** False negatives — real glycosites sitting in the control set — make
the two groups *more alike*. That pulls any measured difference **toward zero**.

- If you expect and find a null, contamination is a genuine worry: it could
  manufacture the null you expected.
- If you find a difference, contamination cannot have created it — the difference
  survived dilution, so it is if anything understated.

This became directly relevant when the eukaryotic secretory set produced a
positive paired retention difference. Recorded in `docs/negative_controls.md`
*before* that result was computed.

---

## Part 2 — questions and observations

### "I'm confused about what the clusters are referring to here, and what cluster awareness is."

An **ortholog cluster** is the same protein across species — human, mouse and
zebrafish versions of one gene — inherited from the ortholog database that
supplies the sites.

Sites in one cluster are not independent observations. If human and mouse serum
albumin both have an occupied sequon at the equivalent position, that is close to
one fact observed twice.

**Cluster-aware** statistics respect that. Treating 16 correlated observations as
16 independent ones gives an interval that is too narrow and a false sense of
precision. Here the resampling unit is coarser still: a single control protein
can be matched to several occupied cases, so the unit is the *connected
component* of the graph joining occupied clusters to shared control proteins.

That is why the counts step down: **16 contrasts → 13 clusters → 12 units.**

### Your figures for the primary comparison

> occupied protein/ortholog clusters/resampling units: 13/13/12

**Small correction: it is 16 / 13 / 12.** Sixteen distinct occupied proteins (one
per contrast — no protein contributes twice), collapsing to 13 ortholog clusters
and 12 resampling units. Everything else you wrote is right: +0.61 log-odds,
+0.458 SD, CI [−0.227, +1.098], occupied higher in 9 of 16.

> on average ProteinMPNN assigns higher conditional sequon scores to occupied
> sites than to matched internal controls, but the interval is wide and includes
> zero.

**Correct.** And the permutation test agrees: p = 0.251. The direction is real in
the sense that it is positive under every matching we tried; the magnitude is
entirely undetermined.

### "There are six figures…"

**Nine now.** Figures 7–9 were added after you wrote this: retention by class,
the control-set provenance diagram you asked for, and the combined
score-versus-retention panel. They are also not all design-branch — figures 2, 3,
6 and 8 are scoring/design-of-experiment, 4, 5 and 7 are retention, 9 is both.
All explained in `docs/figures.md`.

### "Explain this pairing process… what are the 32 internal controls, what does internal mean here… what is a match here?"

**"Internal"** means they come from *inside the same corpus of structures* as the
occupied sites — same organisms, same secretory compartment, same kind of
experiment. Nothing about them is external to the positive set, which is what
makes them the cleanest comparison available.

**A match** pairs one occupied site with one control site sitting in a similar
local structural environment, so a score difference cannot be blamed on the
environment. Three metrics, all computed from the structure around the
asparagine:

| Metric | What it is |
|---|---|
| Relative solvent accessibility | exposed surface area, divided by the maximum possible for that residue type. 0 = buried, 1 = fully exposed. Shrake–Rupley rolling ball, glycans stripped first |
| Neighbour count within 8 Å | how many residues have an atom within 8 Å — a packing-density measure |
| Hydrophobic fraction within 8 Å | of those neighbours, what proportion are hydrophobic — the chemical character of the pocket |

A pair only forms if the two sites are within a **caliper** of 0.25 pooled
standard deviations across those three dimensions combined. Beyond that, no match
is made rather than a bad one — which is why 28 scoreable controls yield only 16
pairs. Twelve controls have no admissible partner.

### "The bacterial (n=278) and the cytosolic (n=270)"

**Now 280 and 273.** Those sets were still using seeded greedy matching when the
primary switched to deterministic assignment; I propagated the change afterwards,
which shifted both slightly.

### "(diagnostic, what does this mean??)"

Informative about how the measurement behaves, but not a valid answer to the
question. Both sets are confounded on purpose and in opposite directions —
cytosolic differs in compartment, bacterial in kingdom. They were built to be
read against each other, not quoted alone.

The eukaryotic secretory set you commissioned is labelled **parallel** rather
than diagnostic, because it carries neither confound.

### "The range is really large for the internal control, likely because the small sample size?"

**Yes, and that is the whole story of that figure.** 16 pairs, 12 independent
resampling units. Compare the eukaryotic secretory comparison: 262 pairs, 72
units, interval four times narrower.

### "It's unclear to me what the data these are based on actually is — ProteinMPNN scores?"

**Yes.** All four comparisons use the identical conditional sequon score,
computed the same way on every site. The only thing that changes between them is
which control group the occupied sites are compared against.

### "How was the equivalence margin chosen?"

Somewhat arbitrarily, and it is labelled that way throughout — *an exploratory
statistical threshold, not a biologically validated one*. 0.2 standard deviations
is a conventional "small effect" from the effect-size literature; nothing in
glycobiology tells us what a meaningful shift in ProteinMPNN log-odds would be.

What matters more than the number is that it was written into
`config/scoring_frozen.toml` **before any comparison was computed**. Choosing a
margin after seeing the differences is not a test.

### "I don't understand why an ordinary significance test cannot show us the results."

Because the expected answer is a **null**, and an ordinary test cannot deliver
one. It can only reject the null or fail to reject it, and failing to reject is
compatible with both "there is no difference" and "there might be a big
difference but 16 pairs cannot see it". At this sample size the second is very
much live, so a non-significant p-value would say almost nothing.

**Equivalence testing** inverts it: state in advance how small a difference counts
as nothing, then ask whether the confidence interval fits entirely inside that
band. If it does, that is positive evidence of no meaningful effect. If it does
not, you are honestly inconclusive — which is where we are.

### "I don't understand the matching sensitivity — greedy matching, seeds, deterministic."

**Greedy matching** walks the occupied sites one at a time and gives each its
nearest unused control. Greedy because it takes the best option at each step
without considering later steps.

**The seed** sets the walk order. With only 28 controls this matters enormously:
if control C is the only admissible partner for site B, but site A is also near C
and gets processed first, A takes C and B is stranded. Change the order and both
match. So the seed changes both how many pairs form and which.

**Deterministic** matching does not walk at all. It solves the whole assignment
at once — the arrangement giving the most pairs and, among those, the least total
distance. No seed, no ordering, one answer.

The 200-seed sweep exists to show how much the old answer depended on that
arbitrary choice: **all 200 seeds gave a positive estimate, but only 75 produced
an interval excluding zero.** Direction robust, significance decided by a random
number. That is the argument for the deterministic version.

### "I don't understand the quintiles… why this graph matters… where that score comes from."

Take all scoreable sites, sort by conditional score, cut into five equal groups.
For each group, ask what fraction of the 32 designs kept the sequon.

The score is **the same site-level conditional score** as everywhere else — three
residues, not the protein — and the same run, not old data.

**Why it matters:** the two branches measure different things and could easily
have been unrelated. The score is what the model believes while reading the
native sequence; retention is what it writes when generating fresh, which
involves hundreds of interacting decisions. The figure shows they track each
other monotonically (0% retention in the lowest quintile, 30% in the highest,
Spearman +0.547). That means the abstract number is a valid proxy for behaviour,
and the two halves of the project describe one underlying quantity.

### "32 unconstrained designs each (so like 77,000???)"

**Yes — about 77,500 designed sequences.** Designs are generated per chain, and
every sequon on that chain is read from the same 32, so the cost is per chain
rather than per site.

### "What determined the scorability?"

Whether ProteinMPNN will process all three sequon residues, which comes down to
whether the backbone atoms (N, CA, C, O) are all present. Crucially this is
decided **from the coordinates alone, before matching** — see the figure 6 answer
below for why that ordering matters.

### "It looks like ProteinMPNN removes the experimentally validated site as a matter of course — is this true?"

**True for occupied sites, but the 81.6% headline is not only occupied sites.**
That figure is across all 2,423 scoreable sites, which breaks down as:

| Class | n |
|---|---|
| bacterial control | 1,096 |
| cytosolic control | 1,021 |
| **occupied** | **285** |
| eukaryotic secretory control | 255 |
| internal control | 21 |

For occupied sites specifically it is **78.6%** lost in every one of 32 designs,
mean retention 0.080. So yes — hand ProteinMPNN a glycoprotein backbone and the
glycosylation sequons usually disappear.

### "Do we have comparisons for how the model removes the non-glycosylated motifs, and is this different from the background rate of sequence mutation for non-motifs?"

**First question: yes, and it is now the most interesting result.** Paired, each
occupied site against the control matched to it:

| Comparison | Pairs | Occupied | Control | Difference | p (cluster permutation) |
|---|---|---|---|---|---|
| internal control | 16 | 0.160 | 0.008 | +0.152 | 0.125 |
| eukaryotic secretory | 245 | 0.081 | 0.037 | +0.043 | 0.030 |
| bacterial | 254 | 0.079 | 0.076 | +0.002 | 0.816 |
| cytosolic | 251 | 0.078 | 0.075 | +0.003 | 0.821 |

Occupied sequons are retained *more* in both comparisons that hold eukaryotic
secretory context constant, and not at all in the two confounded sets. **But
nothing survives correction for the eight tests run** — smallest corrected p is
0.120. It is a lead, not a finding.

**Second question: we cannot answer it.** The design runner stored only the
classification at the three sequon positions, not the full sequences, so overall
sequence recovery is not recoverable without regenerating everything.
ProteinMPNN's published recovery is roughly 40–50%, which would suggest sequons
are lost far more often than residues in general, but that is a literature
comparison rather than a matched internal one. Storing the sequences would fix it
and costs nothing but disk — it is on the open list.

### "I don't really understand the point of the 6th graph at all… why? How was it fixed?"

**Why the model refuses:** it builds its representation of a residue from the
geometry of its four backbone atoms. If one is missing from the deposited
structure — common at flexible loops and chain ends — the geometry is undefined
and it marks the position unusable.

**The twenty-one ones:** internally the model makes a table with one row per
residue and 21 columns (20 amino acids plus an unknown token). It creates the
table full of **zeros** and writes real values only into rows it decoded. Skipped
rows stay zero. Those stored values are *log* probabilities, so the final step
exponentiates them — and e⁰ = 1. A skipped row therefore comes out as twenty-one
1s. Read as probabilities that says P(asparagine) = 1 and P(serine) +
P(threonine) = 2. A probability of 2 is impossible; that is the tell.

**How it was fixed — two separate things, and you are right that the ordering was
one of them.**

1. **Guards.** The scorer now refuses to score a site if the model reports it did
   not decode that position, *and* separately refuses any row whose
   probabilities do not sum to 1. Two independent checks, so neither is a single
   point of failure. This stops bad values existing.
2. **Ordering.** Scoreability is now settled *before* matching. This is not about
   the correctness of individual scores — it is study design. Previously sites
   were matched into balanced pairs and *then* some were dropped for being
   unscoreable, which unbalanced the sets matching had just carefully balanced.

Conveniently, scoreability depends only on coordinates, so it can be checked
without running the model at all.

### "I do not understand the small paragraph about the corrections meaning that the logic behind the controls did not survive?"

Only that the triangulation argument lost its evidence. Before the fix:
cytosolic −0.237, bacterial −0.145, internal −0.057 — all negative, shrinking as
matching improved, which looked like meaningful support for "the apparent
differences come from compartment and taxonomy". After the fix: +0.062, −0.174,
+0.458 — three directions, no ordering. The reasoning was fine; the numbers it
rested on were wrong.

### "I don't understand what is meant by 'since the corrections, NXS and NXT identical requirement'… I don't have clarity on why."

NXT sequons are glycosylated more often than NXS in nature, and serine and
threonine are chemically different, so the model may score them differently for
reasons having nothing to do with occupancy.

Before this rule, about **45% of pairs matched an occupied NXS against an
unoccupied NXT**. Any difference measured could then have been a subtype
difference wearing occupancy's clothes. Requiring them identical removes that, at
the cost of pairs.

### "Does deterministic matching mean we allow a member of a pair to have multiple matches?"

**No** — still one control per occupied site, each control used once.
"Deterministic" refers only to removing the randomness in *how* pairs are chosen,
not to how many are allowed. A variant permitting up to five controls per case
was run as a sensitivity check and gives a similar answer with fewer independent
cases.

### "I would be interested to see what the results are before the surface exposure matching."

Run:

| | Occupied | Control | Difference |
|---|---|---|---|
| **Unmatched**, all scoreable sites | n=314 | n=28 | +0.213 SD, p = 0.25 |
| **Matched** pairs | n=16 | n=16 | +0.458 SD, CI includes 0 |

Both positive, neither significant. Worth noting the two groups were *already*
similar structurally before matching (SMD 0.11–0.15), so matching bought little
here while costing 314 sites — quite unlike the bacterial and cytosolic sets
where the imbalance was 0.5+. That observation is part of what made the
eukaryotic secretory set worth building.

### "The score averages the log prob of the N at the first position and of the S or T at the third — what does this mean?"

Take a sequon, say **N-I-T**. Show the model the backbone and the rest of the
native sequence. At the first position it outputs a probability for each amino
acid — perhaps 12% asparagine. Convert to **log odds**: log(p/(1−p)), which makes
the scale symmetric around zero so p = 0.5 → 0, p = 0.12 → −2.0, p = 0.88 → +2.0.
Without this, probabilities bunch near zero and averages distort.

Do the same at the third position, but add P(serine) + P(threonine) since either
satisfies the motif. Average the two log-odds values. That is the site's score.

### "If there is a proline there though, the presence is recorded as the sequon no longer being functional."

**Correct in both branches.** In scoring, P(proline at +1) is recorded as a
diagnostic but excluded from the score. In retention, a design that puts proline
at the middle position counts as motif loss even if N and S/T are both intact.

### "I don't understand the cluster bootstrap. Nor the blinded convergence check. Nor the 200-seed matching sweep."

**Cluster bootstrap** — a bootstrap estimates uncertainty by resampling the data
many times and seeing how much the answer moves. A naive one resamples individual
sites; a cluster bootstrap resamples whole groups, so all of a group is in a
given resample or none of it is. That respects the fact that group members move
together.

**Blinded convergence check** — the model's conditional probabilities depend on a
randomly sampled decoding order, so one order is a single draw. We average over
several. The check asked how many are enough: on 50 sites chosen *without looking
at their labels*, the median difference between 8 and 16 orders was 0.0022 SD.
Eight was adopted on that basis. "Blinded" means the subset was picked before
knowing which sites were occupied, so the choice could not be tuned to a result.

**200-seed sweep** — answered above.

### "The diagnostics disagree… which I read as there being no glycosylation awareness, or the data being too messy."

**Both readings are defensible and the significance analysis now discriminates
between them slightly.** Across all eight tests — four control sets, two outcomes
— nothing survives multiple-comparison correction. The best-matched, best-powered
comparison sits essentially on zero (+0.073 SD, interval [−0.056, +0.346], point
estimate *inside* the equivalence margin).

So "no detectable glycosylation awareness" is the supported reading. "Too messy"
is the honest caveat on top of it: no comparison is precise enough to establish
equivalence formally, only to fail to find a difference.

### "Weird that occupied NXT sites score higher than matched controls, while NXS all lower in the three negative. Not in the bacterial set though."

**You have spotted a real and replicated pattern.** Exact figures:

| Comparison | NXS | NXT |
|---|---|---|
| internal control | +0.542 (n=7) | +0.394 (n=9) |
| eukaryotic secretory | **−0.155** (n=126) | **+0.285** (n=136) |
| cytosolic | **−0.094** (n=137) | **+0.230** (n=136) |
| bacterial | −0.135 (n=133) | −0.177 (n=147) |

The split appears in **both eukaryotic** control sets, at n > 125 per cell, and is
**absent from the bacterial** one — exactly as you noticed. Subtype is matched
exactly within pairs, so it is not a matching artefact, and the NXS and NXT
contrasts come from disjoint pairs.

**The thing you may be remembering:** NXT sequons are occupied more efficiently
than NXS — this is well established biochemically, and usually attributed to the
threonine methyl group stabilising the interaction with
oligosaccharyltransferase. So NXT is the "more canonical" glycosylation motif.

Whether that explains our observation is **speculation**, and it should be
labelled as such: one could argue occupied NXT sites are more archetypal and
therefore more recognisable to a structure-based model, while occupied NXS sites
are the unusual ones. The absence in the bacterial comparison would then follow
from bacterial sequons being under no such selection at all. This has not been
tested and is not in any result document as a claim.

### "I am not going to grow the internal control class… [the new set] basically as a parallel analysis."

**Done, and it was the right call.** Full account in
`docs/negative_controls.md` § Set 3 and `docs/primary_result.md`.

> the new negative set inclusion criteria is that it's eukaryotic, has PDB,
> secreted/TM/signal and exclusion is KW-0325

**Exactly as implemented.** 3,619 proteins → 4,418 sequons in 1,543 → 2,296 with
coordinates → 1,946 scoreable → **262 matched pairs**.

> also have we checked that the quality filter is excluding the proteins that have
> UniProt's automated glyco annotation? … I think this is encompassed in KW-0325.

**You were right, and it is verified rather than assumed:** 0 of 1,543 control
proteins carry *any* CARBOHYD feature, manual or automated. UniProt's rule-based
and ML annotations propagate the keyword, so excluding it removes them too.
`pipeline/13_name_audit.py`.

> one check I would like to do is a name check…

**Run.** 15 proteins share a gene symbol with a known glycoprotein — **RNASE1** is
the clearest case, since RNase B is the glycosylated form of the same gene
product, exactly the database-artefact case you described. CALR, TTR, CST3,
PNLIP, GZMB are documented glycoproteins elsewhere. A further 17 sit in families
where N-glycosylation is the rule. 32 suspect proteins in total.

Lectins were deliberately **excluded** from the suspect list: galectins bind
glycans but are cytosolic and unglycosylated, so flagging them would inflate the
count with proteins that correctly belong in the control set.

Only **8 of 262 matched pairs (3.1%)** come from suspect proteins, and removing
them does not move the result: +0.073 → +0.074 SD.

> the matching process led to 262 pairs, pre-matching the imbalance in the RSA was
> +0.554, which is significant meaning that the controls are more buried. After
> matching the RSA imbalance is just +0.004.

**Correct on all counts.** Median RSA 0.314 for the new controls against 0.433 for
occupied. Matching did real work here, unlike for the internal controls.

> wondering what the design pipeline looks like for these right now. Is it doing
> it for a selection of all 5 sets?

At the time you wrote this: **four of five — the new set was not covered.** It has
since been run (262 sites, 233 chains, 100 minutes) and all five classes are now
in `results/designs/`.

> Results came in: for the eukaryotic secretory, 262 sites, the interval is much
> narrower, and the result is still inconclusive, meaning the range is large.

**Correct, with one refinement worth having.** +0.073 SD, CI [−0.056, +0.346].
The interval is four times narrower than the primary's and the **point estimate
sits inside the equivalence margin**; only the upper bound escapes it. So it is
inconclusive by a narrow margin rather than a wide one — roughly twice the pairs
would likely settle it.

---

## Diagrams you asked for

| Diagram | Status |
|---|---|
| ProteinMPNN redesigning an exposed patch — fine normally, quietly aglycosylating a real site | **not built.** The clearest statement of the problem; deferred at your request |
| Dataset basic statistics — protein type, taxa, key stats per set | **not built.** Suggest: taxa distribution, secreted/TM/cytosolic split, sequon density per 1,000 residues, RSA distribution, sequons per protein |
| Scoring process — what scoring a sequon actually involves | **not built** |
| Design process + retention analysis | **not built** |
| Glycosylation awareness spectrum | **not built.** Needs the axis defining first — worth a conversation before drawing |
| Control-set provenance and filtering | **built** — figure 8 |

---

## Corrections to your account, collected

1. **16 / 13 / 12**, not 13/13/12 — sixteen distinct occupied proteins.
2. **Nine figures**, not six, and they are not all design-branch.
3. **Bacterial 280, cytosolic 273**, not 278/270 — deterministic matching shifted both.
4. **The 81.6% headline covers all 2,423 sites**, not only experimentally validated ones; for occupied sites specifically it is 78.6%.
5. **Three comparisons is now four** — the eukaryotic secretory set is labelled *parallel*, not diagnostic, because it carries neither confound.
