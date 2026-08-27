# Figures and captions

Every figure in the project, what it shows and what it does not. Merges the
former `figures.md` and `glyco_context/docs/figure_captions.md`.

Context figures carry titles, axes, legends and significance markers only —
everything explanatory lives here, so a figure reads quickly and is understood
properly. Terms are defined in [`glossary.md`](glossary.md); current results in
[`OVERVIEW.md`](OVERVIEW.md).

> **⚠ Some benchmark figures below quote numbers that predate the 2026-08-20
> alphabet and 2026-08-25 sequon-indexing corrections.** What each figure *shows*
> is unchanged; the quantities in the surrounding prose may not be. OVERVIEW is
> the maintained source.

---

# Part 1 — the occupancy benchmark


> **⚠ Numbers below predate the 2026-08-20 alphabet correction.**
> `mpnn_scoring.ALPHABET` held the wrong string, so `p_asn_at_n` was reading
> P(aspartate). Every ProteinMPNN score and every retention figure produced
> before that date is superseded. Scores have since been regenerated; retention
> has not. **The argument and the method here still stand — the specific
> quantities do not.** See
> [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) for
> what changed and [`OVERVIEW.md`](OVERVIEW.md) for current numbers.

Each one answers a single question. Regenerate them all with
`python pipeline/20_figures_summary.py`.

---

## Figure 1 — Where the sites go

![attrition](../results/figures/fig1_attrition.png)

**What it shows:** how many sites survive each stage of the pipeline, for the
two groups being compared.

Start with 922 occupied sites and 32 internal controls. Ask which have a solved
structure: the occupied group drops hard, to 332, but the controls do not drop at
all — they only exist *because* they came from structures. Ask which have a
backbone complete enough for ProteinMPNN to read: 314 and 28. Finally, ask which
can be paired with a partner of similar local environment: 16 each.

**The point:** the occupied group is never the bottleneck. Every constraint in
this study comes from the 32 internal controls, and by the end only half of them
are usable. If you want a better answer, this is the figure that tells you where
to spend effort.

---

## Figure 2 — Four comparisons, and the best-powered sits on zero

![three comparisons](../results/figures/fig2_three_comparisons.png)

**What it shows:** occupied sites compared against each of the four control
sets. The dot is the estimate, the bar is the 95% confidence interval, and the
grey band is the difference agreed in advance to be too small to care about.

The primary comparison (top) sits at +0.458 but on only 16 pairs, so its interval
is enormous. The eukaryotic secretory comparison — the one that matches the
occupied sites on both taxonomy and compartment — supplies 262 pairs and lands at
+0.073, essentially on zero, with an interval four times narrower.

**The point:** the better-powered and better-matched the comparison, the closer
to zero it sits. An earlier version showed all comparisons negative and shrinking
neatly as matching improved, which looked like evidence; that pattern came from a
scoring bug and is withdrawn.

---

## Figure 3 — The seed was deciding the answer

![matching sensitivity](../results/figures/fig3_matching_sensitivity.png)

**What it shows:** what happens when you run the old matching procedure 200
times with different random starting points. Each bar counts how many runs gave
that answer.

Every single one is positive — the whole distribution sits to the right of zero,
between +0.29 and +0.70. But only 75 of the 200 produced a confidence interval
that excluded zero.

**The point:** two things are true at once, and they matter separately. The
*direction* is robust: nothing we tried made occupied sites score lower. The
*statistical significance* is not: whether the result "worked" depended on an
arbitrary random seed inside the matching code. That is why the matching is now
deterministic — the vertical teal line is the single answer it gives, with no
seed involved.

---

## Figure 4 — The score predicts the behaviour

![retention bridge](../results/figures/fig4_retention_bridge.png)

**What it shows:** sites sorted into five equal groups by the probability score
the model assigns them, and then, for each group, how often the model actually
keeps the sequon when asked to redesign the protein.

The lowest-scoring fifth keeps the sequon in 0% of designs. The highest-scoring
fifth keeps it 30% of the time. It rises smoothly all the way across.

**The point:** the abstract number we measure and the concrete thing the model
does are two views of the same quantity. This matters because it links this
analysis to the earlier scoping analysis, which looked at generated sequences
rather than internal probabilities.

---

## Figure 5 — ProteinMPNN destroys sequons routinely

![retention distribution](../results/figures/fig5_retention_distribution.png)

**What it shows:** for each of 2,423 sites, the fraction of 32 unconstrained
designs that kept the sequon intact. Note the vertical axis is a log scale — the
first bar is far taller than it looks.

Over 80% of sites lose the sequon in **every single design**. Only 50 sites, about
2%, keep it every time.

**The point:** this confirms, at the level of individual sites and at scale, what
the earlier scoping analysis suggested from a handful of proteins. If you hand a glycoprotein to ProteinMPNN and ask it to
redesign the surface, you should expect the glycosylation sites to disappear.
Nothing in the model's training told it they mattered.

---

## Figure 6 — The bug that inverted the first answer

![scorer defect](../results/figures/fig6_scorer_defect.png)

**What it shows:** the distribution of every score computed, with the invalid
ones marked in brown. Real scores run from about −5 to +1. The brown cluster
sits at +13.8, which is not a possible value.

The cause: ProteinMPNN silently refuses to process any residue with a missing
backbone atom, and returns zeros instead. Exponentiating zeros gives a row of
ones, which reads as "probability 1 that this is asparagine, and probability 2
that it's serine or threonine" — arithmetic nonsense that nothing was checking
for.

**The point:** only 105 of 2,564 sites were affected, about 4%. But every result
in the study is expressed in units of the spread of these scores, and those few
absurd values nearly doubled that spread, from 1.33 to 2.62. A small data problem
became a large error in the measuring stick, and the first version of the primary
result came out with the wrong sign. Two independent checks now reject anything
that is not a genuine probability distribution.

---

## Figure 7 — Retention differs between control sets, so a class average is not the test

![retention by class](../results/figures/fig7_retention_by_class.png)

**What it shows:** how often ProteinMPNN keeps the sequon, split by whether the
site is a real glycosylation site or a control. Error bars are bootstrapped over
whole proteins rather than sites, because several sequons on one protein share
one set of designs and are not independent.

Occupied sites keep the sequon 8.0% of the time. Bacterial 7.7%, cytosolic 6.6%
— but the eukaryotic secretory set only 3.6%.

**The point:** the control sets differ from each other, so a class average cannot
answer the question. Whatever separates eukaryotic secretory proteins from
bacterial and cytosolic ones is riding along in this comparison.

Because every control was matched site-by-site to an occupied site, the paired
contrast is available and is the right test — see Figure 9 and
[`significance.md`](archive/significance_SUPERSEDED_2026-08-25.md). An earlier version of this figure was
titled "occupied sequons are destroyed at the same rate as unoccupied ones",
which was drawn from these class averages before a taxonomy- and
compartment-matched control set existed. That claim has been withdrawn.

---

## Figure 8 — How the control sets are built, and what each costs

![control provenance](../results/figures/fig8_control_provenance.png)

**What it shows:** left, every filter applied to build the eukaryotic secretory
control set, from 7,423 candidate proteins down to 262 usable matched pairs.
Right, all four control sets placed against the two things that matter — how well
the population matches the occupied sites, and how defensible the negative label
is. Marker area is the number of usable pairs.

**The point:** the internal controls sit in the corner you would want, strong on
both axes — they are not a compromise, they are simply tiny. Every other set buys
two orders of magnitude more pairs by giving up one axis: the cytosolic and
bacterial sets change the population, the eukaryotic secretory set weakens the
label. Nothing large occupies the good corner, and that is the whole problem.

The highlighted bar on the left is the step that constitutes the entire negative
label for the new set — the absence of a UniProt glycoprotein keyword.

---

## Figure 9 — Both measurements, all classes

![all classes](../results/figures/fig9_all_classes.png)

**What it shows:** left, what the model *writes* — design retention for all five
classes. Right, what the model *believes* — the conditional-score difference
against each control set.

**The point:** two different quantities, measured on different scales, and they
agree. Neither shows the model distinguishing occupied sequons by a margin this
design can resolve.

The paired versions of the retention contrast — each occupied site against the
control matched to it, rather than class averages — are in
[`significance.md`](archive/significance_SUPERSEDED_2026-08-25.md). They show occupied sequons retained more
often in both comparisons that hold eukaryotic secretory context constant, and
not at all in the two confounded sets, but nothing that survives correction for
the number of comparisons run.

---

## Draft schematics (10–14)

Five drafts made on request. Only figure 11 carries data; the rest explain what
the pipeline does and what is at stake. **All are drafts — the shapes and
wording are meant to be argued with.**

| Figure | What it shows |
|---|---|
| `fig10_mechanism.png` | Why any of this matters. Three panels: what is really there (an exposed asparagine carrying a glycan), what ProteinMPNN sees (the same backbone, no glycan — they are identical to it), and what it designs (a hydrophobic residue that packs better and deletes the site). Usually a sound design move; here it aglycosylates the protein. Schematic, not a real structure. |
| `fig11_dataset_summary.png` | What is actually in each dataset — size before and after matching, taxonomic spread, and how exposed the sequons are. The taxonomy panel is a **caveat**: the eukaryotic secretory set matches on kingdom and compartment but not species (28% human against occupied's 43%), and species composition is not matched anywhere in the study. |
| `fig12_scoring_process.png` | What scoring a sequon involves, end to end: structure → scoreability check → model reads backbone plus native sequence → probabilities at the three positions → log-odds → one number. |
| `fig13_design_process.png` | What the design branch involves: structure → 32 unconstrained sequences → read the three original positions → classify → retention fraction. Nothing is fixed or biased, so the sequon is free to vanish. |
| `fig14_awareness_spectrum.png` | A proposed scale from glycan-blind to glycan-aware, with what the evidence says about where ProteinMPNN sits and what would move a model along it. **The axis itself is the open question** — worth settling before this is used in a write-up. |

Regenerate with `python pipeline/24_figures_schematics.py`.

---

## Figure 15 — Three models, two measurements

![model comparison](../results/figures/fig_model_comparison.png)

*Written 2026-08-27. **Unlike most of Part 1, every number here is post-correction**
— ProteinMPNN index-corrected, ESM-IF index-corrected, ESMC as run.*

**What it shows:** the whole benchmark in four panels, each carrying one claim.

**A — the scoring effect.** Occupied minus matched control, in each model's own
reference SD. ProteinMPNN +0.282 [+0.114, +0.426]; ESM-IF +0.431 [+0.267, +0.726];
ESMC +0.261 [+0.063, +0.598]. All three intervals exclude zero, which is why all
three carry a marker.

**B — what survives hiding the motif.** Filled dot is the motif visible, open dot
hidden, arrow is the change. ESMC falls +0.261 → −0.113 and crosses zero: its
apparent preference *was* motif recognition. The two structure-conditioned models
barely move — ESM-IF +0.431 → +0.398, ProteinMPNN +0.282 → +0.315.

**C — the design result.** Sequon retention, occupied against matched control.
ProteinMPNN 0.121 vs 0.051; ESM-IF 0.151 vs 0.059. ESMC's blank row is a fact
about the model, not missing data: it has no `SequenceDesigner`.

**D — the per-site contrasts panel A averages.** Standardised so the three are
comparable at all.

**The point:** structure-conditioned models retain their preference when the
motif is hidden and a sequence-only model does not. That is the strongest single
result in the benchmark, because both arms of every pair contain a sequon — so
whatever separates them once the motif is hidden cannot be motif recognition.

**What it does not show, and four cautions:**

- **Panel D undercuts panel A, and should.** After standardisation the three
  distributions are nearly identical, and **35–40% of sites have negative
  contrasts**. These are shifts in broad overlapping distributions, not
  separation. Nothing here would classify a site.
- **Panel B's markers refer to the change, not the effect** — and all three
  changes are significant, including ProteinMPNN's, which moves the *wrong way*:
  hiding the motif slightly *improved* its contrast. There is no account of why.
  "Robust to masking" and "a small significant effect in an unexpected direction"
  are different statements, and this is the second.
- **Panel B mixes two scales.** The dots are SD-standardised; the significance
  comes from `15_masking_comparison`, which works in log odds. The change values
  in log odds are ProteinMPNN −0.035 [−0.088, −0.008], ESM-IF +0.043
  [+0.010, +0.073], ESMC +0.558 [+0.447, +0.724].
- **ProteinMPNN's rows are not the same site set as the others** — 232 contrasts
  against 262, and 216 retention pairs against 245. Sound within each model;
  do not rank the models off it.

**Panel C's missing ProteinMPNN marker is a genuine disagreement, not an
oversight.** Its cluster-aware interval is [−0.004, +0.119] and includes zero,
while Wilcoxon gives p=0.0033. The bootstrap respects protein and ortholog
clustering; Wilcoxon does not. Significance is marked from the interval, which is
the conservative and repository-standard reading. Note also that **129 of its 216
pairs are tied** — both sites lose the sequon — so the effect rests on 87
informative pairs; ESM-IF's rests on 93 of 245.

CARBonAra is absent by design: with no motif-visible condition it has no
within-model masking contrast, so it cannot appear in panel B at all.

---

## Figure 16 — What the models actually output

![raw score distributions](../results/figures/fig_score_distributions.png)

*Written 2026-08-27.*

**What it shows:** the raw conditional sequon score — log odds, not standardised
— for experimentally occupied sites against their matched partners with no
annotated glycan. Dashed lines are medians. One panel per model.

| Model | occupied median | no-annotation median | occupied sites below the other group's median |
|---|---|---|---|
| ProteinMPNN | −1.355 | −1.763 | 34.1% |
| ESM-IF | −0.907 | −1.691 | 29.0% |
| ESMC | −0.741 | −1.552 | 29.4% |

**The point:** it is the same result as Figure 15A without the standardisation,
and it makes the size of the effect legible in a way an SD cannot. If the two
groups were identical, half of the occupied sites would fall below the
no-annotation median. About 30% do. That is a real shift and a modest one.

**The three panels must not be read against each other.** Raw log odds are
comparable *within* a model and not between: ProteinMPNN averages eight decoding
orders, ESM-IF runs a single causal pass, ESMC never sees the backbone. The
x-axes differ deliberately, and the panels are three separate comparisons that
happen to share a figure.

**These are the matched pairs, not the populations.** Each case has exactly one
control, so nothing is averaged and no spread is compressed. The unmatched
version of this comparison is the one the archived context analysis showed
collapsing once composition was controlled — occupancy was confounded with
protein identity — so the matched form is the only one worth plotting.

**"No annotated glycan" is not "unoccupied".** These sites are
`control_secretory_eukaryotic_unannotated`: absence of annotation, not evidence
of absence. A separate and much smaller group, `observed_unmodified` (31 sites),
carries real evidence of absence and is not shown here. Reading the grey
distribution as unglycosylated sites is the single easiest mistake to make with
this figure.

---

## What the results figures say together

The dataset is sound and the measurement works — Figure 4 shows the score tracks
real model behaviour, and Figure 5 shows the behaviour is dramatic and worth
studying.

**The answer is a null that does not quite close.** Across four control sets and
two outcomes, no comparison shows the model distinguishing occupied sequons by a
margin the design can resolve, and no test survives correction for the number of
comparisons made (Figures 2 and 9, `significance.md`).

**The best-matched, best-powered comparison sits on zero.** The eukaryotic
secretory set matches taxonomy and compartment and gives 262 pairs: +0.073 SD,
interval [−0.056, +0.346]. The point estimate is inside the equivalence margin;
only the upper bound escapes it. Perhaps twice the pairs would settle it.

**One pattern is worth testing properly.** Occupied sequons are retained more
often than their matched partners in both comparisons that hold eukaryotic
secretory context constant, and not at all in the two confounded sets — the
opposite of what confounding would produce. It does not reach significance and
it was found by looking, so it is a hypothesis for new data, not a result.

**The original 16-pair comparison is not wrong, just uninformative.** Figure 1
explains why there are only 16, and Figure 8 why no larger set has an equally
firm negative label.

Figure 6 is the reason to trust the current numbers more than the earlier ones.

---

# Part 2 — the glyco-site context analysis

Figures carry titles, axes, legends and significance markers only. Everything
that explains, qualifies or interprets lives here, so a figure can be read
quickly and understood properly.

Terms are defined in [`glossary.md`](glossary.md).

---

## fig2 — the natural reference

**`glyco_context/results/figures/fig2_occupied_context.png`**

> Local environment of 318 experimentally supported N-linked glycosylation sites
> with usable structures (the `triplet_core` view). **Left:** relative solvent
> accessibility at the asparagine and at the +1 and +2 positions of the sequon.
> **Centre:** DSSP secondary structure, coarsened to loop, sheet and helix.
> **Right:** backbone region at the asparagine, from phi and psi.
>
> The left-handed alpha region holds 18% of sites. This is expected rather than
> anomalous: asparagine and glycine occupy left-handed conformations far more
> readily than other residues, so the category is real and not a rounding
> artefact.
>
> No comparison group appears here. This is the reference picture that the
> design distances in fig4 and fig6 are measured against.

## fig4 — the context-retention result

**`glyco_context/results/figures/fig4_context_retention.png`**

> Fifty occupied sites across 38 proteins, each redesigned 32 times by
> ProteinMPNN with the sequon held fixed. **Top:** distance from natural
> occupied context for the wild type, a mutation-count-matched random control,
> and the designs. **Bottom:** the paired change from wild type. Bars are 95%
> intervals from resampling proteins; markers give the two-sided bootstrap p.
>
> Every design preserves N-X-S/T, verified by reading the triplet back out of
> the designed sequences.
>
> The random control is what makes the top row readable. A wild type is itself a
> natural occupied site, so it sits inside the reference by construction and any
> perturbation moves it outward. The control changes the same number of residues
> at random and so measures how much of the drift is simply that — about half.
>
> Designs move further than the control (+0.039), but that difference is not
> significant at this sample size.

## fig5 — which features move

**`glyco_context/results/figures/fig5_context_features.png`**

> Per-feature shift from wild type to design, in standard deviations of the
> natural reference, for the five features of fifteen whose intervals separate
> from zero. Bars are 95% intervals from resampling proteins; q values are
> Benjamini-Hochberg across all fifteen. Faint bars are weaker evidence.
>
> **Nothing survives correction at q < 0.05.** The three strongest — less
> aromatic content near ND2, more glycine in both the flanking window and the
> shell — are chemically coherent with one another, which is a reason to check
> them in a larger run and not a reason to believe them yet.

## fig6 — feature distributions

**`glyco_context/results/figures/fig6_feature_distributions.png`**

> Empirical cumulative distributions for each of the fifteen panel features.
> **Grey:** all 285 natural occupied sites, the reference. **Black:** the 50
> wild-type sites that were redesigned — a subset of the grey, drawn separately
> as a check that the tested sites are representative. **Blue:** 1,600 designs,
> 32 from each of the 50 sites. Panel titles give the Benjamini-Hochberg q for
> the paired wild-type-to-design shift.
>
> A shift between distributions appears as horizontal separation; a mass of
> identical values appears as a vertical jump.
>
> Cumulative distributions rather than histograms because these features are
> heavily zero-inflated — 91% of natural sites have no cysteine in the ND2 shell
> — so a histogram's zero bar dominates and a kernel density would spread mass
> across values that never occur. An eCDF uses every observation exactly and
> requires no binning choice.
>
> The 1,600 designs carry roughly 50 sites' worth of independent information, so
> the smoothness of the blue curve should not be read as precision.

---

## Archived — the comparative analysis

These belong to a question that was asked, answered negatively and set aside.
See [`../archive/comparative_analysis/README.md`](../glyco_context/archive/comparative_analysis/README.md).

### fig1 — effects collapse under matching

**`glyco_context/archive/comparative_analysis/results/fig1_effects_collapse.png`**

> **Left:** standardised difference between occupied sites and secretory-
> unannotated comparison sites, before and after controlling composition by
> matching. Ringed points are the three matching variables, which are balanced
> by construction and cannot be tested here — a null on them means the matching
> worked. **Right:** the proportion of each comparison set's proteins that also
> carry an occupied site.
>
> The occupied and secretory-unannotated sets share no proteins and no chains,
> so occupancy is confounded with protein identity and every difference is also
> a difference between two sets of proteins. The internal controls largely are
> not, which is why their near-zero estimates carry weight despite n=31.

### fig3 — the same features under four framings

**`glyco_context/archive/comparative_analysis/results/fig3_framings.png`**

> Five features under each of four analyses: population-level and matched-pair,
> against internal controls and against secretory-unannotated sites. Asterisks
> mark survival of Benjamini-Hochberg within that comparison.
>
> The two matched-pair framings control composition; the two population framings
> do not, and their larger values are the confounded ones. Only beta-sheet at
> the +2 position holds its direction throughout.
