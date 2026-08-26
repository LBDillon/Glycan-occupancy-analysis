# The figures, in plain language


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
