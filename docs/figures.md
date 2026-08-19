# The six figures, in plain language

Each one answers a single question. Regenerate them all with
`python runners/summary_figures.py`.

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

## Figure 2 — The three comparisons disagree

![three comparisons](../results/figures/fig2_three_comparisons.png)

**What it shows:** the result of comparing occupied sites against each of the
three control sets. The dot is the estimate, the bar is the 95% confidence
interval, and the grey band is the difference we agreed in advance would be too
small to care about.

The primary comparison (top, teal) sits at +0.458 — occupied sites score higher —
but its interval crosses zero, so we cannot rule out no difference. The two
diagnostic comparisons point elsewhere: bacterial slightly negative, cytosolic
essentially nothing.

**The point:** these three do not tell a consistent story. An earlier version of
this analysis showed all three negative and shrinking neatly as the matching got
better, which looked like meaningful evidence. That pattern came from a bug. It
is gone, and nothing has replaced it — the diagnostics simply do not corroborate
the primary result.

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
analysis to the preprint, which was about generated sequences rather than
internal probabilities.

---

## Figure 5 — ProteinMPNN destroys sequons routinely

![retention distribution](../results/figures/fig5_retention_distribution.png)

**What it shows:** for each of 2,423 sites, the fraction of 32 unconstrained
designs that kept the sequon intact. Note the vertical axis is a log scale — the
first bar is far taller than it looks.

Over 80% of sites lose the sequon in **every single design**. Only 50 sites, about
2%, keep it every time.

**The point:** this replicates the preprint's central finding at the level of
individual sites. If you hand a glycoprotein to ProteinMPNN and ask it to
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

## What the six say together

The dataset is sound and the measurement works — Figure 4 shows the score tracks
real model behaviour, and Figure 5 shows the behaviour is dramatic and worth
studying.

The answer to the actual question is weak. Occupied sites do tend to score higher
(Figures 2 and 3), and that direction survives everything we threw at it. But 16
pairs cannot make it precise, and Figure 1 explains exactly why there are only 16.

Figure 6 is the reason to trust the current numbers more than the earlier ones.
