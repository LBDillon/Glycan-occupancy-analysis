# Why we are building a glyco-site context analysis

*Written 2026-08-23, at the point where the occupancy benchmark finished and the
next piece began. It explains where this work sits in the wider project, what
the benchmark could and could not answer, and why the answer to that gap is a
description of natural glycosylation contexts rather than another model.*

---

## The wider question

Protein design models are now good enough to be used. They are not obviously
good enough to be trusted with glycoproteins, and nobody has established which
of those two statements is closer to the truth.

The concern is specific. N-linked glycosylation is not an optional decoration:
it governs folding, secretion, stability, half-life and immune recognition. A
design tool that silently removes a glycan, or adds a sequon that never gets
used, has changed the protein in ways its scoring function never accounted for.
Whether current models avoid that — and whether they could be made to — is the
programme this work belongs to.

That question is too big to attack directly, so it decomposes into four stages:

1. **Observe what the models have learned.** Without correcting them, without
   assuming they are wrong, measure what they currently distinguish.
2. **Describe the biology they ought to represent.** Build the reference against
   which model behaviour can be judged.
3. **Determine whether the missing biology matters.** A model can fail a
   benchmark without that failure causing design errors.
4. **Test whether grounding the design process helps.** Interventions, then
   wet-lab.

Stage 1 has just finished. This document is about why stage 2 follows from it,
and why it is not merely dataset construction.

---

## What stage 1 asked, and why the design is unusual

The obvious experiment — "do models know about glycosylation?" — is not
answerable, because it conflates several very different capabilities. A model
might recognise that N-X-S/T is a motif worth preserving while knowing nothing
about whether any particular instance is used. Those are different claims and
they need different tests.

So the benchmark asked the narrowest useful version:

> Does a model treat an experimentally **occupied** sequon differently from a
> structurally matched sequon that carries **no glycan**?

The design has one feature that makes it hard and one that makes it honest.

**Both arms carry the motif.** Every site in the comparison is an N-X-S/T
sequon. The occupied ones have experimental evidence of a glycan; the controls
are matched sequons without it. This means motif recognition cannot produce a
difference — a model that has merely learned "N-X-S/T is special" scores both
arms identically. Anything that separates them has to come from somewhere else.

**Absence of annotation is never treated as absence of glycan.** A site with no
evidence is `unknown`, not negative, because well-studied proteins accumulate
annotations and obscure ones do not. Treating unannotated sites as negatives
would measure curation effort and report it as glycobiology. The only sites
called controls are those where a structure models glycans elsewhere but not
here — so sugars demonstrably survived preparation and this depositor
demonstrably modelled them.

That honesty is expensive: it leaves only 32 such internal controls, which is why
a larger, weaker-labelled eukaryotic-secretory set carries most of the
statistical weight.

## What stage 1 found

Three models, two outcomes, four comparisons, corrected across the whole family:

- **ESM-IF distinguishes occupied sequons from matched controls** on both what it
  scores and what it writes when redesigning. Both survive correction.
- **ProteinMPNN distinguishes them too, on the conditional score**, about 1.5x
  more weakly than ESM-IF. On design retention it does not survive correction.
  *(Figures are maintained in [`../../docs/OVERVIEW.md`](../../docs/OVERVIEW.md)
  and deliberately not repeated here.)* *This replaces the earlier reading that ProteinMPNN did
  not discriminate on either outcome, which rested on a sequon-indexing defect
  corrected on 25/08 — see
  [`correction_2026-08-25_sequon_indexing.md`](../../docs/correction_2026-08-25_sequon_indexing.md).*
- **A sequence-only model (ESMC) also distinguishes them**, which means the
  signal is not purely structural.
- **The two diagnostic control sets — bacterial and cytosolic — produce large,
  significant effects in contradictory directions across models.** They are
  measuring composition, not glycosylation, and cannot support a biological
  claim.

Then a control sharpened it. If you hide the whole sequon and ask the model to
judge the site from its surroundings alone, the sequence-only model's
discrimination **vanishes** — it was judging how well the motif fits its context,
not reading the context by itself. The structure-conditioned model's
discrimination **does not change**, because it barely had any to lose.

## The question stage 1 cannot answer

The benchmark measures *whether* models discriminate. It cannot say *what they
discriminate on*.

This matters more than it sounds. We now know ESM-IF separates occupied from
unoccupied sequons that share an identical motif, so it must be responding to the
surrounding sequence and structure. But "the surroundings" is not a measurement.
Until we can say which properties of a site's environment differ between occupied
and unoccupied sequons, three things stay out of reach:

**We cannot explain the model difference.** ESM-IF and ProteinMPNN both read a
backbone. One separates the classes and the other does not. Without a description
of the environments involved, that is an observation with no mechanism.

**We cannot evaluate a designed site.** The eventual aim is to ask whether a
designed glycosylation site sits in a plausible environment. Plausible compared
to what? There is currently no distribution to compare against.

**We cannot tell a real signal from a confound.** The bacterial control set fires
hard in every model, which we interpret as composition rather than glycosylation.
That interpretation is currently an inference from the pattern of results. With a
description of the contexts, it becomes something measurable.

## What the context analysis is

A description of the sequence and structural environments in which glycosylation
actually occurs, built from sites with experimental evidence, and expressed so
that an arbitrary site — natural or designed — can be located within it.

**What was originally planned** was a per-feature profile plus a density score,
across three populations, so that any site could be located within the natural
distribution and so that a feature could be shown to be *discriminative* rather
than merely typical.

**What that turned into.** The discriminative half was attempted and archived.
Comparing occupied sites against unoccupied ones needs a negative label that does
not exist, and the two populations share no proteins and no chains — so occupancy
is confounded with protein identity, and once composition is matched away nearly
every difference disappears. The account is in
[`../archive/comparative_analysis/`](../archive/comparative_analysis/README.md);
it is a real negative result and nothing is built on it.

**What is being done instead** keeps the descriptive half, which was never in
doubt, and pairs it with a question that does not need a negative label at all:
take a natural occupied site, redesign the protein around it with the sequon
held fixed, and ask whether the environment the model builds still resembles the
ones glycosylation actually occurs in. Same protein, same backbone, same
position — the comparison is paired by construction, so the confound that
stopped the comparative analysis cannot arise. Result:
[`findings_2026-08-26_context_retention.md`](findings_2026-08-26_context_retention.md).

## Why now, and why in this order

The first analysis is deliberately **explanatory, not predictive**.

The tempting move is to train a classifier: occupied versus not, report an AUROC,
call it an occupancy predictor. That would be premature and probably wrong. The
control sites are unannotated rather than known-negative, and roughly half of
eukaryotic secretory proteins with structures carry a glycoprotein keyword — so
the negative class is contaminated by construction. A classifier would learn some
mixture of real biology and annotation patterns, and report the total as
accuracy.

So the order is:

1. **Describe** the natural distribution, stratified by evidence strength and by
   NXS/NXT, and report where structural data is missing.
2. **Explain** — on the frozen matched pairs, test whether within-pair feature
   differences account for the model score differences. This needs no new labels,
   reuses matching that is already fixed, and attacks the sharpest open question:
   what did ESM-IF learn that ProteinMPNN did not?
3. **Discriminate, cautiously** — with cluster-held-out validation, ask whether
   simple interpretable features separate occupied from unannotated-secretory
   sites, and describe it as exactly that rather than as occupancy prediction.

Step 2 is the one that makes this more than dataset construction. It turns the
work from a description into an instrument for interrogating the models.

## What this cannot establish

The honest limit is the label. Until there is a resource recording, per site and
per experiment, whether a glycopeptide was observed in modified or unmodified
form, "unoccupied" means "not reported as occupied". That blocks any claim to
predict occupancy.

It does not block this work, and it does not block asking what models respond to
— both of those are questions about the contexts of sites that *are* known to be
occupied, and about which features explain model behaviour. Those are answerable
now. The occupancy predictor is not, and should not be claimed.

## Where it leads

ESM-IF's advantage over ProteinMPNN is now much smaller than it looked — both
discriminate, ESM-IF roughly 1.5x more strongly — so "explain the architectural
gap" is no longer the obvious next question. What remains formulable, and more useful, is
the perturbation test: change specific context features while preserving the
sequon and see whether a model's preference moves. That is a causal test of what
a model uses rather than an associational one, and it needs the feature
vocabulary to exist first.

Further out, the comparison that matters for the programme is the same model
before and after glycosylation-specific training, evaluated on the same frozen
pairs. Comparing different architectures confounds training data with
architecture and scale; comparing one model with itself does not. Stage 1 set the
bar — ESM-IF's effect is what a glycan-aware model has to beat. Stage 2 supplies
the vocabulary to say *how* it beats it, rather than merely *that* it does.

---

### Related

- [`OVERVIEW.md`](../../docs/OVERVIEW.md) — current results and what they support
- [`adding_models_explainer.md`](../../docs/models.md) — how the benchmark
  came to hold three models, and what broke on the way
- [`negative_controls.md`](../../docs/control_sets.md) — what evidence stands behind
  each control set
