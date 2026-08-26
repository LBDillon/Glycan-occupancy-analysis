# Findings — is a sequon lost more often than anything else?

*2026-08-26. Needs no occupancy labels, so it is the one question here that does
not depend on the weakest part of the dataset.*

## The question

> When unconstrained design is applied to proteins with known occupied sequons,
> how often are those sequons lost **relative to** the overall mutation rate and
> to comparable non-glycan three-residue motifs?

"ProteinMPNN destroys sequons" is true and, on its own, uninformative:
ProteinMPNN changes most residues. The question is whether sequons are lost
*more* than the chain around them.

## The answer

**No detectable excess loss. ProteinMPNN shows no evidence of selectively
protecting occupied sequons during unconstrained redesign, and their high loss
rate is broadly consistent with its overall sequence-recovery rate.**

318 occupied sites, 220 proteins, 32 unconstrained designs each, nothing held
fixed. Intervals resample proteins.

| Quantity | Mean | 95% CI |
|---|---|---|
| Sequon retained, **exact triplet** | 7.6% | [4.9, 10.4] |
| Control triplet retained, exact | 9.5% | [9.0, 10.0] |
| **control − sequon** | **1.8 pp** | **[−1.0, +4.6]** |
| Sequon retained, N-X-S/T pattern | 13.6% | [10.2, 16.8] |
| Background mutation rate | 55.7% | [54.7, 56.8] |

The control is every three-residue window in the same chains that does not touch
a sequon, measured on the same designs. The comparison that matters is exact
against exact: **7.6% against 9.5%, a difference of 1.8 percentage points whose
interval includes zero.**

## Why the control is the whole experiment

A sequon is three residues. With 44% sequence recovery, three consecutive
residues surviving is unlikely for reasons that have nothing to do with
glycosylation — the control shows it happens 9.5% of the time. Comparing sequon
loss against the *per-residue* rate would have made an ordinary outcome look
alarming.

Both readings are given because they differ. **Pattern** retention (13.6%) is how
the benchmark defines retention: asparagine, not proline, then serine or
threonine — so a sequon can survive while residues change. An arbitrary triplet
has no such latitude, so pattern-against-exact is not a fair comparison and only
exact-against-exact is.

## Two checks

**Against an independent run.** The ARC unconstrained retention run, different
hardware and code path, gives 13.0% pattern retention against this run's 13.6%,
correlated 0.993 across 283 shared sites.

**Against residue composition.** Sequon residues are not unusually vulnerable:
N 44.6%, S 35.3%, T 45.2%, mean 41.7% against 45.5% across all residues — well
inside the spread from proline at 82% to glutamine at 11%. Multiplying the
per-residue rates as though independent predicts about 8% exact-triplet
retention, close to the observed 7.6%. There is no motif-level effect beyond the
residues the motif happens to contain.

## What this settles, and a correction to how it was first framed

This was initially written up as showing that sequon loss reflects "general
redesign rather than glycan blindness". **That is a false dichotomy and the
conclusion does not follow.** A glycan-blind model would be *expected* to treat
the sequon as ordinary mutable sequence and lose it at about the ordinary rate.
Finding the ordinary rate is therefore consistent with blindness, not evidence
against it. The two are not alternatives.

What the result does establish is narrower and still useful:

> ProteinMPNN shows no evidence of **selectively protecting** occupied sequons.
> It treats a known biological requirement as ordinary mutable sequence.

That is precisely the concern the project started from, now measured rather than
assumed. Without an explicit constraint the model has no reason to keep the
motif, and it does not.

It also agrees with the fixed-sequon composition control, which found the
compositional changes around a *protected* sequon to be global to the chain
rather than local to the site. Neither analysis finds glycosylation sites being
treated differently from anything else — and neither can distinguish "the model
knows nothing about glycosylation" from "the model knows something but is not
acting on it during unconstrained design".

The practical consequence is unchanged. Sequons are lost in roughly 92% of
designs, so keeping them requires instructing the model to.

## What it does not establish

- **Not that the designs would be glycosylated** if the sequon survived. This
  measures motif survival, nothing downstream.
- **Not that other models behave the same.** Only ProteinMPNN was run here.
- **Not that the control is composition-matched.** Control triplets are arbitrary
  windows, not matched on residue composition. The per-residue check above makes
  a strong confound unlikely, and a matched control would sharpen a null rather
  than overturn it.
- **Not that the rates are equal.** The interval on control minus sequon runs to
  **+4.6 percentage points**, so a real excess loss of that size is entirely
  compatible with these data. "No detectable excess" is the claim; "the same" is
  not.
- **Not a complete answer about the functional motif.** Exact-against-exact is
  the fair *sequence-recovery* comparison, and it is the one reported, but the
  biological requirement is the permissive N-X-S/T pattern. An arbitrary triplet
  has no equivalent latitude, so no clean control exists for the pattern reading,
  and whether the functional motif receives special treatment is not fully
  settled by this design.
