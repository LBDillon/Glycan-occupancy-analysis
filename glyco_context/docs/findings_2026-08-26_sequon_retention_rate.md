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

**No. Sequons are lost at about the rate any three-residue motif is lost.**

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

## What this settles

It was an open question whether sequon loss reflects **glycan blindness** or
**general sequence redesign**. On this evidence it is general redesign. A model
that specifically disregarded glycosylation would lose sequons *faster* than
comparable motifs, and this one does not.

This agrees with the fixed-sequon composition control, which found the
compositional changes around a protected sequon to be global to the chain rather
than local to the site. Two independent analyses, the same conclusion: **nothing
about ProteinMPNN's behaviour at glycosylation sites is special.**

The practical consequence is unchanged and arguably strengthened. Sequons are
lost in roughly 92% of designs, and if you want them kept you must say so — no
model is going to preserve them incidentally.

## What it does not establish

- **Not that the designs would be glycosylated** if the sequon survived. This
  measures motif survival, nothing downstream.
- **Not that other models behave the same.** Only ProteinMPNN was run here.
- **Not that the control is composition-matched.** Control triplets are arbitrary
  windows, not matched on residue composition. The per-residue check above makes
  a strong confound unlikely, and a matched control would sharpen a null rather
  than overturn it.
