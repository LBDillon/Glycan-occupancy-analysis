# Pre-specification — the fixed-sequon context-retention test

*Written and committed 2026-08-24, before any design was generated or scored.*

## Amendment, 2026-08-25 — the motivating contrast changed, the question did not

This was written while ProteinMPNN's conditional-score result stood at +0.090 SD
(BH 0.30, inconclusive) against ESM-IF's +0.431 SD. Correcting ProteinMPNN's
sequon indexing moved it to **+0.282 SD (BH 0.021)**, so the sharp architectural
contrast that made "ProteinMPNN preserves the motif but perhaps not its context"
an attractive hypothesis is much weaker: both models discriminate, one about
1.5x more strongly.

**Nothing below changes.** The question — does fixing the motif protect the
biology around it, or only the three letters — does not depend on how the two
models compare. It is a question about what fixed-backbone redesign does to a
site's environment, and it would be worth answering if only one model existed.

What does change is what a result would license. A drift away from natural
context can no longer be read as *the* explanation for a ProteinMPNN-specific
deficit, because there is no longer much of a ProteinMPNN-specific deficit to
explain. It would be a statement about protected-sequon redesign, full stop.

## The question

> When ProteinMPNN is forced to preserve a naturally occupied N-X-S/T sequon,
> does it redesign the surrounding residues into a local environment that is
> unusual for natural occupied sites?

In one sentence: **does fixing the motif protect the biology around it, or only
the three letters?**

This tests one specific failure mode — a model preserving the visible motif
while removing the context that normally accompanies its use. It is not an
occupancy predictor, and no result here shows that ProteinMPNN predicts
occupancy.

## Why this design avoids the problem that stopped the comparative analysis

The comparison is **paired within a site**: same protein, same backbone, same
position, only the sequence changes. There is no comparison group drawn from
other proteins, so occupancy cannot be confounded with protein identity. The
reference distribution is used as a yardstick, not as a control arm.

## What can and cannot move

ProteinMPNN is fixed-backbone. Solvent accessibility, secondary structure,
distance to the C-terminus and disulfide *geometry* are inherited from the input
structure and cannot change. **A null on them would be arithmetic, not
evidence**, so they are excluded from the outcome and used to condition and
stratify the reference instead.

The primary outcome is therefore the local chemical environment, which is what
fixed-backbone sequence design actually determines.

## The panel

**Primary outcome — sequence window (±5 residues around the Asn, the three fixed
sequon positions excluded):** fractions of hydrophobic, aromatic, charged, polar,
glycine, proline and cysteine residues.

**Primary outcome — three-dimensional shell (same-chain residues with any heavy
atom within 8 Å of ND2, the sequon itself excluded):** the same seven fractions,
plus net formal charge (D/E = −1, K/R = +1, histidine neutral).

Fifteen features. Deliberately small and interpretable; no multivariate density
model at this stage.

**Conditioning and stratification, never outcomes:** NXS versus NXT, RSA at the
Asn, secondary structure at the Asn, distance to the C-terminus, disulfide
proximity.

## Sample and design protocol

- **Reference distribution:** the 318 occupied sites in `context_triplet_core`.
- **Test set:** 30–50 proteins carrying supported occupied sites, drawn to
  balance NXS and NXT, and to spread across buried/exposed and loop/sheet/helix
  settings rather than concentrating in one structural class.
- **Design:** ProteinMPNN, temperature 0.1, 32 designs per chain, matching the
  frozen `STANDARD_CONDITION` so this is commensurate with the retention work.
- **Held fixed during design:** the three sequon positions (N, X, S/T); every
  cysteine participating in a detected disulfide.

## Scoring

Each site's environment is summarised as a distance from the natural occupied
distribution:

    D = median over the panel of | (x - mu_ref) / sigma_ref |

`mu_ref` and `sigma_ref` come from occupied sites **in other proteins**
(leave-one-protein-out), so no site contributes to the reference it is scored
against.

**Primary result:** the paired change

    dD = D_design - D_wild_type

Positive means the design has moved away from natural occupied context.

**Secondary:** the proportion of designs whose D exceeds the 95th percentile of
held-out natural sites; and the per-feature contributions, to say *which*
features move rather than only that something did.

**Control:** a mutation-count-matched random redesign — the same positions
changed, replacement residues drawn from the wild-type chain's own amino-acid
composition. Drawing uniformly over the twenty would break composition
trivially and make ProteinMPNN look good for no reason. The control answers
whether ProteinMPNN preserves context better than arbitrary alteration of the
same size.

## Uncertainty and multiplicity

Cluster bootstrap over **proteins**, 2,000 replicates, percentile intervals.
Benjamini–Hochberg across the fifteen per-feature tests. Designs within a chain
are replicates of one draw, not independent observations, and are averaged
within site before any protein-level resampling.

## What each outcome would mean

| Result | Reading |
|---|---|
| dD ≈ 0 and below the random control | ProteinMPNN preserves coarse local context when the sequon is protected |
| dD > 0 | protecting N-X-S/T is not sufficient; contextual constraints are needed |
| only some features move | those become specific candidates for constraints and for experimental test |
| sites differ by structural class | the result is which settings are vulnerable, not a single global verdict |

## Bounded conclusion

Whatever the outcome, the claim is about whether a natural-like local chemical
environment survives protected-sequon redesign. It is not a claim about
occupancy, about whether the designed protein folds, or about whether any
feature causes glycosylation.
