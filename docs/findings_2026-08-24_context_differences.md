# Findings — occupancy-associated context differences

*2026-08-24. Covers Step 2 (description), Step 3 (population comparison) and
Amendment 1 (matched pairs). Supersedes any earlier reading of the Step 3
output, including the interpretation embedded in commit 8bf84fb.*

## Headline

Once composition is controlled by matching, **almost every context difference
between occupied and unoccupied sequons disappears.** In the well-powered
comparison the only feature surviving correction is where the site sits in the
chain, not any property of its local environment.

This is a negative result about the features measured, not about glycosylation.
It says these seventeen properties do not distinguish an occupied sequon from a
structurally matched unoccupied one — which is informative, because several of
them looked like strong effects before matching.

## What occupied sites look like (Step 2, descriptive)

318 sites across 233 proteins. The asparagine is solvent-exposed (median RSA
0.43), in a loop 64% of the time, with about 41 protein atoms within 8 Å of ND2.
Backbone regions are β 48%, α_R 33%, α_L 18% — the last is high but expected,
since asparagine occupies left-handed conformations far more readily than most
residues. NXT and NXS are near-balanced at 52/48.

## Population comparison (Step 3) and why it is not evidence

Against secretory-unannotated sites, sixteen features cleared correction, several
strongly: exposure +0.57 SD, ND2 crowding −0.59, nearest aromatic +0.39.

**These are not occupancy effects.** Occupied and secretory-unannotated sites
share **zero proteins and zero chains**, so occupancy is entirely confounded with
protein identity. The comparison measures the difference between two sets of
proteins.

Note this is not the label-contamination problem. Contamination of the
unannotated set biases toward the null and cannot manufacture an effect.
Between-protein composition can, and did.

## Matched pairs (Amendment 1) — the designed test

Each occupied site against the control it was matched to on RSA, neighbour count
and hydrophobic fraction. 261 secretory pairs over 72 resample units; 16 internal
pairs over 12.

Balance is good in the secretory pairs: RSA differs by 0.000, neighbour count by
−0.004. **Hydrophobic fraction is imbalanced in the internal pairs** (+0.076,
95% CI [0.027, 0.136]), which is worth knowing independently of anything here.

In the secretory comparison, one feature survives correction:

| Feature | Mean difference | Standardised | 95% CI | q |
|---|---|---|---|---|
| `distance_to_n_terminus_resolved` | −64.7 residues | −0.27 | [−107, −40] | **0.012** |

Occupied sequons sit closer to the N-terminus of the resolved chain than their
matched controls. Everything else is null, including every feature that looked
strong before matching — ND2 crowding falls from −0.59 to −0.08, nearest aromatic
from +0.39 to +0.10, helix occupancy at the asparagine from −0.16 to −0.10.

The internal comparison produces several apparently strong effects — sequence
position, aromatic presence, loop length, all with |standardised| ≥ 0.5 — but at
16 pairs over 12 dependency units, with a known imbalance in one matching
variable. **They are not reliable and are recorded as hypothesis-generating.**

## The one feature that has held up across every framing

`plus2_ss_coarse == sheet` — β-sheet at the +2 position — is positive in all
four analyses: population secretory (+0.15, q < 0.01), population internal
(+0.21, q = 0.027), matched secretory (+0.10, q = 0.32), matched internal
(+0.31, q = 0.16). It was the only feature meeting the pre-specified
"significant in both comparisons, same direction" bar, and matching weakens it
without reversing it. It is the strongest candidate the data currently offer and
it is not yet a result.

## What cannot be tested this way

Matching controls composition by removing it, so the matching variables and
anything strongly collinear with them cannot be tested here. RSA, neighbour count
and hydrophobic fraction are balanced by construction. ND2 packing, residue
counts within 8 Å and side-chain crowding are all near-inverses of exposure, so
their collapse under matching is partly mechanical rather than evidential.

**This analysis therefore cannot say whether exposure governs occupancy.**
Exposure has been matched away. Answering that needs a within-protein comparison
— occupied and unoccupied sequons in the same protein — which currently has 31
sites.

## A recorded error

On first reading the Step 3 output I argued that because the bacterial and
cytosolic diagnostics reproduced the secretory effect sizes, the secretory
effects must be compositional. **The conclusion happened to survive the proper
test; the reasoning was invalid.** Bacterial and cytosolic sequons cannot be
occupied in any compartment sense, so they are unoccupied regardless of local
structure, and the difference between them and occupied sites carries no
information about whether local structure governs occupancy. Equal effect sizes
can have unequal causes.

The correct argument is the one above and it never needed those sets: the
populations share no proteins, and matching removes the effects.

## What follows

- The within-protein comparison is the analysis that would answer the actual
  question, and it is currently limited to 31 internal-control sites. Expanding
  that set is worth more than any further work on the existing populations.
- Step 4 (do models respond to context differences?) should use the matched-pair
  feature contrasts, not the population differences, for the same reason.
- No feature here is ready to act as a design filter.
