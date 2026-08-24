# Amendment 1 to the context-difference pre-specification

*2026-08-24, written after the population-level Step 3 ran and before the
amended analysis was computed. It records a design flaw in the original
specification, a wrong interpretation I drew from the first results, and the
analysis that replaces it.*

## The flaw in the original specification

Step 3 as specified compares occupied sites against comparison sites at the
level of whole populations. The populations do not overlap at all:

    internal-control proteins that also carry an occupied site:   15 of 24
    internal-control sites sharing a chain with an occupied site: 13 of 31

    secretory proteins that also carry an occupied site:           0 of 1098
    secretory chains sharing a chain with an occupied site:        0

The occupied-versus-secretory comparison therefore shares **no proteins and no
chains**. Occupancy is completely confounded with protein identity — fold class,
size, expression system, depositor, resolution. Any difference found there is a
difference between two sets of proteins, and cannot be attributed to occupancy.

This is not the label-contamination problem. Contamination of the
secretory-unannotated set biases toward the null and cannot manufacture an
effect. Between-protein composition can.

The internal-control comparison does not have this flaw, which is why its
near-zero estimates are informative rather than merely underpowered: it is
largely a within-protein comparison, and for `n_rsa`, `nd2_atoms_8a_same_chain`,
`nd2_residues_8a_same_chain` and `nearest_aromatic_sidechain_nd2` its 95%
interval excludes the secretory point estimate.

## A wrong interpretation, recorded so it is not repeated

On first reading the Step 3 output I argued that because the bacterial and
cytosolic diagnostics reproduced the secretory effect sizes almost exactly, the
secretory effects must therefore be compositional rather than
occupancy-related.

**That inference is invalid.** Bacterial and cytosolic sequons cannot be
occupied at all — wrong compartment, no access to the oligosaccharyltransferase
machinery. They are unoccupied irrespective of their local structure, so the
difference between occupied sites and those sets cannot be evidence about
whether local structure governs occupancy. Two comparisons can produce the same
effect size for entirely different reasons, and matching magnitudes are
suggestive at best.

What the diagnostics legitimately provide is the scale of a purely compositional
difference in these features, which is what they were extracted for: interpreting
the bacterial and cytosolic *model* effects in Step 4. They are not a calibration
for the occupancy contrast.

The commit message on 8bf84fb states that the diagnostics "are what make the
result readable". That is the wrong interpretation above, and this document
supersedes it. History is not rewritten; the error is recorded instead.

## The amended analysis

Occupied sites are compared to their **frozen matched controls, within pair**,
using the matching the occupancy benchmark already rests on:

- `matched_pairs_secretory.csv` — 262 pairs, one control per case, 72 resample
  units;
- `matched_pairs_optimal.csv` — 16 internal-control pairs.

Matching was built from RSA, neighbour counts and hydrophobic fraction, and
never from model output, so composition is controlled by construction rather
than by assumption. Using the same pairs also keeps this branch commensurate
with the benchmark it is meant to explain.

For each feature the quantity is the within-pair difference

    dx_i = x_case_i - x_control_i

with the null that its mean is zero. Uncertainty is a cluster bootstrap over
`resample_unit` (72 units, not 262 pairs), 2,000 replicates, percentile
interval, BH across the family within each comparison.

## Features that cannot be tested here, and why

`n_rsa`, `n_neighbours_8a` and `neighbour_hydrophobic_fraction_8a` are the
**matching variables**. They are balanced between arms by construction, so a
null result on them means the matching worked, not that occupancy is unrelated
to exposure. They are reported as a balance check and excluded from the tested
family.

This is the central limitation of the amended analysis: it can say whether
occupied sites differ from *structurally matched* unoccupied sequons in features
beyond the three that were matched. It cannot say whether exposure itself
governs occupancy, because exposure has been matched away. That question needs
the within-protein comparison, which at present has 31 sites.

## Status of the original Step 3

Retained and reported, reclassified as **descriptive of population differences**
rather than confirmatory about occupancy. Its secretory column is not evidence
of occupancy-associated biology and is not to be quoted as such.
