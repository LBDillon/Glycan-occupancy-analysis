# Findings — the fixed-sequon context-retention test

*2026-08-26. Pre-specified in
[`prespecification_fixed_sequon_context_retention.md`](prespecification_fixed_sequon_context_retention.md),
including the amendment recorded when the motivating model contrast changed.*

## The question

> When ProteinMPNN is forced to preserve a naturally occupied N-X-S/T sequon,
> does it redesign the surrounding residues into a local environment that is
> unusual for natural occupied sites?

Or: does fixing the motif protect the biology around it, or only the three
letters?

## The answer

Only the three letters — but the margin over arbitrary change is small.

50 sites across 38 proteins, 32 designs each, sequon held fixed and verified
present in every design. Scored against 285 natural occupied sites, each site's
own protein held out of its own reference.

| Quantity | Mean | 95% CI | p |
|---|---|---|---|
| D(wild type) | 0.633 | [0.584, 0.674] | |
| D(design) | 0.703 | [0.659, 0.754] | |
| D(random control) | 0.665 | [0.644, 0.686] | |
| **design − wild type** | **+0.071** | **[0.025, 0.122]** | **0.002** |
| random − wild type | +0.032 | [−0.002, 0.069] | 0.069 |
| design − random | +0.039 | [−0.001, 0.080] | 0.056 |

Designs drift away from natural occupied context with the motif fully protected.
**ProteinMPNN is not measurably better at preserving that context than changing
the same number of residues at random** — numerically it is slightly worse, and
that difference sits on the significance boundary.

## Why the random control decides the reading

The wild types *are* natural occupied sites, so they sit inside the reference
distribution by construction and any perturbation moves outward. Without a
control, "+0.071, p = 0.002" would read as a finding about ProteinMPNN when much
of it is the arithmetic of disturbing a point that starts near the centre. The
mutation-count-matched control puts roughly half the drift down to that.

29 of the 50 sites move outward and 21 move inward. The effect is a tendency
across sites, not something that happens to each one.

## No single feature carries it

Nothing survives correction. Three reach q < 0.10, and they are at least
chemically coherent with one another:

| Feature | Shift (reference SDs) | q |
|---|---|---|
| aromatic content in the ND2 shell | −0.251 | 0.075 |
| glycine in the flanking window | +0.233 | 0.100 |
| glycine in the ND2 shell | +0.201 | 0.100 |

More glycine and fewer aromatics near the attachment point. That is a lead worth
checking in a larger run, not a result.

## What this licenses

**It says**: preserving the sequon is not sufficient to preserve its
environment. For SugarFix's preserve mode that is directly useful — retention as
a metric does not capture what redesign does *around* a protected site, so
contextual constraints are a real candidate rather than a hypothetical one.

**It does not say** the designs would be unglycosylated. The panel measures
resemblance to natural occupied context, not occupancy, and nothing here has
been tested experimentally.

**It is not about ProteinMPNN specifically.** Per the amendment, the
architectural contrast that originally motivated this — ESM-IF discriminating
where ProteinMPNN did not — largely dissolved when ProteinMPNN's sequon indexing
was corrected on 25/08. This is a statement about protected-sequon redesign.

## Scope and limits

- 50 sites, 38 proteins. Small. `design − random` at [−0.001, +0.080] would be
  settled either way by a larger run, and that is the number that matters.
- Fixed-backbone design cannot change accessibility, secondary structure or
  distance to the termini, so those are excluded from the outcome by
  construction and condition the reference instead. This says nothing about
  whether the *structural* setting is preserved, because it cannot move.
- The reference is 285 of 318 occupied sites; the rest are chains whose two
  parses cannot be reconciled and are refused rather than guessed.
- The distance is a median absolute standardised departure over 15 features. It
  is a summary, deliberately not a density model, and not calibrated against any
  experimental occupancy measurement.

## Figures

Captions in [`figure_captions.md`](figure_captions.md); the figures themselves
carry only titles, axes, legends and significance markers.

- `fig4_context_retention.png` — the paired result and the random control
- `fig5_context_features.png` — per-feature shifts with intervals
- `fig6_feature_distributions.png` — natural, wild type and design as empirical
  cumulative distributions, per feature
- `fig2_occupied_context.png` — the natural reference the distances are measured
  against

Terms are defined in [`glossary.md`](glossary.md).
