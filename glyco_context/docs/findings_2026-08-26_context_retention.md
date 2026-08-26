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

Only the three letters — and ProteinMPNN is measurably *worse* than arbitrary
change of the same size.

**285 sites across 207 proteins**, 32 designs each, sequon held fixed and
verified present in every design. Every site that passes the mapping gates is
included. Scored against the same 285 natural occupied sites, each site's own
protein held out of its own reference.

| Quantity | Mean | 95% CI | p |
|---|---|---|---|
| D(wild type) | 0.656 | [0.639, 0.674] | |
| D(design) | 0.706 | [0.691, 0.721] | |
| D(random control) | 0.675 | [0.667, 0.683] | |
| **design − wild type** | **+0.050** | **[0.033, 0.067]** | **0.0005** |
| random − wild type | +0.018 | [0.006, 0.031] | 0.002 |
| **design − random** | **+0.032** | **[0.019, 0.045]** | **0.0005** |

Designs drift away from natural occupied context with the motif fully protected,
and they drift **further than changing the same number of residues at random**.
That last row was the open question at the first pass, where it read
+0.039 [−0.001, 0.080] and sat on the boundary. At full scale it resolves.

*A first pass over 50 sites gave +0.071 for design and +0.039 for design minus
random. Both estimates shrank at full scale — the direction held, the magnitude
was somewhat inflated by the smaller sample, which is the usual direction for
that error and worth recording.*

## Why the random control decides the reading

The wild types *are* natural occupied sites, so they sit inside the reference
distribution by construction and any perturbation moves outward. Without a
control, "+0.071, p = 0.002" would read as a finding about ProteinMPNN when much
of it is the arithmetic of disturbing a point that starts near the centre. The
mutation-count-matched control puts roughly half the drift down to that.

The effect is a tendency across sites rather than something that happens to each
one: a substantial minority move inward.

## What moves: proline and glycine in, aromatics out

Five of fifteen features survive correction, and they describe one coherent
change:

| Feature | Shift (reference SDs) | q |
|---|---|---|
| proline in the ND2 shell | **+0.213** | 0.0025 |
| aromatic content in the ND2 shell | **−0.209** | 0.0025 |
| proline in the flanking window | **+0.186** | 0.0025 |
| aromatic content in the flanking window | **−0.118** | 0.015 |
| glycine in the flanking window | **+0.106** | 0.033 |
| glycine in the ND2 shell | +0.124 | 0.093 |

**Designs put proline and glycine near the site and take aromatics away**, in
both the sequence window and the three-dimensional shell.

At 50 sites the leading feature appeared to be glycine, with proline at q = 0.21.
At full scale proline is the strongest of the set. Which feature leads was not
stable at the smaller sample; the direction of the whole pattern was.

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

- 285 sites, 207 proteins — every site that passes the mapping gates, so this
  cannot be enlarged without relaxing them or extending to other populations.
- **The random control preserves composition; ProteinMPNN does not.** Its
  replacements are drawn from the wild-type chain's own amino-acid frequencies,
  so it changes *which* residues sit where without changing the overall mix.
  ProteinMPNN has its own composition preferences, so part of `design − random`
  may be a global bias rather than anything specific to sequon surroundings. The
  control that would separate these is ProteinMPNN's shift at non-sequon
  positions in the same proteins, which has not been measured.
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

Captions in [`figure_captions.md`](../../docs/figures_and_captions.md); the figures themselves
carry only titles, axes, legends and significance markers.

- `fig4_context_retention.png` — the paired result and the random control
- `fig5_context_features.png` — per-feature shifts with intervals
- `fig6_feature_distributions.png` — natural, wild type and design as empirical
  cumulative distributions, per feature
- `fig2_occupied_context.png` — the natural reference the distances are measured
  against

Terms are defined in [`glossary.md`](../../docs/glossary.md).
