# Figure captions

Figures carry titles, axes, legends and significance markers only. Everything
that explains, qualifies or interprets lives here, so a figure can be read
quickly and understood properly.

Terms are defined in [`glossary.md`](glossary.md).

---

## fig2 — the natural reference

**`results/figures/fig2_occupied_context.png`**

> Local environment of 318 experimentally supported N-linked glycosylation sites
> with usable structures (the `triplet_core` view). **Left:** relative solvent
> accessibility at the asparagine and at the +1 and +2 positions of the sequon.
> **Centre:** DSSP secondary structure, coarsened to loop, sheet and helix.
> **Right:** backbone region at the asparagine, from phi and psi.
>
> The left-handed alpha region holds 18% of sites. This is expected rather than
> anomalous: asparagine and glycine occupy left-handed conformations far more
> readily than other residues, so the category is real and not a rounding
> artefact.
>
> No comparison group appears here. This is the reference picture that the
> design distances in fig4 and fig6 are measured against.

## fig4 — the context-retention result

**`results/figures/fig4_context_retention.png`**

> Fifty occupied sites across 38 proteins, each redesigned 32 times by
> ProteinMPNN with the sequon held fixed. **Top:** distance from natural
> occupied context for the wild type, a mutation-count-matched random control,
> and the designs. **Bottom:** the paired change from wild type. Bars are 95%
> intervals from resampling proteins; markers give the two-sided bootstrap p.
>
> Every design preserves N-X-S/T, verified by reading the triplet back out of
> the designed sequences.
>
> The random control is what makes the top row readable. A wild type is itself a
> natural occupied site, so it sits inside the reference by construction and any
> perturbation moves it outward. The control changes the same number of residues
> at random and so measures how much of the drift is simply that — about half.
>
> Designs move further than the control (+0.039), but that difference is not
> significant at this sample size.

## fig5 — which features move

**`results/figures/fig5_context_features.png`**

> Per-feature shift from wild type to design, in standard deviations of the
> natural reference, for the five features of fifteen whose intervals separate
> from zero. Bars are 95% intervals from resampling proteins; q values are
> Benjamini-Hochberg across all fifteen. Faint bars are weaker evidence.
>
> **Nothing survives correction at q < 0.05.** The three strongest — less
> aromatic content near ND2, more glycine in both the flanking window and the
> shell — are chemically coherent with one another, which is a reason to check
> them in a larger run and not a reason to believe them yet.

## fig6 — feature distributions

**`results/figures/fig6_feature_distributions.png`**

> Empirical cumulative distributions for each of the fifteen panel features.
> **Grey:** all 285 natural occupied sites, the reference. **Black:** the 50
> wild-type sites that were redesigned — a subset of the grey, drawn separately
> as a check that the tested sites are representative. **Blue:** 1,600 designs,
> 32 from each of the 50 sites. Panel titles give the Benjamini-Hochberg q for
> the paired wild-type-to-design shift.
>
> A shift between distributions appears as horizontal separation; a mass of
> identical values appears as a vertical jump.
>
> Cumulative distributions rather than histograms because these features are
> heavily zero-inflated — 91% of natural sites have no cysteine in the ND2 shell
> — so a histogram's zero bar dominates and a kernel density would spread mass
> across values that never occur. An eCDF uses every observation exactly and
> requires no binning choice.
>
> The 1,600 designs carry roughly 50 sites' worth of independent information, so
> the smoothness of the blue curve should not be read as precision.

---

## Archived — the comparative analysis

These belong to a question that was asked, answered negatively and set aside.
See [`../archive/comparative_analysis/README.md`](../archive/comparative_analysis/README.md).

### fig1 — effects collapse under matching

**`archive/comparative_analysis/results/fig1_effects_collapse.png`**

> **Left:** standardised difference between occupied sites and secretory-
> unannotated comparison sites, before and after controlling composition by
> matching. Ringed points are the three matching variables, which are balanced
> by construction and cannot be tested here — a null on them means the matching
> worked. **Right:** the proportion of each comparison set's proteins that also
> carry an occupied site.
>
> The occupied and secretory-unannotated sets share no proteins and no chains,
> so occupancy is confounded with protein identity and every difference is also
> a difference between two sets of proteins. The internal controls largely are
> not, which is why their near-zero estimates carry weight despite n=31.

### fig3 — the same features under four framings

**`archive/comparative_analysis/results/fig3_framings.png`**

> Five features under each of four analyses: population-level and matched-pair,
> against internal controls and against secretory-unannotated sites. Asterisks
> mark survival of Benjamini-Hochberg within that comparison.
>
> The two matched-pair framings control composition; the two population framings
> do not, and their larger values are the confounded ones. Only beta-sheet at
> the +2 position holds its direction throughout.
