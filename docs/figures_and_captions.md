# Figures and captions

One entry per figure: the file, and what it shows. Figures carry titles, axes,
legends and significance markers only — everything explanatory lives here, so a
figure reads quickly and is understood properly.

Terms are defined in [`glossary.md`](glossary.md). Current results are in
[`OVERVIEW.md`](OVERVIEW.md). Exact numbers for the model-comparison figures are
in the `*_values.json` written beside them, which is regenerated with the figure
and does not go stale.

---

## Part 1 — the occupancy benchmark

Regenerate figures 1–7 with `python pipeline/20_figures_summary.py`, figure 8
with `22_figures_control_provenance.py`, figure 9 with `21_figures_all_classes.py`.

### `fig1_attrition.png`

![attrition](../results/figures/fig1_attrition.png)

How many sites survive each stage of the pipeline, for the two groups being
compared. The occupied group falls steeply once a solved structure is required;
the internal controls do not fall at all, because they exist only by virtue of
having come from structures. Both groups are small by the final matching step.

The occupied group is never the bottleneck. Every constraint in this study comes
from the internal controls, and by the end only half of them are usable — so this
is the figure that says where more effort would actually buy something.

### `fig2_three_comparisons.png`

![three comparisons](../results/figures/fig2_three_comparisons.png)

Occupied sites against each of the four control sets. The dot is the estimate,
the bar the 95% confidence interval, the grey band the difference agreed in
advance to be too small to care about.

The primary comparison has the largest point estimate and by far the widest
interval, resting on very few pairs. The eukaryotic secretory comparison — the
one matching occupied sites on both taxonomy and compartment — supplies the most
pairs and lands essentially on zero. The better-powered and better-matched the
comparison, the closer to zero it sits.

### `fig3_matching_sensitivity.png`

![matching sensitivity](../results/figures/fig3_matching_sensitivity.png)

The old matching procedure run 200 times from different random starting points;
each bar counts how many runs gave that answer. The whole distribution sits to
the right of zero, but only a minority of runs produced a confidence interval
excluding zero.

Two things are true at once and matter separately: the *direction* is robust,
and the *statistical significance* was not — it depended on an arbitrary seed
inside the matching code. The matching is deterministic now, and the vertical
line is the single answer it gives.

### `fig4_retention_bridge.png`

![retention bridge](../results/figures/fig4_retention_bridge.png)

Sites sorted into five equal groups by the score the model assigns them, against
how often the model actually keeps the sequon when asked to redesign the
protein. Retention rises smoothly from the lowest fifth to the highest.

The abstract number measured and the concrete thing the model does are two views
of one quantity. This is what connects the score-based analysis to the earlier
scoping work, which looked at generated sequences rather than internal
probabilities.

### `fig5_retention_distribution.png`

![retention distribution](../results/figures/fig5_retention_distribution.png)

For each site, the fraction of 32 unconstrained designs that kept the sequon
intact. The vertical axis is a log scale, so the first bar is far taller than it
looks. Most sites lose the sequon in *every* design; a small tail keeps it every
time.

Hand a glycoprotein to ProteinMPNN and ask it to redesign the surface, and the
glycosylation sites should be expected to disappear. Nothing in the model's
training told it they mattered.

### `fig6_scorer_defect.png`

![scorer defect](../results/figures/fig6_scorer_defect.png)

The distribution of every score computed, with invalid ones marked. Real scores
occupy a plausible range; the flagged cluster sits at a value the quantity
cannot take. ProteinMPNN silently returns zeros for any residue with a missing
backbone atom, and exponentiating zeros yields "probability 1 that this is
asparagine, probability 2 that it is serine or threonine".

Only a few percent of sites were affected, but every result in the study is
expressed in units of the spread of these scores, and those few absurd values
nearly doubled that spread. A small data problem became a large error in the
measuring stick. Two independent checks now reject anything that is not a
genuine probability distribution.

### `fig7_retention_by_class.png`

![retention by class](../results/figures/fig7_retention_by_class.png)

Sequon retention split by whether the site is a real glycosylation site or a
control. Error bars bootstrap over whole proteins rather than sites, because
several sequons on one protein share one set of designs and are not independent.

The control sets differ from *each other*, so a class average cannot answer the
question — whatever separates eukaryotic secretory proteins from bacterial and
cytosolic ones rides along in the comparison. Because every control was matched
site by site to an occupied site, the paired contrast is available and is the
right test.

### `fig8_control_provenance.png`

![control provenance](../results/figures/fig8_control_provenance.png)

Left, every filter applied to build the eukaryotic secretory control set, from
candidate proteins down to usable matched pairs; the highlighted bar is the step
that constitutes the entire negative label for that set — the absence of a
UniProt glycoprotein keyword. Right, all four control sets placed against the two
things that matter: how well the population matches the occupied sites, and how
defensible the negative label is. Marker area is the number of usable pairs.

The internal controls sit in the corner you would want, strong on both axes —
they are not a compromise, they are simply tiny. Every other set buys two orders
of magnitude more pairs by giving up one axis. Nothing large occupies the good
corner, and that is the whole problem.

### `fig9_all_classes.png`

![all classes](../results/figures/fig9_all_classes.png)

Left, what the model *writes* — design retention for all five classes. Right,
what the model *believes* — the conditional-score difference against each control
set. Two different quantities on different scales, and they agree: neither shows
the model distinguishing occupied sequons by a margin this design can resolve.

---

## Part 2 — schematics

Five drafts. Only figure 11 carries data; the rest explain what the pipeline does
and what is at stake. **All are drafts — the shapes and wording are meant to be
argued with.** Regenerate with `python pipeline/24_figures_schematics.py`.

| File | What it shows |
|---|---|
| `fig10_mechanism.png` | Why any of this matters. What is really there (an exposed asparagine carrying a glycan), what ProteinMPNN sees (the same backbone, no glycan — identical to it), and what it designs (a hydrophobic residue that packs better and deletes the site). Usually a sound design move; here it aglycosylates the protein. Schematic, not a real structure. |
| `fig11_dataset_summary.png` | What is in each dataset — size before and after matching, taxonomic spread, and sequon exposure. The taxonomy panel is a **caveat**: the eukaryotic secretory set matches on kingdom and compartment but not species, and species composition is not matched anywhere in the study. |
| `fig12_scoring_process.png` | Scoring a sequon end to end: structure → scoreability check → model reads backbone plus native sequence → probabilities at the three positions → log odds → one number. |
| `fig13_design_process.png` | The design branch: structure → 32 unconstrained sequences → read the three original positions → classify → retention fraction. Nothing is fixed or biased, so the sequon is free to vanish. |
| `fig14_awareness_spectrum.png` | A proposed scale from glycan-blind to glycan-aware, and what the evidence says about where a model sits on it. **The axis itself is the open question** — worth settling before this is used in a write-up. |

---

## Part 3 — the model comparison

Regenerate the first four with `python pipeline/25_figures_model_comparison.py`,
the heatmaps with `26_figures_sequon_heatmap.py`. Exact values are in
`fig_model_comparison_values.json` and `fig_sequon_heatmap_values.json`.

### `fig_scoring.png`

![scoring](../results/figures/fig_scoring.png)

**A:** occupied minus matched control for each model, in that model's own
reference SD; a marker means the interval excludes zero. **B:** the per-site
contrasts those averages are computed from, as cumulative distributions.

B is the corrective to A and should be read with it. The distributions overlap
heavily and a substantial minority of sites have negative contrasts. These are
shifts in broad overlapping distributions, not separation — nothing here would
classify a site.

Raw magnitudes are not comparable between models — eight averaged decoding
orders, one causal pass, one masked position, one shot — so each model's
contrasts are divided by its own reference SD.

### `fig_masking.png`

![masking](../results/figures/fig_masking.png)

What survives hiding the motif. Filled dot is the motif visible, open dot hidden,
arrow the change. Both arms of every comparison contain a sequon, so motif
recognition alone cannot separate them; whatever remains once the motif is hidden
is coming from the surroundings.

Three things to hold while reading it. The markers refer to **the change**, not
the effect. The dots are SD-standardised but the significance is computed in log
odds, so the two scales in the figure are not the same scale. And a model with no
row here has no masked condition at all — CARBonAra is one-shot and has no
motif-visible arm — which is a property of the model, not missing data.

### `fig_retention.png`

![retention](../results/figures/fig_retention.png)

Sequon retention, occupied against matched control, for the models that can
redesign a chain. A blank row is a property of the model — no `SequenceDesigner`
— not a gap in the data.

Significance is marked from the cluster-aware interval, which respects protein
and ortholog clustering. A Wilcoxon test on the same pairs can disagree with it,
and where they disagree the interval is the reading used, being the conservative
one. Many pairs are **tied** — both sites lose the sequon — so the effect rests
on considerably fewer informative pairs than the pair count suggests.

### `fig_score_distributions.png`

![raw score distributions](../results/figures/fig_score_distributions.png)

The raw conditional sequon score — log odds, not standardised — for
experimentally occupied sites against their matched partners with no annotated
glycan. Dashed lines are medians, one panel per model, laid out as a grid.

It makes the size of the effect legible in a way an SD cannot: if the two groups
were identical, half the occupied sites would fall below the other group's
median. Rather fewer do, which is a real shift and a modest one.

**The panels must not be read against each other.** Raw log odds are comparable
within a model and not between, so the x-axes differ deliberately and these are
separate comparisons that happen to share a figure.

**"No annotated glycan" is not "unoccupied".** These are
`control_secretory_eukaryotic_unannotated`: absence of annotation, not evidence
of absence. A separate and much smaller group carries real evidence of absence
and is not shown. Reading the grey distribution as unglycosylated sites is the
single easiest mistake to make with this figure.

### `fig_sequon_heatmap.png`

![sequon heatmap](../results/figures/fig_sequon_heatmap.png)

The full predicted distribution over amino acids at each of the three sequon
positions — the `probs_n`, `probs_plus1` and `probs_plus2` the score files
already hold. Rows are models, the two columns are the two arms of the matched
comparison, each panel three positions down by twenty amino acids across. An
interactive version with hover values is alongside as `.html`.

Amino acids are grouped by side-chain chemistry rather than alphabetically. **N**,
**S** and **T** are bold because they are the residues the score actually reads;
that the three sit adjacent is a property of the polar-uncharged group, not an
arrangement chosen to flatter the result.

Rings mark an amino acid whose **paired** difference between the arms clears zero
on a bootstrap over resampling units — the same clustering the main analysis
uses, not a site-level test.

Three cautions. There are sixty cells per model and no multiplicity correction,
so a few rings per model are expected by chance and individual rings on
low-probability cells are not findings; this decomposes an effect established
elsewhere, it is not a screen. **Nothing is confidently predicted anywhere** —
most cells sit near chance for twenty residues, so the models express mild
preferences rather than calls. And the rows are not all the same site set, so
they are sound within a model but should not be used to rank models. The two ESM3
arms *are* the same sites as one another and may be compared row to row.

#### The aspartate column

Asparagine and aspartate are isosteric, and the structure-conditioned models
shift **both** between the arms — for some of them the aspartate shift is the
larger of the two. What a structure model registers at occupied sites is
therefore not specifically a preference for asparagine but a preference for an
amide-shaped residue, of which asparagine is one and aspartate the other. The
sequence-only models behave differently, carrying much more asparagine mass and
almost no aspartate shift: a sequence model reads the motif in context and can
name the residue, a structure model sees a shape.

This does not undermine the score — the paired design subtracts whatever is
common to both arms, and pooling Asn with Asp barely changes which sites rank
highly. It does mean the asparagine term should be described as *amide-shape
preference* rather than *asparagine recognition* when a structure model is the
subject.

### `fig_sequon_heatmap_diff.png`

![sequon heatmap difference](../results/figures/fig_sequon_heatmap_diff.png)

The paired difference of the previous figure — occupied minus its matched partner
— plotted directly rather than left to the reader to subtract by eye. Both arms
contain a sequon by construction, so side by side they look nearly identical and
the separation is invisible at that scale. Diverging scale centred on zero,
symmetric so equal shifts in either direction read the same size; rings are the
same bootstrap test.

It makes the per-residue structure of each model's effect legible on one page,
and it is where a model with an inconclusive scalar score becomes interpretable —
such a model appears here as a nearly flat row.

**Read it with `fig_sequon_heatmap.png` open.** A difference of +0.06 means
something quite different where the underlying probability is 0.25 than where it
is 0.02, and this figure cannot show that. The absolute levels are in the other.

---

## Part 4 — the glyco-site context analysis

### `glyco_context/results/figures/fig2_occupied_context.png`

Local environment of the experimentally supported N-linked glycosylation sites
with usable structures (the `triplet_core` view). **Left:** relative solvent
accessibility at the asparagine and at the +1 and +2 positions. **Centre:** DSSP
secondary structure, coarsened to loop, sheet and helix. **Right:** backbone
region at the asparagine, from phi and psi.

The left-handed alpha region holds a substantial minority of sites. This is
expected rather than anomalous: asparagine and glycine occupy left-handed
conformations far more readily than other residues, so the category is real and
not a rounding artefact.

No comparison group appears here. This is the reference picture that the design
distances in the next two figures are measured against.

### `glyco_context/results/figures/fig4_context_retention.png`

Occupied sites across a few dozen proteins, each redesigned 32 times by
ProteinMPNN with the sequon held fixed. **Top:** distance from natural occupied
context for the wild type, a mutation-count-matched random control, and the
designs. **Bottom:** the paired change from wild type. Bars are 95% intervals
from resampling proteins.

Every design preserves N-X-S/T, verified by reading the triplet back out of the
designed sequences.

The random control is what makes the top row readable. A wild type is itself a
natural occupied site, so it sits inside the reference by construction and any
perturbation moves it outward; the control changes the same number of residues at
random and so measures how much of the drift is simply that. Designs move further
than the control, but not significantly at this sample size.

### `glyco_context/results/figures/fig5_context_features.png`

Per-feature shift from wild type to design, in standard deviations of the natural
reference, for the features whose intervals separate from zero. Bars are 95%
intervals from resampling proteins; q values are Benjamini-Hochberg across all
fifteen features. Faint bars are weaker evidence.

**Nothing survives correction at q < 0.05.** The strongest few are chemically
coherent with one another, which is a reason to check them in a larger run and
not a reason to believe them yet.

### `glyco_context/results/figures/fig6_feature_distributions.png`

Empirical cumulative distributions for each of the fifteen panel features.
**Grey:** all natural occupied sites, the reference. **Black:** the wild-type
sites that were redesigned — a subset of the grey, drawn separately as a check
that the tested sites are representative. **Blue:** the designs. Panel titles
give the Benjamini-Hochberg q for the paired wild-type-to-design shift.

A shift between distributions appears as horizontal separation; a mass of
identical values appears as a vertical jump.

Cumulative distributions rather than histograms because these features are
heavily zero-inflated, so a histogram's zero bar would dominate and a kernel
density would spread mass across values that never occur. An eCDF uses every
observation exactly and requires no binning choice. The designs carry roughly as
much independent information as there are sites, so the smoothness of the blue
curve should not be read as precision.

---

## Archived — the comparative analysis

These belong to a question that was asked, answered negatively and set aside. See
[`comparative_analysis/README.md`](../glyco_context/archive/comparative_analysis/README.md).

### `glyco_context/archive/comparative_analysis/results/fig1_effects_collapse.png`

**Left:** standardised difference between occupied sites and secretory-unannotated
comparison sites, before and after controlling composition by matching. Ringed
points are the three matching variables, which are balanced by construction and
cannot be tested here — a null on them means the matching worked. **Right:** the
proportion of each comparison set's proteins that also carry an occupied site.

The occupied and secretory-unannotated sets share no proteins and no chains, so
occupancy is confounded with protein identity and every difference is also a
difference between two sets of proteins. The internal controls largely are not,
which is why their near-zero estimates carry weight despite being few.

### `glyco_context/archive/comparative_analysis/results/fig3_framings.png`

Five features under each of four analyses: population-level and matched-pair,
against internal controls and against secretory-unannotated sites. Asterisks mark
survival of Benjamini-Hochberg within that comparison.

The two matched-pair framings control composition; the two population framings do
not, and their larger values are the confounded ones. Only beta-sheet at the +2
position holds its direction throughout.
