# Glossary — what the terms mean

Every term this project uses, for someone arriving fresh. Merges the former
`concepts.md` and `glyco_context/docs/glossary.md`, which split the same job
between two files.

**No current results appear here.** Explanations of *how* a measurement works are
stable; the numbers they produced are not, and were the reason the old
`concepts.md` needed a staleness banner. For any figure, see
[`OVERVIEW.md`](OVERVIEW.md).

---

## 1. The biology

**N-linked glycosylation** — attachment of a sugar chain to the side-chain
nitrogen (ND2) of an asparagine. Not decoration: it governs folding, secretion,
stability, half-life and immune recognition.

**Sequon** — the motif N-X-S/T that N-linked glycosylation requires: asparagine,
then almost any residue except proline, then serine or threonine. **Necessary
but not sufficient** — most sequons in a proteome are never glycosylated.

**Occupied** — a sequon with experimental evidence that a glycan is attached.
The opposite is not "unoccupied" but *unknown*, which is the central difficulty
of this project: absence of evidence accumulates in understudied proteins.

**ND2** — asparagine's side-chain nitrogen, the atom the glycan bonds to.
Features named "ND2 shell" are centred there rather than on the backbone,
because that is where the sugar sits.

**NXS / NXT** — the two subtypes, by whether the third position is serine or
threonine. NXT is glycosylated more often in nature, and the two residues are
chemically different, so a model may score them differently for reasons having
nothing to do with occupancy. They are never pooled without stratifying, and
matching pairs them exactly.

**Oligosaccharyltransferase (OST)** — the enzyme that attaches the glycan.
Present in the secretory pathway, absent from the cytosol and from bacteria,
which is why sequons in those compartments cannot be occupied at all.

---

## 2. Site sets and populations

**`occupied`** — sites with positive occupancy evidence.

**`internal_control`** — a sequon with no modelled glycan in a structure that
*does* model glycans elsewhere on the same protein, from an organism that can
glycosylate. Normally a bare asparagine proves nothing: the glycan may have been
trimmed before crystallisation, or the sugar too mobile to see. But where the
crystallographer demonstrably *could* see glycans, a bare asparagine is a
decision rather than a silence. "Internal" means from inside the same set of
structures as the occupied sites — same organisms, same compartment, same kind
of experiment. The cleanest comparison available, and there are only 32 of them.
They are **not** proven negatives.

**`secretory_unannotated`** — eukaryotic secretory sequons with no occupancy
annotation. Numerous, weakly labelled: roughly half of these proteins carry a
glycoprotein keyword, so the set is contaminated by construction. Contamination
biases toward the null, so it cannot manufacture an effect.

**`bacterial` / `cytosolic`** — diagnostics, not tests. Informative about how the
measurement behaves; not a valid answer to the question. Both are confounded on
purpose and in opposite directions — cytosolic proteins are eukaryotic but never
meet the machinery (confound: compartment); bacterial secreted proteins are the
right kind of protein from organisms with no N-glycosylation (confound:
kingdom). The original hope was that a signal surviving both confounds would be
more credible. The corrected numbers did not support that reasoning.

**`triplet_core`** — the primary context view: every feature in the row describes
the sequon it names. Requires the observed triplet to match the expected one,
all three residues located, and the mapping continuous.

**`asn_core`** — wider, requiring only that the asparagine itself was measured
correctly. Valid for features centred on the asparagine; **not** for +1 or +2
exposure, structure or geometry.

**`construct_review`** — everything excluded from `triplet_core`, each row
carrying its reason. Inspected, never tested.

---

## 3. How the benchmark measures things

### The conditional sequon score

Show the model a backbone and the native sequence, and ask what probability it
holds at each sequon residue. Take P(asparagine) at the first position and
P(serine) + P(threonine) at the third — either satisfies the motif — convert each
to log odds, and average the two.

**Log odds** converts a probability to a symmetric scale: p/(1−p) is the odds,
and the log makes it symmetric about zero. p = 0.5 → 0; p = 0.12 → −2.0;
p = 0.88 → +2.0. Without it, probabilities bunch near zero and averages distort.

**Why the middle residue is excluded.** Any residue except proline satisfies the
motif there, so a preference at that position is not a preference for the motif.
Proline is recorded separately, because proline there abolishes the sequon.

**Sequon level, not protein level** — three residues. A whole-protein score would
swamp three residues among several hundred.

### Retention

What the model actually *writes* when redesigning, rather than what it believes
while reading. Generation involves hundreds of interacting decisions, so a sequon
could survive or vanish for reasons unrelated to the model's opinion about one
position. The two outcomes are measured separately because they can disagree.

### Matching

Pairing an occupied site with a control site in a similar local structural
environment, so a score difference cannot be blamed on the environment. Three
metrics, all computed from structure and **never from model output**:

- **RSA** — how exposed the residue is, 0 (buried) to 1 (exposed)
- **Neighbour count within 8 Å** — a crude density measure
- **Hydrophobic fraction within 8 Å** — the chemical character of the pocket

A pair is allowed only within a **caliper** — a maximum combined distance in
those three dimensions. Beyond it, no match is made rather than a bad one. NXS
pairs only with NXS and NXT with NXT: before that rule, about 45% of pairs put an
occupied NXS against an unoccupied NXT, so any difference could have been a
subtype difference wearing occupancy's clothes.

**Greedy matching** walks the occupied sites in some order, giving each its
nearest unused control — "greedy" because it takes the best option at each step
regardless of later consequences. **The seed** sets that order. With few controls
this matters enormously: if control C is the only admissible partner for site B
but site A is processed first and takes C, B goes unmatched. **Deterministic**
matching solves the whole assignment at once, maximising pairs and then
minimising total distance. No seed, one answer. A 200-seed sweep showed the
direction was robust across seeds while significance was being decided by the
random number — which is the argument for the deterministic version.

### Clusters and uncertainty

**Ortholog cluster** — the same protein across species. Sites in one cluster are
not independent: human and mouse albumin with an occupied sequon at equivalent
positions is close to one fact observed twice.

**Resample unit** — the connected group formed by linking occupied clusters to
the control proteins they share, since one control protein can serve several
occupied sites. Intervals resample these, not individual sites.

**Cluster bootstrap** — resampling whole groups rather than individual
observations, so members that move together stay together. Treating correlated
observations as independent produces intervals that are too narrow.

**Equivalence testing, and the ±0.2 SD margin.** An ordinary test can only reject
a null or fail to; failing to reject is compatible with both "no difference" and
"too small a sample to see one". Equivalence testing states in advance how small
a difference would count as none, then asks whether the interval fits inside that
band. The ±0.2 SD margin is **an exploratory statistical threshold, not a
biologically validated one** — 0.2 SD is a conventional "small effect" and
nothing in glycobiology says what a meaningful shift in model log-odds would be.
What matters is that it was fixed before any comparison was computed.

**Permutation rather than Wilcoxon.** Wilcoxon and the sign test assume pairs
are independent. They are not: occupied sites in one ortholog cluster are near
copies, and one control protein can serve several occupied cases. The test used
here flips the sign of every contrast within a whole resample unit at once, so
the null respects the dependency. Wilcoxon p-values are reported alongside and
are systematically smaller — that gap is the dependency, not extra evidence.

**Correcting across all eight tests.** Four control sets times two outcomes were
run, and correcting within one outcome only would understate how many chances
the result had to appear. Benjamini-Hochberg and Holm are both reported: BH
controls the expected false-discovery proportion, Holm the chance of any false
positive at all.

### Scoreability

Whether the model will process all three sequon residues, which comes down to
whether the backbone atoms N, CA, C and O are present. The model builds a
residue's representation from that geometry; if one atom is missing — common at
flexible loops and chain ends — the position is unusable. Scoreability depends
only on coordinates, so it is determined **before** matching: deciding it
afterwards drops sites from pairs that matching had just balanced.

---

## 4. The context features

Each says what it measures, why it is worth measuring, and how it is obtained.

### Sequence context

**Flanking composition** (`flank_*_fraction`) — proportion of residues in each
chemical class among the ±5 residues either side of the sequon, *excluding the
sequon*. **Why:** the immediate neighbourhood affects enzyme access and local
folding. **How:** counted from sequence; the window is clipped at chain ends, so
the denominator varies between sites.

**`uniprot_residues_after_asn` / `_after_sequon`** — residues following the site
in the full-length protein. **Why:** position along the chain relates to
co-translational timing — a site near the C-terminus meets the enzyme later.
**How:** from the UniProt sequence, not the structure, so a truncated construct
cannot shorten it.

### Local structure

**RSA** — as above. **How:** Shrake-Rupley surface area over a per-residue
maximum, with glycans excluded, so an occupied site's own sugar cannot make it
look buried.

**Secondary structure** — helix, sheet or loop at each sequon position, from
DSSP's eight classes coarsened to three. **Why:** glycosylation favours
flexible, accessible regions.

**Phi / psi and backbone region** — the two backbone dihedrals and a coarse
region label (alpha-right, beta, alpha-left). **How:** from backbone atoms, never
across a gap in the model, where neighbouring atoms belong to a different stretch
of chain. Angles are circular, so they enter comparisons as a region rather than
a mean.

**Loop-run length, and censoring** — how many consecutive loop residues contain
the asparagine, and whether that run reaches an unresolved boundary. A run
reaching the end of the model is a **lower bound** and is flagged; unflagged it
would bias loop lengths downward exactly where density is poor.

### Environment around the attachment point

**ND2 shell** (`shell_*_fraction`, `nd2_atoms_8a_*`, `nd2_residues_8a_*`) —
residues with any heavy atom within 8 Å of ND2. **Why:** this is the space the
glycan must occupy. **How:** measured from ND2 rather than the alpha carbon, and
counted separately for the same chain and other chains, so an oligomer interface
or crystal contact cannot read as local sequence context.

**`shell_net_charge`** — formal charge in that shell: D/E = −1, K/R = +1,
histidine neutral.

**`sidechain_neighbour_residues_5a`** — count of *residues* with a heavy atom
within 5 Å of the asparagine side chain. Named for what it counts; an earlier
name implied atoms and the difference mattered.

**`nearest_aromatic_sidechain_nd2`, `nearest_disulfide_sg_nd2`** — distance from
ND2 to the nearest aromatic side chain, and to the nearest sulfur in a disulfide.
**Note:** disulfide coverage differs sharply by compartment, because disulfides
form in the oxidising secretory pathway — so comparing arms on it compares
missingness rather than biology.

---

## 5. Names this project uses

**Contrast** — one occupied site minus its matched control, on either outcome.
The unit of the benchmark.

**Variant** — a tag naming which run produced a set of numbers (`esm_if`,
`proteinmpnn_index_corrected`). Reading the wrong score file is silent, so stages
refuse to guess which is meant.

**`model_index`** — a residue's ordinal position in a chain as this project's
parser reads it. **Not** interchangeable with a model's own indexing; see
[`methods_sequon_indexing.md`](methods_sequon_indexing.md).

**Fixed-sequon context-retention test** — the experiment asking whether
protecting the motif during redesign also protects its surroundings. Short form:
*does fixing the motif protect the biology around it, or only the three letters?*

**D, and ΔD** — D is a site's distance from natural occupied context: the median
absolute standardised departure across the fifteen-feature panel, with the site's
own protein excluded from the reference. ΔD is the paired change from wild type
to design. Positive means further from natural.

**Panels-only** — a mode computing features without designing anything, used to
build the natural reference, which needs every site but none of the designs.
