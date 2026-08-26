# Glossary

For someone arriving fresh. Terms are grouped by what they describe: the
biology, the site sets, the measured features, and the names this project has
invented for things.

---

## The biology

**N-linked glycosylation** — attachment of a sugar chain to the side-chain
nitrogen (ND2) of an asparagine. Not decoration: it governs folding, secretion,
stability, half-life and immune recognition.

**Sequon** — the sequence motif N-X-S/T that N-linked glycosylation requires:
asparagine, then almost any residue except proline, then serine or threonine.
**Necessary but not sufficient** — most sequons in a proteome are never
glycosylated.

**Occupied** — a sequon with experimental evidence that a glycan is actually
attached there. The opposite is not "unoccupied" but *unknown*, which is the
central difficulty of this project: absence of evidence accumulates in
understudied proteins and is not evidence of absence.

**ND2** — the side-chain nitrogen of asparagine, the atom the glycan bonds to.
Features named "ND2 shell" are centred on it rather than on the backbone,
because that is where the sugar actually sits.

**NXS / NXT** — the two sequon subtypes, by whether the third position is serine
or threonine. Threonine sequons are glycosylated more often, so the two are
never pooled without stratifying.

**Oligosaccharyltransferase (OST)** — the enzyme that performs the attachment.
Present in the secretory pathway, absent from the cytosol and from bacteria,
which is why sequons in those compartments cannot be occupied at all.

---

## Site sets and populations

**`occupied`** — sites with positive occupancy evidence. 332 have usable
structures.

**`internal_control`** — a sequon in a structure where the depositor modelled
glycans at *other* sites but not this one. The most informative negative label
available, because the sugars demonstrably survived preparation and the
depositor demonstrably modelled them. Only 32 sites exist.

**`secretory_unannotated`** — eukaryotic secretory sequons with no occupancy
annotation. Plenty of them (2,296), weak label: roughly half of these proteins
carry a glycoprotein keyword, so the set is contaminated by construction.
Contamination biases toward the null, so it cannot manufacture an effect.

**`bacterial` / `cytosolic`** — diagnostic sets, not tests. These sequons cannot
be occupied for compartment reasons, so they measure what a purely compositional
difference looks like in these features.

**`triplet_core`** — the primary analysis view: every feature in the row
describes the sequon it names. Requires the observed triplet to match the
expected one, all three residues to be structurally located, and the mapping to
be continuous. 2,556 sites.

**`asn_core`** — a wider view requiring only that the asparagine itself was
measured correctly. Valid for features centred on the asparagine; **not** for
+1 or +2 exposure, structure or geometry. 2,624 sites.

**`construct_review`** — everything excluded from `triplet_core`, each row
carrying the reason. Inspected, never tested. 104 sites.

---

## The measured features

Each feature says what it measures, why it is worth measuring, and how it is
obtained.

### Sequence context

**Flanking composition** (`flank_*_fraction`) — the proportion of residues in
each chemical class among the ±5 residues either side of the sequon, *excluding
the sequon itself*. **Why:** the immediate sequence neighbourhood affects how
accessible the site is to the enzyme and how the region folds. **How:** counted
from the sequence; the window is clipped at chain ends, so the denominator
varies between sites.

**`uniprot_residues_after_asn` / `_after_sequon`** — how many residues follow the
site in the full-length protein. **Why:** position along the chain relates to
co-translational timing — a site near the C-terminus meets the enzyme later,
after more of the protein has folded. **How:** from the UniProt sequence length,
not the structure, so a construct that stops short cannot shorten it.

### Local structure

**RSA (relative solvent accessibility)** — how exposed a residue is, from 0
(buried) to about 1 (fully exposed). **Why:** the enzyme has to reach the
asparagine. **How:** Shrake-Rupley accessible surface area divided by a
reference maximum for that residue type, with glycans excluded from the
calculation so that an occupied site's own sugar cannot make it look buried.

**Secondary structure** — helix, sheet or loop at each sequon position. **Why:**
glycosylation favours flexible, accessible regions. **How:** DSSP, coarsened
from its eight classes.

**Phi / psi and backbone region** — the two backbone dihedral angles, and a
coarse region label (alpha-right, beta, alpha-left). **Why:** local geometry
constrains how a glycan can be accommodated. **How:** computed from backbone
atoms, and never across a gap in the model, where the neighbouring atoms belong
to a different stretch of chain. Angles are circular, so they enter comparisons
as a region rather than as a mean.

**Loop-run length, and censoring** — how many consecutive loop residues contain
the asparagine, and whether that run reaches an unresolved boundary. **Why:**
loop length indicates local flexibility. **How:** walked along the DSSP
assignment. A run that reaches the end of the model is a *lower bound* and is
flagged, because unflagged it would bias loop lengths downward exactly where
density is poor.

### Environment around the attachment point

**ND2 shell** (`shell_*_fraction`, `nd2_atoms_8a_*`, `nd2_residues_8a_*`) —
residues with any heavy atom within 8 Å of ND2. **Why:** this is the space the
glycan must occupy, so it is the neighbourhood that matters. **How:** measured
from ND2 rather than from the alpha carbon, and counted separately for the same
chain and other chains, so an oligomer interface or a crystal contact cannot be
read as local sequence context.

**`shell_net_charge`** — formal charge in that shell: D/E count −1, K/R count
+1, histidine is treated as neutral. **Why:** electrostatics affect enzyme
approach.

**`sidechain_neighbour_residues_5a`** — the count of *residues* with a heavy atom
within 5 Å of the asparagine side chain. Named for what it counts: an earlier
name implied atoms and the difference mattered.

**`nearest_aromatic_sidechain_nd2`, `nearest_disulfide_sg_nd2`** — distance from
ND2 to the nearest aromatic side chain, and to the nearest sulfur participating
in a disulfide. **Why:** aromatic stacking and disulfide proximity both
constrain the local environment. **Note:** disulfide coverage differs sharply by
compartment (52% of occupied sites, 1.6% of cytosolic), because disulfides form
in the oxidising secretory pathway — so comparing arms on it compares
missingness.

---

## Names this project uses

**Conditional sequon score** — the probability a model assigns to the sequon
residues when shown the native backbone and surrounding sequence. What the model
*believes*, without generating anything.

**Retention** — whether the sequon survives when a model actually redesigns the
protein. What the model *does*. The two can disagree.

**Contrast** — one occupied site minus its matched control, on whichever
outcome. The unit of the benchmark.

**Resample unit** — a connected group of sites that are not independent, because
they share a protein or an ortholog cluster. Intervals resample these rather
than individual sites; treating 19 sequons in one protein as 19 observations
would narrow every interval for reasons that have nothing to do with the
biology.

**Matched pair** — an occupied site paired with an unoccupied sequon of similar
RSA, neighbour count and hydrophobic fraction, matched exactly on NXS/NXT.
Matching never uses model output, so every model is scored on identical pairs.

**Variant** — a tag naming which run produced a set of numbers
(`esm_if`, `proteinmpnn_index_corrected`). Reading the wrong score file is
silent, so stages refuse to guess which one is meant.

**`model_index`** — a residue's ordinal position in a chain as this project's
parser reads it. **Not** interchangeable with a model's own indexing: see
[`../../docs/methods_sequon_indexing.md`](../../docs/methods_sequon_indexing.md).

**Fixed-sequon context-retention test** — the experiment asking whether
protecting the motif during redesign also protects its surroundings. Short form:
*does fixing the motif protect the biology around it, or only the three letters?*

**D, and ΔD** — D is a site's distance from natural occupied context: the median
absolute standardised departure across the fifteen-feature panel, with the
site's own protein excluded from the reference. ΔD is the paired change from
wild type to design. Positive means further from natural.

**Panels-only** — a mode that computes features without designing anything, used
to build the natural reference, which needs every site but none of the designs.
