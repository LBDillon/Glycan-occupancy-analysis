# Negative control sets

## Why these exist

The resource holds 922 sites with experimental evidence of a glycan, 32 observed
unmodified under an informative experimental context, and 3,353 whose occupancy is
unknown. That is enough to ask whether a protein model has learned the sequon motif or
its biological use — but not enough to answer it convincingly, because 32 negatives
cannot support a null result, and because any comparison between occupied and
unoccupied sites is open to the objection that the two groups differ in something other
than occupancy.

These control sets exist to close that gap. They are **not** part of the dataset. They
are never pooled with the 4,307 candidate sites, never counted in any site total, and
always carry a `control_set` label recording where they came from.

## The problem they solve

Occupied sites are not a random sample of sequons. They sit in secreted and membrane
proteins, in solvent-exposed loops, in well-studied organisms. Any of those properties
can be read by a structure-based model, so a model that separates occupied from
unoccupied sites may be detecting one of them rather than occupancy itself.

The defence is not a better single control. It is several controls whose confounds do
not overlap, so that no one alternative explanation can account for a result that holds
across all of them.

| Control set | Matches the positives on | Confounded by | Approx. sequons |
|---|---|---|---|
| Structural internal controls (the 32) | organism, compartment, experiment | very small; structurally biased | 32 |
| Cytosolic eukaryotic | taxonomy, curation depth | subcellular compartment | 19,337 |
| Bacterial extracytoplasmic | compartment, membrane translocation | taxonomy, fold repertoire | 5,865 |

Read them together. If a model separates occupied sites from the cytosolic set but not
from the 32, it has learned where proteins live. If it separates them from the bacterial
set but not from the cytosolic set, it has learned taxonomy. **The gap between the
comparisons is the measurement**, not any single comparison.

## Set 1 — cytosolic eukaryotic

Proteins localised to the cytosol, with no signal peptide and no transmembrane region.

N-linked glycosylation happens in the ER lumen, catalysed by oligosaccharyltransferase.
A protein that never enters the secretory pathway never meets OST, so its sequons are
not merely unobserved — they are unreachable. This is positive knowledge from the
biology rather than an inference from missing annotation, which makes these stronger
negatives than anything the annotation databases can supply.

The cost is that the comparison now differs in compartment. Secreted and cytosolic
proteins differ in composition, disulfide content, hydrophobic patterning and fold
family, all visible to a structural model.

**A precision point:** cytosolic and nuclear proteins *are* glycosylated, by O-GlcNAc
transferase, on serine and threonine. That is different chemistry and does not
contaminate an N-linked negative set, but the claim must be written as "never
N-glycosylated", not "not glycosylated".

## Set 2 — bacterial extracytoplasmic

Periplasmic, outer-membrane and secreted bacterial proteins, from clades with no known
N-glycosylation machinery.

This fixes what is wrong with set 1. These proteins cross a membrane, carry signal
peptides, and fold in an oxidising, disulfide-forming compartment — the periplasm is the
closest bacterial analogue to the ER lumen. What differs from the positives is the
presence of an OST, which is much closer to isolating the variable of interest.

The cost is taxonomic: bacterial and eukaryotic secreted proteins differ in fold
repertoire, composition, and evolutionary history with respect to glycans.

### Exclusions, and why each one

Bacterial N-glycosylation is rare but real, and excluding it by annotation alone would
repeat the absence-of-evidence error this project exists to avoid. Clades are therefore
excluded by **known machinery**:

| Excluded | Reason |
|---|---|
| Archaea (all) | Glycosylate via AglB; a genuine OST, not an edge case |
| *Campylobacter* | PglB, the best-characterised bacterial OST |
| *Helicobacter* | PglB-family machinery |
| *Haemophilus* | HMW1C-type cytoplasmic N-glycosyltransferase, acts on N-X-S/T |
| *Actinobacillus* | HMW1C homologue |
| *Yersinia* | HMW1C homologue |
| *Kingella* | HMW1C homologue |

The annotation filter is kept as a second line of defence, not as the primary one.

## A second use, independent of the modelling

In an organism with no OST, sequons are under no glycosylation-related selection. The
bacterial set therefore describes the **neutral sequon distribution** — what sequon
density and local context look like with the selective pressure removed.

That makes it a null model for the evolutionary half of the project: whether sequons in
eukaryotic secreted proteins are enriched or depleted relative to neutral expectation, at
matched amino acid composition. This question needs no model scoring at all.

## What the built sets look like

| Set | Proteins | Sequons | Per 1000 residues |
|---|---|---|---|
| Cytosolic eukaryotic | 9,561 | 19,337 | 4.54 |
| Bacterial extracytoplasmic | 1,410 | 5,865 | 8.75 |

Bacterial extracytoplasmic proteins carry nearly twice the sequon density, and that is
composition rather than biology: they are richer in asparagine (0.060 against 0.039) and
threonine (0.069 against 0.052), and poorer in the proline that blocks a sequon.

The useful check is against the density expected from composition alone, treating each
position independently:

| Set | Expected | Observed | Ratio |
|---|---|---|---|
| Cytosolic eukaryotic | 4.55 | 4.54 | 1.00 |
| Bacterial extracytoplasmic | 8.26 | 8.75 | 1.06 |

Both sit essentially at neutral expectation, which is what proteins under no
glycosylation-related selection should do. That is a direct check on the premise these
sets rest on, and it passed. It also sets up the comparison that matters for the
evolutionary question: whether eukaryotic *secreted* proteins deviate from the same
expectation, which would be selection on sequons rather than chance.

## Rules these sets must obey

1. **Never pooled.** Each carries a `control_set` label and lives in its own table.
   They are not sites in the resource and are excluded from every site count.
2. **Never counted as evidence.** They say nothing about the 4,307 candidates.
3. **Experimental structures only.** If positives used deposited structures and controls
   used predicted ones, structure source would confound compartment and taxonomy on top
   of everything else.
4. **Composition reported.** Amino acid usage differs between the sets, so sequon density
   differs for reasons unrelated to glycosylation. Any comparison that depends on it must
   match on it.
5. **Outside the ortholog framework.** These proteins have no ortholog-pair context and
   cannot participate in the sequon-loss analysis.

## Considered and rejected

*Plasmodium falciparum* would control taxonomy and compartment at once, being a eukaryote
with almost no N-glycosylation. Its proteome is extremely AT-rich and low-complexity and
few structures exist, so it trades a clean confound for a messy one.

## How to build them

Stage 1 — sequence inventory and sequon scan. Fast, no structure downloads:

```bash
python -m experimental_glycosylation_sites fetch-controls
```

Writes `results/negative_control_sites.csv` (one row per sequon, labelled by
`control_set`) and `results/negative_control_summary.json` (counts, sequon density,
amino acid composition, and the excluded taxa with their reasons). Protein sequences go
to `data/cache/negative_control_proteins.csv.gz`, since they are bulky and only the
feature stage needs them.

**Stage 2 — structural features — is not built yet.** It needs deposited structures for
the control proteins, which is a download on the order of thousands of entries, so it is
deliberately separated from stage 1. The feature code itself already exists and is
shared with the main dataset (`features.residue_features`), so the work is fetching the
structures and pointing the existing pipeline at this table. Until that runs, the control
sets support sequence-level and composition analyses but not the structure-based scoring
comparison.
