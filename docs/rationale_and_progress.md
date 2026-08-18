# Why this module exists, and what has been built

*A plain-language account of the work, kept current as the module changes. Last updated
2026-08-14, after the observed-unmodified class and structural features were added. Technical detail appears only where the argument depends
on it; the README covers how to run things and `evidence_sources.md` covers what each
database can and cannot establish.*

---

## The problem in cataloguing sequons by motif

The ortholog sequon-conservation database finds pairs of orthologous proteins where one member carries an N-linked sequon and the other has lost it. Each such comparison is one row.

Storing the data as one row per comparison suits the question the main ortholog database was built to answer of how often sequons are lost across evolution and if we can learn the context changes other than the motif shift. However we need to adapt it for the occupancy analysis because a sequon is only a motif in the sequence. Whether it is actually used depends on local structure, on membrane topology, and the presence of an oligosaccharyltransferase (OST). The sequence motif being insufficient for a glycan attachment means distinguishing "this pair lost a sequon" and "this pair lost a glycan" is important to make the database biologically relevant.

The main reasons is the inflating the statitical signifigance if we treat pairs as the main comparison. Take human P00709 with an asparagine at position 64. The database compares that protein against its bovine, murine and porcine orthologs, and each comparison is stored as its own row, meaning position 64 P00709 appears three times. If counting rows, there would be three glycosylation observations where biology has only one asparagine. Across the whole dataset that is the difference between 13,816 rows and 4,307 actual sites.

Asking "does this asparagine carry a glycan" is a question about one residue, so the site is counted once. Asking "how often is this sequon lost, and against which orthologs" is a question about comparisons, so all three rows stay available in `site_pair_associations.csv`, carrying cluster, source and homology quality. Separating them means an evidence count cannot be inflated by how many orthologs a protein happens to have been compared against.

For the occupancy analysis data, every candidate site for a protein in the database we ask if there is experimental evidence that the asparagine specifically carries a glycan. To check if it does we use data for the protein from UniProt and GlyGen, as well as looking at deposited structures and GlyConnect for supporting evidence.

**UniProt** is the primary source, because it is the only one with site-level annotations for per-feature evidence codes. Its glycosylation features are read at exact residue positions only and a range or an uncertain position is rejected.

**GlyGen** aggregates glycosylation data from many labs and databases. It is independent of UniProt where its records cite mass-spectrometry repositories and published papers, and not independent where it is re-exporting UniProt's own predictions. The code separates those cases by reading the sources GlyGen itself cites for each site. Where the only citation is UniProt, the record is not counted; where the citation is a PubMed paper or a mass-spectrometry repository, it counts. Of the sites GlyGen holds, 1,415 are UniProt data for the motif identification.

**Deposited structures**. When a crystallographer models a sugar covalently bonded to a specific asparagine, that is a direct physical observation of occupancy. These bonds are recorded in the structure files already cached, and the code reads them.

**GlyConnect** contributes supporting detail, but its coverage is thin and GlyGen already has most of it.

## Evidence handling

One of the realities of the data we have available, is that an absence of annotation is likely to mean that glycosylation was not screened for. Well-studied proteins accumulate annotations, thus we do not want the database to primarily reflect the current state of glycoprotein curation biased towards such proteins

So each site carries two separate facts: whether it qualifies for the experimentally-supported
set, and what is actually known about its occupancy. For most sites the second answer is
`unknown`, and that is the honest answer rather than a gap to be filled.

A residue visible in a structure with no sugar attached is not, on its own, an unoccupied
site. Glycans are routinely trimmed before crystallisation, proteins are often expressed in
bacteria that cannot glycosylate at all, and sugars are frequently too mobile to appear in
the density. So by default we record only that the residue is resolved and no glycan is
modelled, which is a statement about the experiment rather than about the biology.

There is one circumstance where that changes, and it is the basis of the third class. If the
same structure models a glycan at some *other* residue, then sugars demonstrably survived
sample preparation and this depositor demonstrably modelled them. If the protein was also
expressed in a host that can glycosylate, a bare asparagine in that structure is a decision
rather than a silence. Sites meeting both conditions are classified `observed_unmodified`:
**32 of them, across 25 proteins**, in structures modelling a median of three glycans
elsewhere. Nothing else can enter that class, and the criteria are pinned by tests so it
cannot quietly widen.

**A note on what to call these.** `observed_unmodified` is the stored class
value, and it claims more than the evidence supports. These sites are not shown
to be chemically unmodified; they are sites with **no modelled glycan under
internal-control conditions** — conditions unusually informative about absence,
but still a statement about the model deposited rather than about the molecule.
They are the best available internal controls and they are not definitive
biochemical negatives. Prose and figures use the longer phrase; the short value
survives in the data files as a legacy identifier.

Mass spectrometry could in principle supply many more. A PNGase F digest in H2-18O converts
occupied asparagines to labelled aspartate, so a sequon peptide detected with the asparagine
intact is a genuine, quantified negative. That is a data-acquisition route rather than a code
one, and it is the obvious way to grow this class beyond 32.
 
Note : The evidence code `ECO:0007744` was not serving the purpose that i had originally thought.

Applying the strictest reading of UniProt evidence to the 4,307 candidate sites gives 505 sites across 401 proteins with direct experimental support. Restricting to the most confident ortholog comparisons — the pairs the homology QC labels
`strict_ortholog_like`, where orthology is supported by more than one line of evidence and
the alignment is high quality — narrows this to 321 sites in 278 proteins.

Reading glycan bonds out of the deposited structures found 184 sites with direct physical evidence of a sugar attached. Forty-four of those fail UniProt's evidence bar, as they are sites where no curator has recorded experimental glycosylation, but have been modelled with crystallisation. Each was checked individually to confirm it maps to a genuine asparagine.

1,714 of the candidate accessions carry a GlyGen cross-reference, and only those were asked
for over the API. That restriction matters practically: GlyGen answers HTTP 500 rather than
404 for an accession it has no entry for, so querying all 2,878 spent the run retrying
failures. Filtering first took the fetch from an estimated 35 hours to 40 minutes. The enriched total is 922 sites across 703 proteins. Of the 417 sites added on top of the UniProt baseline, 383 have GlyGen support and 44 have a structural glycan linkage, with 10 supported by both. Restricting to the most confident ortholog comparisons leaves 396 sites in 333 proteins.

Total: 922 sites across 703 proteins. The most confident ortholog comparisons: 396 sites in 333 proteins.

Occupancy prediction now is testable. We can also ask: of the sequons that get lost between orthologs, are the ones that in fact have a glycan attached lost at a different rate from the ones that merely matched the motif? Yet the comparison will also need controls for how well-studied each protein is.

 
## Limitations
Only 484 of
the 2,878 candidate proteins have any experimental structure, so 3,546 sites cannot be
assessed structurally at all. Within those 484 every available structure is examined: 611
further PDB entries were downloaded for the 123 proteins cross-referencing more than one, and
each site is scored against all of them. That widening found 12 additional glycan-bearing
sites and resolved 21 residues invisible in the first structure, but it also largely exhausts
the deposited structures for this protein set.

Chain matching is conservative. When a structure contains several similar chains, we prioritise a close sequence match before crediting a glycan to one of them. This was because of an issue of crediting sugars to unrelated chains. The current setting can, in principle, miss a glycan
that sits on an engineered mutant when an unmutated copy is also present. This is a deliberate favouring of false negatives over false positives. 

Taxonomic reach is uneven: GlyGen and GlyConnect are strongly biased toward human and mouse, so enrichment will be much stronger for well-studied model organisms than for the rest of the tree

## Where the work stands

The module is built, reviewed and tested: 205 tests, covering configuration, UniProt parsing
and evidence grading, the site universe, the evidence join, the frozen baseline, all four
evidence layers, the combination logic, provenance, the command line, structural feature
extraction, ProteinMPNN scoring, matching and the contrast statistics. Both the UniProt
baseline and the enriched totals are frozen as regression fixtures, so any drift in the
inputs surfaces as a failing test rather than a quietly different number.

The dataset now separates the three classes the occupancy analysis needs:

| Class | Sites |
|---|---|
| `occupied_supported` | 922 |
| `observed_unmodified` | 32 |
| `unknown` | 3,353 |

Structural context features — relative solvent accessibility, exposure bin, neighbourhood
composition, distance to the chain terminus — are computed for the 364 sites with mapped
coordinates (all 32 unmodified, 332 occupied) in `results/site_structural_features.csv`.
Sites without coordinates are retained and flagged rather than filtered out.

One correction is worth recording, because it would have inverted a conclusion. Solvent
accessibility initially included non-protein atoms, so each occupied asparagine was being
shaded by its own glycan. Occupied sites appeared buried 106 times out of 332; with glycans
stripped it is 6 out of 332, and median accessibility moves from 0.185 to 0.433. Burial was
encoding the label. Corrected, occupied sites are overwhelmingly solvent-exposed, which is
what oligosaccharyltransferase access requires.

## The occupancy analysis: done for ProteinMPNN, and frozen

The point of the resource is to ask whether protein models have learned the sequon motif or
its biological use. That question has now been put to one model, and the answer is weak but
honest. Full account in [`primary_result.md`](primary_result.md); the plain-language version
of what changed along the way is in [`correction_2026-08-18.md`](correction_2026-08-18.md).

**What was run.** ProteinMPNN v_48_020, scored zero-shot on the native backbone and native
surrounding sequence, with no training and no alteration of the sequon. The score is the mean
log-odds of asparagine at the first position and of serine-or-threonine at the third,
averaged over eight seeded decoding orders. Occupied sites were compared with internal
controls matched on relative accessibility, neighbourhood packing and hydrophobic fraction,
with NXS/NXT required to be identical.

**What it found.** Occupied sites tend to score higher — every point estimate under every
matching we tried is positive — but 16 matched pairs do not establish a precise or
statistically robust difference. The primary estimate is +0.458 SD with a 95% interval of
[−0.227, +1.098], which includes zero.

**Why the number of pairs is so small.** 32 internal controls; 28 survive the requirement
that ProteinMPNN can actually decode all three sequon residues; 16 find a partner within the
matching caliper. That attrition, not the modelling, is the constraint on this analysis.

**Two corrections worth carrying forward.** ProteinMPNN silently returns a row of ones rather
than a probability distribution for any residue with an incomplete backbone, which inflated
the reference scale and reversed the sign of the first result. And greedy matching made the
significance of the answer depend on a random seed. Both are fixed, both are recorded as
formal amendments, and the analysis is frozen at that point.

**The 32 remain a calibration probe, not a training class.** Positive-unlabelled learning, as
in the KinoPlex kinase atlas, is still the intended route for the 3,353 unlabelled sites:
roughly 116,000 known phosphosites against 1.7 million uncharacterised residues there, with
an Elkan-Noto correction to recover true performance. This resource is the same shape — 922
positives against 3,353 unlabelled. The 32 earn their place by testing whether the unlabelled
pool behaves like a mixture containing negatives, which is the assumption the method rests on.
That work has not begun, deliberately: it should not start while the zero-shot result is this
imprecise.

Two cautions carry into it. The 32 skew exposed, 20 of 32, so matching is done pair by pair
rather than assumed from group averages — which is what the analysis does. And they exist
only where someone solved a glycoprotein structure, so they represent well-studied secreted
and membrane proteins rather than a random draw from the unlabelled pool, which matters
precisely because that pool is much broader.

## What would actually move this forward

Not more models. The binding constraint is 16 pairs, and no amount of ESM-IF or TriFlow
changes that — running a panel now would produce several imprecise answers instead of one.

The thing that would change it is more internal controls, and that is a data-acquisition
problem: a PNGase F digest in H₂¹⁸O converts occupied asparagines to labelled aspartate, so a
sequon peptide detected with its asparagine intact is a genuine, quantified negative. That is
the route from "tends to score higher" to a result worth defending.
