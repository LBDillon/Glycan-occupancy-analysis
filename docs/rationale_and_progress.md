# Why this module exists, and what has been built

*A plain-language account of the work, kept current as the module changes. Last updated
2026-08-14, after the full enrichment run (Task 12 of 12). Technical detail appears only where the argument depends
on it; the README covers how to run things and `evidence_sources.md` covers what each
database can and cannot establish.*

---

## The problem in cataloguing sequons by motif

The ortholog sequon-conservation database finds pairs of orthologous proteins where one member carries an N-linked sequon and the other has lost it. Each such comparison is one row. 

Storing the data as one row per comparison suits the question the database was built to
answer — how often sequons are lost across evolution — but it makes a second question hard
to ask, because a sequon is only a motif in the sequence. Whether it is actually used depends on local structure, on membrane topology, and the presence of an oligosaccharyltransferase (OST). The sequence motif being insufficient for a glycan attachment means distinguishing "this pair lost a sequon" and "this pair lost a glycan" is important to make the database biologically relevant. 

The second problem is arithmetic. Take human P00709 with an asparagine at position 64. The database compares that protein against its bovine, murine and porcine orthologs, and each comparison is stored as its own row, meaning position 64 P00709 appears three times. If counting rows, there would be three glycosylation observations where biology has only one asparagine. Across the whole dataset that is the difference between 13,816 rows and 4,307 actual sites.

Those three ortholog relationships do matter, and nothing is thrown away — the intuition in
the second half of that question is the right one. It is the same data indexed twice for two
different questions. Asking "does this asparagine carry a glycan" is a question about one
residue, so the site is counted once. Asking "how often is this sequon lost, and against
which orthologs" is a question about comparisons, so all three rows stay available in
`site_pair_associations.csv`, carrying cluster, source and homology quality. Separating them
means an evidence count cannot be inflated by how many orthologs a protein happens to have
been compared against, while the comparisons themselves remain fully in play.

There is a diagram of this collapse: (https://claude.ai/code/artifact/d0248a3d-8f3e-4ebc-950e-31cd245dc835). 

For every candidate site for a protein in the database we ask if there is experimental evidence that the asparagine specifically carries a glycan. To check if it does we use data for the protein from UniProt and GlyGen, as well as looking at deposited structures and GlyConnect for supporting evidence.

**UniProt** is the primary source, because it is the only one with site-level annotations for per-feature evidence codes. Its glycosylation features are read at exact residue positions only and a range or an uncertain position is rejected.

**GlyGen** aggregates glycosylation data from many labs and databases. It is independent of UniProt where its records cite mass-spectrometry repositories and published papers, and not independent where it is re-exporting UniProt's own predictions. The code separates those cases by reading the sources GlyGen itself cites for each site. Where the only citation is UniProt, the record is not counted; where the citation is a PubMed paper or a mass-spectrometry repository, it counts. Of the sites GlyGen holds, 1,415 are UniProt data for the motif identification.

**Deposited structures**. When a crystallographer models a sugar covalently bonded to a specific asparagine, that is a direct physical observation of occupancy. These bonds are recorded in the structure files already cached, and the code reads them.

**GlyConnect** contributes supporting detail, but its coverage is thin and GlyGen already has most of it.

## Evidence handling

One of the realities of the data we have available, is that an absence of annotation is likely to mean that glycosylation was not screened for. Well-studied proteins accumulate annotations, thus we do not want the database to primarily reflect the current state of glycoprotein curation biased towards such proteins

So each site carries two separate facts: whether it qualifies for the experimentally-supported
set, and what is actually known about its occupancy. For most sites the second answer is
`unknown`. There is a classification for "looked at and found no glycan", but it is
deliberately never populated, because no source available today can establish it. Mass
spectrometry reports the peptides it detected, not the ones it searched for and failed to
find, so a site missing from a glycoproteomics run may have been unoccupied, or may simply
never have produced an observable peptide. Distinguishing those two requires an experiment
designed to test that specific site, and no database records that.

The same applies to structures. A residue visible in a crystal structure with no sugar attached is not automatically classified as an unoccupied site as glycans are routinely trimmed off before crystallisation, proteins are often expressed in bacteria that cannot glycosylate at all, and sugars are frequently too floppy to appear in the density. Instead we record is that the residue is resolved and no glycan is modelled.
 
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

Structural coverage is now bounded by what exists rather than by how we use it. Only 484 of
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

The module is built and tested. Eleven of twelve planned pieces of work are complete and
have passed independent review: configuration, UniProt parsing and evidence grading, the
site universe, the evidence join, the frozen baseline, all four evidence layers, the
combination logic, provenance recording, the command-line interface, and the documentation.

The last step is the full enrichment run — populating the GlyGen cache and freezing a second
regression fixture for the enriched totals. Those enriched numbers are deliberately not
predicted in advance; they will be generated, checked, and only then locked in.

Development turned up several real defects, all of which were caught by review rather than
by luck: output files that would have been written into the wrong directory, a residue-
mapping bug that attributed glycans to the wrong protein chain, a retry loop that would have
crashed on a particular configuration, and an API access pattern that would have taken
thirty-five hours instead of forty minutes. Each is fixed and has a test pinning the
behaviour so it cannot return quietly.
