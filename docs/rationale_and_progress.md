# Why this module exists, and what has been built

*A plain-language account of the work, kept current as the module changes. Last updated
2026-08-14, after the full enrichment run (Task 12 of 12). Technical detail appears only where the argument depends
on it; the README covers how to run things and `evidence_sources.md` covers what each
database can and cannot establish.*

---

## The problem in cataloguing sequons by motif

The ortholog sequon-conservation database asks an evolutionary question. It finds pairs of orthologous proteins where one member carries an N-linked sequon, which is the N-X-S/T motif that is necessary (but not sufficient) for N-linked glycosylation, and the other has lost it. Each such comparison is one row.

Storing the data as one row per comparison suits the question the database was built to
answer — how often sequons are lost across evolution — but it causes two problems as soon as
you want to say anything about glycans.

The first is that a sequon is only a motif in the sequence. Whether it is actually used
depends on local structure, on membrane topology, and on whether the protein ever meets an
oligosaccharyltransferase (OST). So "this pair lost a sequon" and "this pair lost a glycan"
are different claims, and only the first is supported by the motif alone.

The second problem is arithmetic, and it is easiest to see with an example. Take human
P00709 with an asparagine at position 64. The database compares that protein against its
bovine, murine and porcine orthologs, and each comparison is stored as its own row — so
position 64 appears three times. Nothing about the protein has changed; only the number of
comparisons has. If you counted rows you would record three glycosylation observations where
biology has only one asparagine. Across the whole dataset that is the difference between
13,816 rows and 4,307 actual sites, and any statistic built on the row count would be
roughly three times more confident than the evidence allows.

There is a diagram of exactly this collapse, along with the full provenance of every source,
at the [evidence provenance
page](https://claude.ai/code/artifact/d0248a3d-8f3e-4ebc-950e-31cd245dc835).

The module therefore reorganises the evidence around one protein and one residue position.
Ortholog pairs are still kept, but in a separate table, so they provide context without
inflating any count.

For every candidate site for a protein in the database we ask if there is
experimental evidence that the asparagine specifically carries a glycan. To check if it does we use data for the protein from UniProt and GlyGen, as well as looking at deposited structures and GlyConnect for supporting evidence.

**UniProt** is the primary source, because it is the only one with site-level annotations for per-feature evidence codes. Its glycosylation features are read at exact residue positions only and a range or an uncertain position is rejected.

**GlyGen** aggregates glycosylation data from many labs and databases. It is independent of
UniProt where its records cite mass-spectrometry repositories and published papers, and not
independent where it is re-exporting UniProt's own predictions. The code separates those
cases by reading the sources GlyGen itself cites for each site. Where the only citation is
UniProt, the record is treated as an echo and carries no weight; where the citation is a
PubMed paper or a mass-spectrometry repository, it counts. This matters more than it sounds:
of the sites GlyGen holds, 1,415 turn out to be UniProt's own sequence rule coming back
round, and admitting them would have quietly put predictions into a dataset whose entire
value is that it contains none.

**Deposited structures**. When a crystallographer models a sugar covalently bonded to a specific asparagine, that is a direct physical observation of occupancy. These bonds are recorded in the structure files already cached, and the code reads them.

**GlyConnect** contributes supporting detail, but its coverage is thin and GlyGen already has most of it.

## Why the evidence handling is so fussy

There is one mistake this module is built to avoid, and it is an easy one to make.

It is tempting to treat "annotated in UniProt" as glycosylated and "not annotated" as not
glycosylated, then compare the two groups. But absence of annotation almost always means
*nobody looked*. Well-studied proteins accumulate annotations; obscure ones do not. A study
built on that comparison would be measuring how much attention each protein has received
from curators, and reporting the result as glycobiology.

So the module keeps two separate facts about every site: whether it qualifies for the
experimentally-supported set, and what is actually known about its occupancy. For most
sites the honest answer to the second question is *unknown*. There is a slot for "looked at
and found bare", but it is deliberately empty — no source available today can establish
that, and the code has a test asserting the module never claims otherwise.

The same care applies to structures. A residue visible in a crystal structure with no sugar
attached is not an unoccupied site: glycans are routinely trimmed off before
crystallisation, proteins are often expressed in bacteria that cannot glycosylate at all,
and sugars are frequently too floppy to appear in the density. So the module records "the
residue is resolved and no glycan is modelled" as exactly that, and never as evidence of
absence.

One older convention in this repository had to be corrected along the way. The evidence code
`ECO:0007744` had been labelled as though a deposited structure stood behind it, implying
the site had been seen in one. It does not mean that — it means a curator combined
computational and experimental evidence. Seven sites rest on it alone, and they are flagged so
they can be set aside in a sensitivity check. Nothing in the module now describes them as
structurally observed.

## What has been found so far

Applying the strictest reading of UniProt evidence to the 4,307 candidate sites yields **505
sites across 401 proteins** with direct experimental support. Restricting to the most
confident ortholog comparisons narrows this to 321 sites in 278 proteins. These figures are
frozen as a regression fixture: if a future run disagrees, the tests fail loudly and point at
which input changed, rather than silently adopting the new number.

Reading glycan bonds out of the cached structures then found **172 sites** with direct
physical evidence of a sugar attached. Thirty-two of those fail UniProt's evidence bar
entirely — sites where no curator has recorded experimental glycosylation, but a
crystallographer has modelled the sugar sitting on the residue. Each was checked
individually to confirm it maps to a genuine asparagine. That is the multi-layer design
earning its place: a third of a hundred sites recovered that any single-source approach
would have missed.

The GlyGen layer has since been populated — 1,714 of the candidate accessions carry a GlyGen
cross-reference, and only those were requested — bringing the enriched total to **912 sites
across 697 proteins**. Of the 407 sites added on top of the UniProt baseline, 383 have GlyGen
support and 32 have a structural glycan linkage, with 8 supported by both. Restricting to the
most confident ortholog comparisons leaves 396 sites in 333 proteins. These enriched figures
are frozen as their own regression fixture beside the UniProt baseline, so a change in the
GlyGen cache or the cached structures shows up as a failing test rather than a quietly
different number.

## Why this matters beyond the immediate count

Two things follow from having this dataset.

The first is that occupancy prediction becomes testable. Models that predict whether a
sequon is used need ground truth to be measured against, and the scarce ingredient has never
been sequences — it has been sites whose glycosylation status is known and whose provenance
can be audited. This gives a positive set where every entry can be traced to the evidence
behind it.

The second, and the reason the ortholog database was built in the first place, is that the
evolutionary question becomes answerable in a sharper form. Of the sequons that get lost
between orthologs, are the ones that were actually *used* lost at a different rate from the
ones that merely matched the motif? That question is only askable once occupancy is a
site-level fact with a stated evidence standard rather than an assumption. The honest
caveat is that the annotation bias described above does not disappear just because it is
acknowledged — any such comparison will need controls for how well-studied each protein is,
which is why the analysis guide treats that as the first problem to solve rather than a
footnote.

The dataset is also a standalone artefact. It references the ortholog database rather than
copying it, records the exact provenance of every input, and regenerates deterministically
from cache, so it can travel with a paper without dragging the whole pipeline behind it.

## What is known to be imperfect

Three limitations are worth stating plainly.

Structural coverage is thin. Only a minority of these proteins have a cached structure at
all, and the module currently examines one structure per protein even where many are
available. Widening that is the clearest route to more structural evidence.

Chain matching is deliberately conservative. When a structure contains several similar
chains, the module insists on a close sequence match before crediting a glycan to one of
them. This was tightened three times during development, because the first version credited
sugars to entirely unrelated chains. The current setting can, in principle, miss a glycan
that sits on an engineered mutant when an unmutated copy is also present. That trade is
intentional: the module prefers to miss evidence than to invent it.

Taxonomic reach is uneven. GlyGen and GlyConnect are strongly biased toward human and mouse,
so enrichment will be much stronger for well-studied model organisms than for the rest of
the tree — which is itself a form of the annotation bias the whole design is guarding
against.

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
