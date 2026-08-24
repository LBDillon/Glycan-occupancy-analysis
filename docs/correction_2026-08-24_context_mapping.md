# Correction, 2026-08-24 — sequence-to-structure mapping in the context extractor

The v2 context extractor measured the wrong residue in six different ways. Every
context feature produced before this date has been regenerated. No biological
conclusion changes, because none had been drawn yet — the extraction had
finished three days earlier and nothing downstream had consumed it.

## What the extractor is for

The occupancy benchmark asked whether a model treats an experimentally occupied
sequon differently from a structurally matched sequon carrying no glycan. It
answered that: ProteinMPNN does not, ESM-IF does. It could not say *what* about
those sites a model is responding to, because it never described the sites.

The context atlas is that description. For each of 2,660 sites it records what
the sequon's neighbourhood actually looks like — solvent exposure, secondary
structure, backbone geometry, what sits within 8 Å of the attachment point. The
whole point is that a feature attributed to a site must belong to that site. A
description that is subtly about the wrong residue is worse than no description,
because it looks like evidence.

## What was wrong

**Insertion codes never reached the extractor.** Some depositions number
residues 36, 36A, 36B, 36C — chymotrypsin-numbered proteases do it routinely.
The evidence table recorded which one a site sat on; the manifest dropped the
column, and the extractor could not pass what it never received, so it asked for
"residue 36" and got whichever one had a blank insertion code. For Q99895 at
position 52 in 4H4F the sequon is at 36B and the extractor read 36, reporting
the triplet `LKN` where the site is `NDT`. One site in the whole set is affected.
It was still measuring a different residue than the one it named.

**The DSSP table was keyed on residue number with the insertion code thrown
away.** Independent of the above, and it matters more than it sounds: an
insertion block collapsed onto a single entry, so whichever residue DSSP emitted
last supplied the secondary structure for all four. Plumbing insertion codes
through the manifest would not have fixed this.

**+1 and +2 were taken as the next *resolved* residues.** This was deliberate,
and the docstring explained why: adding 1 and 2 to a residue number is wrong when
insertion codes exist. But walking to the next resolved residue is wrong for a
different reason — when the deposition never observed a residue, the walk steps
over the gap and lands somewhere else entirely. P17936 at position 116 sits at
residue 89, and the next residue present in the model is 182. That 93-residue
jump was recorded as the sequon's +1.

Twenty-five sites were affected. Nine of them passed the triplet check, because
the residue on the far side of the gap happened to have the same identity as the
one expected. Those nine are the reason the triplet check cannot be the only
guard: it compares letters, and letters collide.

**Backbone dihedrals used those same non-adjacent neighbours.** Phi needs the
preceding residue's carbonyl carbon and psi the following residue's nitrogen.
Across a gap those atoms belong to a different stretch of chain, and the angle
computed from them is not a torsion of anything.

**Terminal distances were author-number arithmetic.** `chain_length_resolved`
counted residues; `distance_to_n_terminus_resolved` subtracted author residue
numbers. The two are only commensurate when a chain is numbered without gaps,
which most are not. 1,218 of 2,660 rows failed d_N + d_C + 1 = chain length.

**DSSP failed outright on every multi-character chain identifier.** Forty-three
sites, all of them in large assemblies — precisely where the control arms draw
most heavily, so the loss was not evenly spread. The cause is that the legacy
DSSP output format has a one-column chain field and Biopython always requests
that format. Converting the input to mmCIF is not enough; mkdssp refuses to write
a chain called `AB` into a format that cannot hold it. The fix is to relabel the
single-chain extract before DSSP sees it, since secondary structure comes from
geometry rather than from the label.

## How it was found, and what I got wrong

Not by inspection. An external audit of the finished table found four of the six,
and I found the other two while verifying its claims.

I had reported the opposite conclusion three days earlier. Looking at 96 triplet
mismatches out of 2,660, I classified them and concluded there was no systematic
mapping bug — the signature of one would be a coherent positional shift, and I
counted only three. That reasoning was sound and the conclusion was wrong, for a
reason worth recording: **I was looking for the defect in the rows the QC column
had already flagged.** The nine most damaging cases were not in that set. They
passed. A check that can only see disagreements is blind to a defect that
produces agreement, and I had implicitly trusted it to be exhaustive.

I also asserted that 44 mismatches with a correct asparagine were "isoform or
construct differences". I had not demonstrated that, and about half of them were
gap-jumps. The audit was right to say so.

Two numbers in the audit needed correcting in turn. It reported 182 inconsistent
terminal distances, counting rows where a distance was negative or exceeded the
chain; the invariant d_N + d_C + 1 = chain length catches 1,218. And its count of
N→Q substitutions differed from mine because it counted rows whose first position
is Q while I counted rows where *only* the first position differs — both correct
under their own definition, neither an error.

## What changed, and what did not

Every one of the 1,259 changed rows is attributable to a named correction, with
none left over:

| Correction | Sites |
|---|---|
| Terminal distances now count residues | 1,167 |
| DSSP recovered on multi-character chains | 43 |
| Dihedrals withdrawn across gaps | 23 |
| +1/+2 no longer step over gaps | 16 |
| Previously invisible gap-jumps | 9 |
| Insertion code propagated | 1 |
| **Unexplained** | **0** |

That last row is the one that matters. An audit that merely summarises a diff
cannot tell you whether something else moved while you were not looking; this one
assigns each changed row to a cause and reports anything it cannot explain.

The populations are unchanged: 332 occupied sites, 32 internal controls, 2,296
secretory-unannotated. The frozen model benchmark was not touched and the
original outputs were not discarded.

## What else was built, and why

Fixing the six defects was not sufficient, because each had reached a finished
table without anything objecting.

**The feature panel was drifting from its own specification.** The specification
existed only as a sentence in a chat log — "something like 15 to 20 numbers per
site" — and had never been written down. One feature had drifted outright:
`sidechain_contacts_5a` counted *residues* with a heavy atom near the asparagine
side chain, while its name read as a count of atom contacts, and the measurement
the name implied had never been implemented. It is now
`sidechain_neighbour_residues_5a`, and the count of non-hydrogen atoms within 8 Å
of ND2 has been added rather than substituted. Same-chain and other-chain
contributions are counted separately, so an oligomer interface cannot be read as
local sequence context. Loop-run length is recorded with a censoring flag,
because a loop that runs off the end of the model has no measured boundary and
reporting its length as observed biases loop lengths downward exactly where
density is poor.

**Technical quality was mixed in with biology.** Resolution, experimental method,
per-position DSSP availability, residue numbers, insertion codes and mapping
continuity now sit in a separate block from the biological panel. P17936 sits in
a 3.6 Å cryo-EM map, which is *why* its +1 and +2 were never resolved — an
analysis has to be able to see that rather than inferring a short loop.

**Inputs are resolved by configuration.** The repository holds code, tests and
generated tables; the UniProt release and structure cache are large and live
outside it. Stages reached them with literal `../../data/...` paths, which is why
the corrected rerun could not be started from this repository at all.

Making that change immediately introduced a regression, which is the most
instructive part of this whole exercise. The control sequences load from a
relative `data/cache/` path guarded by `if Path(cache).exists()`. In the new
location it did not exist, the guard skipped it silently, and the manifest came
out with no sequence context for its largest population — meaning the triplet
check was inoperative for 2,296 of 2,660 sites, with nothing in the output to say
so. The same class of failure as the original six: not a crash, a quiet wrong
answer. Stage 41 now names every cache it resolves and refuses to write a
manifest where any population has no sequence context, and the sequence checks
run per row, so a cache that loads *half* is detectable too.

## The gates

Nothing reaches an analysis without passing these, and each one exists because
its absence let something through:

- Shards merge on the full key including `population`; missing shards, duplicate
  keys, recorded failures and short coverage are all fatal rather than warnings.
- Every row must have a UniProt sequence, a complete N-X-S/T triplet with X ≠ P,
  and coordinates inside the sequence. A site that is not a sequon has been
  mapped wrongly upstream and cannot be an N-linked site.
- The invariants are asserted, not described, and the report exits non-zero when
  one breaks. A report that only prints cannot stop a bad table being used.

## The three views

One table cannot serve the atlas, because the reasons a site is imperfect are not
interchangeable. A crystallographer's N→Q knockout, a +1 differing between
isoform and construct, and a +2 that was never resolved are three different
facts.

- **`triplet_core`** (2,556 sites) — triplet agrees, all three residues located,
  mapping continuous. It requires continuity *as well as* agreement, because
  requiring agreement alone readmits the nine sites whose +1 came from across a
  gap.
- **`asn_core`** (2,624) — the asparagine was measured correctly. Valid for
  features centred on it; not for +1 or +2 exposure, structure or geometry.
- **`construct_review`** (104) — everything excluded, each row carrying the
  reason: 36 substituted asparagines, 30 sequence substitutions, 38 unresolved
  positions.

## Known loose end

The `discontinuous_mapping` exclusion category cannot fire. Any discontinuity
produces an unresolved position, so `unresolved_position` always claims the row
first — verified as zero rows across the whole table. It is dead code. The
distinction it was meant to capture is a property of the *change* rather than of
the new row, and lives correctly in the change audit.

## What this does not settle

The extractor is now trustworthy; the atlas is not yet written. Nothing here says
anything biological. The first real use of this table is to ask which contextual
features distinguish occupied sites from controls, and whether those are the
features ESM-IF responds to and ProteinMPNN misses.

Two panel items remain underspecified rather than implemented: residues after the
asparagine and after the sequon are joined from the manifest at view-build time
rather than being extractor outputs, which is correct but means they exist in the
views and not in the feature table. And `nearest_disulfide_sg_nd2` has coverage
of 51% in occupied sites against 15% in secretory ones. That is real biology
rather than a defect, but the gap between arms is wide enough that it should not
enter a comparison without being stated.
