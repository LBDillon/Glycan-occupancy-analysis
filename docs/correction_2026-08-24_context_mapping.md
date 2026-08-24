# Correction, 2026-08-24 — sequence-to-structure mapping in the context extractor

The v2 context extractor measured the wrong residue in six different ways. Every
context feature produced before this date has been regenerated. No biological
conclusion changes, because none had been drawn yet — the first extraction
finished on 23/08 and nothing downstream had consumed it.

This note assumes no familiarity with the project, because the reason the defect
mattered is inseparable from what the table was going to be used for.

---

## Where this sits

Protein design models are now good enough to be used. They are not obviously good
enough to be trusted with glycoproteins.

N-linked glycosylation is not decoration. A glycan attached at an N-X-S/T sequon
governs folding, secretion, stability, half-life and immune recognition. A design
tool that silently removes one, or adds a sequon that never gets used, has changed
the protein in ways its scoring function never accounted for. And they do remove
them: unconstrained ProteinMPNN retains the sequon in roughly **4%** of designs
across the 84 glycoproteins in the glyco-bias benchmark. SugarFix, the wider
project, exists to make redesign glycosylation-aware.

That programme is too large to attack directly, so it decomposes into four stages:

1. **Observe what the models have learned** — without correcting them, measure
   what they currently distinguish.
2. **Describe the biology they ought to represent** — build the reference that
   model behaviour can be judged against.
3. **Determine whether the missing biology matters** — a model can fail a
   benchmark without that failure causing design errors.
4. **Test whether grounding the design process helps** — interventions, then
   wet-lab.

This repository is stages 1 and 2. Stage 1 finished on 23/08. Stage 2 is the
glyco-site context analysis, and this correction is about the machinery that
produces it.

## What came before: the occupancy benchmark

Stage 1 asked the narrowest useful question:

> Does a design model treat an experimentally **occupied** N-X-S/T sequon
> differently from a structurally matched sequon carrying **no glycan**?

Two features of that design make it hard and honest. **Both arms carry the
motif**, so a model that has merely learned "N-X-S/T is special" scores them
identically and anything separating them comes from elsewhere. And **absence of
annotation is never treated as absence of glycan** — a site with no evidence is
`unknown`, not negative, because well-studied proteins accumulate annotations and
obscure ones do not. Only sites where a structure models glycans *elsewhere* but
not here are called controls. That honesty is expensive: it leaves 32 internal
controls (16 pairs), so a larger, weaker-labelled secretory set carries most of
the statistical weight.

The answer was architecture-dependent. In the secretory comparison, **ESM-IF1**
scores occupied sequons **+0.431 SD** (BH 0.0008) and retains the motif **+0.0925**
more often — 15.1% against 5.9% (BH 0.015). **ProteinMPNN** gives +0.090 SD and
+0.0423, neither surviving correction (BH 0.30 for both). The direction agrees
everywhere; ESM-IF resolves something ProteinMPNN also weakly shows.

A masking control sharpened it. Hide the whole sequon and ask what the
surroundings alone say: sequence-only **ESMC loses its effect entirely** (+0.12
falls to +0.007, a change of +0.113, 95% CI [+0.071, +0.147], p < 0.0001). It was
judging how well the motif fits its context, not reading context by itself.

Two corrections preceded this one. On 18/08, two rounds in a single day: a
defect in the scorer, and an arbitrary dependence removed from the matching. On
20/08 a larger one: `mpnn_scoring.ALPHABET` held a
three-letter lookup table from inside `parse_PDB_biounits` rather than
ProteinMPNN's token alphabet, so `p_asn_at_n` was reading **P(aspartate)**. Every
score and retention figure was regenerated. It changed magnitudes, not
conclusions.

That is three silent wrong answers in seven days, none of which crashed. It is
the recurring failure mode of this codebase and the reason this note is as long
as it is.

## Why the context analysis, and why it has to be exact

The benchmark measures *whether* models discriminate. It cannot say *what they
discriminate on*. ESM-IF separates occupied from unoccupied sequons sharing an
identical motif, so it is responding to the surrounding sequence and structure —
but "the surroundings" is not a measurement. Until we can say which properties of
a site's environment differ between occupied and unoccupied sequons, three things
stay out of reach: explaining the difference between the two models, evaluating
whether a *designed* site sits in a plausible environment (plausible compared to
what?), and distinguishing a real signal from a compositional confound.

The context analysis is that description. For each of 2,660 sites — 332 occupied, 32
internal controls, 2,296 secretory-unannotated — it records what the sequon's
neighbourhood actually looks like: solvent exposure, secondary structure,
backbone geometry, what sits within 8 Å of the attachment point.

Which is why the defect mattered. **A feature attributed to a site must belong to
that site.** It is going to be the reference against which model behaviour
and designed sites are judged, so a description that is subtly about the wrong
residue is worse than no description at all — it looks like evidence.

---

## What was wrong

**Insertion codes never reached the extractor.** Some depositions number residues
36, 36A, 36B, 36C — chymotrypsin-numbered proteases do it routinely. The evidence
table recorded which one a site sat on; the manifest dropped the column, and the
extractor could not pass what it never received, so it asked for "residue 36" and
got whichever one had a blank insertion code. For Q99895 at position 52 in 4H4F
the sequon is at 36B and the extractor read 36, reporting the triplet `LKN` where
the site is `NDT`. One site in the whole set is affected. It was still measuring a
different residue than the one it named.

**The DSSP table was keyed on residue number with the insertion code thrown
away.** Independent of the above, and it matters more than it sounds: an
insertion block collapsed onto a single entry, so whichever residue DSSP emitted
last supplied the secondary structure for all four. Plumbing insertion codes
through the manifest would not have fixed this.

**+1 and +2 were taken as the next *resolved* residues.** This was deliberate,
and the docstring explained why: adding 1 and 2 to a residue number is wrong when
insertion codes exist. But walking to the next resolved residue is wrong for a
different reason — when the deposition never observed a residue, the walk steps
over the gap and lands somewhere else. P17936 at position 116 sits at residue 89,
and the next residue present in the model is 182. That 93-residue jump was
recorded as the sequon's +1.

Twenty-five sites were affected. Nine passed the triplet check, because the
residue on the far side of the gap happened to have the same identity as the one
expected. Those nine are why the triplet check cannot be the only guard: it
compares letters, and letters collide.

**Backbone dihedrals used those same non-adjacent neighbours.** Phi needs the
preceding residue's carbonyl carbon and psi the following residue's nitrogen.
Across a gap those atoms belong to a different stretch of chain, and the angle
computed from them is not a torsion of anything.

**Terminal distances were author-number arithmetic.** `chain_length_resolved`
counted residues; `distance_to_n_terminus_resolved` subtracted author residue
numbers. The two are only commensurate when a chain is numbered without gaps,
which most are not. 1,218 of 2,660 rows failed d_N + d_C + 1 = chain length.

**DSSP failed outright on every multi-character chain identifier.** Forty-three
sites, all in large assemblies — precisely where the control arms draw most
heavily, so the loss was not evenly spread. The legacy DSSP output format has a
one-column chain field and Biopython always requests that format. Converting the
input to mmCIF is not enough; mkdssp refuses to write a chain called `AB` into a
format that cannot hold it. The fix is to relabel the single-chain extract before
DSSP sees it, since secondary structure comes from geometry, not from the label.

## How it was found, and what I got wrong

Not by inspection. An external audit of the finished table found four of the six,
and I found the other two while verifying its claims.

I had reported the opposite conclusion the day before. Looking at 96 triplet
mismatches out of 2,660, I classified them and concluded there was no systematic
mapping bug — the signature of one would be a coherent positional shift, and I
counted three. The reasoning was sound and the conclusion was wrong, for a reason
worth recording: **I looked for the defect only among rows the QC column had
already flagged.** The nine most damaging cases were not in that set. They passed.
A check that can only see disagreements is blind to a defect that produces
agreement, and I had implicitly trusted it to be exhaustive.

I also asserted that 44 mismatches with a correct asparagine were "isoform or
construct differences". I had not demonstrated that, and about half were
gap-jumps. The audit was right to say so.

Two of the audit's numbers needed correcting in turn. It reported 182 inconsistent
terminal distances, counting rows where a distance was negative or exceeded the
chain length; the invariant d_N + d_C + 1 = chain length catches 1,218. And its
count of N→Q substitutions differed from mine because it counted rows whose first
position is Q while I counted rows where *only* the first position differs — both
correct under their own definition, neither an error.

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

The populations are unchanged. The frozen model benchmark was not touched and the
original outputs were not discarded.

## What else was built, and why

Fixing the six defects was not sufficient, because each had reached a finished
table without anything objecting.

**The feature panel was drifting from its own specification.** That specification
existed only as a sentence in a chat log — "something like 15 to 20 numbers per
site" — and had never been written into the repository. One feature had drifted
outright: `sidechain_contacts_5a` counted *residues* with a heavy atom near the
asparagine side chain, while its name read as a count of atom contacts, and the
measurement the name implied had never been implemented. It is now
`sidechain_neighbour_residues_5a`, and the count of non-hydrogen atoms within 8 Å
of ND2 has been added rather than substituted. Same-chain and other-chain
contributions are counted separately, so an oligomer interface cannot be read as
local sequence context. Loop-run length is recorded with a censoring flag, because
a loop running off the end of the model has no measured boundary, and reporting
its length as observed biases loop lengths downward exactly where density is poor.

**Technical quality was mixed in with biology.** Resolution, experimental method,
per-position DSSP availability, residue numbers, insertion codes and mapping
continuity now sit in a separate block from the biological panel. P17936 sits in a
3.6 Å cryo-EM map, which is *why* its +1 and +2 were never resolved — an analysis
has to be able to see that rather than inferring a short loop.

**Inputs are resolved by configuration.** The repository holds code, tests and
generated tables; the UniProt release and structure cache are large and live
outside it. Stages reached them with literal `../../data/...` paths, which is why
the corrected rerun could not be started from this repository at all.

Making that change immediately introduced a regression, which is the most
instructive part of the exercise. The control sequences load from a relative
`data/cache/` path guarded by `if Path(cache).exists()`. In the new location it
did not exist, the guard skipped it silently, and the manifest came out with no
sequence context for its largest population — meaning the triplet check was
inoperative for 2,296 of 2,660 sites, with nothing in the output to say so. The
same class of failure as the original six, introduced by a change explicitly
meant to prevent that class of failure. Stage 41 now names every cache it resolves
and refuses to write a manifest where any population has no sequence context, and
the sequence checks run per row, so a cache that loads *half* is detectable too.

## The gates

Nothing reaches an analysis without passing these, and each exists because its
absence let something through:

- Shards merge on the full key including `population`; missing shards, duplicate
  keys, recorded failures and short coverage are fatal rather than warnings.
- Every row must have a UniProt sequence, a complete N-X-S/T triplet with X ≠ P,
  and coordinates inside the sequence. A site that is not a sequon has been mapped
  wrongly upstream and cannot be an N-linked site.
- The invariants are asserted, not described, and the report exits non-zero when
  one breaks. A report that only prints cannot stop a bad table being used.

## The three views

One table cannot serve the analysis, because the reasons a site is imperfect are not
interchangeable. A crystallographer's N→Q knockout, a +1 differing between isoform
and construct, and a +2 that was never resolved are three different facts.

- **`triplet_core`** (2,556 sites) — triplet agrees, all three residues located,
  mapping continuous. It requires continuity *as well as* agreement, because
  requiring agreement alone readmits the nine sites whose +1 came from across a
  gap.
- **`asn_core`** (2,624) — the asparagine was measured correctly. Valid for
  features centred on it; not for +1 or +2 exposure, structure or geometry.
- **`construct_review`** (104) — everything excluded, each row carrying its
  reason: 36 substituted asparagines, 30 sequence substitutions, 38 unresolved
  positions.

## Known loose end

The `discontinuous_mapping` exclusion category cannot fire. Any discontinuity
produces an unresolved position, so `unresolved_position` always claims the row
first — verified as zero rows across the whole table. It is dead code. The
distinction it was meant to capture is a property of the *change* rather than of
the new row, and lives correctly in the change audit.

## What this does not settle

The extractor is now trustworthy; the analysis is not yet written. Nothing here says
anything biological. The first real use of this table is stage 2's actual
question: which contextual features distinguish occupied sites from controls, and
are those the features ESM-IF responds to and ProteinMPNN misses.

Two panel items remain underspecified rather than implemented: residues after the
asparagine and after the sequon are joined from the manifest at view-build time
rather than being extractor outputs, which is correct but means they exist in the
views and not in the feature table. And `nearest_disulfide_sg_nd2` has coverage of
51% in occupied sites against 15% in secretory ones. That is real biology rather
than a defect, but the gap between arms is wide enough that it should not enter a
comparison without being stated.
