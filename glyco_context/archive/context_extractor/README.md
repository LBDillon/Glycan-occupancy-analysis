# Archived — the structural context extractor

*Archived 2026-08-26. The code is archived; its output is retained and still
used.*

## What it was for

Stages 41 through 47 built a 97-column description of every sequon's structural
environment: solvent accessibility, secondary structure, backbone geometry,
loop runs, neighbour composition, the shell around ND2, and the QC fields that
make those interpretable. It was built for the comparative analysis — occupied
sites against unoccupied ones — which is archived separately in
`../comparative_analysis/`.

## Why it is archived

**The live experiment uses none of it.** The fixed-sequon context-retention test
computes its own fifteen composition features in `local_chemistry.py`, from the
sequence plus a single structure pass. Of the 97 columns this chain produces,
the live analysis reads zero.

That is not a criticism of the extractor, which was correct and heavily tested.
It answered a question that turned out to be confounded. And the reason it
cannot serve the successor is structural: ProteinMPNN is fixed-backbone, so
accessibility, secondary structure and geometry **cannot change** between a wild
type and its design. The features this chain is best at are exactly the ones the
successor experiment cannot move.

## What is retained, and still used

`glyco_context/results/datasets/context_triplet_core.csv` — the site list and its
features. Stage 51 reads it for identifiers; stage 50 reads it for the
descriptive reference figure. Both remain live. The table is treated as a
retained input, in the way an expensive upstream step usually is.

Regenerating it means restoring these stages, and the modules in `src/` alongside
them, to the paths their imports expect.

## What was learned here and kept

The mapping correction of 2026-08-24 — insertion codes, gap-aware +1/+2 walking,
dihedrals not crossing gaps, DSSP on multi-character chains, terminal distances
as residue counts — lives in `docs/correction_2026-08-24_context_mapping.md`
here. The equivalent defect on the model-facing side was found because of it,
and that fix is live: see `docs/methods_sequon_indexing.md` in the repository
root.
