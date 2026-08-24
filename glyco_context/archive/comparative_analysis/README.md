# Archived — the comparative context analysis

*Archived 2026-08-24, the day it was completed. Not deleted: it is a real result
and the reasoning is worth keeping. Not built on either.*

## What it asked

Whether the local environments of experimentally occupied N-X-S/T sequons differ
from those of unoccupied ones — against internal controls, and against a large
secretory-unannotated set.

## What it found

Once composition is controlled by matching, nearly every apparent difference
disappears. In the 261-pair secretory comparison only distance to the N-terminus
survives correction; ND2 crowding falls from −0.59 to −0.08, nearest aromatic
from +0.39 to +0.10. One feature holds its direction under every framing without
reaching significance in the matched ones: β-sheet at the +2 position.

Full account in `docs/findings_2026-08-24_context_differences.md`, including a
recorded error in my own first interpretation of it.

## Why it is archived rather than continued

The question needs a negative label that does not exist, and a comparison group
drawn from other proteins. Occupied and secretory-unannotated sites share **zero
proteins and zero chains**, so occupancy is confounded with protein identity;
the alternative — comparing within a protein — has 31 sites. Neither route is
worth more machinery.

The successor question avoids this entirely by being paired within a site: the
same protein, the same backbone, the same position, with only the sequence
changing. See `glyco_context/docs/` for the fixed-sequon context-retention test.

## What was kept from it

The extractor and its correctness work, which any context analysis needs, and
the descriptive reference distribution of natural occupied sites — which is
exactly the reference the successor experiment scores against.
