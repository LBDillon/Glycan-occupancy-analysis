# Findings — occupied against motif-only sequons in the same protein

*2026-08-31. A rebuild of the old `glycosylation-bias-analysis` experiment on
`gca`'s evidence standard and corrected ProteinMPNN code path. The analysis
design is the old one; the data, the site labels and the model path are not.*

## The question

> When a protein is redesigned with nothing held fixed, is an experimentally
> occupied sequon treated differently from a sequon in the **same chain** that
> carries no occupancy evidence?

The matched-secretory comparison pairs an occupied site in one protein against a
control in another, so protein identity is never fully controlled. This design
is within-protein by construction: both classes sit on one chain, share one
fold, one depositor, and — because designs are generated per chain — **the same
32 designs**.

## The answer

**On the readings that matter, no.** Occupied sequons are not retained
differently from motif-only sequons in the same protein. The old experiment's
null survives the stricter evidence standard rather than being overturned by it.

134 chains carry both classes. Intervals resample proteins; the chain is the
unit of pairing, so each chain contributes one difference however many sequons
it carries.

| Reading | occupied | motif-only | difference | 95% CI |
|---|---|---|---|---|
| `n_retained` | 0.2498 | 0.2262 | +0.0236 | [−0.0436, +0.0905] |
| `pattern_retained` | 0.1178 | 0.0985 | +0.0193 | [−0.0275, +0.0678] |
| `exact_retained` | 0.0782 | 0.0384 | **+0.0398** | **[+0.0064, +0.0763]** |

`pattern_retained` is the biologically meaningful reading — asparagine, not
proline, then serine or threonine — and it is **null**. So is bare asparagine
retention. Only the strictest reading, all three residues identical, separates
the classes, and that is discussed as a lead below rather than as the result.

## N retention by class

Site-weighted means, protein-resampled intervals, matching
`10_analyse_retention_by_class.py`.

| Class | sites | proteins | Asn retained | 95% CI |
|---|---|---|---|---|
| occupied sequon | 290 | 211 | 0.2481 | [0.2057, 0.2897] |
| motif-only sequon | 410 | 138 | 0.2799 | [0.2387, 0.3183] |
| non-sequon asparagine | — | 220 | 0.3738 | [0.3487, 0.3996] |

Two things worth separating here.

**An asparagine inside a sequon is retained less often than one outside it**
(24.8% and 28.0% against 37.4%). That is not a glycosylation effect: sequon
asparagines sit in exposed loops and turns, which is where redesign changes most.

**The unpaired comparison points the opposite way from the paired one** — motif-only
looks *higher* (28.0% against 24.8%) while the within-protein contrast is
slightly positive (+0.024). This is exactly the confound the paired design
removes. The two classes are not drawn from the same proteins: 290 occupied
sites spread over 211 proteins, 410 motif-only sites over 138. Comparing the
pooled means compares two different sets of chains.

The old experiment reported 25.0% against 25.6% and read it as a null. The
occupied arm reproduces closely (24.8%), and its null conclusion stands — but it
reached it through the unpaired comparison, which cannot support it either way.

## Is a sequon lost faster than any other triplet?

No, and this reproduces the existing result rather than revising it.

| Quantity | Mean | 95% CI |
|---|---|---|
| Occupied sequon retained, exact triplet | 7.35% | — |
| Control triplet retained, exact | 9.49% | — |
| **control − sequon** | **+2.14 pp** | **[−0.64, +4.73]** |

Exact against exact on both sides. A control triplet has no pattern latitude, so
comparing it against `pattern_retained` would overstate the case; that comparison
appears nowhere here. Against `56_sequon_retention_rate.py`'s +1.8 pp
[−1.0, +4.6] this is the same answer.

## The exact-retention lead, and why it is not the headline

`exact_retained` is the one reading whose interval excludes zero, and three
things bound it.

**It is one of three readings with no multiplicity correction.** At a 95%
interval, one of three nominally clearing is unremarkable on its own.

**The biologically meaningful reading is null.** A sequon is functional if it
matches N-X-S/T; it does not have to be the *same* N-X-S/T. That reading gives
+0.019 with an interval covering zero.

**It is not a composition artefact, which makes it worth keeping.** The obvious
explanation would be a different NXS/NXT mix between the classes, and within the
134 paired chains there is not one — occupied are 51.1% NXS, motif-only 52.9%.
The effect sits entirely in NXT: occupied NXT retain exactly 10.97% against
motif-only NXT at 6.39%, while the NXS arms are indistinguishable (6.68% against
6.54%).

So: a real pattern, concentrated in one sequon subtype, on the strictest reading
only. That is a lead for a pre-registered test, not a finding to report.

## Validation

The rebuild reproduces the existing ProteinMPNN result closely, which is the
cheapest evidence it is wired correctly.

| Check | This run | Reference |
|---|---|---|
| Pattern retention, `56_`-equivalent site set | 13.49% | 13.58% (`56_`), 13.0% (ARC) |
| Pattern retention, occupied only | 13.18% | 13.06% (`retention_by_class`) |
| Control triplet, exact | 9.49% | 9.49% (`56_`) |
| Background mutation rate | 56.13% | 55.7% (`56_`) |
| Site-by-site correlation with `56_` | **r = 0.989**, 241/290 identical | r = 0.993 (`56_` vs ARC) |

Two comparisons need their populations stated or they look like discrepancies.

**`56_` did not filter on occupancy status.** It verifies every scoreable
manifest site, so its 13.6% covers occupied sites *and* the 28 internal
controls — 318 sites, not 314 occupied ones. The like-for-like check is the
`56_`-equivalent row above; the occupied-only figure is a different population.

**`fig_model_comparison_values.json`'s 0.121 is the matched-secretory subset**,
216 pairs rather than the dataset-wide set. The right dataset-wide reference is
`retention_by_class`'s 0.1306, which this run reproduces at 0.1318.

## Two corrections, and what they cost

**Sequon indexing (2026-08-25).** `verify_sequon_index` passes for **290/290**
occupied sites analysed, with zero refusals. That is not the guard being
vacuous: the manifest's `n_model_index` is already post-correction, so the
refusals were spent when the manifest was built. The guard bites at design time
instead, where **20 chains (24 sites) were rejected** as unmappable onto
ProteinMPNN's parse and are recorded in
`retention_dataset_unconstrained_proteinmpnn_failures.csv`. No new index mapping
was written; `chain_mapping`, `to_manifest_space` and `verify_sequon_index` are
`gca`'s.

**Token alphabet (2026-08-20).** Not reached. Nothing here indexes into a
probability vector — every quantity is read off generated sequences. The old
Figure 1 panel C, which reports MPNN confidence scores, is the part that would
have been exposed, and it is deferred.

## Attrition

| Gate | Chains | Sites |
|---|---|---|
| Scoreable manifest sites | 248 | 342 |
| Structure cached | 248 | 342 |
| Chain mappable onto ProteinMPNN's parse | 228 | 318 |
| `verify_sequon_index` passes | 228 | 318 |
| Design length matches native | 228 | 318 |
| — of which `occupied_supported` | — | 290 |
| — of which `observed_unmodified` | — | 28 |
| Motif-only sequons found in the same chains | 228 | 410 |
| **Chains carrying both classes** | **134** | — |

One protein count and one site count throughout: **228 chains, 220 accessions,
290 occupied sites, 410 motif-only sites, 134 paired chains.**

## What this does not cover

**Three of the four models did not run.** Only ProteinMPNN is runnable on this
machine. ESM3 works but needs 8+ hours on CPU. ESM-IF and CARBonAra have no
environment here at all — ESM-IF hits the `fair-esm` / EvolutionaryScale `esm`
name collision, and no CARBonAra checkout exists. This is environment, not code:
`57_cross_protein_sequon_retention.py` takes `--sequences`, so each model slots
in with one command once its designs exist. See
[`running_locally.md`](../../docs/running_locally.md).

**Deferred by design**, and free once sequences are on disk: secondary-structure
and Ramachandran stratification (needs DSSP), the MPNN confidence panel, and
matching occupied to motif-only on RSA and packing. The `asn_only` and
`full_sequon` design conditions need new runs and are not covered.

**The motif-only class includes the 28 `observed_unmodified` sites**, following
the brief's definition of "not `occupied_supported`". Those have evidence of
*absence* rather than absence of evidence, which is a different thing; they are
flagged `manifest_unsupported` in the site table so the stricter reading is
available without a re-run.

## Reproducing

```bash
export GCA_STRUCTURE_DIRS=/path/to/pdb/cache:/path/to/ortholog/structures
export PROTEINMPNN_DIR=/path/to/ProteinMPNN

python pipeline/08_design.py \
  results/manifests/candidate_manifest_dataset.csv \
  results/designs/retention_dataset_unconstrained_proteinmpnn.csv \
  --model proteinmpnn --save-sequences

python glyco_context/pipeline/57_cross_protein_sequon_retention.py \
  --sequences results/designs/retention_dataset_unconstrained_proteinmpnn_sequences.csv \
  --label cross_protein_proteinmpnn
```

ProteinMPNN `v_48_020`, temperature 0.1, 32 designs per chain, seed 0. Intervals
are percentile bootstraps over proteins, `N_BOOT=4000`, `BOOT_SEED=11`. Use a
fresh `--out` path: resume keys on manifest sites, so resuming onto an existing
table skips chains whose sites are done and would leave their sequences
unwritten. The stage warns when this applies.
