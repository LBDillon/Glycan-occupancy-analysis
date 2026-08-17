# Phase 1 audit — control data verification

Gate for model scoring. Nothing in Phase 2 onward may proceed on data that fails
these checks.

## Verdict

**Pass, after two defects were found and repaired.** Both were real and both
would have invalidated the comparisons had scoring started on the previous
state.

## Defects found

### 1. Provenance scrambled — 43.5% of rows mislabelled

`build_features` re-sorts its output by `(accession, position)`. The caller
attached `control_set` afterwards positionally, pairing labels from input order
against rows in sorted order.

| | bacterial (source) | cytosolic (source) |
|---|---|---|
| **labelled bacterial** | 3,258 | 2,608 |
| **labelled cytosolic** | 2,607 | 3,505 |

5,215 of 11,978 rows carried the wrong set, so the bacterial and cytosolic
matched comparisons were thoroughly mixed. `occupancy_status`, which travelled
with the row rather than being reattached, was correct — which is what localised
the fault.

**Repair.** `build_features` now takes `carry_columns` and copies provenance per
row, so no caller can reattach positionally. The joining script uses a key-based
merge with `validate="one_to_one"`. A regression test supplies input in
reverse-sorted order, which fails against the old behaviour.

### 2. The cytosolic set was not eukaryotic

`go:0005829` carries no taxonomic restriction, and none was applied. The set
contained at least **1,929 bacterial** and **302 archaeal** proteins; *E. coli*
was its third most common organism.

The archaeal contamination is the serious part. Archaea N-glycosylate via AglB,
so those were potentially genuine glycoproteins inside a set whose definition is
that its members cannot be glycosylated. It also explains the 13 accessions that
appeared in *both* control sets: bacterial proteins annotated cytosolic and
periplasmic satisfied both queries.

**Repair.** `taxonomy_id:2759` added to the cytosolic query; inventory and all
downstream files rebuilt.

## Assertions

| Check | Result |
|---|---|
| `control_set` agrees with source inventory | **0** disagreeing rows |
| `occupancy_status` equals `control_<control_set>` | **0** disagreeing rows |
| Each matched comparison contains only its intended controls | **0** violations, across 3 comparisons × 5 seeds |
| No control protein in both biological control sets | **0** accessions |
| Cytosolic set is eukaryotic | **0** bacterial or archaeal organisms |

## Final counts

**Control inventories**

| Set | Proteins | Sequons | With structural features |
|---|---|---|---|
| Cytosolic eukaryotic | 4,764 | 16,386 | 3,543 sequons / 1,719 proteins |
| Bacterial extracytoplasmic | 1,141 | 5,865 | 3,280 sequons / 1,033 proteins |

Cytosolic fell from 19,337 to 16,386 sequons when the non-eukaryotic proteins
were removed. Bacterial was unaffected, as expected.

**Dataset sites with structural features**

| Class | Sites | Proteins |
|---|---|---|
| `occupied_supported` | 332 | 245 |
| `observed_unmodified` | 32 | — |

**Matched sets** (seed 0; 1:1 for the scarce comparison, 1:5 otherwise)

| Comparison | Controls available | Occupied sites matched | Pairs |
|---|---|---|---|
| vs observed-unmodified | 32 | 24 | 24 |
| vs bacterial control | 3,280 | 280 | 1,180 |
| vs cytosolic control | 3,543 | 275 | 1,196 |

### Reading the observed-unmodified comparison correctly

It contains **24 occupied sites matched to 24 observed-unmodified controls**.
Those 24 are drawn from the **332 structurally scoreable occupied sites** — they
are not "24 of 32". The 32 is the size of the available control pool, of which 24
found a partner inside the caliper. The statistical unit is the occupied site,
and there are 24 of them in this comparison.

## Weighted balance

Each occupied site carries total weight one; the controls matched to it share
weight one between them. Without this, a case matched to five controls would
contribute five times as much to the control mean as a case matched to one, so
the statistic would describe the matching's bookkeeping rather than the
comparison. Verified: total control weight equals the number of matched cases.

Standardised mean differences, seed 0, with the range across seeds 0–4:

| Comparison | RSA | Neighbours | Hydrophobic fraction |
|---|---|---|---|
| vs observed-unmodified | +0.017 (−0.049..+0.017) | 0.000 | 0.000 |
| vs bacterial | +0.011 (−0.005..+0.021) | 0.000 | −0.000 |
| vs cytosolic | +0.020 (+0.020..+0.030) | 0.000 | −0.000 |

Worst absolute value across all comparisons, features and seeds: **0.030**,
comfortably inside the conventional 0.1 threshold. Seed-to-seed spread is at most
0.066 and mostly far smaller, so the matching is not seed-sensitive.

**The exact zeros are real, not a bug.** `n_neighbours_8a` is integer-valued and
`hydrophobic_fraction_8a` is a ratio over 5–6 neighbours, so both are effectively
discrete and the matcher finds exact partners for 100% of pairs. Confirmed at
pair level; RSA, being continuous, matches to within about 0.01. It was checked
precisely because a reported SMD of exactly zero usually indicates an error.

## Known limitations carried forward

- Structures above 20 MB are skipped at parsing, mapping and linkage detection,
  so ribosomal and proteasomal control proteins largely lack features. This is
  reported as `structure_too_large`, not silently dropped, but it biases the
  control subset against large-complex components.
- Fewer than half of control sequons have features (6,823 of 15,107), mostly
  through unmapped residues and the size guard.
- The 32 observed-unmodified sites exist only where someone solved a glycoprotein
  structure, so they represent well-studied secreted and membrane proteins rather
  than a random draw from the unlabelled pool.
- 24 pairs cannot support a null result on their own. Any equivalence claim needs
  a smallest effect size of interest fixed before the labelled contrasts are seen.
