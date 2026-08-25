# Correction, 2026-08-25 — locating the sequon in ProteinMPNN's parse

ProteinMPNN was read at the wrong residue for about a quarter of sites. Every
ProteinMPNN score and retention figure produced before this date needs
regenerating. **ESM-IF and ESMC are unaffected**, and their numbers are
unchanged.

The headline moves. ProteinMPNN's conditional-score result goes from **+0.090 SD
(BH 0.30, inconclusive)** to **+0.282 SD (BH 0.021)**. "ProteinMPNN does not
distinguish occupied sequons" does not survive.

## What was wrong

Two parsers index a chain differently, and the manifest's index belongs to one
of them.

`structures._parse_chains` lists the residues actually observed, in order;
`model_index` is an ordinal into that list. ProteinMPNN's `parse_PDB_biounits`
instead walks the author numbering and emits a slot per number, inserting a
placeholder wherever a number is absent:

```python
for resn in range(min_resn,max_resn+1):
    if resn in seq:
        for k in sorted(seq[resn]): seq_.append(aa_3_N.get(seq[resn][k],20))
    else: seq_.append(20)
```

The two coincide only for chains numbered without gaps. Most depositions have
gaps, so the index drifts by the number of absent residue numbers before the
site. In 9G3Q chain A — 402 observed residues spanning 422 numbers — the drift is
20, and `model_index` 181 reads `LKN` where the site is `NES`.

**ESM-IF was never affected**, because it already reconciles the two parses:

```python
if native.sequence == esm_seq:
    to_esm = {i: i for i in range(len(esm_seq))}
else:
    to_esm = {a - 1: b - 1 for a, b in _alignment_pairs(native.sequence, esm_seq)}
```

and then refuses anything whose residues do not reproduce the manifest's
triplet. ESMC checks the triplet too. The ProteinMPNN adapter accepted the same
`expected_triplet` argument and never used it.

The scoring path and the design path failed separately. `07_score` reads
probabilities at `model_index`; `08_design` reads designed sequences at the same
indices, and the ESM-IF adapter's `design()` docstring states the contract that
makes that safe — *"Unconstrained designs, returned in the manifest's index
space"* — which ProteinMPNN's did not honour.

## How it was found

Not by inspection. While building an unrelated experiment that holds a sequon
fixed during redesign, the fixed positions came back as residues that were not
the sequon. The check that exposed it is the cheapest one available: hold three
positions fixed, then read them back.

```
5H5Y:A 257   manifest 'NRS'   design at model_index 'VNI'
9G3Q:A 225   manifest 'NES'   design at model_index 'LKN'
4EBY:A  52   manifest 'NSS'   design at model_index 'NSS'
```

## What it affected

Across the 2,640 sites of the scoring manifest:

| | Sites | Share |
|---|---:|---|
| Index already correct | 1,740 | 65.9% |
| Read at the wrong residue | 668 | 25.3% |
| Not reconcilable, now dropped | 232 | 8.8% |

The 232 are chains where ProteinMPNN's parse is *shorter* than the observed
residue list, which gap-filling cannot produce; the two parsers disagree about
which residues exist, and no arithmetic settles it. They are refused rather than
guessed.

Rescoring confirms the attribution is exact. Comparing corrected against frozen
scores, differences appear in **100% of the 69 dataset sites the mapping moved**
and **0% of the 214 it did not** (median |difference| 2.25 against 2.7e-6).

## What changed in the results

Secretory comparison, cluster-level sign-flip permutation, BH across eight tests:

| Outcome | Frozen | Corrected |
|---|---|---|
| Conditional score | +0.090 SD, BH 0.30 | **+0.282 SD, BH 0.021** (Holm 0.047) |
| Design retention | +0.0423, BH 0.30 | +0.0700, BH 0.225 |

The score result clears correction; retention moves the same way but does not.
Its permutation p (0.112) and Wilcoxon p (0.003) differ by more than an order of
magnitude, which is the ortholog clustering doing real work — the permutation
test is the pre-specified one.

Under the pre-specified equivalence framework the score is **directional,
magnitude undetermined**: the interval [+0.114, +0.426] SD excludes zero and
straddles the ±0.2 SD margin.

Unchanged: the internal-control comparison is null on both outcomes (n=15), the
bacterial diagnostic fires hard negative (−0.528 SD, BH 0.0004), and the
cytosolic diagnostic is flat (−0.002 SD). The compositional-confound reading of
the diagnostics stands.

## What did not change, and how that is known

ESM-IF was rerun through the same code on the same machine as a control. Its
scores differ from the frozen ones by a median of 3.2e-05 against a
between-site standard deviation of 1.513 — unbiased, 47.7% positive — which is
thread-count-dependent floating-point reduction, not indexing. Rerun through the
analysis it reproduces **+0.431 SD exactly**.

That control is what makes the ProteinMPNN change attributable. Both models were
rerun; only the one with the defect moved.

## The claim, restated

Previously:

> ESM-IF distinguishes occupied sequons from matched controls on both outcomes.
> ProteinMPNN does not, on either.

Now:

> Both models distinguish occupied sequons on the conditional score, ESM-IF
> roughly 1.5x more strongly. On design retention only ESM-IF survives
> correction.

The architecture contrast that motivated much of the downstream framing is
substantially smaller than it appeared, and "current models are
glycosylation-blind" is weaker still as a description of the baseline.

## What now guards it

- `mpnn_scoring.build_index_map` reconstructs ProteinMPNN's enumeration and
  returns nothing unless every mapped position matches that parser's own
  sequence residue for residue.
- `ProteinMPNNAdapter.score_from` uses the `expected_triplet` it always
  accepted, refusing residues that do not reproduce the manifest's.
- `ProteinMPNNAdapter.design` projects designs back into the manifest's index
  space, so the contract every caller assumes is met once rather than
  remembered separately.
- Method and worked examples: [`methods_sequon_indexing.md`](methods_sequon_indexing.md).
