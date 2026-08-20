# The second model: ESM-IF1

Added 2026-08-20. `esm_if1_gvp4_t16_142M_UR50`, implementing both benchmark
protocols — `SequonScorer` and `SequenceDesigner`.

The point of a second model is that a null result from one model is a statement
about that model's inductive bias until a second, differently-built model agrees.
ProteinMPNN is a message-passing graph network scoring bidirectionally; ESM-IF is
an autoregressive GVP-transformer. On a 29-site pilot their scores correlate at
Spearman 0.13, so ESM-IF is a genuine second opinion rather than a restatement.

## What ESM-IF's conditional is, and why it is not ProteinMPNN's

ProteinMPNN gives

    P(residue at i | backbone, ALL other native residues)

averaged over eight sampled decoding orders. ESM-IF is trained strictly left to
right, so the only conditional it can honestly produce is the prefix one, from a
single teacher-forced pass:

    P(residue at i | backbone, native residues 1..i-1)

Three consequences that belong in any write-up putting the two side by side:

1. **Asymmetric context.** Scoring the +2 residue, both models see the native
   asparagine upstream; only ProteinMPNN also sees the sequence C-terminal to the
   site. The two numbers answer neighbouring questions.
2. **Raw magnitudes are not comparable.** Compare the *matched-pair contrast*
   within each model — which is what the analysis rests on anyway — never the
   score means.
3. **No decoding-order spread.** A teacher-forced pass is deterministic, so
   `conditional_sequon_score_sd` is structurally 0 and `n_decoding_orders` is 1.
   Those columns exist so both models share a schema. They are not evidence that
   ESM-IF is the more precise model.

`conditioning` is recorded as `autoregressive_prefix` rather than `conditional`,
so the two are never silently pooled.

## The index hazard, and why it is closed

The manifest's `model_index` is an ordinal into the chain as
`structures._parse_chains` reads it: Biopython, keeping residues that satisfy
`is_aa(standard=False)` and carry a CA. ESM-IF reads structures through biotite,
whose residue set is not guaranteed to be the same. Trusting the index would
score the wrong residue and still return a plausible number.

Measured on the four matched sets before any mitigation: **~5% of sites**
disagreed. Two causes, two fixes:

- **biotite raises `KeyError`** converting a non-standard residue name (PCA and
  friends) to one letter, aborting a whole chain over one residue. The residue
  *count* is unaffected, because `coords` is built by
  `get_atom_coords_residuewise` independently of the sequence string — so mapping
  unknown names to `X` recovers the chain **without shifting any index**.
- **Where the sequences still differ**, the module's own `_alignment_pairs`
  decides the correspondence.

A site is then scored only if its three residues map to positions whose
identities reproduce the manifest's triplet (`ChainMapping.check_triplet`);
anything else raises. After both fixes, attrition is **zero** across all four
matched sets, so the frozen matching is preserved exactly:

| Comparison | Pairs | Sites | ESM-IF usable |
|---|---|---|---|
| optimal | 16 | 32 | 32 (100%) |
| secretory | 262 | 525 | 525 (100%) |
| bacterial | 280 | 560 | 560 (100%) |
| cytosolic | 273 | 547 | 547 (100%) |

## Generation

`design_sequences` decodes all 32 designs as **one batch**. ESM-IF is
autoregressive, so a chain of length L costs L sequential decoder steps however
many sequences are wanted; sampling one at a time pays that latency 32 times for
no benefit. Measured against ESM-IF's own `model.sample()` on 1A2W/A, 12 designs
at T=0.1:

| | native recovery | speed |
|---|---|---|
| batched (this adapter) | 0.421 ± 0.010 | 2.63 s/design |
| `model.sample()` | 0.415 ± 0.010 | 3.73 s/design |

Residue-composition total variation distance 0.022 — the two agree in
distribution. The batching speedup is only ~1.4x on CPU, which is compute-bound;
the large win is on GPU, where each decoder step is latency-bound.

Sampling is restricted to the twenty standard amino acids, matching
ProteinMPNN's design pass omitting `X`. Without it the decoder can emit a special
token at a sequon position, which `classify_retention` would score as a lost
motif — a fact about the vocabulary, not about the model.

Sequences are returned **in the manifest's index space**, with `X` where ESM-IF
has no counterpart residue, so `design[n_model_index]` reads correctly for either
model. Those `X` positions are unscoreable by construction, so no scored site
ever reads one.

**One deliberate divergence from the earlier ESM-IF work in
`decoding-design-bias`.** That code seeds each design individually
(`torch.manual_seed(seeds[i])`, then one `model.sample`). Batching cannot
preserve per-design seeds: the whole batch is drawn from one seeded generator, so
run *i* of a batch is not reproducible as an individually-seeded sample. The run
as a whole is reproducible from its seed, and retention is a fraction over
designs rather than a claim about any single one, so nothing downstream depends
on per-design identity. It is recorded here because the two repositories would
otherwise look like they disagree.

## Environment

fair-esm 2.0.0 does not import against biotite >= 1.0 without help; both patches
are applied by `esmif_scoring.patch_biotite()` and are idempotent.

- `filter_backbone` was renamed `filter_peptide_backbone`.
- `ProteinSequence.convert_letter_3to1` raises on unknown residue names; it is
  wrapped to return `X`.

Requires `torch-geometric` **and** `torch-scatter` (the GVP encoder imports both)
but **not** `torch-sparse`, which is absent here and unused — worth knowing,
because it is the one that reliably fails to build.

Carried over from `decoding-design-bias/design/score_esmif_cohort.py`: **GVP and
torch_scatter are most reliable on CPU on macOS**, which is why `--device`
defaults to `cpu` and `resolve_device` downgrades rather than failing. Use
`--device cuda` on Colab.

## Running it

```bash
python pipeline/05_scoreability.py <manifest> <out> --model esm_if
python pipeline/07_score.py        <manifest> <out> --model esm_if [--device cuda]
python pipeline/08_design.py       <manifest> <out> --model esm_if [--device cuda]
```

Defaults are unchanged (`--model proteinmpnn`, `--device cpu`), and the
refactored ProteinMPNN path reproduces the previous scores bit-for-bit.

For GPU: build the bundle with `pipeline/30_package_for_colab.py`, then run
`notebooks/esm_if_and_mpnn_gpu.ipynb`.

## Not done

- A protein-level ESM-IF likelihood (`score_sequence`'s `ll_fullseq`), which is
  the score `decoding-design-bias` reports. It answers a chain-level question and
  cannot resolve a three-residue one, but it would make a reasonable **covariate**
  — does a site-level effect survive conditioning on overall chain likelihood?
  Not built, because it was not asked for.
- ESM-IF has no soluble-weights variant to pair with ProteinMPNN's.
