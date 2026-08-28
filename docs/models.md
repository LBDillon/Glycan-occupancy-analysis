# The models

How the benchmark grew from one model to six, what each one is, and how to
run it. Merges the former `adding_models_explainer.md`, `second_model_esm_if.md`
and `third_model_esmc.md`.

*Written 2026-08-20, covering the day ProteinMPNN stopped being the only model.*

## The shape of the problem

The benchmark asks one question: does a model treat an experimentally occupied
N-X-S/T sequon differently from a structurally matched sequon carrying no glycan?
Answering it for a single model tells you about that model. Answering it for one
model and calling it a result about *models* is a mistake, because a null can
just as easily be a property of one architecture's inductive bias as a fact about
what protein models know.

So the aim was never "add ESM-IF". It was to make the benchmark **model-shaped**:
to get to a state where the question is asked of an interchangeable thing, and
adding the next one is a day's work rather than a rewrite.

## The one decision everything else follows from

Downstream of scoring, nothing is allowed to know which model produced a number.

Matching, contrasts, the cluster bootstrap, significance testing and the figures
all operate on tables keyed by `(accession, position)`. A model contributes two
columns — a conditional sequon score and a retention fraction — and contributes
nothing else. That constraint is what makes a second model an adapter instead of
a fork.

It also has a consequence worth stating plainly: **matching is model-independent
and is never recomputed.** Pairs are built from RSA, neighbour counts and
hydrophobic fraction, never from model output. Every model is therefore scored on
exactly the same pairs, and a disagreement between models is a disagreement about
the same sixteen, or two hundred and sixty-two, comparisons. If matching moved
per model, no cross-model comparison would mean anything.

## The interface, and the two invariants it enforces

`adapters/base.py` declares two protocols. A model implements either or both.

| Protocol | Question | Feeds |
|---|---|---|
| `SequonScorer` | what probability does the model hold at the three sequon residues? | `07_score.py` |
| `SequenceDesigner` | what does it write when redesigning the chain? | `08_design.py` |

Two rules are non-negotiable, and both were learned by being burned.

**Never score a residue the model did not evaluate.** ProteinMPNN returns a row
of zeros for residues with incomplete backbones, which exponentiates to
twenty-one ones and scores about +13.8. That defect inverted the sign of the
first result this project produced. Adapters raise rather than return.

**`decodable_positions` must not need a model pass.** Scoreability has to be
answerable before matching. Settle it afterwards and matched sets quietly lose
members, unbalancing exactly what matching had just balanced.

One addition was made this week. `SequonScorer` gained
`prepare_chain` / `score_from` alongside `score_site`, because the models waste
effort in opposite directions: ProteinMPNN runs a decoder pass *per position*, so
it wants every position on a chain at once; ESM-IF decodes the whole chain in one
pass, so its second sequon should cost nothing. Splitting the once-per-chain work
from the per-sequon read serves both without either model's quirk leaking into a
pipeline stage.

## What each model turned out to be

### ProteinMPNN — the incumbent

Backbone plus every other native residue, averaged over eight sampled decoding
orders. Bidirectional and symmetric.

### ESM-IF1 — a second structure-conditioned opinion

The closest conceptual comparison: also reads a backbone, also generates
sequence. The point is to ask whether ProteinMPNN's null is specific to
ProteinMPNN or general to inverse folding.

Its scoring is genuinely not equivalent, and pretending otherwise would have been
the easy mistake. ESM-IF is autoregressive, so the only conditional it can
honestly produce is prefix-only:

    P(residue at i | backbone, native residues 1..i-1)

There is no decoding-order distribution to average, because permuting a causal
decoder's order is off-distribution. So `conditional_sequon_score_sd` is
structurally zero and `n_decoding_orders` is one, and the `conditioning` column
records `autoregressive_prefix` rather than `conditional` so the two can never be
silently pooled. Scoring the +2 residue, both models see the native asparagine
upstream; only ProteinMPNN also sees C-terminal context. Raw magnitudes are not
comparable between models. The SD-standardised matched-pair contrast is, and that
is what the analysis rests on.

The alternative — falling back to a whole-sequence likelihood — was rejected. A
per-residue mean over several hundred positions cannot resolve a three-residue
question, which is the objection the module already raises against protein-level
scores.

#### What ESM-IF's conditional is, and why it is not ProteinMPNN's

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

#### The index hazard, and why it is closed

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

#### Generation

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

#### Environment

fair-esm 2.0.0 does not import against biotite >= 1.0 without help; both patches
are applied by `esmif_scoring.patch_biotite()` and are idempotent.

- `filter_backbone` was renamed `filter_peptide_backbone`.
- `ProteinSequence.convert_letter_3to1` raises on unknown residue names; it is
  wrapped to return `X`.

Requires `torch-geometric`. It does **not** require `torch-sparse` (absent here
and unused), and `torch-scatter` is optional.

`torch-scatter` is a compiled extension needing a wheel matched to the exact
torch build. PyG stops publishing them for older torch, and building from source
needs `nvcc` — so on a current torch it is simply unavailable, which is what
blocked the first ARC setup. ESM-IF imports two names from it, calls one
(`scatter_add`, to count edges into each node) and never uses the other, so
`_torch_scatter_shim.py` supplies both in native torch and registers itself as
`torch_scatter` when the real package is missing. Verified `torch.equal` against
the real package, and ESM-IF scores through the shim are bit-identical on the
same device.

**Cross-device reproducibility, separately:** the same sites scored on CPU and on
a GPU differ by up to ~2e-3 in the sequon score. That is ordinary float32
arithmetic, not a defect, but it means a score table should record which device
produced it if runs are ever to be compared at that precision.

Carried over from `decoding-design-bias/design/score_esmif_cohort.py`: **GVP and
torch_scatter are most reliable on CPU on macOS**, which is why `--device`
defaults to `cpu` and `resolve_device` downgrades rather than failing. Use
`--device cuda` on Colab.

#### Running it

```bash
python pipeline/05_scoreability.py <manifest> <out> --model esm_if
python pipeline/07_score.py        <manifest> <out> --model esm_if [--device cuda]
python pipeline/08_design.py       <manifest> <out> --model esm_if [--device cuda]
```

Defaults are unchanged (`--model proteinmpnn`, `--device cpu`), and the
refactored ProteinMPNN path reproduces the previous scores bit-for-bit.

For GPU: build the bundle with `pipeline/30_package_for_colab.py`, then run
`notebooks/esm_if_and_mpnn_gpu.ipynb`.

#### Not done

- A protein-level ESM-IF likelihood (`score_sequence`'s `ll_fullseq`), which is
  the score `decoding-design-bias` reports. It answers a chain-level question and
  cannot resolve a three-residue one, but it would make a reasonable **covariate**
  — does a site-level effect survive conditioning on overall chain likelihood?
  Not built, because it was not asked for.
- ESM-IF has no soluble-weights variant to pair with ProteinMPNN's.


### ESMC 300M — the sequence-only control

The first model here that sees no structure at all. It exists to ask whether
sequence context alone reproduces whatever the structure-conditioned models
report. If it does, the effect need not be structural. If it does not, they are
doing work sequence cannot.

It is scored on the chain sequence as the module's own parser reads it, rather
than the full UniProt sequence, so indices, scoreability and matched pairs stay
identical and "sequence alone versus sequence plus structure" is like-for-like
rather than also changing the context window.

It implements `SequonScorer` only. Sampling from a masked language model would
condition on sequence rather than a backbone, so "retention" would not mean what
it means for the other two.


#### Environment: ESMC and ESM-IF cannot coexist

`fair-esm` (ESM-IF1) and EvolutionaryScale's `esm` (ESMC, ESM3) both install a
top-level package named **`esm`**. Installing one shadows the other. This is not
a version conflict that can be pinned around — it is the same import name.

So each needs its own environment. The registry lazy-imports, so the absent model
is simply unavailable rather than breaking the package, and `tests/` skips its
model-dependent test when the SDK is missing.

```bash
python -m venv --system-site-packages esmcenv
esmcenv/bin/pip install --no-deps esm==3.2.2
esmcenv/bin/pip install --no-deps huggingface_hub'<1.0' 'tokenizers>=0.21,<0.22' \
    'transformers<4.48.2' tenacity httpx zstd msgpack-numpy cloudpathlib \
    brotli attrs einops regex safetensors
```

`--no-deps` throughout is deliberate: the SDK declares `torchtext`, which is
dead against modern torch and is not imported by ESMC.

#### What is scored, and on which sequence

The chain sequence as `structures._parse_chains` reads it — the same string the
manifest's `model_index` is an ordinal into. Chosen over the full UniProt
sequence so that model indices, scoreability and the matched pairs stay identical
to the structure-conditioned models: "sequence alone versus sequence plus
structure" is then a like-for-like comparison rather than one that also changes
the context window.

The cost is that ESMC sees only the resolved chain. Unresolved loops and
truncated termini are absent — exactly as they are for the other two models. The
full-UniProt variant remains a reasonable sensitivity analysis; it is not built.

Because a sequence model has no backbone requirement, `decodable_positions`
returns all True. ESMC's scoreable set is a superset of the others', which never
widens a matched set: the pairs were frozen on ProteinMPNN's scoreability.

#### The token offset, verified rather than assumed

The tokenizer prepends `<cls>`, so sequence position *i* is token index *i + 1*.
`_assert_token_offset` round-trips a probe sequence through the tokenizer at load
and raises if it does not reproduce it.

This is the check the ProteinMPNN alphabet defect went undetected for. An
assumption about how a model indexes or encodes its own input is not a fact until
something reproduces the input from it.

#### Masking: two estimands, not two estimates

| Mode | Reads | Use |
|---|---|---|
| `single` (default) | `P(residue at i \| every other native residue)` | primary — closest sequence-only analogue of ProteinMPNN's conditional |
| `joint` | all three sequon positions masked together | sensitivity |

`joint` exists because `single` leaves a confound. Masking only the +2 residue
still shows the model a native asparagine two positions upstream, and N-X-S/T is
a heavily learned motif — so the model can infer S/T from the N rather than from
anything about this site.

Measured on 13 dataset sites:

| | single | joint | delta |
|---|---|---|---|
| P(N) | 0.2244 | 0.1767 | −0.0477 |
| P(S)+P(T) | 0.6022 | 0.5404 | −0.0618 |
| score | −0.3795 | −0.7195 | **−0.3400** |

Correlation between modes r = 0.837. So roughly **0.34 log-odds of the score is
the motif reinforcing itself** rather than site-specific information.

This does not invalidate `single`. Occupied and control sites both carry a
sequon, so the inflation applies to both arms and largely cancels in the paired
contrast — but it compresses dynamic range and could hide a real difference,
which is why the sensitivity is worth running.

**The same confound applies to ProteinMPNN and ESM-IF**, which also condition on
the rest of the native sequon. **Both now have a joint-masking variant.**

ProteinMPNN's was straightforward — `conditional_probs` already takes a position
set — and was run on 26/08. Its preference *survives* masking, unlike ESMC's:
see the masking section of [`OVERVIEW.md`](OVERVIEW.md).

ESM-IF's is harder and is done differently. A causal decoder cannot hide an
upstream residue without hiding it from everything downstream, and substituting
a `<mask>` token is off-distribution — it sends most of the probability mass
onto aromatics. So `marginalised_probabilities` marginalises over the hidden
residues instead, sampling them from the model's own belief rather than
substituting a token that means nothing to it. `MASK_MODES = ("single",
"joint")` and `DEFAULT_MARGINAL_SAMPLES = 16`, a count fixed by measurement
rather than choice (the check is archived in `archive/settled/`).

**It has never been run.** The capability exists and no ESM-IF joint scores are
on disk. It is the more informative comparison than ESMC, because ESM-IF and
ProteinMPNN are both structure-conditioned and so are blindable to the same
degree — where ESMC, seeing only sequence, loses more when its motif is hidden
than a structure-conditioned model can.

#### Why there is no SequenceDesigner

ESMC is a masked language model. Sampling from it would condition on sequence
rather than on a backbone, so "retention" would not mean what it means for the
inverse-folding models. Scoring only.

#### Running it

```bash
esmcenv/bin/python pipeline/07_score.py <manifest> <out> --model esmc
```

Cost on CPU: ~1.6 s/site under `single` (three masked variants per sequon,
batched per chain), ~0.6 s/site under `joint` (one). Cheaper than either
inverse-folding model.

#### Not done

- ESM3-open, which would give sequence-only / structure-only / sequence+structure
  from one set of weights. Gated behind a licence acceptance and a token, ~1.4B
  parameters. The function and annotation tracks must be left empty, or they can
  leak the glycosylation label.
- ESM-2, which occupies the same position in the benchmark as ESMC.
- The full-UniProt-sequence sensitivity described above.

### CARBonAra — a fourth parser, scorer and designer

*Added 2026-08-27.* CARBonAra (Krapp et al., Nat Commun 2024) is a geometric
transformer over atomic coordinates and element names, built on PeSTo. It has no
amino-acid-specific parameterisation, which is what lets it read arbitrary
molecular context — and is exactly the capability this integration switches off.

**The input is protein only.** One chain, MSE converted to MET, every glycan,
ligand, ion, water, hydrogen and other chain removed, heteroatom and water
removal set on the loader as well for good measure. That is the same view
ProteinMPNN gets. Scoring CARBonAra on its default context-aware input would not
be a fourth opinion on the same question; it would be a different question, and a
partly circular one, since a NAG modelled onto ND2 is covalently bonded to the
asparagine whose probability is being read.

**The conditional is `conditional_all_other_native`, `n_orders = 1`.**
CARBonAra is not autoregressive and has no decoding order. Residue identity
enters as an "imprint": a per-residue one-hot that is zero wherever identity is
withheld. `apply_model` sees element types plus that imprint and nothing else —
side chains are stripped and an idealised C-beta placed at every residue,
glycine included — so a withheld position genuinely leaks no identity, from the
sequence or from a side chain. Scoring position *i* means imprinting every other
mappable canonical residue natively and zeroing row *i*:

    P(residue at i | backbone, all other native residues)

That is nearer ProteinMPNN's bidirectional conditional than ESM-IF's causal one,
but it is not the same: there is no decoding-order distribution to average, so
`conditional_sequon_score_sd` is structurally 0 and `n_decoding_orders` is 1, as
they are for ESM-IF. One forward pass per position, three per site.

**The probabilities are calibrated before use.** The network emits independent
sigmoids that do not sum to one. Scoring them directly would put a logit on
something that is not a probability. The checkpoint's own empirical confidence
map — `CARBonAra.conf`, a smoothed CDF shipped beside the weights — runs first,
and every row is then checked for twenty finite entries in `[0, 1]` summing to
one before it is scored.

**The alphabet is sorted by abundance, not alphabetically.** Upstream
`src/data_encoding.py`:

```python
std_aminoacids = np.array([
    'LEU', 'GLU', 'ARG', 'LYS', 'VAL', 'ILE', 'PHE', 'ASP', 'TYR',
    'ALA', 'THR', 'SER', 'GLN', 'ASN', 'PRO', 'GLY', 'HIS', 'TRP',
    'MET', 'CYS'])
```

Asparagine is column **13**, serine 11, threonine 10, proline 14. Every
`probs_*` column this model writes is twenty entries in that order — not
twenty-one like ProteinMPNN's, and not alphabetical. Assuming otherwise is the
2026-08-20 defect exactly, so `verify_alphabet` runs against the checkout at
model load and stops the run if upstream ever reorders.

**MSE has no special case upstream, and that is an indexing hazard.**
ProteinMPNN maps `HETATM MSE` to methionine. CARBonAra's `clean_structure` does
not: dropped as a heteroatom the residue vanishes and every index after it
shifts; kept as MSE it is not in `std_aminoacids` and becomes a ligand. Structure
preparation converts it, and emits every record as `ATOM` so that nothing else
can be removed by the heteroatom flag either.

**The residue enumeration is friendlier than ProteinMPNN's.** `clean_structure`
renumbers observed residues consecutively from 1 in file order, so there is no
gap-filling and no placeholder token, and a numbering gap shifts nothing. The
mapping is therefore an identity — but it is *verified* against what CARBonAra
parsed rather than assumed, residue by residue, and a chain that disagrees is
dropped. An unverified identity mapping is precisely what the 25.3% ProteinMPNN
misindexing looked like before anyone checked.

#### Generation, and why it is not upstream's procedure

*Added 2026-08-28, reversing the earlier decision to omit generation entirely.*

Upstream's `imprint_sampling` samples from RAW confidences at a chosen
`imprint_ratio`, injecting a sampled sequence back as a prior. Neither the
uncalibrated values nor that free parameter has a counterpart in the other
models' protocols, so `design()` does not use it. Instead:

- nothing is imprinted (`yt = 0`), so the model sees the backbone and no residue
  identity at all — the interface requires this, since fixing anything would
  measure the constraint rather than the model;
- the calibrated distribution is sharpened by `p**(1/T)` at the frozen T = 0.1
  and sampled independently per position.

Design rows record `generation = independent_calibrated_sampling` with
`native_procedure = False`, so a retention table can never be mistaken for one
produced by CARBonAra's own sampler.

**One forward pass per chain, not per position** — the reverse of scoring. All 32
designs come from a single unconditioned pass, so retention is the cheap stage
here and the expensive one for ProteinMPNN.

#### What its retention rate does and does not support

CARBonAra is one-shot, so positions are sampled independently and the sequence
carries no correlation between them. Its retention therefore supports a
**within-CARBonAra occupancy-associated difference in independent-marginal
retention**, and nothing stronger. It does not license claiming the architecture
difference cancels in the paired contrast: that would require the correlation
loss to be equal in both arms, which has not been shown.

**The control arm, run 2026-08-28.** 262 of 262 matched secretory sites, zero
failures, 36 minutes — against 244 for scoring, because generation needs one
unconditioned pass per chain and scoring needs one per position.

    exact (closed form)      0.0963      NXS 0.080   NXT 0.111
    sampled, 32 designs      0.0949      mean |difference| 0.0090

The two agree to 0.0014 in the mean, which is the compatibility check passing:
the closed form is what sampling estimates, so a gap between them would mean the
two paths had diverged. **The occupied arm has not been run**, so there is no
paired result yet and the figure below is a control-arm rate, not an effect.

**Whether scoring marginals predict design retention depends on the model.**
An earlier draft claimed flatly that they do not, for any model, and that the
autoregressive models retain the sequon 1.4–1.8× more often than their own
marginals imply. Both statements were built on comparing marginals at T = 1
against retention measured at T = 0.1. Corrected, the picture splits:

| Model | independence-implied at T = 0.1 | measured |
|---|---:|---:|
| ProteinMPNN | 0.141 | 0.121 |
| ESM-IF | 0.192 | 0.151 |
| CARBonAra (control arm) | 0.103 | **0.096** |

For the autoregressive models the prediction overshoots by 15–20%, because
scoring conditions on native context and generation conditions on the model's
own output — different distributions, so one does not substitute for the other.

For CARBonAra the prediction lands within 7%, and its two distributions agree
closely in the mean at T = 0.1 (P(Asn) 0.248 design against 0.250 scoring,
P(Ser|Thr) 0.400 against 0.419). That is not a licence to skip the design pass:
the per-site correlation between the two P(Asn) estimates is only 0.667, so the
agreement is in aggregate and not site by site, and `expected_retention` is
computed from a design-time pass for that reason.

#### Environment and running it

CARBonAra is a script layout rather than an installable package, so the checkout
directory itself goes on `sys.path`. Clone it — the weights ship in the
repository under `model/save/` — and set `CARBONARA_DIR`, or leave it beside this
one. It needs `gemmi` and `blosum`, which are not core dependencies here.
Discovery is deferred to first model use, so nothing else notices its absence.

```bash
python pipeline/05_scoreability.py results/manifests/scoring_manifest_secretory.csv \
  results/manifests/scoreability_secretory_carbonara.csv --model carbonara

python pipeline/07_score.py results/manifests/scoring_manifest_secretory.csv \
  results/scores/scores_secretory_carbonara.csv --model carbonara --device cpu
```

#### Smoke test against the real checkpoint

Run 2026-08-27 on `s_v6_4_2022-09-16_11-51`, against the three chains
[`methods_sequon_indexing.md`](methods_sequon_indexing.md) already verifies
end-to-end for ProteinMPNN.

| Chain | Residues | Manifest idx | Mapped idx | Reads | Score |
|---|---|---|---|---|---|
| 4EBY:A | 200 | 27 | 27 | `NSS` | −1.2495 |
| 5H5Y:A | 286 | 226 | 226 | `NRS` | −0.8520 |
| 9G3Q:A | 402 | 181 | 181 | `NES` | −1.3440 |

Every residue count matches the documented figure, and **the mapping is the
identity on all three** — including the chains with 10 and 20 numbering gaps,
where ProteinMPNN needs +10 and +20 shifts. That is the predicted consequence of
consecutive renumbering, confirmed rather than assumed: the triplet check is what
establishes it. All nine probability rows carry twenty finite entries in [0, 1]
summing to one. All three guards fire on real data (wrong triplet, out-of-range
index, unevaluated position).

Protein-only preparation was exercised on real glycoproteins: 4EBY carries `NAG`
and `BMA`, 9G3Q carries `NAG` and `MLT`, and the prepared input for both contains
zero `HETATM` records and no non-amino-acid residue.

**Sequence recovery, as an independent check on the alphabet.** Nothing
imprinted, so this is plain inverse folding:

```
4EBY:A  62.0%      5H5Y:A  52.8%      9G3Q:A  59.5%
```

against upstream's reported ~51% on CATH. A scrambled column order would give
about 5%. This settles the alphabet against the model's own behaviour rather than
against a constant transcribed from a source file — the check the 2026-08-20
defect went a month without.

#### The Asn/Asp ambiguity, and what it actually turned out to be

*Revised 2026-08-28, after the diagnostic. The first version of this section said
the confusion was CARBonAra-specific and would compress its effect relative to
the other models. Both halves were wrong.*

At two of three smoke-test sites CARBonAra's most likely residue at the sequon
asparagine was aspartate. The two are isosteric — an oxygen where an amide sits —
and CARBonAra is shown no side chains, so the worry was that half of
`conditional_sequon_score` was a coin flip for this model in particular.

Measured across all four models on the same 262 sequons, it is not
CARBonAra-specific. Asparagine beats aspartate at 55% of sites for CARBonAra and
ESM-IF and 64% for ProteinMPNN, against 50% for a coin flip, and asparagine
recovery is roughly half of X-position recovery in **every** model, including the
two that see side chains elsewhere in the chain. It is an intrinsically hard call
from a backbone, not a defect of one architecture.

Nor is the term noise. The two halves of the score are near-uncorrelated
(r = −0.02 to 0.14) and both track the total, and pooling Asn with Asp barely
changes which sites rank highly (Spearman 0.88–0.95). The analysis is a paired
difference, so the near-constant shift pooling produces cancels.

What the diagnostic did establish is worth more than what it ruled out. The
paired heatmap ([Figure 17](figures_and_captions.md)) shows aspartate shifting
between the arms alongside asparagine:

| Model | ΔP(Asn) at position 1 | ΔP(Asp) at position 1 |
|---|---:|---:|
| ProteinMPNN | +0.038 | **+0.061** |
| ESM-IF | +0.052 | +0.046 |
| ESMC | **+0.045** | +0.012 |
| CARBonAra | +0.065 | **+0.072** |

For ProteinMPNN and CARBonAra the aspartate shift is the larger of the two. So
what a structure-conditioned model registers at an occupied site is better
described as a preference for an **amide-shaped residue** than as recognition of
asparagine. ESMC is the exception — large asparagine shift, negligible aspartate
shift — which is what a model reading the motif in the sequence rather than
inferring a shape should look like.

The ordering among the three structure models is not explained by how much
sequence each sees: ProteinMPNN conditions on every other native residue and has
the largest aspartate excess. That is unexplained, and should not be narrated
into a story.

#### The result

*Run on ARC 2026-08-27/28. All three sets, 100% coverage, zero failures.*

    occupied vs matched control, secretory:  +0.288 SD  [+0.130, +0.543]
    262 contrasts, 72 resampling units       directional, magnitude undetermined

That places it beside ProteinMPNN (+0.282) rather than ESM-IF (+0.431). Its
interval excludes zero but crosses the ±0.2 SD equivalence margin, which is the
same verdict ProteinMPNN and ESMC receive.

**Scores are not bit-reproducible across machines.** The same 262 secretory sites
run on a laptop and on ARC agree to a mean |Δ| of 0.009 and a maximum of 0.076
log-odds, and *no* site reproduces exactly. The cause is `extract_topology`,
which takes a hard `topk(D, 64)` nearest-neighbour cut: near-ties at that
boundary resolve differently across architectures, so the message-passing graph
itself differs. Being a discrete change, it moves the output by ~1e-2 rather than
the ~1e-6 float accumulation alone would give. Against a reference SD of 1.28
that is about 0.01 SD, negligible beside effects of 0.3 SD — but every arm of a
comparison must come from one environment, and CARBonAra numbers should not be
quoted beyond two decimal places.

#### Not done

No context-conditioned scoring: no glycan-present arm, no ligands, ions or water.
No generation, retention or stage 08, so CARBonAra appears in neither the
retention panel nor the masking panel of [Figure 15](figures_and_captions.md) —
with no motif-visible condition it has no within-model masking contrast at all.
Matching was not rerun; it is model-independent and shared by every variant.

### ProGen2 — causal, sequence-only

*Added 2026-08-28, from `score_proteins_progen2_colab.ipynb`.*

The notebook's forward pass is reused; its statistic is not. It reports
`total_logp / L`, a whole-chain mean, which is the protein-level score this
document already argues cannot resolve a three-residue question.

Fills the empty cell in the conditioning grid:

| | masked / bidirectional | causal |
|---|---|---|
| structure + sequence | ProteinMPNN, CARBonAra | ESM-IF |
| sequence only | ESMC | **ProGen2** |

so **ESM-IF is the like-for-like comparison, not ESMC**: both are causal and
prefix-only, and the difference between them isolates what the backbone adds
under identical conditioning. `conditioning` is ESM-IF's string.

#### The result, and how to read it

    secretory:  -0.004 SD  [-0.259, +0.327]   261 contrasts   inconclusive

Essentially zero, against ESMC's +0.261 [+0.063, +0.598]. Both are
sequence-only, so "sequence carries nothing" cannot be the explanation.

The difference is causality, and it is mechanical. ProGen2's P(Asn at i)
conditions on residues 1..i-1 and **has not seen the serine or threonine two
positions downstream**. ESMC, masked and bidirectional, sees the whole rest of
the chain including that hydroxyl, so its asparagine term can use the downstream
motif as evidence. ProGen2's structurally cannot.

Which places ProGen2 at about zero and ESMC's *motif-hidden* arm at -0.113 —
two different routes to removing motif visibility, landing in the same place.
ProGen2 arrives there by architecture rather than by masking. Stated carefully:
its interval is wide and includes ESMC's +0.261, so this is a consistent picture
rather than a demonstration.

#### Joint masking is a smaller manipulation here

`mask_mode="joint"` integrates the residue at the asparagine position out of the
+1 and +2 terms, weighted by the model's own distribution there rather than
uniformly. The asparagine row is returned unchanged, because a causal prefix
ending at i-1 has nothing of the motif to hide.

So two of three terms change, against three of three for ESMC and ESM-IF. **A
near-zero masking change for ProGen2 therefore means less than it would for
them.** Recorded as `autoregressive_prefix_marginalised`, ESM-IF's string for
the same operation.

#### Two limits worth carrying

**P(Asn) is 0.064-0.092** on the smoke chains, against 0.16-0.25 for every other
model. A causal model predicting from a truncated crystal chain is off its
training distribution: it expects a protein from the N-terminus and receives a
fragment starting at residue 21 or 25. That is a sharper version of the caveat
already noted for ESMC, because for a causal model the prefix *is* the
conditioning.

**A 2048-token context.** 3JAV:A is 2328 residues and cannot be scored; it fails
closed with a named reason. Truncating the prefix was rejected — it would answer
a different question silently.

### ESM3 — the structure track, switched on and off

*Added 2026-08-28, from `score_proteins_esm3_colab.ipynb`. Its masking scheme is
kept, its whole-chain statistic is not.*

Every other structure-versus-sequence comparison here is **between models** and
so confounds the question with architecture, training data and tokenisation.
ESM3 carries a structure track that can simply be withheld:

    struct_cond   VQ-VAE structure tokens from the backbone, intact
    seq_only      the same model, same tokeniser, same masking, no structure

Crossed with `mask_mode`, that is a 2x2 inside one model — structure on/off by
motif visible/hidden — and the four conditioning strings are distinct so no two
arms can be pooled.

Scoring is masked and bidirectional, as for ESMC. Only the three sequon
positions are masked, not all L as the notebook does, so a site costs one or
three forward passes rather than one per residue.

#### The structure track functions

Smoke test on the three chains `methods_sequon_indexing` verifies:

| Chain | struct_cond | seq_only | difference |
|---|---:|---:|---:|
| 4EBY:A | -1.874 | -3.265 | 1.390 |
| 5H5Y:A | -1.905 | -2.786 | 0.881 |
| 9G3Q:A | -1.524 | -2.020 | 0.495 |

Withholding structure lowers the score every time, and most of the loss is in
P(Asn) — 0.089 to 0.009 on 4EBY. **Three occupied sites, no controls**: this
shows the track works and can genuinely be switched off, and says nothing yet
about occupancy discrimination.

Its parse agrees exactly with the manifest on all three chains, unlike ESM-IF's,
which disagreed on about 5%. Checked per chain regardless.

## The recurring bug, which is really one bug

Three separate failures this week were the same failure wearing different
clothes: **an assumption about how a model represents its own input, believed
rather than checked.**

**ProteinMPNN's alphabet.** `mpnn_scoring.ALPHABET` held
`ARNDCQEGHILKMFPSTWYVX`, which is a local three-letter lookup table from inside
`parse_PDB_biounits` — not the model's token alphabet, which is
`ACDEFGHIKLMNPQRSTVWYX`. `p_asn_at_n` was reading P(aspartate). P(Ser) and P(Thr)
were correct by coincidence, because `S` and `T` happen to be two of only four
fixed points between the two orderings. The test asserted the constant against a
copy of itself, so it locked the defect in rather than catching it.

**ESM-IF's residue indices.** The manifest's `model_index` is a Biopython
ordinal; ESM-IF reads structures through biotite. They disagreed on about 5% of
matched-set sites. Trusting the index would have scored a neighbouring residue
and returned a perfectly plausible number.

**ESMC's token offset.** The tokenizer prepends `<cls>`, so sequence position *i*
is token index *i+1* — obvious, easy to get right, and exactly the kind of thing
that is silently wrong when a tokenizer changes.

The fix in all three cases is the same and it is cheap: **round-trip the model's
own representation back to something you already know.** Decode ProteinMPNN's `S`
tensor and check it reproduces the native sequence (19.97% with the wrong
alphabet, 99.53% with the right one). Map manifest indices through ESM-IF's
parser and check the residues reproduce the manifest's triplet. Encode a probe
sequence with ESMC's tokenizer and read it back.

All three checks now run in code — `_assert_token_offset` at load,
`check_triplet` per site, and an assertion in the Colab preflight that refuses to
score if alphabet agreement drops below 95%. None of them is clever. All of them
would have caught their defect on day one.

## Where the models genuinely differ, and why that is the point

The six models are not six measurements of one quantity. They condition on
different things, and the grid is now complete in both directions:

| Model | Sees | Conditional |
|---|---|---|
| ProteinMPNN | backbone + all other residues | bidirectional, 8 orders averaged |
| ESM-IF1 | backbone + native prefix | autoregressive, single pass |
| ESMC | sequence only | masked position, single pass |
| CARBonAra | backbone + all other residues | one-shot, single pass per position |
| ProGen2 | **sequence only** | autoregressive, single pass |
| ESM3 | backbone **or not**, switchable | masked position, single pass |

|  | masked / bidirectional | causal | one-shot |
|---|---|---|---|
| structure + sequence | ProteinMPNN, ESM3 (struct) | ESM-IF | CARBonAra |
| sequence only | ESMC, ESM3 (seq) | ProGen2 | — |

ESM3 appears twice on purpose: it is the only model that can occupy both rows,
which is what makes it the one within-model test of whether structure matters.

That is a **conditioning spectrum**, and it is more informative than three
attempts at the same number would be. Sequence-only versus structure-conditioned
is the comparison that says whether structure is doing any work. It is also why
the `conditioning` column exists and why raw scores are never pooled across
models.

A confound they share was found while building ESMC's masking. Masking only the
+2 residue still shows the model a native asparagine two positions upstream, and
N-X-S/T is a heavily learned motif — so a model can infer S/T from the N rather
than from anything about the site. Measured on ESMC, that is worth about 0.34
log-odds of the score. It inflates both arms of a matched pair and therefore
largely cancels in the paired contrast, but it compresses dynamic range, which
matters when the effects under discussion are 0.1–0.4 SD. All three models now have a joint-masking
variant: ESMC and ProteinMPNN by masking the positions, ESM-IF by marginalising
over them, since a causal decoder cannot hide an upstream residue from what
follows it. ESMC's and ProteinMPNN's have been run; ESM-IF's has not.

## The operational lessons

**Two packages can own the same name.** `fair-esm` (ESM-IF) and
EvolutionaryScale's `esm` (ESMC, ESM3) both install a top-level module called
`esm`. Not a pinnable version conflict — the same import name. They need separate
environments. The registry lazy-imports, so an absent model is unavailable rather
than fatal, and the test suite skips model-dependent tests rather than failing.

**Filenames are an interface.** The analysis stages hard-coded every input path.
With one model that is tidy; with three it is a silent-wrong-answer generator,
because corrected and per-model outputs land under new names and the stages went
on reading the old ones without complaint. Every stage now takes `--variant`,
the empty variant reproduces the original filenames byte-for-byte, and a variant
whose score file is missing **stops** rather than falling back.

**Provenance belongs to the data.** The analysis used to write
`"model": "ProteinMPNN v_48_020"` into its output JSON as a literal. It now reads
`model` / `conditioning` / `n_orders` from the score file. A label restated by
hand is a label that eventually lies.

## Adding the seventh

1. Write `adapters/<name>.py` implementing one or both protocols, plus
   `describe()` for the provenance columns.
2. Register it in `adapters/__init__.py`.
3. Add a round-trip check that the model's own representation reproduces
   something you already know.
4. Run stages 05, 07 and 08 with `--model <name>`, then the analysis with
   `--variant <name>`.

Nothing else should need touching. If something does, that is worth
investigating before working around it — it usually means a model-specific quirk
has escaped its adapter.

### Read next

- [`second_model_esm_if.md`](models.md) — ESM-IF's conditional, its
  index mapping, batched generation
- [`third_model_esmc.md`](models.md) — ESMC's masking schemes and the
  environment split
- [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) — the
  alphabet defect in full
