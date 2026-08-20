# How the benchmark grew from one model to three

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

The three models are not three measurements of one quantity. They condition on
different things:

| Model | Sees | Conditional |
|---|---|---|
| ProteinMPNN | backbone + all other residues | bidirectional, 8 orders averaged |
| ESM-IF1 | backbone + native prefix | autoregressive, single pass |
| ESMC | sequence only | masked position, single pass |

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
matters when the effects under discussion are 0.1–0.4 SD. ESMC now has a
joint-masking variant as a sensitivity; ProteinMPNN could gain one easily, and
ESM-IF cannot without a causal decoder hiding the residue from everything
downstream too.

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

## Adding the fourth

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

- [`second_model_esm_if.md`](second_model_esm_if.md) — ESM-IF's conditional, its
  index mapping, batched generation
- [`third_model_esmc.md`](third_model_esmc.md) — ESMC's masking schemes and the
  environment split
- [`correction_2026-08-20_alphabet.md`](correction_2026-08-20_alphabet.md) — the
  alphabet defect in full
