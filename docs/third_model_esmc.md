# The third model: ESMC 300M

Added 2026-08-20. `esmc_300m`, implementing `SequonScorer` only.

The first model in the benchmark that sees **no structure**. ProteinMPNN and
ESM-IF both condition on a backbone; ESMC conditions only on surrounding
sequence, so it answers a question neither can:

> Does sequence context alone distinguish occupied sequons from structurally
> matched sequons carrying no glycan?

That makes it the control the structure-conditioned results need. If ESMC
reproduces their effect, the effect need not be structural. If it does not, the
structure-conditioned models are doing work sequence alone cannot.

## Environment: ESMC and ESM-IF cannot coexist

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

## What is scored, and on which sequence

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

## The token offset, verified rather than assumed

The tokenizer prepends `<cls>`, so sequence position *i* is token index *i + 1*.
`_assert_token_offset` round-trips a probe sequence through the tokenizer at load
and raises if it does not reproduce it.

This is the check the ProteinMPNN alphabet defect went undetected for. An
assumption about how a model indexes or encodes its own input is not a fact until
something reproduces the input from it.

## Masking: two estimands, not two estimates

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
the rest of the native sequon. Neither currently has a joint-masking variant.
Adding one to ProteinMPNN is straightforward (its `conditional_probs` already
takes a position set); ESM-IF is harder, since a causal decoder cannot hide an
upstream residue without also hiding it from everything downstream.

## Why there is no SequenceDesigner

ESMC is a masked language model. Sampling from it would condition on sequence
rather than on a backbone, so "retention" would not mean what it means for the
inverse-folding models. Scoring only.

## Running it

```bash
esmcenv/bin/python pipeline/07_score.py <manifest> <out> --model esmc
```

Cost on CPU: ~1.6 s/site under `single` (three masked variants per sequon,
batched per chain), ~0.6 s/site under `joint` (one). Cheaper than either
inverse-folding model.

## Not done

- ESM3-open, which would give sequence-only / structure-only / sequence+structure
  from one set of weights. Gated behind a licence acceptance and a token, ~1.4B
  parameters. The function and annotation tracks must be left empty, or they can
  leak the glycosylation label.
- ESM-2, which occupies the same position in the benchmark as ESMC.
- The full-UniProt-sequence sensitivity described above.
