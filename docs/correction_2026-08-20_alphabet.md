# Correction, 2026-08-20 — ProteinMPNN's token alphabet

`mpnn_scoring.ALPHABET` held the wrong string. Every ProteinMPNN score and every
retention figure produced before this date needs regenerating.

## What was wrong

```python
ALPHABET = "ARNDCQEGHILKMFPSTWYVX"      # wrong
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"      # ProteinMPNN's actual token alphabet
```

The old string is real, but it is not the model's alphabet. It is `alpha_1`, a
local three-letter-to-one-letter lookup table defined *inside*
`parse_PDB_biounits` (`protein_mpnn_utils.py:61`), where it is paired with
`alpha_3` and ends in `-`. It never leaves that function.

The alphabet the model speaks is set in `tied_featurize`
(`protein_mpnn_utils.py:193`) and read back in `_S_to_seq`
(`protein_mpnn_utils.py:50`). Both use `ACDEFGHIKLMNPQRSTVWYX`.

The comment above the constant said "index into probability vectors with this,
never with a hand-written ordering." The intent was right; the string was a
hand-written ordering.

## How it was found

Not by inspection — by an index audit while adding ESM-IF. Decoding
ProteinMPNN's own `S` tensor back to text with each candidate alphabet, over
2,559 residues of the dataset manifest:

| Alphabet | Reproduces the native sequence |
|---|---|
| `ARNDCQEGHILKMFPSTWYVX` (old) | **19.97%** |
| `ACDEFGHIKLMNPQRSTVWYX` (correct) | **99.53%** |

The residual 0.5% is residues `parse_PDB` maps to `X`, not disagreement. On
1A2W/A the correct alphabet returns `KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRC`,
which is the chain's actual sequence; the old one returns
`HDTAAAHCDPFELNSSTSAASSSKVRKFLLHSPKITHNPR`.

## What it affected, and what it did not

Only four of the 21 tokens are fixed points between the two orderings — `A`,
`S`, `T` and `X` — which is why the damage is uneven rather than total.

### Scores

| Column | Index used | What ProteinMPNN actually returned |
|---|---|---|
| `p_asn_at_n` | 2 | **P(aspartate)** — wrong |
| `p_ser_at_plus2` | 15 | P(serine) — correct, by coincidence |
| `p_thr_at_plus2` | 16 | P(threonine) — correct, by coincidence |
| `p_pro_at_plus1` | 14 | **P(arginine)** — wrong; reported, not scored |

`conditional_sequon_score` is `0.5 * (logit(p_asn) + logit(p_ser_or_thr))`. The
hydroxyl term is sound. **The asparagine term was measuring aspartate**, so
every conditional sequon score is half right and half a different residue.

Aspartate is asparagine's closest structural neighbour — same shape, carboxylate
instead of amide — so the contaminated term is correlated with the intended one
rather than random. That makes the error *less* likely to be obvious in a
distribution and no less fatal to the claim.

### Retention

Worse, because `retention.design_sequences` decoded sampled sequences with the
same constant. A designed asparagine came back as `K`; a decoded `N` was really
`D`. So `classify_retention`'s asparagine test could not fire correctly, while
its serine/threonine test was fine. Every retention number is affected.

`omit_AAs_np` is unaffected — `X` is at index 20 in both orderings.

### Not affected

- **Scoreability.** It depends on ProteinMPNN's backbone mask, not on residue
  identity, so `scoreability_*.csv` stands.
- **Matching.** Built from RSA, neighbour counts and hydrophobic fraction —
  never from model output.
- **The manifests**, the evidence layers, the control-set construction.
- **ESM-IF.** Its adapter indexes through `alphabet.get_idx(...)`, ESM-IF's own
  dictionary, so it never had a hand-written ordering to get wrong.

## Why the tests did not catch it

`tests/test_mpnn_scoring.py` asserted:

```python
assert ALPHABET == "ARNDCQEGHILKMFPSTWYVX"
```

That is a copy of the value under test, not an independent check, so it locked
the defect in place. It now asserts the verified alphabet plus its fixed points,
and the Colab notebook's preflight decodes ProteinMPNN's `S` tensor and refuses
to score if agreement falls below 95% — a check that would have failed loudly on
day one.

**The general lesson, and it is the same one as the other three corrections:** a
constant copied from a model's source is an assumption, not a fact, until
something round-trips it against that model's own behaviour.

## What has to be regenerated

Everything downstream of a ProteinMPNN score:

- `results/scores/scores_{dataset,controls,secretory}.csv`
- `results/designs/*retention*.csv`
- `results/analysis/analysis_*.json`, `contrasts_*.csv`, `significance.csv`
- every figure in `results/figures/`
- the numbers quoted in `docs/primary_result.md`, `docs/significance.md`,
  `docs/figures.md`, and the reference SD in `config/scoring_frozen.toml`

The reference SD is pooled over all scoreable dataset sites, so it moves too:
every effect size expressed in SD units rescales even where a raw contrast
happens not to.
