# ProGen2 Masked Rerun Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correct ProGen2 sequence framing and whole-sequon marginalisation, verify the corrected estimator, and launch isolated full-corpus single and joint runs on ARC.

**Architecture:** Keep ProGen2 behind the existing `SequonScorer` adapter and the model-agnostic stage-07 runner. Use the official forward-direction token `1`, and implement the joint arm as the same 16-sample sequential marginalisation used by ESM-IF: sample N from the native prefix, sample X conditional on N, then average the +1 and +2 distributions. Add a dedicated ARC environment and write corrected runs under new variants so the invalid first pass cannot be overwritten or mistaken for the rerun.

**Tech Stack:** Python, PyTorch, Hugging Face Transformers, pandas, pytest, SLURM.

**Spec:** `docs/models.md` ProGen2 section plus the user request in this task; official framing reference is Salesforce ProGen2 `progen2/likelihood.py`.

## Global Constraints

- The benchmark outcome remains an occupancy-associated statistical preference, not evidence of model understanding.
- The asparagine term is unchanged in a causal model because downstream sequon residues cannot enter its prefix.
- Whole-sequon marginalisation must integrate both N and X out of the +2 prefix.
- Corrected outputs must use new variant names and must not overwrite today's invalid `progen2` files.
- The frozen marginal sample count is 16 with seed 0, matching ESM-IF.

---

### Task 1: Correct sequence framing

**Files:**
- Modify: `src/experimental_glycosylation_sites/progen2_scoring.py`
- Test: `tests/test_progen2_scoring.py`

**Interfaces:**
- Consumes: a Hugging Face tokenizer with one token for the literal `1`.
- Produces: `direction_token_id(tokenizer) -> int`; `load_model(...)` returns that token as its third element.

- [x] **Step 1: Write a failing regression test**

```python
def test_direction_token_is_used_when_gpt2_bos_is_different():
    tokenizer = FakeTokenizer()
    tokenizer.bos_token_id = VOCAB.index("<|eos|>")
    assert pg.direction_token_id(tokenizer) == VOCAB.index("1")
```

- [x] **Step 2: Verify the test fails because `direction_token_id` does not exist**

Run: `pytest -q tests/test_progen2_scoring.py::test_direction_token_is_used_when_gpt2_bos_is_different`

- [x] **Step 3: Implement and round-trip the literal direction token**

```python
def direction_token_id(tokenizer) -> int:
    ids = tokenizer("1", add_special_tokens=False)["input_ids"]
    if len(ids) != 1 or tokenizer.convert_ids_to_tokens(ids)[0] != "1":
        raise TokenisationError("forward direction marker '1' is not one token")
    return int(ids[0])
```

- [x] **Step 4: Replace tokenizer BOS use with this token and rerun the focused test**

Run: `pytest -q tests/test_progen2_scoring.py`

### Task 2: Marginalise both visible prefix residues

**Files:**
- Modify: `src/experimental_glycosylation_sites/progen2_scoring.py`
- Modify: `src/experimental_glycosylation_sites/adapters/progen2.py`
- Test: `tests/test_progen2_scoring.py`

**Interfaces:**
- Consumes: `marginalised_probabilities(..., n_samples=16, seed=0, max_batch=8)`.
- Produces: a probability matrix where N is native-prefix, +1 averages over sampled N, and +2 averages over sampled N and X.

- [x] **Step 1: Write a failing test whose +2 prediction echoes X**

The fake model predicts W at X after sampling N and echoes the token at X when predicting +2. Assert the joint +2 argmax is W rather than the native K.

- [x] **Step 2: Verify the existing N-only implementation fails that assertion**

Run: `pytest -q tests/test_progen2_scoring.py::test_joint_marginalisation_integrates_out_x_position`

- [x] **Step 3: Implement sequential 16-draw marginalisation**

Draw N from the restricted standard-amino-acid distribution at the N row, forward the batch, draw X from each +1 row, forward again, and average the relevant rows. Preserve batching and OOM backoff.

- [x] **Step 4: Record estimator provenance in the adapter**

For joint mode, `describe()` adds `marginal_samples: 16`; pass `marginal_samples` and `seed` into the scorer.

- [x] **Step 5: Run the focused and complete unit suites**

Run: `pytest -q tests/test_progen2_scoring.py tests/test_adapters.py`

### Task 3: Make ARC capable of running ProGen2

**Files:**
- Modify: `scripts/arc/glyco_setup.sh`
- Modify: `scripts/arc/glyco_score.slurm`
- Modify: `docs/running_on_arc.md`

**Interfaces:**
- Consumes: `sbatch scripts/arc/glyco_score.slurm progen2 joint`.
- Produces: `venv-progen2`, an offline cached `hugohrban/progen2-base`, and sharded `scores_*_<variant>.csv` outputs.

- [x] **Step 1: Add `venv-progen2` with Transformers 4.40-4.48 and required scientific packages**

- [x] **Step 2: Prefetch the checkpoint on the login node before enabling offline mode**

- [x] **Step 3: Route `progen2` jobs to the new environment and update usage text**

- [x] **Step 4: Verify both shell scripts parse**

Run: `bash -n scripts/arc/glyco_setup.sh scripts/arc/glyco_score.slurm`

### Task 4: Withdraw invalid outputs from current interpretation

**Files:**
- Modify: `docs/models.md`
- Modify: `README.md`
- Modify: `pipeline/25_figures_model_comparison.py`

**Interfaces:**
- Consumes: corrected variants `progen2_direction1` and `progen2_joint_direction1`.
- Produces: documentation and figures that cannot silently reuse the wrong-start-token run.

- [x] **Step 1: Mark today's `progen2` result superseded by the token-framing audit**

- [x] **Step 2: Point comparison figures at corrected variant names**

- [x] **Step 3: Document that ProGen2 joint masking is a causal marginalisation, not literal mask-token substitution**

### Task 5: Verify and launch

**Files:**
- Verify only: all files above and the remote ARC checkout.

**Interfaces:**
- Consumes: clean tests, cached model, structures, manifests, and an updated ARC checkout.
- Produces: SLURM job identifiers for isolated corrected single and joint runs.

- [x] **Step 1: Run the full local test suite**

Run: `pytest -q`

- [x] **Step 2: Run one real-checkpoint single/joint smoke test with the literal `1` token**

- [ ] **Step 3: Update the ARC checkout without overwriting its result directory**

- [ ] **Step 4: Verify remote environment, checkpoint, manifests, and empty destination variants**

- [ ] **Step 5: Submit the corrected single and whole-sequon-marginalised arrays**

Use distinct variants and matching shard counts; capture both job IDs immediately.
