# Running the model stages on ARC

Colab was the wrong tool for retention. Sessions die at 12 hours, ProteinMPNN
projected to ~32 hours before batching, and a dropped runtime loses the work.
ARC removes all three problems and adds the one that matters: **job arrays**, so
1,725 chains run in parallel rather than in series.

**These jobs run CPU-only.** ARC's `short`/`medium` partitions have no GPUs and
reject `--gres` outright with `Invalid feature specification`. Throughput comes
from array width instead. `--device auto` picks CUDA if it is ever present, so
pointing this at a GPU partition needs no edit:

```bash
sbatch --partition=<gpu-partition> --gres=gpu:1 --array=0-15 \
       scripts/arc/glyco_retention.slurm esm_if
```

Scripts live in `scripts/arc/`, inside this module, so they travel with the standalone repository.

## Once, on the login node

Compute nodes have no outbound network, so everything needing the internet — pip,
git, model weights — happens here and is cached on disk.

```bash
bash scripts/arc/glyco_setup.sh
```

It creates `/data/chem-proteindesign/sjoh5764/glyco_occupancy/` containing:

| Path | What |
|---|---|
| `module/` | this repository |
| `ProteinMPNN/` | checkpoints ship inside the repo |
| `venv-if/` | torch, fair-esm, torch-geometric, torch-scatter |
| `venv-esmc/` | torch, EvolutionaryScale `esm` |
| `cache/torch`, `cache/hf` | pre-fetched weights |
| `module/data/cache/pdb/` | 1,824 structures from the release bundle |
| `env.sh` | sourced by every job |

**Two environments, and it is not optional.** `fair-esm` (ESM-IF) and
EvolutionaryScale's `esm` (ESMC) both install a top-level package named `esm`.
Installing one shadows the other.

`env.sh` sets `HF_HUB_OFFLINE=1` and `TORCH_HOME`, so a compute node never tries
to reach the network and fail halfway through a chain.

## Submitting

```bash
sbatch scripts/arc/glyco_score.slurm     esm_if
sbatch scripts/arc/glyco_score.slurm     esmc single
sbatch scripts/arc/glyco_score.slurm     esmc joint
sbatch scripts/arc/glyco_score.slurm     proteinmpnn

sbatch scripts/arc/glyco_retention.slurm esm_if
sbatch scripts/arc/glyco_retention.slurm proteinmpnn
```

Scoring uses `--array=0-7` on `short` (4 h); retention `--array=0-63` on
`medium` (12 h, 64 GB, 8 cores). Widen the array to trade queue time against wall
time — everything resumes, so a task that runs out of time costs only the chain
in flight.

## How sharding works, and why it is safe

`--shard K/N` keeps chain groups where `index % N == K`. Sharding is **by chain
group, not by site**, so a chain's expensive per-chain work happens exactly once
in exactly one task. Interleaving rather than contiguous blocks matters: chains
are ordered by PDB id, which correlates with nothing, so contiguous blocks would
hand one task all the long chains.

Verified as a genuine partition — every chain in exactly one shard, none twice —
and a three-shard run merges to a table **bit-identical** to the unsharded one
(max absolute score difference 0.000e+00).

Every task is resumable: each reads its own output first and skips finished work.
A timed-out task costs the chain in flight and nothing else. Resubmit the same
array indices.

## Merging — do not skip the checks

```bash
python pipeline/31_merge_shards.py \
    'results/retention/retention_scoring_manifest_esm_if.shard*.csv' \
    results/designs/retention_scoring_manifest_esm_if.csv \
    --expect-manifest results/manifests/scoring_manifest.csv --shards 16
```

Quote the glob so the shell does not expand it. The merge reports three things
the analysis cannot work out for itself:

- **every shard present** — a missing one is a silently short table, not an error;
- **no duplicate sites** — duplicates mean two tasks ran the same chain;
- **coverage against the manifest**, warning below 95%.

That last check exists because the first Colab retention run reported success
while dropping 21% of its sites: `08_design.py` records failures to a side file
and carries on. A short table that looks complete is the most dangerous output
this pipeline can produce.

Always read the merged `*_failures.csv` before treating a run as done.

## Then the analysis

Copy the merged tables back, and run the model-agnostic stages with the matching
variant:

```bash
python pipeline/09_analyse_scores.py optimal --variant esm_if
python pipeline/10_analyse_retention_by_class.py --variant esm_if
python pipeline/10b_analyse_retention_paired.py --variant esm_if
python pipeline/11_significance.py --variant esm_if
```

Filenames must follow the convention or the stages will not find them:
`scores_{set}_{variant}.csv` and `retention_{tag}_{variant}.csv`. A variant whose
score file is missing stops with an error rather than falling back — see
[`../README.md`](../README.md).

## Expected cost

| Stage | Model | Serial | 16-wide array |
|---|---|---|---|
| scoring | ESM-IF | ~35 min | minutes |
| scoring | ESMC | ~30 min/mode | minutes |
| retention | ESM-IF | hours | ~1 h |
| retention | ProteinMPNN | ~32 h before batching | ~2 h |

**Read the first ETA a task prints rather than trusting that table** — it is an
extrapolation, not a measurement.

## Memory, and the mistake worth not repeating

The first ARC submission lost **63 of 64 tasks to OUT_OF_MEMORY**. Retention
decodes several designs at once, and activation memory scales with
*batch x chain length*. These chains vary six-fold — 299 residues at the median,
1,287 at the longest — so a fixed batch of 32 is 4,800 residue-slots for a short
chain and 41,000 for a long one.

The adaptive halve-and-retry guard did not help, and could not: a **host** OOM is
delivered by the kernel's OOM killer or the SLURM cgroup, so the process dies
without Python seeing an exception. Only CUDA OOM is catchable. The batch has to
be bounded *before* allocating.

So the batch is now chosen from the chain length against a budget of 6,000
residue-slots — 32 designs at 150 residues, 20 at 299, 8 at 729, 4 at 1,287.
`--max-batch N` overrides it. If tasks still OOM, lower the budget or raise
`--mem` rather than narrowing the array.

Check exit states before believing a run finished:

```bash
sacct -X --name=glyco-ret --format=State -Pn | sort | uniq -c
```

Expect `64 COMPLETED`. `OUT_OF_MEMORY` there is why the merge coverage check
exists.
