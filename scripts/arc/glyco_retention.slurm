#!/bin/bash
# Sequon retention: the expensive stage. 32 unconstrained designs per chain,
# over 1725 chains, for one model.
#
#   sbatch scripts/arc/glyco_retention.slurm esm_if
#   sbatch scripts/arc/glyco_retention.slurm proteinmpnn
#   sbatch scripts/arc/glyco_retention.slurm carbonara
#
# CARBonAra is the CHEAP one here, unlike scoring: one unconditioned forward pass
# per chain gives every position's distribution at once, so 32 designs cost one
# pass rather than 32. Its memory is the risk instead -- extract_topology builds
# an atom-by-atom distance matrix, which killed a scoring shard at 32G.
#
# The array splits chains across tasks, so wall time is roughly the serial cost
# divided by the array width. Each task writes its own shard; merge afterwards
# with glyco_merge.py. Every task is resumable -- resubmitting skips work
# already on disk, so a timed-out task costs only the chain in flight.
#
#
# Runs CPU-only by default: ARC's short/medium partitions have no GPUs and
# reject --gres with "Invalid feature specification". Retention is the expensive stage, so the array is wide.
# --device auto picks CUDA when it is actually present, so adding a GPU is just:
#
#   sbatch --partition=<gpu-partition> --gres=gpu:1 --array=0-15 <this script> <model>
#
# Command-line options override the directives below.
#SBATCH --job-name=glyco-ret
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G          # carbonara wants more: see the note above
#SBATCH --time=12:00:00
#SBATCH --array=0-63
#SBATCH --output=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.out
#SBATCH --error=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/chem-proteindesign/sjoh5764/glyco_occupancy}"
source "${PROJECT_ROOT}/env.sh"
# env.sh pins OMP_NUM_THREADS=1 for reproducibility on a laptop; here the
# allocation is the unit, so let torch use the cores it was given.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

MODEL="${1:?usage: sbatch glyco_retention.slurm <esm_if|proteinmpnn|carbonara>}"
# ESM-IF and ProteinMPNN share venv-if; CARBonAra needs gemmi, blosum,
# scikit-learn and h5py. ESMC does no retention -- it has no SequenceDesigner.
case "${MODEL}" in
  carbonara) VENV="${PROJECT_ROOT}/venv-carbonara" ;;
  *)         VENV="${PROJECT_ROOT}/venv-if" ;;
esac
# The denominator is FIXED, not SLURM_ARRAY_TASK_COUNT. Resubmitting a few
# failed indices (--array=36,41) would otherwise set the count to 2 and compute
# "shard 36/2", which is invalid -- and silently a different partition of the
# chains even where it is not. It must match the width of the original array.
SHARDS="${SHARDS:-64}"
# See glyco_score.slurm: override to write beside a frozen set, never over it.
VARIANT="${VARIANT:-${MODEL}}"
SHARD="${SLURM_ARRAY_TASK_ID}/${SHARDS}"

module purge
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" "${CUDA_MODULE:-CUDA/12.6.0}"
cd "${PROJECT_ROOT}/module"

OUT="${PROJECT_ROOT}/results/retention"
mkdir -p "${OUT}"

# Both sources stage 10 reads. scoring_manifest carries all five classes;
# manifest_matched_secretory is the parallel comparison.
for pair in "scoring_manifest:scoring_manifest" "secretory:manifest_matched_secretory"; do
  tag="${pair%%:*}"; manifest="${pair##*:}"
  echo "=== ${MODEL} / ${tag} / shard ${SHARD} ==="
  "${VENV}/bin/python" pipeline/08_design.py \
      "results/manifests/${manifest}.csv" \
      "${OUT}/retention_${tag}_${VARIANT}.shard${SLURM_ARRAY_TASK_ID}.csv" \
      --model "${MODEL}" --device auto --shard "${SHARD}"
done

echo "shard ${SHARD} done"
