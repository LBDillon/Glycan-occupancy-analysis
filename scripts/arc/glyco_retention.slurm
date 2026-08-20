#!/bin/bash
# Sequon retention: the expensive stage. 32 unconstrained designs per chain,
# over 1725 chains, for one model.
#
#   sbatch scripts/arc/glyco_retention.slurm esm_if
#   sbatch scripts/arc/glyco_retention.slurm proteinmpnn
#
# The array splits chains across tasks, so wall time is roughly the serial cost
# divided by the array width. Each task writes its own shard; merge afterwards
# with glyco_merge.py. Every task is resumable -- resubmitting skips work
# already on disk, so a timed-out task costs only the chain in flight.
#
#SBATCH --job-name=glyco-ret
#SBATCH --partition=medium
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-15
#SBATCH --output=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.out
#SBATCH --error=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/chem-proteindesign/sjoh5764/glyco_occupancy}"
source "${PROJECT_ROOT}/env.sh"

MODEL="${1:?usage: sbatch glyco_retention.slurm <esm_if|proteinmpnn>}"
VENV="${PROJECT_ROOT}/venv-if"     # both live here; ESMC does not do retention
SHARD="${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}"

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
      "${OUT}/retention_${tag}_${MODEL}.shard${SLURM_ARRAY_TASK_ID}.csv" \
      --model "${MODEL}" --device cuda --shard "${SHARD}"
done

echo "shard ${SHARD} done"
