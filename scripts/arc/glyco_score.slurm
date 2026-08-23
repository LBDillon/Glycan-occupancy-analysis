#!/bin/bash
# Conditional sequon scores. Far cheaper than retention -- one pass per chain
# for ESM-IF and ESMC, one decoder pass per position for ProteinMPNN -- so a
# narrow array on a short partition is enough.
#
#   sbatch scripts/arc/glyco_score.slurm esm_if
#   sbatch scripts/arc/glyco_score.slurm esmc single
#   sbatch scripts/arc/glyco_score.slurm proteinmpnn
#
#
# Runs CPU-only by default: ARC's short/medium partitions have no GPUs and
# reject --gres with "Invalid feature specification". Scoring is cheap; a narrow array is plenty.
# --device auto picks CUDA when it is actually present, so adding a GPU is just:
#
#   sbatch --partition=<gpu-partition> --gres=gpu:1 --array=0-15 <this script> <model>
#
# Command-line options override the directives below.
#SBATCH --job-name=glyco-score
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-7
#SBATCH --output=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.out
#SBATCH --error=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/chem-proteindesign/sjoh5764/glyco_occupancy}"
source "${PROJECT_ROOT}/env.sh"
# env.sh pins OMP_NUM_THREADS=1 for reproducibility on a laptop; here the
# allocation is the unit, so let torch use the cores it was given.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

MODEL="${1:?usage: sbatch glyco_score.slurm <proteinmpnn|esm_if|esmc> [mask_mode]}"
MASK_MODE="${2:-}"
# The denominator is FIXED, not SLURM_ARRAY_TASK_COUNT. Resubmitting a few
# failed indices (--array=36,41) would otherwise set the count to 2 and compute
# "shard 36/2", which is invalid -- and silently a different partition of the
# chains even where it is not. It must match the width of the original array.
SHARDS="${SHARDS:-8}"
SHARD="${SLURM_ARRAY_TASK_ID}/${SHARDS}"

# ESMC needs the other environment: fair-esm and EvolutionaryScale's esm both
# claim the import name `esm` and cannot be installed together.
if [[ "${MODEL}" == "esmc" ]]; then VENV="${PROJECT_ROOT}/venv-esmc"; else VENV="${PROJECT_ROOT}/venv-if"; fi

VARIANT="${MODEL}${MASK_MODE:+_${MASK_MODE}}"
EXTRA=(); [[ -n "${MASK_MODE}" ]] && EXTRA=(--mask-mode "${MASK_MODE}")

module purge
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" "${CUDA_MODULE:-CUDA/12.6.0}"
cd "${PROJECT_ROOT}/module"

OUT="${PROJECT_ROOT}/results/scores"
mkdir -p "${OUT}"

for pair in "dataset:candidate_manifest_dataset" "controls:manifest_matched_controls" "secretory:manifest_matched_secretory"; do
  tag="${pair%%:*}"; manifest="${pair##*:}"
  echo "=== ${VARIANT} / ${tag} / shard ${SHARD} ==="
  "${VENV}/bin/python" pipeline/07_score.py \
      "results/manifests/${manifest}.csv" \
      "${OUT}/scores_${tag}_${VARIANT}.shard${SLURM_ARRAY_TASK_ID}.csv" \
      --model "${MODEL}" --device auto --shard "${SHARD}" "${EXTRA[@]}"
done

echo "shard ${SHARD} done"
