#!/bin/bash
# Conditional sequon scores. Far cheaper than retention -- one pass per chain
# for ESM-IF and ESMC, one decoder pass per position for ProteinMPNN -- so a
# narrow array on a short partition is enough.
#
#   sbatch scripts/arc/glyco_score.slurm esm_if
#   sbatch scripts/arc/glyco_score.slurm esmc single
#   sbatch scripts/arc/glyco_score.slurm proteinmpnn
#
#SBATCH --job-name=glyco-score
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.out
#SBATCH --error=/data/chem-proteindesign/sjoh5764/glyco_occupancy/logs/%x-%A_%a.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/chem-proteindesign/sjoh5764/glyco_occupancy}"
source "${PROJECT_ROOT}/env.sh"

MODEL="${1:?usage: sbatch glyco_score.slurm <proteinmpnn|esm_if|esmc> [mask_mode]}"
MASK_MODE="${2:-}"
SHARD="${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}"

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
      --model "${MODEL}" --device cuda --shard "${SHARD}" "${EXTRA[@]}"
done

echo "shard ${SHARD} done"
