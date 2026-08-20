#!/bin/bash
# Run ONCE on the ARC login node. Compute nodes generally have no outbound
# network, so everything that needs the internet -- pip, git clone, model
# weights -- has to happen here and be cached on disk.
#
#   bash scripts/arc/glyco_setup.sh
#
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/chem-proteindesign/sjoh5764/glyco_occupancy}"
REPO_URL="${REPO_URL:-https://github.com/LBDillon/Glycan-occupancy-analysis.git}"
BUNDLE_URL="${BUNDLE_URL:-https://github.com/LBDillon/Glycan-occupancy-analysis/releases/download/bundle-2026-08-20/colab_bundle.tar}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
CUDA_MODULE="${CUDA_MODULE:-CUDA/12.6.0}"

module purge
module load "${PY_MODULE}" "${CUDA_MODULE}"

mkdir -p "${PROJECT_ROOT}"/{logs,cache,results}
cd "${PROJECT_ROOT}"

# ARC home directories have a small quota, and pip caches wheels -- torch alone
# is ~2.5 GB, twice over for two venvs. Left at their defaults these land in
# $HOME/.cache and blow the quota partway through, which surfaces later as an
# unrelated-looking "Disk quota exceeded" on the next thing that touches home.
# Point every cache and temp directory at /data before installing anything.
export PIP_CACHE_DIR="${PROJECT_ROOT}/cache/pip"
export XDG_CACHE_HOME="${PROJECT_ROOT}/cache/xdg"
export TMPDIR="${PROJECT_ROOT}/cache/tmp"
mkdir -p "${PIP_CACHE_DIR}" "${XDG_CACHE_HOME}" "${TMPDIR}"

# --- code -----------------------------------------------------------------
[[ -d module ]] || git clone --depth 1 "${REPO_URL}" module
[[ -d ProteinMPNN ]] || git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git ProteinMPNN

# --- environments ---------------------------------------------------------
# Two of them, because fair-esm (ESM-IF) and EvolutionaryScale's esm (ESMC)
# both install a top-level package called `esm` and cannot coexist.
if [[ ! -d venv-if ]]; then
  python -m venv venv-if
  ./venv-if/bin/pip install -q --upgrade pip
  ./venv-if/bin/pip install -q torch numpy pandas scipy biopython "biotite>=1.0" \
      fair-esm==2.0.0 torch-geometric
  # torch-scatter must match the torch build; torch-sparse is NOT needed.
  TV=$(./venv-if/bin/python -c "import torch;print(torch.__version__.split('+')[0])")
  CU=$(./venv-if/bin/python -c "import torch;print('cu'+torch.version.cuda.replace('.','') if torch.cuda.is_available() or torch.version.cuda else 'cpu')")
  ./venv-if/bin/pip install -q torch-scatter -f "https://data.pyg.org/whl/torch-${TV}+${CU}.html" || \
    echo "WARNING: torch-scatter wheel unavailable for torch-${TV}+${CU}; ESM-IF will not run"
fi

if [[ ! -d venv-esmc ]]; then
  python -m venv venv-esmc
  ./venv-esmc/bin/pip install -q --upgrade pip
  ./venv-esmc/bin/pip install -q torch numpy pandas scipy biopython "biotite>=1.0"
  # --no-deps: the SDK declares torchtext, which is dead against modern torch
  # and is never imported by ESMC.
  ./venv-esmc/bin/pip install -q --no-deps esm==3.2.2 "huggingface_hub<1.0" \
      "tokenizers>=0.21,<0.22" "transformers<4.48.2" tenacity httpx zstd \
      msgpack-numpy cloudpathlib brotli attrs einops regex safetensors
fi

# --- structures -----------------------------------------------------------
STRUCT="${PROJECT_ROOT}/module/data/cache/pdb"
if [[ ! -d "${STRUCT}" || -z "$(ls -A "${STRUCT}" 2>/dev/null)" ]]; then
  mkdir -p "${STRUCT}" bundle
  [[ -f bundle.tar ]] || wget -q --show-progress -O bundle.tar "${BUNDLE_URL}"
  tar -xf bundle.tar -C bundle
  find bundle/structures -name '*.gz' -exec sh -c 'gunzip -c "$1" > "'"${STRUCT}"'/$(basename "$1" .gz)"' _ {} \;
  find bundle/structures \( -name '*.pdb' -o -name '*.cif' \) -exec cp {} "${STRUCT}/" \;
  mkdir -p module/results/manifests module/results/matching
  cp bundle/manifests/*.csv module/results/manifests/ 2>/dev/null || true
  cp bundle/matching/*.csv  module/results/matching/  2>/dev/null || true
fi
echo "structures: $(ls -1 "${STRUCT}" | wc -l) (expect 1824)"

# --- model weights, pre-fetched because compute nodes are offline ---------
export TORCH_HOME="${PROJECT_ROOT}/cache/torch"
export HF_HOME="${PROJECT_ROOT}/cache/hf"
mkdir -p "${TORCH_HOME}" "${HF_HOME}"

./venv-if/bin/python - <<'PY' || echo "WARNING: ESM-IF weights not cached"
import esm.pretrained
esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
print("ESM-IF weights cached")
PY

./venv-esmc/bin/python - <<'PY' || echo "WARNING: ESMC weights not cached"
import torch
from esm.models.esmc import ESMC
ESMC.from_pretrained("esmc_300m", device=torch.device("cpu"))
print("ESMC weights cached")
PY

cat > "${PROJECT_ROOT}/env.sh" <<ENVEOF
export PROJECT_ROOT="${PROJECT_ROOT}"
export TORCH_HOME="${PROJECT_ROOT}/cache/torch"
export HF_HOME="${PROJECT_ROOT}/cache/hf"
export PIP_CACHE_DIR="${PROJECT_ROOT}/cache/pip"
export XDG_CACHE_HOME="${PROJECT_ROOT}/cache/xdg"
export TMPDIR="${PROJECT_ROOT}/cache/tmp"
export HF_HUB_OFFLINE=1
export PROTEINMPNN_DIR="${PROJECT_ROOT}/ProteinMPNN"
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
ENVEOF

echo
echo "setup complete: ${PROJECT_ROOT}"
echo "next: sbatch scripts/arc/glyco_retention.slurm esm_if"
