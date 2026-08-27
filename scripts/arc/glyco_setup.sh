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

# CARBonAra ships its own weights, so there is no separate checkpoint to fetch --
# but the repository is 1.1 GB, of which 838 MB is the authors' own result files.
# Sparse-checkout takes only the entry point, src/ and one checkpoint (~7 MB).
CARBONARA_MODEL="${CARBONARA_MODEL:-s_v6_4_2022-09-16_11-51}"
if [[ ! -d CARBonAra ]]; then
  git clone --depth 1 --filter=blob:none --sparse \
      https://github.com/LBM-EPFL/CARBonAra.git CARBonAra \
    && git -C CARBonAra sparse-checkout set src "model/save/${CARBONARA_MODEL}" \
    || { echo "sparse clone failed; falling back to a full one"; rm -rf CARBonAra;
         git clone --depth 1 https://github.com/LBM-EPFL/CARBonAra.git CARBonAra; }
fi
if [[ -s "CARBonAra/model/save/${CARBONARA_MODEL}/model.pt" ]]; then
  echo "CARBonAra weights: $(du -h "CARBonAra/model/save/${CARBONARA_MODEL}/model.pt" | cut -f1)"
else
  echo "WARNING: CARBonAra checkpoint missing; --model carbonara will not run"
fi

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
  # torch-scatter is a compiled extension with no wheel for current torch and no
  # source build without nvcc. It is optional here: the package ships a
  # native-torch shim for the two names ESM-IF imports, verified bit-identical
  # to the real thing. Try the wheel anyway -- if it is there, use it.
  echo "torch-scatter: trying a wheel for torch-${TV}+${CU} (optional)"
  ./venv-if/bin/pip install -q torch-scatter -f "https://data.pyg.org/whl/torch-${TV}+${CU}.html" 2>/dev/null \
    || ./venv-if/bin/pip install -q torch-scatter -f "https://data.pyg.org/whl/torch-${TV}+cpu.html" 2>/dev/null \
    || echo "  no wheel for torch-${TV}; ESM-IF will use the built-in shim instead"
fi

if [[ ! -d venv-esmc ]]; then
  python -m venv venv-esmc
  ./venv-esmc/bin/pip install -q --upgrade pip
  ./venv-esmc/bin/pip install -q torch numpy pandas scipy biopython "biotite>=1.0"
  # --no-deps: the SDK declares torchtext, which is dead against modern torch
  # and is never imported by ESMC.
  ./venv-esmc/bin/pip install -q --no-deps esm==3.2.2 "huggingface_hub<1.0" \
      "tokenizers>=0.21,<0.22" "transformers<4.48.2" tenacity httpx zstd \
      msgpack-numpy cloudpathlib brotli attrs einops regex safetensors \
      tqdm filelock pyyaml requests packaging
fi

# A third environment for CARBonAra. Not merged into venv-if: it needs gemmi,
# blosum, scikit-learn and h5py -- upstream's src/__init__.py imports its whole
# src package, so the scoring and dataset modules' dependencies are needed even
# though this integration calls none of them -- and adding four packages to a
# working ESM-IF environment to save 2.5 GB is not a trade worth making.
if [[ ! -d venv-carbonara ]]; then
  python -m venv venv-carbonara
  ./venv-carbonara/bin/pip install -q --upgrade pip
  ./venv-carbonara/bin/pip install -q torch numpy pandas scipy biopython \
      gemmi blosum scikit-learn h5py tqdm
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

# Fetch weights as FILES rather than by importing torch. `import torch` needs
# ~2 GB of address space, which login nodes cap at 2 GB, so a torch-based
# prefetch can only run on a compute node -- and compute nodes may have no
# outbound network. wget needs neither, so this works in both places.
IF_CKPT="${TORCH_HOME}/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt"
mkdir -p "$(dirname "${IF_CKPT}")"
if [[ ! -s "${IF_CKPT}" ]]; then
  echo "fetching ESM-IF weights (1.7 GB)..."
  wget -q --show-progress -O "${IF_CKPT}.part" \
    https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt \
    && mv "${IF_CKPT}.part" "${IF_CKPT}" \
    || { rm -f "${IF_CKPT}.part"; echo "WARNING: could not fetch ESM-IF weights"; }
fi
[[ -s "${IF_CKPT}" ]] && echo "ESM-IF weights: $(du -h "${IF_CKPT}" | cut -f1)"

if [[ ! -d "${HF_HOME}/hub/models--EvolutionaryScale--esmc-300m-2024-12" ]]; then
  echo "fetching ESMC weights..."
  ./venv-esmc/bin/huggingface-cli download EvolutionaryScale/esmc-300m-2024-12 \
      >/dev/null 2>&1 || echo "WARNING: could not fetch ESMC weights"
fi
[[ -d "${HF_HOME}/hub/models--EvolutionaryScale--esmc-300m-2024-12" ]] && echo "ESMC weights cached"

cat > "${PROJECT_ROOT}/env.sh" <<ENVEOF
export PROJECT_ROOT="${PROJECT_ROOT}"
export TORCH_HOME="${PROJECT_ROOT}/cache/torch"
export HF_HOME="${PROJECT_ROOT}/cache/hf"
export PIP_CACHE_DIR="${PROJECT_ROOT}/cache/pip"
export XDG_CACHE_HOME="${PROJECT_ROOT}/cache/xdg"
export TMPDIR="${PROJECT_ROOT}/cache/tmp"
export HF_HUB_OFFLINE=1
export PROTEINMPNN_DIR="${PROJECT_ROOT}/ProteinMPNN"
export CARBONARA_DIR="${PROJECT_ROOT}/CARBonAra"
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
ENVEOF

echo
echo "setup complete: ${PROJECT_ROOT}"
echo "next: sbatch scripts/arc/glyco_retention.slurm esm_if"
