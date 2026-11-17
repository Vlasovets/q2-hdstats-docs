#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial"
CODEX_ROOT="/Users/cmueller/Documents/Codex/2026-06-26/Atacama-q2-gglasso-q2-classo-tutorial"
OUT="${ROOT}/outputs/fresh-atacama-workflow/q2-gglasso-top300"
QIIME_ENV="${ROOT}/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
QIIME="${QIIME_ENV}/bin/qiime"

export PYTHONPATH="${CODEX_ROOT}/scripts/qiime_no_q2_composition"
export CONDA_PREFIX="${QIIME_ENV}"
export MPLCONFIGDIR="/private/tmp/mpl-cache"
export NUMBA_CACHE_DIR="/private/tmp/numba-cache"

mkdir -p "${OUT}"

"${QIIME}" gglasso solve-problem \
  --i-covariance-matrix "${OUT}/atacama-top300-53-clr-pseudocount1-correlation.qza" \
  --p-n-samples 53 \
  --p-latent \
  --p-lambda1-path 0.65 \
  --p-mu1-path 60 50 40 30 25 20 15 12.5 10 7.5 5 \
  --p-gamma 0.25 \
  --o-solution "${OUT}/atacama-top300-53-clr-pseudocount1-slr-lambda065-mu-sweep-gamma025.qza"
