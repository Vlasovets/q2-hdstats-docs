#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial"
CODEX_ROOT="/Users/cmueller/Documents/Codex/2026-06-26/Atacama-q2-gglasso-q2-classo-tutorial"
OUT="${ROOT}/outputs/fresh-atacama-workflow/q2-gglasso-top300"
CLASSO300="${ROOT}/outputs/fresh-atacama-workflow/q2-classo-top300"
TAXONOMY="${ROOT}/work/atacama-full/top-300-taxonomy.qza"
QIIME_ENV="${ROOT}/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
QIIME="${QIIME_ENV}/bin/qiime"

export PYTHONPATH="${CODEX_ROOT}/scripts/qiime_no_q2_composition"
export CONDA_PREFIX="${QIIME_ENV}"
export MPLCONFIGDIR="/private/tmp/mpl-cache"
export NUMBA_CACHE_DIR="/private/tmp/numba-cache"

mkdir -p "${OUT}"

"${QIIME}" gglasso transform-features \
  --i-table "${CLASSO300}/atacama-top300-prevalence-53-table.qza" \
  --i-taxonomy "${TAXONOMY}" \
  --p-transformation clr \
  --p-pseudo-count 1 \
  --o-transformed-table "${OUT}/atacama-top300-53-clr-pseudocount1.qza"

"${QIIME}" gglasso calculate-covariance \
  --i-table "${OUT}/atacama-top300-53-clr-pseudocount1.qza" \
  --p-method scaled \
  --p-bias \
  --o-covariance-matrix "${OUT}/atacama-top300-53-clr-pseudocount1-correlation.qza"

"${QIIME}" gglasso solve-problem \
  --i-covariance-matrix "${OUT}/atacama-top300-53-clr-pseudocount1-correlation.qza" \
  --p-n-samples 53 \
  --p-no-latent \
  --p-lambda1-path 1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40 0.35 0.30 \
  --p-gamma 0.30 \
  --o-solution "${OUT}/atacama-top300-53-clr-pseudocount1-sgl-ebic-gamma030.qza"

"${QIIME}" gglasso solve-problem \
  --i-covariance-matrix "${OUT}/atacama-top300-53-clr-pseudocount1-correlation.qza" \
  --p-n-samples 53 \
  --p-no-latent \
  --p-lambda1-path 1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40 0.35 0.30 \
  --p-gamma 0.10 \
  --o-solution "${OUT}/atacama-top300-53-clr-pseudocount1-sgl-ebic-gamma010.qza"
