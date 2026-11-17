#!/usr/bin/env bash
# Fit SGL and sparse-plus-low-rank graphical models to all filtered mOTU samples.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
QIIME="${ENV}/bin/qiime"
ART="${ROOT}/artifacts"

export PATH="${ENV}/bin:${PATH}"
export R_HOME="${ENV}/lib/R"
export MPLCONFIGDIR="${ROOT}/.mplconfig"
export NUMBA_CACHE_DIR="${ROOT}/.numba-cache"
# Keep QIIME's deployment cache in this reproducible workspace, not the
# protected per-user Library directory.
export CONDA_PREFIX="${ROOT}/.q2cli-env"
mkdir -p "${CONDA_PREFIX}/var/q2cli"
export PYTHONPATH="/Users/cmueller/Documents/Codex/2026-06-26/Atacama-q2-gglasso-q2-classo-tutorial/scripts/qiime_no_q2_composition"

TABLE="${ART}/motus-prev10pct.qza"
TAXONOMY="${ART}/motus-complete-taxonomy.qza"
LAMBDA=(3.00 2.00 1.50 1.25 1.00 0.80 0.60 0.50 0.40 0.30 0.25 0.20 0.15 0.10 0.075 0.05)
MU=(60 45 30 20 15 10 7.5 5 2.5 1)

# Standard SGL follows the mCLR workflow used for the Atacama network examples.
"${QIIME}" gglasso transform-features \
  --i-table "${TABLE}" \
  --i-taxonomy "${TAXONOMY}" \
  --p-transformation mclr \
  --o-transformed-table "${ART}/motus-60-mclr.qza"

"${QIIME}" gglasso calculate-covariance \
  --i-table "${ART}/motus-60-mclr.qza" \
  --p-method scaled \
  --p-bias \
  --o-covariance-matrix "${ART}/motus-60-mclr-correlation.qza"

"${QIIME}" gglasso solve-problem \
  --i-covariance-matrix "${ART}/motus-60-mclr-correlation.qza" \
  --p-n-samples 34 \
  --p-no-latent \
  --p-lambda1-path "${LAMBDA[@]}" \
  --p-gamma 0.25 \
  --o-solution "${ART}/motus-60-mclr-sgl-ebic-gamma025.qza"

# The latent model uses CLR with a unit zero replacement, as in the S2 robust-PCA workflow.
"${QIIME}" gglasso transform-features \
  --i-table "${TABLE}" \
  --i-taxonomy "${TAXONOMY}" \
  --p-transformation clr \
  --p-pseudo-count 1 \
  --o-transformed-table "${ART}/motus-60-clr-pseudocount1.qza"

"${QIIME}" gglasso calculate-covariance \
  --i-table "${ART}/motus-60-clr-pseudocount1.qza" \
  --p-method scaled \
  --p-bias \
  --o-covariance-matrix "${ART}/motus-60-clr-pseudocount1-correlation.qza"

"${QIIME}" gglasso solve-problem \
  --i-covariance-matrix "${ART}/motus-60-clr-pseudocount1-correlation.qza" \
  --p-n-samples 34 \
  --p-latent \
  --p-lambda1-path "${LAMBDA[@]}" \
  --p-mu1-path "${MU[@]}" \
  --p-gamma 0.25 \
  --o-solution "${ART}/motus-60-clr-slr-ebic-gamma025.qza"
