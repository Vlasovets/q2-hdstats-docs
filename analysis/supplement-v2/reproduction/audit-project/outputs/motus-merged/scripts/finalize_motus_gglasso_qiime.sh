#!/usr/bin/env bash
# Create QIIME visualizations from the selected gamma=0.10 mOTU fits.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
QIIME="${ENV}/bin/qiime"
ART="${ROOT}/artifacts"

export PATH="${ENV}/bin:${PATH}"
export R_HOME="${ENV}/lib/R"
export MPLCONFIGDIR="${ROOT}/.mplconfig"
export NUMBA_CACHE_DIR="${ROOT}/.numba-cache"
export CONDA_PREFIX="${ROOT}/.q2cli-env"
mkdir -p "${CONDA_PREFIX}/var/q2cli"
export PYTHONPATH="/Users/cmueller/Documents/Codex/2026-06-26/Atacama-q2-gglasso-q2-classo-tutorial/scripts/qiime_no_q2_composition"

"${QIIME}" gglasso summarize \
  --i-solution "${ART}/motus-60-mclr-sgl-ebic-gamma010.qza" \
  --o-visualization "${ART}/motus-60-mclr-sgl-ebic-gamma010-summary.qzv"

"${QIIME}" gglasso summarize \
  --i-solution "${ART}/motus-60-clr-slr-ebic-gamma010.qza" \
  --o-visualization "${ART}/motus-60-clr-slr-ebic-gamma010-summary.qzv"

# The eBIC-selected low-rank component has rank 1.  q2-gglasso pca plots
# PC1 against PC2 and explicitly requires a rank >= 2 latent subspace, so no
# PCA visualization is valid for this selected solution.
