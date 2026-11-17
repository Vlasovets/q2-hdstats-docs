#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial"
WORKFLOW_DIR="${PROJECT_ROOT}/outputs/fresh-atacama-workflow"
OUT="${WORKFLOW_DIR}/q2-gglasso-top300"
QIIME="${PROJECT_ROOT}/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7/bin/qiime"

# Required on this local macOS environment so QIIME plugin discovery can import
# q2-composition/rpy2 without crashing before q2-gglasso is reached.
export R_HOME="${PROJECT_ROOT}/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7/lib/R"

TABLE="${OUT}/atacama-top300-53-clr-pseudocount1.qza"
SOLUTION="${OUT}/atacama-top300-53-clr-pseudocount1-slr-lambda065-mu-sweep-gamma025.qza"
METADATA="${OUT}/atacama-top300-53-robust-pca-metadata-q2.tsv"

for COLOR_BY in \
  elevation \
  ph \
  average-soil-relative-humidity \
  average-soil-temperature \
  percentcover \
  log10-extract-concen
do
  "${QIIME}" gglasso pca \
    --i-table "${TABLE}" \
    --i-solution "${SOLUTION}" \
    --m-sample-metadata-file "${METADATA}" \
    --p-n-components 2 \
    --p-color-by "${COLOR_BY}" \
    --o-visualization "${OUT}/atacama-top300-slr-rank2-robust-pca-color-by-${COLOR_BY}.qzv"
done
