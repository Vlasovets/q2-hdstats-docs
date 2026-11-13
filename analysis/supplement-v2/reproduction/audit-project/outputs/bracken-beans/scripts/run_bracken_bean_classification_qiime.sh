#!/usr/bin/env bash
# One compact QIIME-native bean-type classification workflow for the MOSHPIT Bracken table.
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

# Retain species observed in every sample and supported by >=800,000 reads in total.
"${QIIME}" feature-table filter-features \
  --i-table "${ART}/bracken-species-raw.qza" \
  --p-min-samples 14 \
  --p-min-frequency 800000 \
  --o-filtered-table "${ART}/bracken-species-prev14-total800k.qza"

"${QIIME}" classo transform-features \
  --i-features "${ART}/bracken-species-prev14-total800k.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "${ART}/bracken-species-prev14-total800k-clr05.qza"

"${QIIME}" feature-table filter-samples \
  --i-table "${ART}/bracken-species-prev14-total800k.qza" \
  --m-metadata-file "${ROOT}/metadata.tsv" \
  --p-where "[analysis_set]='train'" \
  --p-no-filter-empty-features \
  --o-filtered-table "${ART}/bracken-species-prev14-total800k-train.qza"

"${QIIME}" feature-table filter-samples \
  --i-table "${ART}/bracken-species-prev14-total800k.qza" \
  --m-metadata-file "${ROOT}/metadata.tsv" \
  --p-where "[analysis_set]='test'" \
  --p-no-filter-empty-features \
  --o-filtered-table "${ART}/bracken-species-prev14-total800k-test.qza"

"${QIIME}" classo transform-features \
  --i-features "${ART}/bracken-species-prev14-total800k-train.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "${ART}/bracken-species-prev14-total800k-train-clr05.qza"

"${QIIME}" classo transform-features \
  --i-features "${ART}/bracken-species-prev14-total800k-test.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "${ART}/bracken-species-prev14-total800k-test-clr05.qza"

"${QIIME}" classo classify \
  --i-features "${ART}/bracken-species-prev14-total800k-train-clr05.qza" \
  --m-y-file "${ROOT}/metadata.tsv" \
  --m-y-column bean_type \
  --p-no-path \
  --p-cv-numerical-method P-PDS \
  --p-cv-subsets 5 \
  --p-cv-seed 20261021 \
  --p-cv-one-se \
  --p-cv-lamin 0.01 \
  --p-stabsel-b 100 \
  --p-stabsel-seed 20261021 \
  --p-stabsel-method lam \
  --p-stabsel-lam 0.40 \
  --p-no-stabsel-true-lam \
  --p-stabsel-numerical-method P-PDS \
  --p-stabsel-threshold 0.4 \
  --p-stabsel-threshold-label 0.4 \
  --o-result "${ART}/bracken-bean-type-train-clr05-classo.qza"

"${QIIME}" classo predict \
  --i-features "${ART}/bracken-species-prev14-total800k-train-clr05.qza" \
  --i-problem "${ART}/bracken-bean-type-train-clr05-classo.qza" \
  --o-predictions "${ART}/bracken-bean-type-train-clr05-classo-predictions.qza"

"${QIIME}" classo predict \
  --i-features "${ART}/bracken-species-prev14-total800k-test-clr05.qza" \
  --i-problem "${ART}/bracken-bean-type-train-clr05-classo.qza" \
  --o-predictions "${ART}/bracken-bean-type-test-clr05-classo-predictions.qza"
