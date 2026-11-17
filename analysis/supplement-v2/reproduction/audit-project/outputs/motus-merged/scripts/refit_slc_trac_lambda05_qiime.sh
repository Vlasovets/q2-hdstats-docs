#!/usr/bin/env bash
# Refit the already prepared QIIME 2 mOTU designs over lambda/lambda_max >= 0.05.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
export PATH="$ENV/bin:$PATH"
export R_HOME="$ENV/lib/R"
export MPLCONFIGDIR="$ROOT/.mplconfig"
export NUMBA_CACHE_DIR="$ROOT/.numba_cache"
ART="$ROOT/artifacts"
META="$ROOT/exports/motus-outcomes-all36.tsv"

qiime classo regress \
  --i-features "$ART/motus-clr-pseudocount0.5-train.qza" \
  --m-y-file "$META" --m-y-column host_age_days \
  --p-concomitant \
  --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.05 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
  --p-cv--nlam 60 --p-cv-lamin 0.05 --p-cv-logscale \
  --p-stabsel --p-stabsel-seed 1 --p-stabsel-b 100 --p-stabsel-q 10 \
  --p-stabsel-threshold 0.4 --p-stabsel-threshold-label 0.2 \
  --p-no-lamfixed \
  --o-result "$ART/motus-age-slc-r3.qza"

qiime classo predict --i-features "$ART/motus-clr-pseudocount0.5-train.qza" --i-problem "$ART/motus-age-slc-r3.qza" --o-predictions "$ART/motus-age-slc-r3-train-predictions.qza"
qiime classo predict --i-features "$ART/motus-clr-pseudocount0.5-test.qza" --i-problem "$ART/motus-age-slc-r3.qza" --o-predictions "$ART/motus-age-slc-r3-test-predictions.qza"

qiime classo regress \
  --i-features "$ART/motus-clr-trac-train.qza" --i-weights "$ART/motus-trac-weights.qza" \
  --m-y-file "$META" --m-y-column host_age_days \
  --p-concomitant \
  --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.05 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
  --p-cv--nlam 60 --p-cv-lamin 0.05 --p-cv-logscale \
  --p-stabsel --p-stabsel-seed 1 --p-stabsel-b 100 --p-stabsel-q 10 \
  --p-stabsel-threshold 0.4 --p-stabsel-threshold-label 0.2 \
  --p-no-lamfixed \
  --o-result "$ART/motus-age-trac-r3.qza"

qiime classo predict --i-features "$ART/motus-clr-trac-train.qza" --i-problem "$ART/motus-age-trac-r3.qza" --o-predictions "$ART/motus-age-trac-r3-train-predictions.qza"
qiime classo predict --i-features "$ART/motus-clr-trac-test.qza" --i-problem "$ART/motus-age-trac-r3.qza" --o-predictions "$ART/motus-age-trac-r3-test-predictions.qza"
