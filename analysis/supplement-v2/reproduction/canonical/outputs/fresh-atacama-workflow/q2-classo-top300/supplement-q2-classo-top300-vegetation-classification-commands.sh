#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial"
OUT="$ROOT/outputs/fresh-atacama-workflow/q2-classo-top300"
QIIME_ENV="$ROOT/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
MM="$ROOT/work/bin/micromamba"
export MPLCONFIGDIR=/private/tmp/mpl-cache
export XDG_CACHE_HOME=/private/tmp/mamba-cache-top300
export NUMBA_CACHE_DIR=/private/tmp/numba-cache

# The vegetation metadata use the same 53-sample train/test split as the
# top-15 q2-classo example and the top-300 elevation examples.
cp "$ROOT/outputs/fresh-atacama-workflow/q2-classo-top15/vegetation-train-53-metadata-q2.tsv" \
   "$OUT/vegetation-train-53-metadata-q2.tsv"
cp "$ROOT/outputs/fresh-atacama-workflow/q2-classo-top15/vegetation-test-53-metadata-q2.tsv" \
   "$OUT/vegetation-test-53-metadata-q2.tsv"

"$MM" run -p "$QIIME_ENV" qiime classo classify \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --m-y-file "$OUT/vegetation-train-53-metadata-q2.tsv" \
  --m-y-column vegetation \
  --p-no-path \
  --p-cv \
  --p-cv-one-se \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 100 \
  --p-cv-lamin 0.001 \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-no-huber \
  --p-intercept \
  --o-result "$OUT/atacama-top300-vegetation-classification-cv.qza"

"$MM" run -p "$QIIME_ENV" qiime classo predict \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-vegetation-classification-cv.qza" \
  --o-predictions "$OUT/atacama-top300-vegetation-classification-cv-train-predictions.qza"

"$MM" run -p "$QIIME_ENV" qiime classo predict \
  --i-features "$OUT/atacama-top300-test-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-vegetation-classification-cv.qza" \
  --o-predictions "$OUT/atacama-top300-vegetation-classification-cv-test-predictions.qza"

# Derived TSV summaries and figures were generated from the QIIME 2 artifacts by:
"$MM" run -p "$QIIME_ENV" python "$OUT/plot_top300_vegetation_classification.py"


# C1 stability selection for the same vegetation classification problem.
"$MM" run -p "$QIIME_ENV" qiime classo classify \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --m-y-file "$OUT/vegetation-train-53-metadata-q2.tsv" \
  --m-y-column vegetation \
  --p-no-path \
  --p-no-cv \
  --p-stabsel \
  --p-stabsel-seed 1 \
  --p-stabsel-b 100 \
  --p-stabsel-q 10 \
  --p-stabsel-threshold 0.4 \
  --p-stabsel-threshold-label 0.2 \
  --p-no-lamfixed \
  --p-no-huber \
  --p-intercept \
  --o-result "$OUT/atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04.qza"

"$MM" run -p "$QIIME_ENV" qiime classo predict \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04.qza" \
  --o-predictions "$OUT/atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-train-predictions.qza"

"$MM" run -p "$QIIME_ENV" qiime classo predict \
  --i-features "$OUT/atacama-top300-test-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04.qza" \
  --o-predictions "$OUT/atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-test-predictions.qza"
