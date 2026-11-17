# Fresh q2-classo top-15 Atacama tutorial workflow
# Input files are generated from the full Atacama feature table and metadata.
# YUN1005.1.1 is excluded first; the top-15 ASVs are then selected by the shared
# prevalence ranking on the remaining 53 analysis samples.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT="$ROOT/outputs/fresh-atacama-workflow/q2-classo-top15"

# Regenerate both canonical prevalence-filtered tables from the latest upstream
# Atacama table before importing this top-15 table into QIIME 2.
python "$ROOT/outputs/fresh-atacama-workflow/q2-classo-top300/prepare_top300_inputs.py"

# Import the prevalence-filtered top-15 count table as a QIIME 2 feature table.
# This BIOM table is generated immediately above from the full source table.
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ${OUT}/top15-count-table-53-samples.biom \
  --output-path ${OUT}/atacama-top15-prevalence-53-table.qza

# Split samples into training and held-out test sets.
qiime feature-table filter-samples \
  --i-table ${OUT}/atacama-top15-prevalence-53-table.qza \
  --m-metadata-file ${OUT}/train-test-split-53-q2.tsv \
  --p-where "[split]='train'" \
  --p-no-filter-empty-features \
  --o-filtered-table ${OUT}/atacama-top15-train-53-table.qza

qiime feature-table filter-samples \
  --i-table ${OUT}/atacama-top15-prevalence-53-table.qza \
  --m-metadata-file ${OUT}/train-test-split-53-q2.tsv \
  --p-where "[split]='test'" \
  --p-no-filter-empty-features \
  --o-filtered-table ${OUT}/atacama-top15-test-53-table.qza

# CLR transformation with q2-classo zero replacement coef=0.5.
qiime classo transform-features \
  --i-features ${OUT}/atacama-top15-train-53-table.qza \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x ${OUT}/atacama-top15-train-53-clr.qza

qiime classo transform-features \
  --i-features ${OUT}/atacama-top15-test-53-table.qza \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x ${OUT}/atacama-top15-test-53-clr.qza

# R1: non-concomitant least-squares c-lasso with 5-fold CV and one-SE selection.
qiime classo regress \
  --i-features ${OUT}/atacama-top15-train-53-clr.qza \
  --m-y-file ${OUT}/elevation-train-53-metadata-q2.tsv \
  --m-y-column elevation \
  --p-no-path \
  --p-cv \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 100 \
  --p-cv-lamin 0.001 \
  --p-cv-one-se \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-no-concomitant \
  --p-no-huber \
  --p-intercept \
  --o-result ${OUT}/atacama-top15-elevation-r1-cv.qza

# R3: concomitant least-squares c-lasso with 5-fold CV and one-SE selection.
qiime classo regress \
  --i-features ${OUT}/atacama-top15-train-53-clr.qza \
  --m-y-file ${OUT}/elevation-train-53-metadata-q2.tsv \
  --m-y-column elevation \
  --p-no-path \
  --p-cv \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 100 \
  --p-cv-lamin 0.001 \
  --p-cv-one-se \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-concomitant \
  --p-no-huber \
  --p-intercept \
  --o-result ${OUT}/atacama-top15-elevation-r3-cv.qza

# Held-out prediction for CV-selected models.
qiime classo predict \
  --i-features ${OUT}/atacama-top15-train-53-clr.qza \
  --i-problem ${OUT}/atacama-top15-elevation-r3-cv.qza \
  --o-predictions ${OUT}/atacama-top15-elevation-r3-cv-train-predictions.qza

qiime classo predict \
  --i-features ${OUT}/atacama-top15-test-53-clr.qza \
  --i-problem ${OUT}/atacama-top15-elevation-r1-cv.qza \
  --o-predictions ${OUT}/atacama-top15-elevation-r1-test-predictions.qza

qiime classo predict \
  --i-features ${OUT}/atacama-top15-test-53-clr.qza \
  --i-problem ${OUT}/atacama-top15-elevation-r3-cv.qza \
  --o-predictions ${OUT}/atacama-top15-elevation-r3-test-predictions.qza

# Native q2-classo visual summaries, including CV and prediction displays.
qiime classo summarize \
  --i-problem ${OUT}/atacama-top15-elevation-r1-cv.qza \
  --i-predictions ${OUT}/atacama-top15-elevation-r1-test-predictions.qza \
  --p-maxplot 50 \
  --o-visualization ${OUT}/atacama-top15-elevation-r1-cv-and-test-summary.qzv

qiime classo summarize \
  --i-problem ${OUT}/atacama-top15-elevation-r3-cv.qza \
  --i-predictions ${OUT}/atacama-top15-elevation-r3-test-predictions.qza \
  --p-maxplot 50 \
  --o-visualization ${OUT}/atacama-top15-elevation-r3-cv-and-test-summary.qzv

# R3 stability selection. With the PR60 c-lasso backend, threshold 0.4 retains
# the same two-ASV contrast selected by CV; q and B remain unchanged.
qiime classo regress \
  --i-features ${OUT}/atacama-top15-train-53-clr.qza \
  --m-y-file ${OUT}/elevation-train-53-metadata-q2.tsv \
  --m-y-column elevation \
  --p-no-path \
  --p-no-cv \
  --p-stabsel \
  --p-stabsel-seed 1 \
  --p-stabsel-b 100 \
  --p-stabsel-q 4 \
  --p-stabsel-percent-ns 0.5 \
  --p-stabsel-threshold 0.4 \
  --p-stabsel-threshold-label 0.3 \
  --p-no-lamfixed \
  --p-concomitant \
  --p-no-huber \
  --p-intercept \
  --o-result ${OUT}/atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza


# In-sample prediction for the R3 stability-selection refit.
qiime classo predict \
  --i-features ${OUT}/atacama-top15-train-53-clr.qza \
  --i-problem ${OUT}/atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza \
  --o-predictions ${OUT}/atacama-top15-elevation-r3-stabsel-threshold04-train-predictions.qza

qiime classo predict \
  --i-features ${OUT}/atacama-top15-test-53-clr.qza \
  --i-problem ${OUT}/atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza \
  --o-predictions ${OUT}/atacama-top15-elevation-r3-stabsel-threshold04-test-predictions.qza

qiime classo summarize \
  --i-problem ${OUT}/atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza \
  --i-predictions ${OUT}/atacama-top15-elevation-r3-stabsel-threshold04-test-predictions.qza \
  --p-maxplot 50 \
  --o-visualization ${OUT}/atacama-top15-elevation-r3-stabsel-threshold04-summary.qzv
