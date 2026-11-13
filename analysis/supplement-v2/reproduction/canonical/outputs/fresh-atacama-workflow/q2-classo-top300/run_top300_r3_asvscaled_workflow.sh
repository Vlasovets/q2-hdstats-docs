#!/usr/bin/env bash
set -euo pipefail

# Top-300 q2-classo elevation workflow, revised for the tutorial.
# Uses R3/concomitant least-squares c-lasso, matching the top-15 regression example.
# Covariates are pre-imputed, centered on the training set, and scaled to the
# median training-set CLR ASV column norm; therefore q2-classo rescaling is off.
# Covariates are already scaled to the median CLR-ASV training column norm.
# They therefore use the same penalty weights as individual ASV coordinates (1.0).

OUT="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/outputs/fresh-atacama-workflow/q2-classo-top300"

python "$OUT/prepare_top300_inputs.py"

qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path "$OUT/top300-count-table-53-samples.biom" \
  --output-path "$OUT/atacama-top300-prevalence-53-table.qza"

# Prediction uses training-supported ASVs; the unchanged 300-feature artifact
# above remains the source for all-sample graphical models.
qiime feature-table filter-features \
  --i-table "$OUT/atacama-top300-prevalence-53-table.qza" \
  --m-metadata-file "$OUT/top300-training-supported-feature-ids.tsv" \
  --o-filtered-table "$OUT/atacama-top300-training-supported-53-table.qza"

qiime feature-table filter-samples \
  --i-table "$OUT/atacama-top300-training-supported-53-table.qza" \
  --m-metadata-file "$OUT/train-test-split-53-q2.tsv" \
  --p-where "[split]='train'" \
  --p-no-filter-empty-features \
  --o-filtered-table "$OUT/atacama-top300-train-53-table.qza"

qiime feature-table filter-samples \
  --i-table "$OUT/atacama-top300-training-supported-53-table.qza" \
  --m-metadata-file "$OUT/train-test-split-53-q2.tsv" \
  --p-where "[split]='test'" \
  --p-no-filter-empty-features \
  --o-filtered-table "$OUT/atacama-top300-test-53-table.qza"

qiime classo transform-features \
  --i-features "$OUT/atacama-top300-train-53-table.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "$OUT/atacama-top300-train-53-clr.qza"

qiime classo transform-features \
  --i-features "$OUT/atacama-top300-test-53-table.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "$OUT/atacama-top300-test-53-clr.qza"

python "$OUT/create_classo_constraint_weights.py"

qiime tools import \
  --type ConstraintMatrix \
  --input-path "$OUT/top300-zero-sum-constraint.zip" \
  --output-path "$OUT/top300-zero-sum-constraint.qza"

qiime tools import \
  --type Weights \
  --input-path "$OUT/top300-unit-weights.zip" \
  --output-path "$OUT/top300-unit-weights.qza"

qiime classo add-covariates \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --i-c "$OUT/top300-zero-sum-constraint.qza" \
  --i-weights "$OUT/top300-unit-weights.qza" \
  --m-covariates-file "$OUT/top300-elevation-augmented-covariates-normalized-from-train.tsv" \
  --p-to-add asvscaled-log10-extract-concen --p-to-add asvscaled-ph --p-to-add asvscaled-ec --p-to-add asvscaled-log10-amplicon-concentration --p-to-add asvscaled-toc --p-to-add asvscaled-depth \
  --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale \
  --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 \
  --o-new-features "$OUT/atacama-top300-train-53-clr-plus-asvscaled-covariates-equalweights.qza" \
  --o-new-c "$OUT/top300-plus-asvscaled-covariates-equalweights-constraint.qza" \
  --o-new-w "$OUT/top300-plus-asvscaled-covariates-equalweights-weights.qza"

qiime classo add-covariates \
  --i-features "$OUT/atacama-top300-test-53-clr.qza" \
  --i-c "$OUT/top300-zero-sum-constraint.qza" \
  --i-weights "$OUT/top300-unit-weights.qza" \
  --m-covariates-file "$OUT/top300-elevation-augmented-covariates-normalized-from-train.tsv" \
  --p-to-add asvscaled-log10-extract-concen --p-to-add asvscaled-ph --p-to-add asvscaled-ec --p-to-add asvscaled-log10-amplicon-concentration --p-to-add asvscaled-toc --p-to-add asvscaled-depth \
  --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale --p-no-rescale \
  --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 --p-w-to-add 1.0 \
  --o-new-features "$OUT/atacama-top300-test-53-clr-plus-asvscaled-covariates-equalweights.qza" \
  --o-new-c "$OUT/top300-plus-asvscaled-covariates-equalweights-test-constraint.qza" \
  --o-new-w "$OUT/top300-plus-asvscaled-covariates-equalweights-test-weights.qza"

qiime classo regress \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --m-y-file "$OUT/elevation-train-53-metadata-q2.tsv" \
  --m-y-column elevation \
  --p-no-path \
  --p-cv \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 120 \
  --p-cv-lamin 0.0001 \
  --p-cv-one-se \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-concomitant \
  --p-no-huber \
  --p-intercept \
  --o-result "$OUT/atacama-top300-elevation-asv-only-r3-cv.qza"

qiime classo regress \
  --i-features "$OUT/atacama-top300-train-53-clr-plus-asvscaled-covariates-equalweights.qza" \
  --i-c "$OUT/top300-plus-asvscaled-covariates-equalweights-constraint.qza" \
  --i-weights "$OUT/top300-plus-asvscaled-covariates-equalweights-weights.qza" \
  --m-y-file "$OUT/elevation-train-53-metadata-q2.tsv" \
  --m-y-column elevation \
  --p-no-path \
  --p-cv \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 120 \
  --p-cv-lamin 0.0001 \
  --p-cv-one-se \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-concomitant \
  --p-no-huber \
  --p-intercept \
  --o-result "$OUT/atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza"

qiime classo predict \
  --i-features "$OUT/atacama-top300-train-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-elevation-asv-only-r3-cv.qza" \
  --o-predictions "$OUT/atacama-top300-elevation-asv-only-r3-train-predictions.qza"

qiime classo predict \
  --i-features "$OUT/atacama-top300-test-53-clr.qza" \
  --i-problem "$OUT/atacama-top300-elevation-asv-only-r3-cv.qza" \
  --o-predictions "$OUT/atacama-top300-elevation-asv-only-r3-test-predictions.qza"

qiime classo predict \
  --i-features "$OUT/atacama-top300-train-53-clr-plus-asvscaled-covariates-equalweights.qza" \
  --i-problem "$OUT/atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza" \
  --o-predictions "$OUT/atacama-top300-elevation-covariate-augmented-equalweights-r3-train-predictions.qza"

qiime classo predict \
  --i-features "$OUT/atacama-top300-test-53-clr-plus-asvscaled-covariates-equalweights.qza" \
  --i-problem "$OUT/atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza" \
  --o-predictions "$OUT/atacama-top300-elevation-covariate-augmented-equalweights-r3-test-predictions.qza"
