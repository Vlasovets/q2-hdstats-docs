# Regression: Predicting Temperature

This tutorial demonstrates how to use log-contrast regression to predict average soil temperature from microbial community data. We'll show two approaches: simple log-contrast regression and log-contrast regression with taxonomic information (TRAC).

## Simple Log-Contrast Regression

### Step 1: Transform Features

Apply CLR transformation to the count data:

```bash
qiime classo transform-features \
    --p-transformation clr \
    --p-coef 0.5 \
    --i-features data/atacama-counts.qza \
    --o-x data/xclr
```

### Step 2: Add Covariates

Include environmental metadata with custom weights for each covariate:

```bash
qiime classo add-covariates \
    --i-features data/xclr.qza \
    --m-covariates-file data/atacama-selected-covariates-veg.tsv \
    --p-to-add ph average-soil-relative-humidity elevation average-soil-temperature vegetation \
    --p-w-to-add 1. 0.1 0.1 0.1 1 \
    --o-new-features data/xcovariates_lc \
    --o-new-c data/ccovariates_lc \
    --o-new-w data/wcovariates_lc
```

### Step 3: Split Data

Create training and test sets for regression analysis:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates_lc.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column average-soil-temperature \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/regress-xtraining_lc.qza \
    --o-test-table data/regress-xtest_lc.qza \
    --o-training-targets data/regress-training-targets_lc.qza \
    --o-test-targets data/regress-test-targets_lc.qza
```

### Step 4: Train the Regression Model

Use log-contrast regression with stability selection to identify the most stable predictive features:

```bash
qiime classo regress \
    --i-features data/regress-xtraining_lc.qza \
    --i-c data/ccovariates_lc.qza \
    --i-weights data/wcovariates_lc.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column average-soil-temperature \
    --p-concomitant False \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 1 \
    --p-no-cv-one-se \
    --o-result data/regresstaxa_lc.qza
```

Key parameters:
- `--p-stabsel`: Enable stability selection for robust feature selection
- `--p-stabsel-threshold 0.5`: Features selected in >50% of subsamples
- `--p-cv`: Use cross-validation for model selection
- `--p-concomitant`: Include concomitant variables in the model

### Step 5: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/regress-xtest_lc.qza \
    --i-problem data/regresstaxa_lc.qza \
    --o-predictions data/regress-predictions_lc.qza
```

### Step 6: Visualize Results

Generate a comprehensive summary of the regression results:

```bash
qiime classo summarize \
    --i-problem data/regresstaxa_lc.qza \
    --i-predictions data/regress-predictions_lc.qza \
    --o-visualization data/regresstaxa_R1_lc.qzv
```

The visualization will show selected taxa, model performance metrics, and prediction accuracy.

**Explanation:**
- Generates an interactive QIIME 2 visualization of the estimated network.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).

**Note:** The TRAC approach leverages the hierarchical structure of taxonomic information to improve feature selection and provides more interpretable results in terms of taxonomic relationships.