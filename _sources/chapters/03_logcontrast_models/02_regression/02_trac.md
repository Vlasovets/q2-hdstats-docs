# Log-Contrast Regression with TRAC (Taxonomic Information)

This tutorial demonstrates how to leverage taxonomic hierarchical information in log-contrast regression to predict average soil temperature. TRAC (Taxonomic Regularized Adaptive Contrasts) improves feature selection by incorporating the taxonomic structure.

## Overview

TRAC regression incorporates taxonomic classifications and computes adaptive weights based on the taxonomic hierarchy. This approach can identify which taxonomic groups (not just individual taxa) are most predictive of your target variable, leading to more interpretable and potentially more robust predictions.

## Step 1: Transform Features

Apply CLR transformation to the count data:

```bash
qiime classo transform-features \
    --p-transformation clr \
    --p-coef 0.5 \
    --i-features data/atacama-counts.qza \
    --o-x data/xclr
```

## Step 2: Add Taxonomic Information

Incorporate taxonomic classifications and compute adaptive weights based on taxonomy:

```bash
qiime classo add-taxa \
    --i-features data/xclr.qza  \
    --i-taxa data/classification.qza \
    --o-x data/xtaxa \
    --o-aweights data/wtaxa
```

This step:
- Uses taxonomic hierarchy to structure feature relationships
- Computes adaptive weights based on taxonomic distances
- Encourages selection of taxa that share taxonomic relatives

## Step 3: Add Covariates

Include environmental metadata with custom weights, combined with taxonomic features:

```bash
qiime classo add-covariates \
    --i-features data/xtaxa.qza \
    --i-weights data/wtaxa.qza \
    --m-covariates-file data/atacama-selected-covariates-veg.tsv \
    --p-to-add ph average-soil-relative-humidity elevation average-soil-temperature vegetation \
    --p-w-to-add 1. 0.1 0.1 0.1 1 \
    --o-new-features data/xcovariates_trac \
    --o-new-c data/ccovariates_trac \
    --o-new-w data/wcovariates_trac
```

## Step 4: Split Data

Create training and test sets for regression analysis:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates_trac.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column average-soil-temperature \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/regress-xtraining_trac.qza \
    --o-test-table data/regress-xtest_trac.qza \
    --o-training-targets data/regress-training-targets_trac.qza \
    --o-test-targets data/regress-test-targets_trac.qza
```

## Step 5: Train the Regression Model

Use log-contrast regression with taxonomic information and stability selection:

```bash
qiime classo regress \
    --i-features data/regress-xtraining_trac.qza \
    --i-c data/ccovariates_trac.qza \
    --i-weights data/wcovariates_trac.qza \
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
    --o-result data/regresstaxa_trac.qza
```

**Key parameters:**
- `--p-stabsel`: Enable stability selection for robust feature selection
- `--p-stabsel-threshold 0.5`: Features selected in >50% of subsamples
- `--p-cv`: Use cross-validation for model selection
- `--p-concomitant False`: Use standard formulation

## Step 6: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/regress-xtest_trac.qza \
    --i-problem data/regresstaxa_trac.qza \
    --o-predictions data/regress-predictions_trac.qza
```

## Step 7: Visualize Results

Generate a comprehensive summary including taxonomic information:

```bash
qiime classo summarize \
    --i-problem data/regresstaxa_trac.qza \
    --i-taxa data/classification.qza \
    --i-predictions data/regress-predictions_trac.qza \
    --o-visualization data/regresstaxa_R1_trac.qzv
```

**Output visualization:**
The `.qzv` file contains:
- Selected taxa with taxonomic context
- Model performance metrics (R², RMSE, etc.)
- Taxonomic group importances
- Cross-validation curves
- Prediction accuracy on test data

View the results at [QIIME 2 View](https://view.qiime2.org/).

## Comparing Simple Log-Contrast vs TRAC

**Simple Log-Contrast:**
- Treats each taxon independently
- Faster computation
- Useful for exploratory analysis

**TRAC:**
- Leverages taxonomic relationships
- More interpretable results
- Better feature selection for hierarchically structured data
- Slightly increased computation time
