# Log-Contrast Classification with TRAC (Taxonomic Information)

This tutorial demonstrates how to leverage taxonomic hierarchical information in log-contrast classification to predict vegetation presence. TRAC (Taxonomic Regularized Adaptive Contrasts) improves feature selection by incorporating the taxonomic structure.

## Overview

TRAC classification incorporates taxonomic classifications and computes adaptive weights based on the taxonomic hierarchy. This approach can identify which taxonomic groups (not just individual taxa) are most predictive of your target variable, leading to more interpretable and potentially more accurate classifications.

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

Incorporate taxonomic classifications and compute adaptive weights:

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

Create training and test sets for classification analysis:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates_trac.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column vegetation \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/classify-xtraining_trac.qza \
    --o-test-table data/classify-xtest_trac.qza \
    --o-training-targets data/classify-training-targets_trac.qza \
    --o-test-targets data/classify-test-targets_trac.qza
```

## Step 5: Train the Classification Model

Train a log-contrast classifier with taxonomic information:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_trac.qza \
    --i-c data/ccovariates_trac.qza \
    --i-weights data/wcovariates_trac.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column vegetation \
    --p-huber False \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 42 \
    --p-no-cv-one-se \
    --o-result data/classifytaxa_trac.qza
```

**Parameters explained:**
- `--i-features`: Training feature table
- `--i-c`: C matrix for log-contrast constraints
- `--i-weights`: Feature weights (including taxonomic weights)
- `--m-y-column vegetation`: Target variable
- `--p-huber False`: Use standard logistic loss
- `--p-stabsel`: Enable stability selection for feature selection
- `--p-cv`: Perform cross-validation
- `--p-stabsel-threshold 0.5`: Stability selection threshold

## Step 6: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/classify-xtest_trac.qza \
    --i-problem data/classifytaxa_trac.qza \
    --o-predictions data/classify-predictions_trac.qza
```

## Step 7: Generate Summary Visualization

Create a comprehensive summary including taxonomic information:

```bash
qiime classo summarize \
    --i-problem data/classifytaxa_trac.qza \
    --i-taxa data/classification.qza \
    --i-predictions data/classify-predictions_trac.qza \
    --o-visualization data/classifytaxa_C1_trac.qzv
```

**Output visualization:**
The `.qzv` file contains:
- Selected taxa with taxonomic context
- Model performance metrics (accuracy, precision, recall, F1-score)
- Taxonomic group importances
- Confusion matrix for test predictions
- Cross-validation curves

View the results at [QIIME 2 View](https://view.qiime2.org/).

## Comparing Log-Contrast vs TRAC

**Log-Contrast:**
- Treats each taxon independently
- Faster computation
- Useful for exploratory analysis

**TRAC:**
- Leverages taxonomic relationships
- More interpretable results
- Better feature selection for hierarchically structured data
- Slightly increased computation time
