# Classification: Predicting Vegetation

This tutorial demonstrates how to predict vegetation presence using log-contrast classification. We'll show two approaches: simple log-contrast classification and log-contrast classification with taxonomic information (TRAC).

## Simple Log-Contrast Classification

### Overview

Log-contrast classification is useful for predicting categorical outcomes (like vegetation presence/absence) from compositional data. This approach uses regularized logistic regression with log-contrast penalties to handle the compositional nature of microbiome data.

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

Include environmental metadata with custom weights:

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

Create training and test sets for classification analysis:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates_lc.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column vegetation \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/classify-xtraining_lc.qza \
    --o-test-table data/classify-xtest_lc.qza \
    --o-training-targets data/classify-training-targets_lc.qza \
    --o-test-targets data/classify-test-targets_lc.qza
```

### Step 4: Train the Classification Model

Train a log-contrast classifier using cross-validation to predict vegetation presence:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_lc.qza \
    --i-c data/ccovariates_lc.qza \
    --i-weights data/wcovariates_lc.qza \
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
    --o-result data/classifytaxa_lc.qza
```

**Parameters explained:**
- `--i-features`: Training feature table
- `--i-c`: C matrix for log-contrast constraints
- `--i-weights`: Feature weights
- `--m-y-column vegetation`: Target variable (vegetation presence/absence)
- `--p-huber`: Use Huber loss for robustness
- `--p-stabsel`: Enable stability selection for feature selection
- `--p-cv`: Perform cross-validation
- `--p-stabsel-threshold 0.5`: Stability selection threshold

### Step 5: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/classify-xtest_lc.qza \
    --i-problem data/classifytaxa_lc.qza \
    --o-predictions data/classify-predictions_lc.qza
```

### Step 6: Generate Summary Visualization

Create a comprehensive summary of the classification results:

```bash
qiime classo summarize \
    --i-problem data/classifytaxa_lc.qza \
    --i-predictions data/classify-predictions_lc.qza \
    --o-visualization data/classifytaxa_C1_lc.qzv
```

---

## Log-Contrast Classification with TRAC (Taxonomic Information)

This approach incorporates taxonomic hierarchical information to identify which taxonomic groups are most predictive of vegetation.

### Step 1: Transform Features

Apply CLR transformation to the count data:

```bash
qiime classo transform-features \
    --p-transformation clr \
    --p-coef 0.5 \
    --i-features data/atacama-counts.qza \
    --o-x data/xclr
```

### Step 2: Add Taxonomic Information

Incorporate taxonomic classifications and compute adaptive weights:

```bash
qiime classo add-taxa \
    --i-features data/xclr.qza  \
    --i-taxa data/classification.qza \
    --o-x data/xtaxa \
    --o-aweights data/wtaxa
```

### Step 3: Add Covariates

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

### Step 4: Split Data

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

### Step 5: Train the Classification Model

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

### Step 6: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/classify-xtest_trac.qza \
    --i-problem data/classifytaxa_trac.qza \
    --o-predictions data/classify-predictions_trac.qza
```

### Step 7: Generate Summary Visualization

Create a comprehensive summary including taxonomic information:

```bash
qiime classo summarize \
    --i-problem data/classifytaxa_trac.qza \
    --i-taxa data/classification.qza \
    --i-predictions data/classify-predictions_trac.qza \
    --o-visualization data/classifytaxa_C1_trac.qzv
```

**Explanation:**
- Generates an interactive QIIME 2 visualization of the estimated network.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).

**Note:** The TRAC approach leverages the hierarchical structure of taxonomic information to improve feature selection and provides more interpretable results in terms of taxonomic relationships.