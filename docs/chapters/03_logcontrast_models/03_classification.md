# Classification: Predicting Vegetation

This tutorial demonstrates how to predict vegetation presence using log-contrast classification.

## Overview

Log-contrast classification is useful for predicting categorical outcomes (like vegetation presence/absence) from compositional data. This approach uses regularized logistic regression with log-contrast penalties to handle the compositional nature of microbiome data.

## Step 1: Train the Classification Model

First, we train a log-contrast classifier using cross-validation to predict vegetation presence:

```bash
# Log-contrast classification with cross-validation
qiime classo classify \
    --i-features data/classify-xtraining.qza \
    --i-c data/ccovariates.qza \
    --i-weights data/wcovariates.qza \
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
    --o-result data/classifytaxa.qza
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

## Step 2: Make Predictions

Apply the trained model to test data:

```bash
qiime classo predict \
    --i-features data/classify-xtest.qza \
    --i-problem data/classifytaxa.qza \
    --o-predictions data/classify-predictions.qza
```

## Step 3: Generate Summary Visualization

Create a comprehensive summary of the classification results:

```bash
qiime classo summarize \
    --i-problem data/classifytaxa.qza \
    --i-taxa data/classification.qza \
    --i-predictions data/classify-predictions.qza \
    --o-visualization data/classifytaxa_C1.qzv
```

**Explanation:**
- Generates an interactive QIIME 2 visualization of the estimated network.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).