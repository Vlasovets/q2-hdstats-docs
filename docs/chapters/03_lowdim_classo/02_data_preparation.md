# Data Preparation

This guide demonstrates how to prepare microbiome data for log-contrast analysis using q2-classo. We'll walk through transforming features, adding taxonomic information and covariates, and splitting data for machine learning tasks.

In

## Data Transformation

Transform count data using CLR transformation:

```bash
qiime classo transform-features \
     --p-transformation clr \
     --p-coef 0.5 \
     --i-features data/atacama-counts.qza \
     --o-x data/xclr
```

## Add Taxonomic Information

Incorporate taxonomic classifications and compute adaptive weights:

```bash
qiime classo add-taxa \
    --i-features data/xclr.qza  \
    --i-taxa data/classification.qza \
    --o-x data/xtaxa \
    --o-aweights data/wtaxa
```

## Add Covariates

Include environmental metadata with custom weights for each covariate:

```bash
qiime classo add-covariates \
    --i-features data/xtaxa.qza \
    --i-weights data/wtaxa.qza \
    --m-covariates-file data/atacama-selected-covariates-veg.tsv \
    --p-to-add ph average-soil-relative-humidity elevation average-soil-temperature vegetation \
    --p-w-to-add 1. 0.1 0.1 0.1 1 \
    --o-new-features data/xcovariates \
    --o-new-c data/ccovariates \
    --o-new-w data/wcovariates
```

## Split Data for Regression Analysis

Create training and test sets for continuous target prediction:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column average-soil-temperature \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/regress-xtraining \
    --o-test-table data/regress-xtest \
    --o-training-targets data/regress-training-targets.qza \
    --o-test-targets data/regress-test-targets.qza
```

## Split Data for Classification Analysis

Create training and test sets for categorical target prediction:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column vegetation \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/classify-xtraining \
    --o-test-table data/classify-xtest \
    --o-training-targets data/classify-training-targets.qza \
    --o-test-targets data/classify-test-targets.qza
```

Your data is now prepared for log-contrast modeling with both regression and classification targets.