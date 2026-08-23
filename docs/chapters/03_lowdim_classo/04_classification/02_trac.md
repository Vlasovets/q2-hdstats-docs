# Tree-Aggregated Classification (trac)

trac (tree-aggregation of compositional data) {cite}`bien2021tree` places coefficients on nodes of the taxonomy, not only on individual ASVs. Weights derived from the hierarchy make an internal node — a genus, an order, a phylum — competitive with its members, so the fitted classifier can name the clades whose aggregated abundance tracks the split between vegetated and bare sites. A named clade is easier to interpret than the ASVs beneath it, and may classify more accurately.

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

Attach the taxonomy to the design and derive an adaptive weight for every node:

```bash
qiime classo add-taxa \
    --i-features data/xclr.qza  \
    --i-taxa data/classification.qza \
    --o-x data/xtaxa \
    --o-aweights data/wtaxa
```

The action uses the taxonomic hierarchy to structure the relationships between features and derives each adaptive weight from the taxonomic distances, which encourages the fit to select taxa that share taxonomic relatives.

## Step 3: Add Covariates

Include environmental metadata with custom weights, alongside the taxonomic features:

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

Fit the classifier on the aggregated design:

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
- `--i-c`: C matrix carrying the log-contrast constraints
- `--i-weights`: Feature weights, taxonomic weights included
- `--m-y-column vegetation`: Target variable
- `--p-huber False`: Use standard logistic loss
- `--p-stabsel`: Run stability selection alongside the path fit
- `--p-cv`: Choose the penalty by cross-validation
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

Render the fit, the taxonomy and the predictions as a single report:

```bash
qiime classo summarize \
    --i-problem data/classifytaxa_trac.qza \
    --i-taxa data/classification.qza \
    --i-predictions data/classify-predictions_trac.qza \
    --o-visualization data/classifytaxa_C1_trac.qzv
```

The `.qzv` file contains:
- the selected taxa with their taxonomic context
- model performance metrics (accuracy, precision, recall, F1-score)
- taxonomic group importances
- a confusion matrix for the test predictions
- cross-validation curves

View the results at [QIIME 2 View](https://view.qiime2.org/).

## Log-contrast and trac side by side

**Log-contrast** treats each taxon as an independent predictor, fits faster, and suits exploratory work.

**trac** uses the taxonomic relationships, reports coefficients on named clades, selects features better when the structure is hierarchical, and takes somewhat longer to fit.
