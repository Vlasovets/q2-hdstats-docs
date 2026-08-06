# Adaptive Graphical Lasso

This tutorial demonstrates how to incorporate **environmental covariates** into sparse inverse covariance estimation using the **Adaptive Graphical Lasso** model. This method allows users to assign custom penalty weights to selected covariates, enabling more flexible and biologically informed regularization.

**Purpose**: Incorporates prior knowledge about environmental covariates

**Key characteristics**:
- Applies differential penalties to different types of associations
- Lower penalties for biologically meaningful connections
- Higher penalties for environment-mediated associations
- Integrates environmental metadata directly into the model

**Interpretation**:
- Edges represent environment-independent microbial interactions
- Environmental covariates explicitly modeled rather than confounded
- Most biologically interpretable results
- Requires domain knowledge for penalty weight selection

## Step 1: Add Metadata

First, apply a modified centered log-ratio (mCLR) transformation and standardize environmental covariates alongside microbiome features:

```bash
# standardise covariates 
qiime gglasso transform-features \
     --p-transformation mclr \
     --p-add-metadata True \
     --i-table data/atacama-counts.qza \
     --i-taxonomy data/classification.qza \
     --m-sample-metadata-file data/selected-atacama-sample-metadata.tsv  \
     --o-transformed-table data/atacama-table-mclr-meta.qza
```

**Explanation:**

- `--p-add-metadata True`: Appends covariates (e.g., environmental features) to the feature table.
- `--p-transformation mclr`: Applies mCLR transformation to address compositionality.
- The output table contains both standardized microbial and environmental data.

## Step 2: Calculate new Covariance Matrix

Next, compute the scaled covariance matrix for the transformed table:

```bash
# calculate correlation including covariates
qiime gglasso calculate-covariance \
     --p-method scaled \
     --i-table data/atacama-table-mclr-meta.qza \
     --o-covariance-matrix data/atacama-table-corr-meta.qza
```

**Explanation:**

- `--p-method scaled`: Produces a Pearson correlation matrix (i.e., scaled covariance).
- Covariates are included in the input table and will be treated equally with other features.

## Step 3: Fit the Adaptive Graphical Lasso Model

Use adaptive penalty weights for selected covariates to guide the network estimation:

```bash
# sparse model with specific weights for covariates
qiime gglasso solve-problem \
     --p-n-samples 50 \
     --p-lambda1-min 0.001 \
     --p-lambda1-max 1 \
     --p-n-lambda1 2 \
     --p-gamma 0.01 \
     --p-latent False \
     --p-weights elevation 0.01 ph 0.01 average-soil-relative-humidity 0.01 average-soil-temperature 0.01  \
     --i-covariance-matrix data/atacama-table-corr-meta.qza \
     --o-solution data/atacama-solution-adapt.qza \
     --verbose
```

**Explanation:**
- `--p-n-samples 50`: Number of individuals used to compute the input covariance matrix.
- `--p-lambda1-min`: Lower-bound for the sparsity penalty (λ₁).
- `--p-lambda1-max`: Upper-bound for the sparsity penalty (λ₁).
- `--p-n-lambda1`: Number of grid points between the min and max lambda values.
- `--p-gamma 0.01`: Controls the model selection criterion (e.g., eBIC).
- `--p-weights`: Assigns individual penalty weights to selected covariates.
- `--p-latent False`: Indicates that this is a standard graphical lasso w/o low-rank.
- `--i-covariance-matrix`: Input covariance matrix in QIIME 2 format.
- `--o-solution`: Output artifact containing the estimated sparse inverse covariance matrix.

## Step 4: Visualize the Estimated Network

Generate a visualization that includes covariates in the network:

```bash
# visualize the results
qiime gglasso summarize \
    --i-solution data/atacama-solution-adapt.qza \
    --p-label-size 25pt \
    --p-n-cov 4 \
    --o-visualization data/adapt-summary.qzv
```

**Explanation:**
- Generates an interactive QIIME 2 visualization of the estimated network.
- `--p-label-size 25pt`: Sets the font size of node labels in the network plot.
- `--p-n-cov 4`: Indicates the number of covariates included, so they are treated as separate nodes.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).
