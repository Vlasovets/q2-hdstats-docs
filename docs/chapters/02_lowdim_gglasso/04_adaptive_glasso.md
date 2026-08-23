# Adaptive Graphical Lasso

The adaptive graphical lasso appends the environmental covariates to the feature table
as additional variables and gives them their own penalty weights. Edges between taxa
then represent associations that survive once the measured environment is modelled
explicitly rather than left to confound them.

The weights carry prior knowledge and have to come from you: penalise the associations you
take to be environment-mediated more heavily, and the connections you take to be
biologically meaningful less. Nothing in the data chooses them.

## Adding the covariates to the table

Apply the modified centred log-ratio (mCLR) transform and standardise the environmental
covariates alongside the microbiome features:

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

- `--p-add-metadata True`: appends the covariates from the sample metadata to the feature table.
- `--p-transformation mclr`: applies the mCLR transform, which addresses compositionality.
- The output table holds the standardised microbial and environmental variables together.

## Covariance over features and covariates

Compute the scaled covariance of the combined table:

```bash
# calculate correlation including covariates
qiime gglasso calculate-covariance \
     --p-method scaled \
     --i-table data/atacama-table-mclr-meta.qza \
     --o-covariance-matrix data/atacama-table-corr-meta.qza
```

**Explanation:**

- `--p-method scaled`: produces a Pearson correlation matrix, the scaled covariance.
- The covariates sit in the input table and are treated like any other feature at this stage.

## Fitting the adaptive model

Give the selected covariates their own penalty weights:

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
- `--p-n-samples 50`: the number of samples the input covariance was computed from.
- `--p-lambda1-min`: lower bound of the sparsity penalty λ₁.
- `--p-lambda1-max`: upper bound of the sparsity penalty λ₁.
- `--p-n-lambda1`: number of grid points between the two bounds.
- `--p-gamma 0.01`: the extended BIC parameter.
- `--p-weights`: assigns an individual penalty weight to each named covariate.
- `--p-latent False`: fits the standard graphical lasso, with no low-rank component.
- `--i-covariance-matrix`: the input covariance, as a QIIME 2 artifact.
- `--o-solution`: the output artifact holding the estimated sparse precision matrix.

## Visualising the network with covariates

Draw the covariates into the network alongside the taxa:

```bash
# visualize the results
qiime gglasso summarize \
    --i-solution data/atacama-solution-adapt.qza \
    --p-label-size 25pt \
    --p-n-cov 4 \
    --o-visualization data/adapt-summary.qzv
```

**Explanation:**
- The action writes an interactive QIIME 2 visualization of the estimated network.
- `--p-label-size 25pt`: font size of the node labels in the network plot.
- `--p-n-cov 4`: the number of covariates in the table, so that they are treated as separate nodes.
- Open the resulting `.qzv` at [QIIME 2 View](https://view.qiime2.org/).
