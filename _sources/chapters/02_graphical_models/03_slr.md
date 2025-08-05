# Single Graphical Lasso with Latent Variables

This tutorial demonstrates how to estimate a sparse inverse covariance matrix **while accounting for latent variables** using the **Sparse + Low-Rank (SLR)** Graphical Lasso model. This approach decomposes the inverse covariance matrix into a sparse component (direct dependencies) and a low-rank component (unobserved latent effects).

**Purpose**: Separates direct sparse interactions from systematic low-rank structure

**Key characteristics**:
- Decomposes precision matrix into sparse + low-rank components
- Sparse component: direct microbial interactions
- Low-rank component: systematic variation (environmental gradients, batch effects)
- Uses both L1 penalty (λ₁) and nuclear norm penalty (μ₁)

**Interpretation**:
- **Sparse component**: True direct microbial associations
- **Low-rank component**: Global patterns affecting multiple taxa simultaneously
- Better separation of biological signal from confounding factors
- More conservative edge detection compared to SGL

## Step 1: Estimate a Sparse + Low-Rank Model

We begin by fitting a model that includes both sparse and low-rank components to account for hidden (latent) structure:

```bash
# slr model
qiime gglasso solve-problem \
    --p-n-samples 50 \
    --p-lambda1-min 0.001 \
    --p-lambda1-max 1 \
    --p-mu1-min 0.50118723 \
    --p-mu1-max 0.79432823 \
    --p-n-lambda1 50 \
    --p-n-mu1 50 \
    --p-gamma 0.01 \
    --p-latent True \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --o-solution data/atacama-solution-slr.qza \
    --verbose
```

**Explanation:**

- `--p-n-samples 50`: Number of individuals used to compute the input covariance matrix.
- `--p-lambda1-min`: Lower-bound for the sparsity penalty (λ₁).
- `--p-lambda1-max`: Upper-bound for the sparsity penalty (λ₁).
- `--p-mu1-min`: Lower-bound for the low-rank penalty (μ₁).
- `--p-mu1-max`: Upper-bound for the low-rank penalty (μ₁).
- `--p-n-lambda1`: Number of grid points between the min and max lambda values.
- `--p-n-mu1`: Number of grid points for each penalty parameter.
- `--p-gamma 0.01`: Controls the model selection criterion (e.g., eBIC).
- `--p-latent True`: Enables the latent variable modeling (SLR).
- `--i-covariance-matrix`: Input covariance matrix in QIIME 2 format.
- `--o-solution`: Output artifact containing the estimated sparse + low-rank decomposition of the inverse covariance matrix.

## Step 2: Visualize the Estimated Network

The estimated network, including effects of latent variables, can be visualized using:

```bash
# visualize the results
qiime gglasso summarize \
    --i-solution data/atacama-solution-slr.qza \
    --p-label-size 25pt \
    --o-visualization data/slr-summary.qzv
```

**Explanation:**

- Generates an interactive QIIME 2 visualization of the inferred network structure.
- `--p-label-size 25pt`: Sets the font size of node labels in the network plot.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).
