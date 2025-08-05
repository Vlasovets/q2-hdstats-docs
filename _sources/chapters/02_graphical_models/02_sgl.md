# Single Graphical Lasso

This tutorial demonstrates how to estimate a sparse inverse covariance matrix using the **Single Graphical Lasso** (SGL) method. This approach identifies conditional dependencies between features (e.g., microbial taxa) by solving an L1-penalized maximum likelihood problem, encouraging sparsity in the precision matrix.

**Purpose**: Estimates sparse direct associations between microbial taxa

**Key characteristics**:
- Assumes all edges have equal importance
- Uses uniform L1 penalty (λ₁) across all potential connections
- Identifies direct conditional dependencies in the network
- Best for: Initial network exploration and identifying core microbial interactions

**Interpretation**:
- Edges represent direct associations after removing indirect effects
- Edge thickness indicates association strength
- Network sparsity controlled by λ₁ parameter
- May include environment-mediated associations

## Step 1: Estimate a Sparse Model

We begin by estimating the inverse covariance (precision) matrix from a precomputed covariance matrix:

```bash
# sparse model
qiime gglasso solve-problem \
     --p-n-samples 50 \
     --p-lambda1-min 0.001 \
     --p-lambda1-max 1 \
     --p-n-lambda1 50 \
     --p-gamma 0.01 \
     --p-latent False \
     --i-covariance-matrix data/atacama-table-corr.qza \
     --o-solution data/atacama-solution-sgl.qza \
     --verbose
```

**Explanation:**

- `--p-n-samples 50`: Number of individuals used to compute the input covariance matrix.
- `--p-lambda1-min`: Lower-bound for the sparsity penalty (λ₁).
- `--p-lambda1-max`: Upper-bound for the sparsity penalty (λ₁).
- `--p-n-lambda1`: Number of grid points between the min and max lambda values.
- `--p-gamma 0.01`: Controls the model selection criterion (e.g., eBIC).
- `--p-latent False`: Indicates that this is a standard graphical lasso (not a latent variable model).
- `--i-covariance-matrix`: Input covariance matrix in QIIME 2 format.
- `--o-solution`: Output artifact containing the estimated sparse inverse covariance matrix.

## Step 2: Visualize the Estimated Network

```bash
# visualize the results
qiime gglasso summarize \
    --i-solution data/atacama-solution-sgl.qza \
    --p-label-size 25pt \
    --o-visualization data/sgl-summary.qzv
```

**Explanation:**

- Generates an interactive QIIME 2 visualization of the estimated network.
- `--p-label-size 25pt`: Sets the font size of node labels in the network plot.
- The output `.qzv` file can be viewed using [QIIME 2 View](https://view.qiime2.org/).
