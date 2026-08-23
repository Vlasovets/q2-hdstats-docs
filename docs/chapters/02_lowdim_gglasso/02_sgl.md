# Single Graphical Lasso

The single graphical lasso (SGL) estimates a sparse inverse covariance matrix from a
precomputed covariance by solving an L1-penalized maximum likelihood problem. The penalty
encourages sparsity in the precision matrix. Each non-zero entry is a conditional
dependence between two features — a direct association between two taxa, with the
indirect paths through the remaining taxa removed — and its magnitude gives the
strength of that association.

The penalty is a single uniform L1 weight, λ₁, applied to every candidate edge, so all
pairs are treated alike. Raising λ₁ removes edges. The environmental covariates are absent
from this table, so an edge may also record two taxa responding to the same gradient
rather than interacting.

Fit SGL first on a new table: it gives an initial view of the network and the core
interactions among taxa, before you add weights or a latent block.

## Fitting the model

Estimate the precision matrix from the correlation matrix you computed earlier:

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

- `--p-n-samples 50`: the number of samples the input covariance was computed from.
- `--p-lambda1-min`: lower bound of the sparsity penalty λ₁.
- `--p-lambda1-max`: upper bound of the sparsity penalty λ₁.
- `--p-n-lambda1`: number of grid points between the two bounds.
- `--p-gamma 0.01`: the extended BIC parameter.
- `--p-latent False`: fits the standard graphical lasso, with no low-rank component.
- `--i-covariance-matrix`: the input covariance, as a QIIME 2 artifact.
- `--o-solution`: the output artifact holding the estimated sparse precision matrix.

## Visualising the network

```bash
# visualize the results
qiime gglasso summarize \
    --i-solution data/atacama-solution-sgl.qza \
    --p-label-size 25pt \
    --o-visualization data/sgl-summary.qzv
```

**Explanation:**

- The action writes an interactive QIIME 2 visualization of the estimated network.
- `--p-label-size 25pt`: font size of the node labels in the network plot.
- Open the resulting `.qzv` at [QIIME 2 View](https://view.qiime2.org/).
