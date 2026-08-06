# Graphical Lasso Models

This is the **reference tier**. Every q2-gglasso action gets its one canonical
demonstration here, on a dataset small enough that you can see the whole
covariance matrix at once: **50 samples × 13 ASVs** from the Atacama soil study,
plus five environmental covariates.

Nothing here is meant to be a scientific result. The point is that each command
runs in seconds, and you can check the output by eye before trusting the same
command on 300 or 30,000 features in
[Tier 2](../04_highdim_atacama/00_index.md) and
[Tier 3](../05_metagenomics/00_index.md).

## What each chapter adds

| Chapter | Action it owns | New parameters |
|---|---|---|
| [Data Preparation](01_data_preparation.md) | `transform-features`, `calculate-covariance` | `transformation`, `pseudo_count`, `keep_original_id`, `method`, `bias` |
| [Single Graphical Lasso](02_sgl.md) | `solve-problem` (sparse) | `n_samples`, `lambda1_min/max`, `n_lambda1`, `gamma` |
| [Sparse + Low-Rank](03_slr.md) | `solve-problem` (latent) | `latent`, `mu1_min/max`, `n_mu1` |
| [Adaptive Graphical Lasso](04_adaptive_glasso.md) | `solve-problem` (weighted) | `weights`, `add_metadata`, `scale_metadata`, `n_cov` |
| [Regularization Paths](05_lambda_paths.md) | model selection | `path_scale`, `lambda1_path`, `mu1_path`, `gamma` |
| [Multiple Graphical Lasso](06_multiple_graphical_lasso.md) | `build-groups`, multi-group `solve-problem` | `reg`, `lambda2_min/max`, `n_lambda2`, `non_conforming`, `group_array` |
| [Latent-Component PCA](07_pca.md) | `pca` | `n_components`, `color_by` |
| [Summarizing a Solution](08_summarize.md) | `summarize` | `width`, `height`, `label_size` |
| [Interpretation](09_interpretation.md) | — | comparing the models |

A machine-checkable version of this mapping lives in the
[Command Coverage Matrix](../90_reference/01_command_coverage.md).

## Before you start

These chapters assume `data/` is populated as described in
[Download the Tutorial Data](../00_getting_started/03_download_data.md), and that
`qiime gglasso --help` lists all six actions — see
[Verifying Your Installation](../01_installation/04_verify.md).

```{tip}
Several q2-gglasso actions have sharp edges that are easy to hit and hard to
diagnose — `pca` silently requires a latent solution, `transform-features`
demands a taxonomy it never reads, and `--p-rank` always raises. Each chapter
flags its own at the point you would hit it, and they are collected in
[Troubleshooting](../90_reference/04_troubleshooting.md).
```
