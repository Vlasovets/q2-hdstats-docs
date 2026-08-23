# Graphical Lasso Models

Every q2-gglasso action gets one canonical demonstration on the same dataset: 50 samples
× 13 ASVs from the Atacama soil study, plus five environmental covariates. The tables
are small enough that you can read a whole covariance matrix at once, and each command
finishes in seconds, so you can check its output by eye before trusting the same command
on 300 features in [the high-dimensional Atacama
chapters](../04_highdim_atacama/00_index.md). None of these fits is meant to be a
scientific result.

## Actions by chapter

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

## Prerequisites

Populate `data/` as described in
[Download the Tutorial Data](../00_getting_started/03_download_data.md), and
confirm that `qiime gglasso --help` lists all six actions — see
[Verifying Your Installation](../01_installation/04_verify.md).

```{tip}
Several q2-gglasso actions fail in ways that are easy to hit and hard to
diagnose: `pca` requires a latent solution without saying so,
`transform-features` demands a taxonomy it never reads, and `--p-rank` always
raises. [Troubleshooting](../90_reference/04_troubleshooting.md) collects them.
```
