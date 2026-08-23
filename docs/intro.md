# High-Dimensional Statistics with QIIME 2

Two plugins extend the [QIIME 2](https://qiime2.org/) microbiome multi-omics platform
{cite}`bolyen2019reproducible` with methods for data where the number of features exceeds
the number of samples:

- [**q2-gglasso**](https://github.com/Vlasovets/q2-gglasso) estimates microbial
  association networks by sparse inverse covariance estimation, including the
  sparse-plus-low-rank and multiple-group formulations
  {cite}`Schaipp2021, friedman2008sparse`.
- [**q2-classo**](https://github.com/Vlasovets/q2-classo) fits sparse log-contrast
  regression and classification models, with taxonomic aggregation through trac
  {cite}`Simpson2021, combettes2021regression`.

<!-- ```{figure} images/png/overview.png
---
name: q2-overview
alt: High-dimensional statistics with QIIME2 overview
width: 600px
align: center
---
High-dimensional statistics with QIIME2.
``` -->

Both read and write QIIME 2 artifacts, so a network or a regression carries the same
provenance record as any other step in a QIIME 2 analysis.

## What the plugins do

- **Estimate association networks** under a penalty, where the sample covariance is
  singular and cannot be inverted directly.
- **Separate direct associations from shared latent drivers**, by decomposing the
  precision matrix into a sparse component and a low-rank one.
- **Fit interpretable regression and classification models** on compositional data, where
  only ratios between features are identified.
- **Aggregate features along a taxonomy**, so that coefficients name clades rather than
  individual sequence variants.

## How this book is organised

Each action is demonstrated once on a dataset small enough to check by eye, then applied
to a problem where the answer is not obvious:

| Part | Data | What it establishes |
|---|---|---|
| [Graphical Lasso Models](chapters/02_lowdim_gglasso/00_index.md) | 50 samples × 13 ASVs, Atacama soil | every q2-gglasso action, output verifiable by inspection |
| [Log-Contrast Models](chapters/03_lowdim_classo/00_overview.md) | synthetic, known ground truth | every q2-classo action, against a known answer |
| [Scale up](chapters/04_highdim_atacama/00_index.md) | 300 ASVs and 100 mOTUs species | model selection when $p > n$ and the penalty decides the model |
| [Reference](chapters/90_reference/01_command_coverage.md) | — | every parameter, and the failure modes worth knowing |

Start with [Prerequisites & Installation](chapters/01_installation/01_prerequisites.md),
or read [Why High-Dimensional Statistics?](chapters/00_getting_started/01_intro.md) for
the problem the plugins address.
