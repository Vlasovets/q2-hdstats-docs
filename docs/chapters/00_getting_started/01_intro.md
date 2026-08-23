# Getting Started

## Introduction

A high-throughput sequencing experiment typically identifies a large number of microbial
features, also called taxa. High-dimensional statistics supplies estimators built for
that regime.
They assess interactions between taxa and reveal patterns that may be indicative of
biological phenomena — shifts in community composition, associations with
environmental factors.

## The two plugins

QIIME2 ships a collection of plugins that preprocess raw sequences into feature tables
ready for downstream analysis. For advanced statistical analysis of those tables, you
often have no out-of-the-box option.

```{figure} ../../images/png/overview.png
---
name: q2-overview
alt: High-dimensional statistics with QIIME2 overview
width: 600px
align: center
---
High-dimensional statistics with QIIME2.
```

The pipeline above is implemented by two QIIME2 plugins: **q2-gglasso** for network
learning and **q2-classo** for regression and classification. Together they close gaps in
QIIME2, adding:

- **Network inference** through sparse inverse covariance estimation, which identifies
  microbial associations, for example which taxa co-occur within a community
- **Sparse log-contrast models** for classification and regression, which predict a
  measured covariate, such as disease status or an environmental factor, from the counts
- **Compositionally-aware estimators** that account for the constraints and the sparsity
  of microbiome counts
- **Interactive visualizations** with publication-ready network graphs and model
  performance plots
- **QIIME2 integration** using standard artifact formats and command-line interfaces

## Scope

These chapters cover:

1. Install and set up the q2-gglasso and q2-classo plugins
2. Apply graphical lasso methods for network inference
3. Use log-contrast models for regression and classification
4. Interpret results and create publication-ready visualizations
5. Integrate these methods into reproducible QIIME2 workflows
