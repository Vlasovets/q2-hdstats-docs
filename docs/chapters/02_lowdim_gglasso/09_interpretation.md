# Network Interpretation and Analysis

## Overview

This chapter explains how to interpret the results from previous chapters and different graphical lasso methods implemented in q2-gglasso. We'll compare three key approaches: Single Graphical Lasso (SGL), Sparse + Low-Rank (SLR), and Adaptive Graphical Lasso.

```{figure} ../../images/png/example_gglasso.png
:name: fig-gglasso-comparison
:width: 100%

The three estimators on the same 13-ASV data, **at a fixed $\lambda_1$**.
**(a)** the empirical covariance $\hat{S}_0$ and the six relationships visible in
it. **(b)** the sparse-only solution $\hat{\Theta}_{SP}$ keeps two:
$e_1 = (\text{ASV-6}, \text{ASV-11})$ and $e_2 = (\text{ASV-1}, \text{ASV-5})$.
**(d)** at the same penalty, the sparse + low-rank solution keeps only $e_2$ —
the low-rank block $L$ absorbs $e_1$. **(c)** the adaptive solution, fit on a
*different* covariance that includes the environmental covariates as nodes.
```

```{important}
**This is a fixed-$\lambda_1$ comparison, and that is the only comparison in
which "SLR gives fewer edges" is meaningful.** Panels (b) and (d) share a
penalty, so the single edge that disappears is attributable to the latent block.

Do not read it as "SLR always yields a sparser network". When each model picks
its own penalty by eBIC, the counts can go the other way — on this same dataset
SGL selects $\lambda_1 = 0.4292$ and keeps 2 edges while SLR selects
$\lambda_1 = 0.2442$ and keeps 4. See
[Adding latent variables](03_slr.md) for the measurement and why the two
comparisons answer different questions.

Panel (c) is not comparable to (b) or (d) at all: it is fit on
`atacama-table-corr-meta.qza`, which appends the covariates as additional
variables, so its node set differs — see
[Adaptive Graphical Lasso](04_adaptive_glasso.md).
```

```{note}
**Reading taxon labels.** Because `transform-features` now keeps the original
feature identifiers (see `--p-keep-original-id`), network nodes and heatmap axes
carry the ASV's real ID, so a selected feature can be traced directly back to its
taxonomy. In the [high-dimensional 300-ASV analysis](../04_highdim_atacama/02_model_selection.md), for
example, the first taxon selected by the log-contrast models is a
*Pseudarthrobacter* ASV — a genus emblematic of the Atacama soil biota, known for
desiccation- and oligotrophy-tolerant lifestyles, pigment production, and
survival in deep, hyperarid subsurface soils {cite}`finger2018pseudarthrobacter,horstmann2025subsurface,neilson2017significant`.
```

## Key Findings 

- **Environmental mediation**: Some apparent microbial correlations are mediated by environmental factors (elevation, pH, soil humidity, temperature)
- **Direct interactions**: Genuine microbial associations remain significant after controlling for environmental variables
- **Latent structure**: The low-rank component captures systematic variation potentially from unmeasured factors or global environmental gradients

### Interpretation guidelines:
- Compare SGL vs. SLR to distinguish direct from latent-mediated associations.
- Use adaptive results to identify environment-independent microbial interactions.
- Consider edge weights and stability across different λ values for robust associations.

### When to Use Each Method
| Method | Use When | Key Benefits |
|--------|----------|-------------|
| **SGL** | Exploratory analysis<br> No environmental data<br> Comprehensive associations<br> Speed priority | Fast, complete network view |
| **SLR** | Suspected confounders<br> Separate direct/indirect effects<br> Batch effects present<br> Core interactions focus | Isolates direct associations |
| **Adaptive model** | Environmental data available<br> Prior knowledge exists<br> Environment-independent focus<br> Hypothesis-driven analysis | Uses prior knowledge, controls confounders |

## Next Steps

Once network associations are identified, use **q2-classo** for regression and classification tasks:

1. **Feature selection**: Use identified microbial associations as candidate features for environmental or phenotype prediction
2. **Regression analysis**: Model continuous outcomes with sparse microbial predictors
3. **Classification tasks**: Predict binary outcome using selected microbial features
4. **Model selection**: Validate predictive performance using classo's built-in model selection methods

