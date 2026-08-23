# Regression: Predicting Continuous Outcomes

Log-contrast regression fits a regularized model to log-transformed compositional data and returns the taxa, or the taxonomic groups, whose balance tracks a continuous response: temperature, pH, or any other environmental or clinical measurement.

## Two workflows

### [Log-Contrast Regression](01_logcontrast.md)
- CLR transformation and covariates, no taxonomy
- Faster to fit
- Suited to exploratory work
- Choose it when taxonomic relationships are not part of the question

### [trac: Regression with taxonomic information](02_trac.md)
- Aggregates predictors along the taxonomic hierarchy through adaptive weights
- Attaches coefficients to named clades
- Groups feature selection phylogenetically
- Choose it for results you intend to publish

## Heterogeneous variance

The concomitant formulation estimates the noise scale jointly with the coefficients. Use it when the residual variance is not constant — see [Concomitant Formulation](../05_advanced/01_concomitant_formulation.md).

## Worked example

Both workflows predict average soil temperature from the Atacama desert microbiome dataset and report which taxa track the temperature gradient.

## Prerequisites

Work through [Data Preparation](../02_data_preparation.md) first.

## Reading order

1. Fit without taxonomy in [Log-Contrast Regression](01_logcontrast.md).
2. Refit with clade-aggregated predictors in [trac](02_trac.md).
3. Read the selected coefficients with the [Interpretation guide](../07_interpretation.md).
