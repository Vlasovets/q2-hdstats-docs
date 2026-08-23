# Classification: Predicting Categorical Outcomes

Log-contrast classification fits a regularized logistic model to log-transformed compositional data and returns the taxa, or the taxonomic groups, whose balance separates two categories: disease status, habitat type, or the presence of vegetation.

## Two workflows

### [Log-Contrast Classification](01_logcontrast.md)
- CLR transformation and covariates, no taxonomy
- Faster to fit
- Suited to exploratory work
- Choose it when taxonomic relationships are not part of the question

### [trac: Classification with taxonomic information](02_trac.md)
- Aggregates predictors along the taxonomic hierarchy through adaptive weights
- Attaches coefficients to named clades
- Groups feature selection phylogenetically
- Choose it for results you intend to publish

## Heterogeneous variance and outliers

The concomitant formulation estimates the noise scale jointly with the coefficients, which suits data whose residual variance is not constant. The classification action does not accept it. [Concomitant Formulation](../05_advanced/01_concomitant_formulation.md) explains that restriction and gives the Huber hinge loss as the robust option for a binary outcome.

## Worked example

Both workflows predict vegetation presence and absence from the Atacama desert microbiome dataset and report which taxa track vegetation status.

## Prerequisites

Work through [Data Preparation](../02_data_preparation.md) first.

## Reading order

1. Fit without taxonomy in [Log-Contrast Classification](01_logcontrast.md).
2. Refit with clade-aggregated predictors in [trac](02_trac.md).
3. Read the selected coefficients with the [Interpretation guide](../07_interpretation.md).
