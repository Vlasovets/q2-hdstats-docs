# Classification: Predicting Categorical Outcomes

This section demonstrates how to use log-contrast classification to predict categorical outcomes (like disease status, habitat type, or vegetation presence) from microbial community composition.

## Overview

Log-contrast classification applies regularized logistic regression to log-transformed compositional data, identifying which microbial taxa or taxonomic groups are most predictive of categorical variables.

## Workflows

We present two approaches:

### [Log-Contrast Classification](01_simple_logcontrast.md)
- Uses CLR transformation and basic covariates
- Faster computation
- Good for exploratory analysis
- Best when taxonomic relationships are not critical

### [TRAC: Classification with Taxonomic Information](02_trac.md)
- Incorporates taxonomic hierarchy via adaptive weights
- More interpretable results with taxonomic context
- Better feature selection through phylogenetic grouping
- Recommended for publication-quality analyses

## Advanced Options

After mastering the basic workflows, explore:
- [Concomitant Formulation](../04_advanced/01_concomitant_formulation.md): For data with heterogeneous variance (applicable to classification with Huber loss)

## Example Use Case

In this tutorial, we predict **vegetation presence/absence** from the Atacama desert microbiome dataset. The models identify which taxa show strong associations with vegetation status.

## Prerequisites

Complete [Data Preparation](../01_data_preparation.md) before starting these tutorials.

## Next Steps

1. Start with [Log-Contrast Classification](01_simple_logcontrast.md)
2. Compare with [TRAC approach](02_trac.md)
3. Interpret your results using the [Interpretation guide](../05_interpretation.md)
