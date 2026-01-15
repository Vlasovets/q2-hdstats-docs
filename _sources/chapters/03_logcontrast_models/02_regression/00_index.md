# Regression: Predicting Continuous Outcomes

This section demonstrates how to use log-contrast regression to predict continuous outcomes (like temperature, pH, or other environmental variables) from microbial community composition.

## Overview

Log-contrast regression applies regularized regression models to log-transformed compositional data, identifying which microbial taxa or taxonomic groups are most predictive of continuous environmental or clinical variables.

## Workflows

We present two approaches:

### [Simple Log-Contrast Regression](01_simple_logcontrast.md)
- Uses CLR transformation and basic covariates
- Faster computation
- Good for exploratory analysis
- Best when taxonomic relationships are not critical

### [TRAC: Regression with Taxonomic Information](02_trac.md)
- Incorporates taxonomic hierarchy via adaptive weights
- More interpretable results with taxonomic context
- Better feature selection through phylogenetic grouping
- Recommended for publication-quality analyses

## Advanced Options

After mastering the basic workflows, explore:
- [Concomitant Formulation](../04_advanced/01_concomitant_formulation.md): For data with heterogeneous variance

## Example Use Case

In this tutorial, we predict **average soil temperature** from the Atacama desert microbiome dataset. The models identify which taxa show strong associations with temperature gradients.

## Prerequisites

Complete [Data Preparation](../01_data_preparation.md) before starting these tutorials.

## Next Steps

1. Start with [Simple Log-Contrast Regression](01_simple_logcontrast.md)
2. Compare with [TRAC approach](02_trac.md)
3. Interpret your results using the [Interpretation guide](../04_interpretation.md)
