# Log-Contrast Models Overview

Log-contrast models are powerful tools for analyzing compositional microbiome data while respecting the constrained nature of relative abundance data. This chapter covers both regression and classification tasks using log-contrast transformations.

## What are Log-Contrast Models?

Compositional data (like microbiome relative abundances) sum to a constant and exist in a constrained space. Log-contrast models:
- Transform data using centered log-ratio (CLR) or other log-ratio transformations
- Apply regularized regression/classification in the transformed space
- Provide interpretable results in terms of relative abundance changes

## Chapter Organization

This chapter is organized into the following sections:

### 1. Data Preparation
Learn how to transform your microbiome data and prepare it for log-contrast modeling.

### 2. Regression Models
Predict continuous outcomes (e.g., temperature, pH) from microbiome composition:
- **Simple Log-Contrast**: Basic approach using CLR transformation
- **TRAC**: Incorporates taxonomic hierarchical information for better feature selection

### 3. Classification Models
Predict categorical outcomes (e.g., disease status, habitat type):
- **Simple Log-Contrast**: Basic approach for classification tasks
- **TRAC**: Uses taxonomic structure to identify predictive taxonomic groups

### 4. Advanced Topics
- **Concomitant Formulation**: Joint estimation of coefficients and noise level for heteroscedastic data

### 5. Interpretation
Understand how to interpret log-contrast model results and extract biological insights.

## Key Concepts

**Log-Ratio Transformations**: Convert compositional data to unrestricted space
- CLR (Centered Log-Ratio): Most common, centers data around geometric mean
- ALR (Additive Log-Ratio): Uses one component as reference

**Regularization**: Prevents overfitting in high-dimensional microbiome data
- L1 penalty (Lasso): Encourages sparsity, selects subset of features
- Stability selection: Identifies robust features across subsamples

**Taxonomic Adaptivity (TRAC)**: Leverages phylogenetic structure
- Computes adaptive weights based on taxonomic hierarchy
- Groups related taxa for more interpretable results

## Prerequisites

Before working through this chapter, ensure you have:
- Completed the [Installation](../01_installation/01_prerequisites.md) section
- Familiarized yourself with the [Atacama dataset](../00_getting_started/02_data.md)
- Basic understanding of regression and classification concepts

## Getting Started

Begin with [Data Preparation](01_data_preparation.md) to learn how to transform and prepare your data, then choose your analysis path:
- For continuous outcomes → [Regression](02_regression/01_simple_logcontrast.md)
- For categorical outcomes → [Classification](03_classification/01_simple_logcontrast.md)
