# Data Preparation

The Atacama soil microbiome dataset {cite}`neilson2017significant` contains:

- $N = 50$ samples from Atacama Desert soil
- $p = 13$ microbial taxa (ASVs)
- $q = 5$ environmental covariates: pH, elevation, temperature, humidity, and vegetation

[Data Overview](../00_getting_started/02_datasets.md) describes the dataset in more detail.

## Transforming the counts

Microbiome counts represent relative abundances constrained to sum to a constant. Transform them
before estimating a covariance.

```bash
# Transform compositional data using mCLR transformation
qiime gglasso transform-features \
     --p-transformation mclr \
     --p-add-metadata False \
     --p-scale-metadata False \
     --i-table data/atacama-counts.qza \
     --i-taxonomy data/classification.qza \
     --m-sample-metadata-file data/selected-atacama-sample-metadata.tsv \
     --o-transformed-table data/atacama-table-mclr.qza
```

The modified centred log-ratio moves the table into an unconstrained space, handles zeros
without adding pseudo-counts, and preserves the relative information between taxa. For the
standard centred log-ratio instead, pass `--p-transformation clr`.

## Building the input correlation

```bash
qiime gglasso calculate-covariance \
     --p-method scaled \
     --i-table data/atacama-table-mclr.qza \
     --o-covariance-matrix data/atacama-table-corr.qza
```

A scaled covariance is the Pearson correlation. For the covariance itself, pass
`--p-method unscaled`.

The input to the graphical lasso problem must be a
[positive semi-definite matrix](https://statproofbook.github.io/P/covmat-psd.html).
