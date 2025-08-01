# Atacama soil microbiome analysis

We'll use the Atacama soil microbiome dataset {cite}`neilson2017significant`, which contains:
- 50 samples from Atacama Desert soil
- 13 microbial taxa (ASVs)
- Environmental covariates: pH, elevation, temperature, humidity

The original data is available through the European Nucleotide Archive under accession [ERP019482](https://www.ebi.ac.uk/ena/browser/view/PRJEB17617).

## Load the Data

```bash
# Create data directory if it doesn't exist
mkdir -p data

# Download the example data files
wget -O data/atacama-counts.qza https://github.com/bio-datascience/q2-gglasso/tree/master/data/atacama-counts.qza
wget -O data/classification.qza https://github.com/bio-datascience/q2-gglasso/tree/master/data/classification.qza
wget -O data/selected-atacama-sample-metadata.tsv https://github.com/bio-datascience/q2-gglasso/tree/master/data/selected-atacama-sample-metadata.tsv
wget -O data/atacama-selected-covariates.tsv https://github.com/bio-datascience/q2-gglasso/tree/master/data/atacama-selected-covariates.tsv
```

## Data Transformation

Microbiome data represents relative abundances constrained to sum to a constant.

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

This transformation:
- Converts compositional data to unconstrained space
- Handles zeros without adding pseudo-counts
- Preserves the relative information between taxa

**Note:** You can also use standard CLR transformation with `--p-transformation clr`.

## Create an input correlation

```bash
qiime gglasso calculate-covariance \
     --p-method scaled \
     --i-table data/atacama-table-mclr.qza \
     --o-covariance-matrix data/atacama-table-corr.qza
```
This method:
- Calculate a scaled covariance, also known as the Pearson correlation.
- One can also use a simple covaraince with `--p-method unscaled`

**Note:** Input for the graphical lasso problem must be a [positive semi-definite matrix](https://statproofbook.github.io/P/covmat-psd.html).