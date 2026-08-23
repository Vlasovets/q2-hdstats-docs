# Data Preparation

A log-contrast model needs three things a count table does not carry: a design
matrix on the log-ratio scale, a constraint matrix that forces the coefficients
to sum to zero, and a vector of penalty weights. Three q2-classo actions build
them from the Atacama counts, each consuming the output of the one before it. A
final step, from sample-classifier, holds out a test set for each of the two
response types.

## Data Transformation

CLR-transform the counts:

```bash
qiime classo transform-features \
     --p-transformation clr \
     --p-coef 0.5 \
     --i-features data/atacama-counts.qza \
     --o-x data/xclr
```

The centered log-ratio divides every feature in a sample by that sample's
geometric mean and takes the logarithm, which discards the arbitrary sequencing
depth and moves the composition off the simplex. The action substitutes a
pseudocount for non-positive entries first, because the logarithm of zero is
undefined. It centres along rows and never transposes, so give it a table with
samples in rows; a feature-major input is transformed along the wrong axis
without an error.

## Add Taxonomic Information

Incorporate taxonomic classifications and compute adaptive weights:

```bash
qiime classo add-taxa \
    --i-features data/xclr.qza  \
    --i-taxa data/classification.qza \
    --o-x data/xtaxa \
    --o-aweights data/wtaxa
```

The action multiplies the transformed design by a matrix encoding the taxonomic tree, so
that the columns become internal nodes plus leaves and a coefficient can attach to a
whole clade rather than only to a single feature. It also divides each
node's weight by the number of leaves beneath it, which puts coarse and fine
ranks on a comparable footing under one penalty. There are no parameters to
set, and the rescaling is not optional.

## Add Covariates

Include environmental metadata with custom weights for each covariate:

```bash
qiime classo add-covariates \
    --i-features data/xtaxa.qza \
    --i-weights data/wtaxa.qza \
    --m-covariates-file data/atacama-selected-covariates-veg.tsv \
    --p-to-add ph average-soil-relative-humidity elevation average-soil-temperature vegetation \
    --p-w-to-add 1. 0.1 0.1 0.1 1 \
    --o-new-features data/xcovariates \
    --o-new-c data/ccovariates \
    --o-new-w data/wcovariates
```

Appended covariates enter the design matrix but not the zero-sum constraint: the
extended constraint row carries zeros in their columns, so their coefficients
are unconstrained. Each covariate also carries its own penalty weight, one of
the two ways to keep a covariate measured in metres from dominating CLR values
of order one under a shared penalty. The other is rescaling — centring a numeric
column and dividing it by its norm — which is off unless you ask for it. A
categorical column such as vegetation is expanded to one column per level, and
each level then appears separately in the coefficient plots.

## Split Data for Regression Analysis

Create training and test sets for continuous target prediction:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column average-soil-temperature \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/regress-xtraining \
    --o-test-table data/regress-xtest \
    --o-training-targets data/regress-training-targets.qza \
    --o-test-targets data/regress-test-targets.qza
```

The metadata column named on the command line supplies the targets, and the
random state fixes the partition, so you get the same split every time you run
this command.

## Split Data for Classification Analysis

Create training and test sets for categorical target prediction:

```bash
qiime sample-classifier split-table \
    --i-table data/xcovariates.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --m-metadata-column vegetation \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-stratify False \
    --o-training-table data/classify-xtraining \
    --o-test-table data/classify-xtest \
    --o-training-targets data/classify-training-targets.qza \
    --o-test-targets data/classify-test-targets.qza
```

The two splits differ only in that metadata column. Neither is stratified, so
the vegetation classes are not balanced across the training and test tables by
construction.

Each split leaves a training table, a test table and their targets. Those,
together with the constraint matrix and the weight vector written by the
preparation steps, are what you pass to the regression and classification
actions.