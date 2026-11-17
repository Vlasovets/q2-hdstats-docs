# Getting Started

Start with the [Atacama soil microbiome](tutorial/s2.md) example. It introduces
the sample metadata and feature tables used for regression, classification, and
network estimation. You can follow the smaller 15-ASV example before moving on
to the high-dimensional analyses.

## Before You Begin

You will need a QIIME 2 environment with **q2-classo** and **q2-gglasso** installed,
along with the feature tables, taxonomy, and metadata for your chosen example.
Check [Reproducibility and source versions](reproducibility.md) for the package
revisions used to generate the results shown here.

```{important}
Prepared inputs are available for the Atacama 15-ASV and top-300 examples and
the MAP mOTU example. Download the bundle linked on the relevant data page;
do not substitute older feature tables. These inputs and the scientific
workflows were checked in the established, pinned environment. A fresh-machine
installation test and portable entry points for all optional figures remain
to be completed. Cocoa input redistribution is not yet confirmed.
```

## Follow a Worked Example

1. [Inspect the Atacama data and metadata](tutorial/s2.md).
2. [Get started with the Atacama top-15 ASV data](tutorial/s2-01.md) to predict elevation, classify vegetation,
   and compare graphical models.
3. [Use the larger Atacama table](tutorial/s2-02.md) and include sample covariates
   in the prediction and network models.

Then explore [infant metagenomes from the MAP study](tutorial/s3.md) or the
[cocoa-bean fermentation data](tutorial/s4.md). These examples introduce
taxonomic aggregation and a time-ordered prediction task.

## Read More About the Methods

Follow the applications at your own pace. When you want more detail, refer to:

- [Compositional data and transformations](tutorial/s1-01.md).
- [Regression, classification, and model selection](tutorial/s1-02.md).
- [Graphical models and latent variables](tutorial/s1-03.md).
- [Notation, assumptions, and limitations](tutorial/s5.md).

The statistical-background command blocks illustrate the model interfaces.
Use the case-study commands with their matching prepared inputs when following
an analysis.
