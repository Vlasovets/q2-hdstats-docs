# High-Dimensional Statistics with QIIME 2

Use **q2-classo** to fit regression and classification models, and
**q2-gglasso** to estimate microbial association networks in QIIME 2.

Begin with [Getting started](getting-started.md), then follow the Atacama
example. You will first inspect the data and sample metadata, work with a small
set of ASVs, and then repeat the analyses with a larger feature table.

| Tutorial | What you will do |
| --- | --- |
| [Atacama soil microbiome](tutorial/s2.md) | Predict elevation, classify vegetation, and estimate networks from soil ASVs |
| [MAP study metagenomics](tutorial/s3.md) | Predict infant age with sparse log-contrast and trac models, then estimate mOTU networks |
| [Cocoa-bean fermentation](tutorial/s4.md) | Train a bean-type classifier and evaluate it on later fermentation time points |

You do not need to read the full statistical background before starting.
Consult [Statistical background](tutorial/s1.md) for the model formulations and
[Summary tables](tutorial/s5.md) for notation, assumptions, and limitations.

These examples use the analyses developed for the paper's supplementary
material. The website puts the applications first; the underlying methods,
results, and figure provenance are unchanged.
