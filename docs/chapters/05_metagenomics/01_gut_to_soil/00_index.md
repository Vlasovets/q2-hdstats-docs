# Gut-to-Soil (16S, 335 ASVs x 99 samples)

The gut-to-soil axis study follows human faecal material through a composting
process and out into soil, sampling the community along the way. The tutorial
subsample used here is an ASV table of **335 features across 99 samples**, with
matching taxonomy and sample metadata, and it is the only dataset in this book
that can be downloaded and run today without waiting for a DOI.

It is used for three things, in this order.

| Page | What it does | Actions |
|---|---|---|
| [The Dataset](01_data.md) | provenance, the download, and an honest account of how shallow the subsample is | `qiime feature-table summarize`, `qiime metadata tabulate` |
| [Network Inference](02_network.md) | one network over all samples, then the gut-versus-soil group comparison | `transform-features`, `calculate-covariance`, `solve-problem`, `summarize`, `build-groups` |
| [Log-Contrast Regression](03_regression.md) | predicting a compost-stage outcome from the community, with taxonomy-derived trac weights | `transform-features`, `add-taxa`, `add-covariates`, `regress`, `classify`, `predict`, `summarize` |

## What this dataset is good for, and what it is not

**Good for**: seeing the whole workflow run on a table you did not curate, at a
size where $p$ exceeds $n$ by a few times rather than by a factor of a hundred;
and, uniquely in this book, for a group comparison where the two groups are
genuinely different communities rather than two halves of one. The joint
two-group fit is the one step that does not currently close through the QIIME 2
CLI — [Network Inference](02_network.md) shows how far it gets, why it stops, and
what to do instead.

**Not good for**: producing a result you would put in a paper. The table is a
10% subsample of the published study and the per-sample depths are low. A
covariance matrix estimated from counts that small is noisy before any penalty is
applied, so what these pages demonstrate is the *procedure*, and the procedure is
the thing worth copying. [The Dataset](01_data.md) quantifies the problem as far
as it can be quantified without a live run, and points at the full study for
anyone who wants to do this properly.

```{important}
Every page in this tier reads its metadata rather than assuming it. Column names
and their levels belong to the upstream study, and this book does not restate
them except where it says explicitly that it checked — the three columns named in
[Log-Contrast Regression](03_regression.md) were read off the served
`sample-metadata.tsv`; nothing else was. Where a command needs a column name it
otherwise takes it from a shell variable that you set once after looking at the
metadata, and the value shown in the assignment makes the command shape
unambiguous rather than being a claim about your file.
```

## Conventions on these pages

Commands assume the layout created in
[Download the Tutorial Data](../../00_getting_started/03_download_data.md): a
project directory containing `data/`, with the tier 3 files in
`data/gut-to-soil/`. That subdirectory is not cosmetic — the gut-to-soil metadata
file is called `sample-metadata.tsv`, exactly like the tier 2 one, and keeping
them apart avoids one silently overwriting the other.

Derived artifacts are written back into `data/gut-to-soil/` alongside the inputs,
with a `gts-` prefix so that a directory listing tells you what came from where.

You are assumed to have worked through
[Single Graphical Lasso](../../02_lowdim_gglasso/02_sgl.md) and the tier 1
q2-classo pages; neither the graphical lasso nor the log-contrast formulation is
re-explained here. The parameter-selection reasoning — eBIC, $\gamma$, the linear
path, sizing the latent block with $\mu_1$ — comes from
[Tier 2](../../04_highdim_atacama/00_index.md).

```{warning}
No page in this tier has been run against QIIME 2 2026.7. There are no captured
outputs, no edge counts and no $R^2$ values anywhere in it. See
[the tier overview](../00_index.md).
```
