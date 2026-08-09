# High-Dimensional Example: q2-classo on 300 Atacama ASVs

This chapter applies q2-classo to the same **300-ASV, 54-sample** Atacama design
used in the [high-dimensional graphical-lasso chapter](02_model_selection.md),
predicting each continuous environmental covariate from the microbiome with
cross-validated log-contrast regression. It reproduces the reference analysis of
Christian L. Müller.

```{important}
**Build the design with `qiime classo transform-features`, not with
`qiime gglasso transform-features`.**

Not because the gglasso action is broken — it is not, any more — but because
this keeps the two chains independent, so a change to the network side cannot
quietly alter the regression results. It is also the shorter path: `classo`'s own
transform takes the raw counts directly.

Both are now valid inputs. Until QIIME 2 2026.7, only one was: `gglasso
transform-features` stored its output with the axes swapped, so the table's
"samples" were feature IDs, `regress` found no overlap with the outcome on the
sample index, and c-lasso failed deep inside with

    IndexError: index 0 is out of bounds for axis 0 with size 0

on a design of shape `(0, 54)`. If you hit that on an **older** artifact, that is
the cause — regenerate it rather than transposing by hand, because
`calculate-covariance` was compensating for the same swap and the two only give
the right answer together. See
[Troubleshooting](../90_reference/04_troubleshooting.md).
```

## Setup

Each of the 15 continuous covariates is used **once as the outcome**; missing
outcome values are mean-imputed. Models use the log-contrast formulation with
an intercept, selected by **5-fold cross-validation** with the
**one-standard-error** rule along a log-spaced $\lambda$ path.

Build the design once:

```bash
qiime classo transform-features \
    --i-features data/atacama-top-300-table.qza \
    --o-x data/atacama-top-300-classo-clr.qza
```

Then fit one model per outcome:

```bash
qiime classo regress \
    --i-features data/atacama-top-300-classo-clr.qza \
    --m-y-file data/atacama-classo-outcomes-mean-imputed.tsv \
    --m-y-column "<outcome>" \
    --p-concomitant \
    --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed \
    --o-result "<outcome>-cv5.qza"
```

`--p-cv-nlam`, single dash before `nlam`. The old `--p-cv--nlam` spelling still
works but is deprecated.

## What the recompute produced

All 15 fits were re-run on QIIME 2 2026.7. This table is **generated** from the
solution artifacts by `slurm/12_classo_summary.sh`:

```{csv-table} Cross-validated log-contrast fits, 300 ASVs x 54 samples
:file: ../../_data/atacama-classo-cv.tsv
:delim: tab
:header-rows: 1
:widths: 34, 16, 16, 17, 17
```

Read the sparsity, not the error scale: the CV error is in the outcome's own
units squared, so it is comparable *across $\lambda$ for one outcome* and not
across outcomes. What is comparable is how many of the 300 ASVs survive — from a
single feature for `depth`, `ec` and `toc` up to thirty-three for
`extract-concen`.

## Prediction from ASVs only (base R1)

```{important}
**The $R^2$ table below has NOT been reproduced.** The recompute reports
cross-validated *error*, which is what the solution artifact stores; turning that
into an out-of-sample $R^2$ needs a `qiime classo predict` pass on held-out
samples, which has not been run. The figures here are carried over from the
reference analysis and are kept because the *ranking* is the chapter's argument —
but treat them as unverified, and do not quote an individual value.

The same applies to the joint and filtered $R^2$ values further down, and to the
named first-selected taxon.
```

The out-of-sample $R^2$ (mean across the 5 folds) shows which environmental
variables are predictable from the 300 ASVs alone:

```{figure} ../../images/png/atacama-full/atacama-top-300-r1-cv5-selected-taxa-heatmap.png
:name: fig-classo-asv-only
:width: 100%

One row per outcome, ordered by cross-validated $R^2$ (left column); one column
per ASV selected by at least one 1-SE model. Colour is the coefficient divided
by the largest absolute coefficient **within that row**, so shades are comparable
along a row and *not* down a column.

Three things to read off it. **Most outcomes are not predictable from the
microbiome alone** — only the top two clear $R^2 = 0.4$. **The bottom four rows
(`toc`, `relative-humidity-soil-high`, `ec`, `depth`) have negative $R^2$**: the
model does worse than predicting the mean, and their rows are correspondingly
almost empty. A negative $R^2$ here is the honest answer, not a bug. Finally,
**`Pseudarthrobacter` (6b780e) is the most widely shared predictor**, carrying
weight in nearly every outcome that is predictable at all — the genus discussed
in [Network interpretation](../02_lowdim_gglasso/09_interpretation.md).
```

```{important}
**The `ph` row is contaminated.** Eight samples in the Atacama metadata carry
`ph = 0`, which is a missing-value sentinel rather than a measurement — pH 0 is
not a soil. Those rows enter this regression as genuine values sitting about 2.8
standard deviations below the mean, so the `ph` result in this figure and the
next one should not be interpreted. See
[q2-classo parameters](../90_reference/03_classo_parameters.md) for the check
that surfaces it.
```

| Outcome | CV $R^2$ |
|---------|----------|
| extract-concen | 0.61 |
| percent-relative-humidity-soil-100 | 0.48 |
| amplicon-concentration | 0.35 |
| average-soil-relative-humidity | 0.30 |
| elevation | 0.26 |
| relative-humidity-soil-low | 0.26 |
| percentcover | 0.26 |
| temperature-soil-high | 0.24 |
| average-soil-temperature | 0.20 |
| ph | 0.17 |
| temperature-soil-low | 0.13 |
| toc / relative-humidity-soil-high / ec / depth | $\le 0$ |

The first selected taxon (largest-magnitude coefficient column, after the
intercept) is a *Pseudarthrobacter* ASV — a genus characteristic of the Atacama
soil community (see the [interpretation notes](../02_lowdim_gglasso/09_interpretation.md)).

## Adjusting for environmental covariates (joint and filtered)

Two variants add the *other* covariates as extra predictors (each rescaled and
L2-normalized, covariate penalty weight `0.1626`):

- **Joint** — add **all** other covariates as predictors.
- **Filtered adjustment** — add only covariates that are **not strongly
  correlated** with the outcome (Pearson $\lvert r\rvert < 0.80$), so a covariate that is
  essentially a proxy for the outcome cannot leak it.

This is a two-step workflow. **First**, build the augmented design by adding the
chosen covariate(s) to the ASV design (this produces new features, constraint,
and weights artifacts):

```bash
qiime classo add-covariates \
    --i-features data/atacama-top-300-classo-clr.qza \
    --m-covariates-file data/atacama-classo-outcomes-mean-imputed.tsv \
    --p-to-add <covariate> --p-rescale --p-w-to-add 0.162565105 \
    --o-new-features <design>.qza --o-new-c <c>.qza --o-new-w <w>.qza
```

**Then** run the same cross-validated regression as before, but on the augmented
design (passing the constraint and weights from the previous step):

```bash
qiime classo regress \
    --i-features <design>.qza --i-c <c>.qza --i-weights <w>.qza \
    --m-y-file atacama-classo-outcomes-mean-imputed.tsv --m-y-column <outcome> \
    --p-do-yshift --p-path --p-path-nlam-log 120 --p-path-lamin-log 0.0001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 120 --p-cv-lamin 0.0001 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed --p-no-concomitant --p-no-huber --p-intercept \
    --o-result <outcome>-joint-r1-cv5.qza
```

Adding covariates sharply increases predictability of the physically-coupled
variables — e.g. soil temperature and humidity — but the **filtered** analysis
removes most of that gain, showing it came largely from outcome-correlated
covariates rather than the microbiome:

```{figure} ../../images/png/atacama-full/atacama-filtered-r1-cv5-selected-predictors-heatmap.png
:name: fig-classo-filtered
:width: 100%

The filtered-adjustment models. Same layout as {numref}`fig-classo-asv-only`,
but the predictor block now begins with the added **covariates** (left of the
vertical rule) before the ASVs, and three score columns replace one: the
filtered $R^2$, its gain over the ASV-only model, and its loss relative to
using *all* covariates.

The two gain columns are the argument. `elevation` reaches 0.84, of which
$+0.58$ comes from covariates and only $-0.03$ is given up by filtering — its
predictors were never proxies. The soil-temperature outcomes behave the
opposite way: `temperature-soil-low` gains $+0.38$ over ASVs alone but loses
$-0.45$ against the unfiltered model, so most of that apparent skill was one
covariate standing in for the outcome. `toc`, `ec` and `depth` stay negative
with $\Delta R^2 \approx 0$ throughout — nothing helps, and the near-empty
rows show the model selecting little beyond the intercept.
```

```{note}
The `×` and `−` marks in the covariate block are not decoded here. This figure
is carried over from the reference analysis and its marker legend was not
preserved with it; there is no script in this repository that regenerates it.
Read the colours and the three score columns, which are labelled, and treat the
markers as unexplained rather than inferring a meaning for them.
```

| Outcome | ASV-only | Joint (all cov.) | Filtered ($\lvert r\rvert<0.8$) |
|---------|----------|------------------|----------------------|
| average-soil-temperature | 0.20 | 0.98 | 0.54 |
| temperature-soil-low | 0.13 | 0.96 | 0.51 |
| temperature-soil-high | 0.24 | 0.94 | 0.56 |
| average-soil-relative-humidity | 0.30 | 0.88 | 0.74 |
| relative-humidity-soil-low | 0.26 | 0.88 | 0.59 |
| elevation | 0.26 | 0.87 | 0.84 |
| percentcover | 0.26 | 0.81 | 0.45 |
| percent-relative-humidity-soil-100 | 0.48 | 0.80 | 0.58 |
| extract-concen | 0.61 | 0.59 | 0.59 |

Elevation stands out: its filtered $R^2$ (0.84) stays close to the joint value
(0.87), i.e. its predictability does **not** depend on outcome-correlated
covariates. For the temperature variables, by contrast, the joint gain
(≈ 0.94–0.98) collapses under filtering (≈ 0.51–0.56), a textbook illustration
of why task-specific covariate filtering matters.
