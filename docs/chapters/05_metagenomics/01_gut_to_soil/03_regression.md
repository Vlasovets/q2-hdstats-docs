# Log-Contrast Regression

[Network Inference](02_network.md) asked which taxa are conditionally dependent on
each other. This page asks the complementary question: **can the community
predict where a sample sits in the composting process, and which clades carry that
signal?**

The answer comes from a sparse log-contrast model {cite}`aitchison1984log`
{cite}`shi2016regression`: a linear model on log-abundances whose coefficients are
constrained to sum to zero, so that the prediction depends only on *ratios*
between taxa and not on sequencing depth. On top of that, this page uses **trac**
{cite}`bien2021tree`, which lets the model select whole clades rather than
individual ASVs — the right move on a shallow table, because a genus supported by
twenty sparse ASVs is estimable when none of the twenty is.

The full chain is `transform-features` → `add-taxa` → `add-covariates` →
`split-table` → `regress` → `predict` → `summarize`. Every step is demonstrated
in [Tier 1](../../03_lowdim_classo/03_regression/02_trac.md); what is new here is
the scale, the shallowness, and a dataset whose metadata you have to read for
yourself.

## Choosing the outcome

Two outcomes are natural for this design, and they take different actions:

- **A continuous compost stage** — elapsed time, or a pile-age variable →
  `qiime classo regress`.
- **A binary sample type** — gut-derived versus soil-derived →
  `qiime classo classify`.

```bash
# These are real column names, read off the downloaded sample-metadata.tsv.
# Note the spaces: keep the quotes.
OUTCOME_COLUMN="Composting Time Point"   # numeric, ranges 1-52
CLASS_COLUMN="SampleType"                # categorical -- see the warning below
COVARIATE_COLUMN="Compost pH"            # numeric, ranges 4-10
```

```{note}
These names were verified against the `sample-metadata.tsv` served by the
gut-to-soil tutorial site. That file describes the **full 1660-sample study**,
while `asv-table-ms2.qza` holds 99 of those samples — so QIIME 2 will warn about
metadata rows with no matching sample, which is expected here. Several columns
are also incomplete across the full study (`Composting Time Point` is present for
1267 rows, `Compost pH` for 1209); check coverage on your 99 before relying on a
column.
```

```{warning}
**`SampleType` has 15 levels, not 2.** The two largest are
`Human Excrement Compost` (799 rows) and `Human Excrement Compost Pre-Roll`
(453). `qiime classo classify` accepts binary outcomes only, so you must reduce
it to a two-level contrast first — see the admonition below.
```

```{important}
`qiime classo classify` accepts **binary** outcomes only. The underlying function
collects the distinct values of the column and raises

    ValueError: Metadata column y is supposed to be binary, but takes more than
    2 different values : ...

if there are more than two. A multi-level sample-type column has to be collapsed
to two levels first — with `qiime feature-table filter-samples --p-where` on the
levels you want to contrast, for example — or modelled as a continuous stage
instead.
```

## Step 1: Transform

```bash
qiime classo transform-features \
    --i-features data/gut-to-soil/asv-table-ms2.qza \
    --p-transformation clr \
    --p-coef 0.5 \
    --o-x data/gut-to-soil/gts-classo-clr.qza
```

`clr` is the only value `--p-transformation` accepts in q2-classo — there is no
`mclr` here, and no `Choices()` on the parameter, so a typo becomes a runtime
`ValueError` rather than a command-line error.

`--p-coef 0.5` is a **zero replacement**, not a global pseudo-count: every
non-positive entry is set to `0.5` before the log, and positive counts are left
alone. On a table as shallow as this one that is a real modelling decision. Most
cells are zero, so most features sit on a common floor within a sample, and the
model sees a large block of features whose values move together with the sample's
geometric mean. Raising `--p-coef` toward 1 makes zeros look more like observed
singletons; lowering it pushes them further below the observed counts. Try both
before trusting a coefficient on a rare clade.

```{note}
This transform is q2-classo's own and is independent of the one used in
[Network Inference](02_network.md); the two plugins do not share transformed
artifacts. `qiime gglasso transform-features` emits `FeatureTable[Frequency]` and
`qiime classo transform-features` emits `FeatureTable[Design]`, and the actions
downstream of each expect their own.
```

## Step 2: Aggregate along the taxonomy (trac)

```bash
qiime classo add-taxa \
    --i-features data/gut-to-soil/gts-classo-clr.qza \
    --i-taxa data/gut-to-soil/taxonomy.qza \
    --o-x data/gut-to-soil/gts-taxa.qza \
    --o-aweights data/gut-to-soil/gts-w-taxa.qza
```

`add-taxa` takes no parameters. It builds an aggregation matrix $A$ from the
taxonomic tree, in which column $j$ marks the leaves belonging to node $j$, and
replaces the design with $XA$ divided by the number of leaves under each node. The
penalty weight it returns for each node is $1/\text{nleaves}$, so a broad clade is
penalized less per unit of aggregated signal than a single ASV — that is the
mechanism by which the model prefers to explain the outcome with a genus or a
family when the ASV-level evidence is individually too thin.

Two consequences that surprise people:

- **The design columns are no longer ASVs.** They are tree nodes: internal clades
  as well as leaves. Coefficient labels in the final visualization are clade
  names, and any feature-level joining you had planned no longer applies.
- **The taxonomy is now part of the model, not decoration.** In the network
  chapter, `--i-taxonomy` was required and unread. Here `--i-taxa` determines the
  hypothesis space. A different classifier or a different reference database
  gives a different set of candidate predictors.

## Step 3: Add covariates

`add-covariates` appends metadata columns to the design as ordinary predictors,
and — importantly — updates the zero-sum constraint so the appended columns are
**excluded** from it. Only the compositional part of the design has to sum to
zero; an environmental covariate is not part of the composition.

```bash
qiime classo add-covariates \
    --i-features data/gut-to-soil/gts-taxa.qza \
    --i-weights data/gut-to-soil/gts-w-taxa.qza \
    --m-covariates-file data/gut-to-soil/sample-metadata.tsv \
    --p-to-add "${COVARIATE_COLUMN}" \
    --p-rescale True \
    --p-w-to-add 0.1 \
    --o-new-features data/gut-to-soil/gts-design.qza \
    --o-new-c data/gut-to-soil/gts-c.qza \
    --o-new-w data/gut-to-soil/gts-w.qza
```

```{important}
**Never put the outcome in the design.** The column added here is
`Compost pH` — a genuine covariate, not the thing being predicted. Adding the
outcome to $X$ and then regressing on it leaks the target: the model recovers it
from itself, $R^2$ looks excellent, and the result says nothing about the
microbiome.

The outcome never needs to be in the design. It reaches `split-table` through
`--m-metadata-column` (Step 4) and `regress` through `--m-y-column` (Step 5),
both read straight from the metadata file. Step 3 is *optional* — skip it
entirely if you want a pure log-contrast model on the taxa alone, and pass
`gts-taxa.qza` / `gts-w-taxa.qza` to Step 4 instead.

The filtered-adjustment analysis in
[Log-Contrast Models at Scale](../../04_highdim_atacama/05_classo_cv.md) works
through why a covariate that is a proxy for the outcome causes the same problem
in a subtler form.
```

**Explanation:**

- `--m-covariates-file` and `--p-to-add` are both required; there is no default
  set of columns.
- `--p-rescale` and `--p-w-to-add` are lists with one entry per added column.
  `--p-rescale True` standardizes the column before it joins a design whose other
  columns are log-ratios; without it, a covariate measured in hundreds swamps a
  scalar penalty. `--p-w-to-add 0.1` penalizes the covariate less than the taxa,
  which is the usual choice when the covariate is a control rather than a
  candidate discovery.
- **Categorical columns are one-hot expanded**, and the resulting design columns
  are labelled `"<name> = <value>"` — with spaces around the equals sign. That is
  the label you will see in the summary visualization, and the string you will
  need if you ever match coefficients programmatically.

## Step 4: Split

Split *after* the design is final, so that training and test tables have
identical columns. `predict` matches the stored coefficient labels against the
columns of whatever table you hand it.

```bash
qiime sample-classifier split-table \
    --i-table data/gut-to-soil/gts-design.qza \
    --m-metadata-file data/gut-to-soil/sample-metadata.tsv \
    --m-metadata-column "${OUTCOME_COLUMN}" \
    --p-test-size 0.2 \
    --p-random-state 42 \
    --p-no-stratify \
    --o-training-table data/gut-to-soil/gts-train.qza \
    --o-test-table data/gut-to-soil/gts-test.qza \
    --o-training-targets data/gut-to-soil/gts-train-targets.qza \
    --o-test-targets data/gut-to-soil/gts-test-targets.qza
```

With 99 samples a 20% test split leaves a test set of around twenty. That is
enough to notice a model that has failed and nowhere near enough to estimate its
performance precisely — treat the held-out score as a sanity check and the
cross-validated score inside `regress` as the number you report.

## Step 5: Fit

```bash
qiime classo regress \
    --i-features data/gut-to-soil/gts-train.qza \
    --i-c data/gut-to-soil/gts-c.qza \
    --i-weights data/gut-to-soil/gts-w.qza \
    --m-y-file data/gut-to-soil/sample-metadata.tsv \
    --m-y-column "${OUTCOME_COLUMN}" \
    --p-do-yshift \
    --p-concomitant \
    --p-path --p-path-nlam-log 100 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 100 --p-cv-lamin 0.001 --p-cv-logscale \
    --p-stabsel --p-stabsel-threshold 0.7 \
    --p-no-lamfixed \
    --o-result data/gut-to-soil/gts-regress.qza \
    --verbose
```

**Why these settings, for this table:**

- `--p-concomitant` (the default) estimates the noise scale jointly with the
  coefficients instead of assuming it. On a shallow table the residual scale is
  genuinely unknown, and a $\lambda$ chosen under a wrong assumed scale is the
  wrong $\lambda$. Turn it off with `--p-no-concomitant` only if you are matching
  a published analysis that did.
- `--p-do-yshift` centres the outcome, which is what makes the intercept
  interpretable alongside a zero-sum design.
- `--p-cv --p-cv-subsets 5 --p-cv-one-se` selects $\lambda$ by 5-fold
  cross-validation with the one-standard-error rule — the sparser model within
  one standard error of the best. On $n \approx 80$ training samples the CV curve
  is noisy enough that the one-SE rule is doing real work; without it you will
  select a denser model than the data support.
- `--p-stabsel --p-stabsel-threshold 0.7` reports how often each clade is selected
  across subsamples {cite}`meinshausen2010stability`. This is the output to
  believe. A coefficient that appears at the CV optimum but has low selection
  frequency is a coefficient you should not put in a table.
- `--p-no-lamfixed` skips the fit at a fixed theoretical $\lambda$, which is one
  less thing to compute when you are not using it.

```{note}
`--p-cv-nlam` is the current spelling. The parameter was originally `cv__nlam`
with a double underscore, rendered by QIIME 2 as `--p-cv--nlam`; that spelling
still works but emits a `DeprecationWarning`, and passing both with different
values is an error. See
[Troubleshooting](../../90_reference/04_troubleshooting.md).
```

### The classification variant

```bash
qiime classo classify \
    --i-features data/gut-to-soil/gts-train.qza \
    --i-c data/gut-to-soil/gts-c.qza \
    --i-weights data/gut-to-soil/gts-w.qza \
    --m-y-file data/gut-to-soil/sample-metadata.tsv \
    --m-y-column "${CLASS_COLUMN}" \
    --p-huber --p-rho 0.0 \
    --p-path --p-path-nlam-log 100 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 100 --p-cv-lamin 0.001 --p-cv-logscale \
    --p-stabsel --p-stabsel-threshold 0.7 \
    --p-no-lamfixed \
    --o-result data/gut-to-soil/gts-classify.qza \
    --verbose
```

```{warning}
`classify` is **not** `regress` with a categorical outcome. Three differences bite:

- There is **no `--p-concomitant`**. The parameter does not exist on this action
  and the solver forces the concomitant formulation off for classification.
  `qiime classo classify --p-concomitant True` fails at the command line.
- There is no `--p-do-yshift` either.
- `--p-rho` defaults to `0.0` here rather than `1.345`. With `--p-huber` it sets
  the robustness threshold of the Huber hinge loss, which is the substitute for
  the concomitant formulation when you need tolerance to outliers.
```

Predicting gut-derived versus soil-derived from the community is, on this design,
close to a sanity check rather than a discovery — the two groups are meant to be
different communities, and a model that cannot separate them is telling you
something has gone wrong upstream. Its value is in the *stability selection
output*: which clades the separation actually rests on.

## Step 6: Predict and summarize

```bash
qiime classo predict \
    --i-features data/gut-to-soil/gts-test.qza \
    --i-problem data/gut-to-soil/gts-regress.qza \
    --o-predictions data/gut-to-soil/gts-predictions.qza

qiime classo summarize \
    --i-problem data/gut-to-soil/gts-regress.qza \
    --i-taxa data/gut-to-soil/taxonomy.qza \
    --i-predictions data/gut-to-soil/gts-predictions.qza \
    --p-maxplot 200 \
    --o-visualization data/gut-to-soil/gts-regress-summary.qzv
```

`predict` takes no parameters at all. It reads the coefficient labels stored in
the problem and rebuilds the design matrix column by column, which is why the test
table has to have come through the same `add-taxa` / `add-covariates` chain.

`--p-maxplot 200` (the registered default) caps the number of coefficients drawn
in the stability-selection profile and the beta plots. After `add-taxa` the
design has one column per tree
node — leaves plus internal clades, so comfortably more than the 335 ASVs — and
the cap matters here in a way it did not at $p = 13$. Raising it produces a taller
and slower figure, not more information, since a sparse model leaves most
coefficients at exactly zero.

Open the `.qzv` at [QIIME 2 View](https://view.qiime2.org/).

```{note}
Cross-validated $R^2$, held-out error, selection frequencies and the identity of
the selected clades are **pending verification against QIIME 2 2026.7**. Nothing
in this tier has been run, and no results are quoted.
```

## Reading the result honestly

Three things to keep in view when you look at the summary.

**A selected clade is not a mechanism.** The log-contrast model finds a set of
ratios that tracks the outcome. On a composting series almost everything tracks
the outcome, because the whole community turns over; the model picks a small
sufficient subset, and a different subset would often have worked nearly as well.
Stability selection frequencies tell you which choices were forced by the data
and which were arbitrary.

**Depth limits what is estimable.** A clade whose aggregated abundance is a
handful of reads per sample has a log-ratio dominated by the zero replacement from
Step 1. If such a clade appears with a large coefficient, check its raw counts in
`gts-table-summary.qzv` before writing it down.

**Compare it with the network.** The two models see the same table through
different lenses, and agreement between them is informative. If the clades that
predict compost stage sit on nodes that the sparse network leaves isolated, one of
the two analyses is describing noise. If the stage signal instead lines up with
the low-rank component from
[Step 5 of the network chapter](02_network.md) — the block that absorbs the
dominant gradient — that is the coherent outcome, and it is the same structure
seen twice.
