# Synthetic Data with Known Truth

Every other chapter in this book asks you to download something. This one does
not. `qiime classo generate-data` synthesizes a design matrix, a zero-sum
constraint matrix and a response vector from a model whose **true support is
known**, which makes it the one end-to-end exercise you can run five minutes
after `conda env create` finishes.

Treat it as a smoke test. It answers three questions that
[Verifying Your Installation](../01_installation/04_verify.md) cannot: does the
c-lasso solver converge inside QIIME 2, does the zarr round-trip through
`CLASSOProblem` survive, and does the visualizer render. Because the true
coefficients are planted rather than inferred, you also get a correctness check
for free — the selected features should be the planted ones.

## What the action actually does

`generate-data` builds a design matrix `X` with `d` columns, draws a coefficient
vector with exactly `d_nonzero` non-zero entries subject to a zero-sum
constraint, and produces `n` responses. It emits two artifacts:

| Output | Type | Contents |
|---|---|---|
| `--o-x` | `FeatureTable[Design]` | the `n × d` design matrix |
| `--o-c` | `ConstraintMatrix` | the `1 × d` zero-sum constraint row |

The response `y` is **not** an artifact. It is written to a file called
`randomy.tsv` in your current working directory — see the warning below.

Column names are generated, not drawn from anything: the first `d // 2` columns
are named `A0`, `A1`, … and the remainder `B0`, `B1`, …. For `--p-d 20` you get
`A0`–`A9` followed by `B0`–`B9`.

```{warning}
**Known bug: `randomy.tsv` is written to the current working directory.**

The response vector is written with a hardcoded relative path, so it lands
wherever you happened to invoke `qiime` from — not next to your `--o-x`
artifact, and not in a QIIME 2 provenance-tracked location. Two consequences:

* Re-running `generate-data` **silently overwrites** the previous `randomy.tsv`.
  The `x.qza` you generated a minute ago is then paired with somebody else's
  `y`, and nothing will warn you about it.
* If you run the action from a shared or read-only directory, it fails there
  rather than at the point where you asked for an output path.

Work around it by running the whole chapter from a dedicated scratch directory
and copying `randomy.tsv` alongside the artifacts it belongs to.
```

## Step 0: a directory to work in

```bash
mkdir -p smoke-test
cd smoke-test
```

Every command in this chapter is run from `smoke-test/`, and every path is
relative to it. This is not tidiness for its own sake — it is the workaround for
the `randomy.tsv` side effect.

## Step 1: generate the problem

```bash
qiime classo generate-data \
    --p-n 100 \
    --p-d 20 \
    --p-d-nonzero 3 \
    --o-x synthetic-x.qza \
    --o-c synthetic-c.qza \
    --verbose
```

**Why these values.** The defaults are `--p-n 100 --p-d 80 --p-d-nonzero 5`,
which is already a mildly high-dimensional problem. For a smoke test you want
the easy regime — comfortably more samples than features — so that a failure to
recover the support points at the software rather than at the statistics. With
`n = 100` and `d = 20`, three planted coefficients should be recoverable.

**`--verbose` is not optional here.** The list of planted non-zero labels is
printed to standard output by the function body; it is not stored in either
output artifact. Without `--verbose`, q2cli swallows it and the ground truth is
gone. Copy that line into a file before you do anything else.

```{important}
**`generate-data` is unseeded.** The underlying `random_data` call is made with
`seed=None` and the action exposes no seed parameter, so every invocation
produces a different design, a different `y`, and a *different planted support*.
You cannot reproduce a run, and you cannot compare two runs. Record the printed
support and the artifact UUIDs if you need to refer back to a result.
```

## Step 2: CLR-transform the design

```bash
qiime classo transform-features \
    --i-features synthetic-x.qza \
    --p-transformation clr \
    --p-coef 0.5 \
    --o-x synthetic-xclr.qza
```

The generator is invoked with `exp=True`, so the emitted design is on a
positive, multiplicative scale — the same shape of data as a count table, and
the same reason you would CLR-transform a real one.

The transform is safe with respect to the planted model because the constraint
is zero-sum. CLR subtracts a per-sample constant from every log-feature, and a
coefficient vector that sums to zero annihilates a constant shift exactly. So
for any feasible β the fitted values from `clr(X)β` and `log(X)β` are identical;
the CLR step changes the parameterization, not the model. This is the same
argument that makes log-contrast regression well defined on compositions in the
first place, and it is worth internalizing here on synthetic data before it
matters on real data in [Data Preparation](02_data_preparation.md).

`--p-coef 0.5` is the pseudo count substituted for non-positive entries. The
generated design has no zeros, so it does nothing here; it is left in the
command so the smoke test exercises the same code path your real analyses will.

```{note}
`--p-transformation` accepts only `clr`. The plugin declares no `Choices()`, so a
typo such as `--p-transformation clrr` is accepted by the CLI and fails inside
the function with `ValueError: Unknown transformation name`. See
[Troubleshooting](../90_reference/04_troubleshooting.md).
```

## Step 3: fit

```bash
qiime classo regress \
    --i-features synthetic-xclr.qza \
    --i-c synthetic-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-path \
    --p-cv \
    --p-cv-seed 1 \
    --p-stabsel \
    --p-stabsel-seed 1 \
    --p-lamfixed \
    --o-result synthetic-result.qza \
    --verbose
```

`randomy.tsv` has exactly two columns: an ID column named `id` holding
`0 … n-1`, and the response in a column named `col`. The design matrix rows carry
the same IDs, so the metadata join is on integers-as-strings and needs no
preparation. `id` is one of the identifier headers QIIME 2 recognises, so the
file is usable as a metadata file as written.

`--i-c` is optional; without it c-lasso applies its own zero-sum default. Pass it
anyway — being explicit costs nothing and it is the constraint the data were
actually generated under.

All four model-selection blocks are switched on so that the smoke test covers
all of them. Each seed is pinned because cross-validation and stability
selection resample: `--p-cv-seed` defaults to `1` already, but
`--p-stabsel-seed` has **no** default and leaves stability selection
non-reproducible unless you set it. What each of these four does, and how to
choose between them, is [Advanced: Choosing a Model](05_advanced/02_model_selection.md).

```{note}
The default formulation for `regress` is the concomitant one (`R3`), because
`--p-concomitant` defaults to `True`. That is a reasonable default and the smoke
test does not override it. See
[Advanced: Concomitant Formulation](05_advanced/01_concomitant_formulation.md).
```

## Step 4: check the recovered support

```bash
qiime classo summarize \
    --i-problem synthetic-result.qza \
    --p-maxplot 25 \
    --o-visualization synthetic-result.qzv
```

Open the `.qzv` at [QIIME 2 View](https://view.qiime2.org/) and go to the
**Stability Selection** tab. The table of selected parameters lists the labels
whose selection probability cleared `--p-stabsel-threshold`, and
`StabSel-prob.csv` gives the full profile for every column. Compare that list
against the labels `generate-data` printed in step 1.

`--p-maxplot 25` is chosen to exceed the number of design columns: 20 features
plus the intercept is 21. Below that the bar plots are silently truncated —
[Predict & Summarize](06_predict_and_summarize.md) explains exactly what gets
dropped.

The **Cross-Validation** and **LAM fixed** tabs give you two more independent
reads on the same question, each with its own downloadable coefficient table
(`CV-beta.csv`, `LAM-beta.csv`).

```{note}
No pass/fail numbers are quoted here. None of the commands in this chapter has
been re-run for this text, and `generate-data` is unseeded, so even a verified
run would not reproduce for you. The check to make is a
*comparison within your own run*: the labels printed in step 1 against the labels
selected in step 4. Recovering all three planted coefficients with no or few
extras is the expected outcome; recovering none of them means something is wrong
with the install, not with the statistics.
```

## Classification variant

Setting `--p-classification` makes the response binary, and the fit then goes
through `qiime classo classify` rather than `regress`:

```bash
qiime classo generate-data \
    --p-n 100 \
    --p-d 20 \
    --p-d-nonzero 3 \
    --p-classification \
    --o-x synthetic-class-x.qza \
    --o-c synthetic-class-c.qza \
    --verbose
```

`randomy.tsv` is overwritten again, this time with a boolean `col`, which QIIME 2
reads as a categorical metadata column — which is what `classify` requires,
since its `--m-y-column` is typed `MetadataColumn[Categorical]`.

```bash
qiime classo transform-features \
    --i-features synthetic-class-x.qza \
    --p-transformation clr \
    --p-coef 0.5 \
    --o-x synthetic-class-xclr.qza

qiime classo classify \
    --i-features synthetic-class-xclr.qza \
    --i-c synthetic-class-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-path \
    --p-cv \
    --p-cv-seed 1 \
    --p-stabsel \
    --p-stabsel-seed 1 \
    --p-lamfixed \
    --o-result synthetic-class-result.qza \
    --verbose
```

```{important}
`classify` has **no** `--p-concomitant` flag — the concomitant formulation is
unavailable for classification and the solver forces it off internally. Adding
`--p-concomitant` to the command above fails with an unrecognised-parameter
error. The robust option for classification is `--p-huber` with `--p-rho`, whose
default is `0.0` for `classify` (against `1.345` for `regress`).
```

## Generating data on a taxonomy

Passing `--i-taxa` makes the synthetic columns carry real feature IDs and
aggregates the planted coefficients onto taxonomic nodes, which is useful when
you want a toy problem shaped like a trac problem. This is the one part of the
chapter that is not download-free — it needs the Atacama taxonomy artifact from
[Download the Tutorial Data](../00_getting_started/03_download_data.md), assumed
here to sit in a `data/` directory next to `smoke-test/`:

```bash
qiime classo generate-data \
    --i-taxa ../data/classification.qza \
    --p-n 100 \
    --p-d 20 \
    --p-d-nonzero 3 \
    --o-x synthetic-taxa-x.qza \
    --o-c synthetic-taxa-c.qza \
    --verbose
```

The `--i-taxa` input is optional (its signature default is `None`). When it is
given, the first `d` tip names of the taxonomy replace the generated `A*`/`B*`
column names, and the support printed under `--verbose` is reported in terms of
*aggregated tree nodes* rather than leaves — so those labels will not always look
like feature IDs.

```{warning}
**Keep `--p-d` at or below the number of tips in your taxonomy.** When the
taxonomy has fewer tips than `d`, the tip names are written into a fixed-width
NumPy string array that was sized for the two- and three-character `A0`/`B0`
placeholders. Long feature IDs are then silently truncated to those few
characters, and the resulting column names are neither the real IDs nor unique.
There is no error and no warning.
```

## Where to go next

* [Data Preparation](02_data_preparation.md) — the same pipeline on the Atacama
  soil data.
* [Advanced: Choosing a Model](05_advanced/02_model_selection.md) — what the four
  model-selection blocks you switched on in step 3 actually compute, using this
  same synthetic problem.
* [Predict & Summarize](06_predict_and_summarize.md) — held-out prediction and
  the anatomy of the `.qzv`.
* [Troubleshooting](../90_reference/04_troubleshooting.md) — if any step above
  raised instead of running.
