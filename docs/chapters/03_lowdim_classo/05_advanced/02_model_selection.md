# Advanced: Choosing a Model

A constrained lasso fit is not one model, it is a family of models indexed by
the penalty λ. `qiime classo regress` and `qiime classo classify` offer four
independent ways to pick a member of that family — or to refuse to pick one and
show you the whole family instead.

The four are **not** alternatives you choose between at the command line. They
are four switches, each with its own parameter block, and each writes its own
tab into the `summarize` visualization. Turning all four on in a single call is
the normal thing to do: they answer different questions, and running them
together costs one solve each rather than four separate `regress` invocations.

This chapter puts all four on the same toy problem so the outputs are directly
comparable.

## The toy problem

The synthetic design from
[Synthetic Data with Known Truth](../01_generate_data.md) is used throughout,
because it is the only Tier-1 problem where the correct answer is known: three
planted non-zero coefficients out of twenty features. Every command below
assumes you are in the `smoke-test/` directory created there and have
`synthetic-xclr.qza`, `synthetic-c.qza` and `randomy.tsv` in it.

Everything transfers unchanged to the Atacama artifacts of
[Data Preparation](../02_data_preparation.md) — swap the three inputs and the
`--m-y-file` / `--m-y-column` pair.

## Which one to use

| Block | Flag | What it returns | Roughly what it costs | Reach for it when |
|---|---|---|---|---|
| Lambda path | `--p-path` | β at every λ on a grid; **no** chosen λ | one path solve | you want to see the whole trajectory, or to pick λ by held-out error via `predict` |
| Cross-validation | `--p-cv` | `lambda_min` and `lambda_1SE`, plus a refit β | `cv_subsets` path solves | prediction accuracy is the goal and you have samples to spare for splitting |
| Stability selection | `--p-stabsel` | a selection probability per feature, plus a refit β | `stabsel_b` subsample fits | *which features* is the question, and you want a frequency-based guard against false positives |
| Fixed lambda | `--p-lamfixed` | β at one λ, before and after refit | one solve | you already know λ, or want the theoretical λ as a fixed reference |

The short version: **cross-validation optimizes prediction, stability selection
optimizes selection, and they routinely disagree.** CV tends to keep a generous
support because extra weakly-informative features rarely hurt held-out error;
stability selection discards anything that is not reproducibly picked across
resamples. If you intend to interpret the coefficients biologically, trust the
stability-selection support and read the CV curve for context. If you intend to
predict, do the opposite.

## Lambda path — `--p-path`

```bash
qiime classo regress \
    --i-features synthetic-xclr.qza \
    --i-c synthetic-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-path \
    --p-path-nlam-log 40 \
    --p-path-lamin-log 0.01 \
    --p-path-n-active 0 \
    --p-no-cv \
    --p-no-stabsel \
    --p-no-lamfixed \
    --o-result ms-path.qza
```

| Parameter | Default | What it controls |
|---|---|---|
| `--p-path-nlam-log` | `40` | number of grid points |
| `--p-path-lamin-log` | `1e-2` | smallest λ, **as a fraction of λ_max** |
| `--p-path-n-active` | `0` | stop once this many coefficients are active; `0` disables early stopping |
| `--p-path-numerical-method` | `not specified` | solver; see [Numerical methods](#numerical-methods) |

The grid is always logarithmic. The plugin sets the path block's `logscale`
unconditionally, which is why there is no `--p-path-logscale` and why the two
grid parameters carry `_log` in their names. `--p-path-lamin-log 0.01` therefore
means "span two decades below λ_max", not "stop at λ = 0.01".

`--p-path-n-active` is the parameter that makes paths affordable on wide
problems. Setting it to, say, 30 tells the solver to stop descending once thirty
coefficients have entered the model — everything below that λ is dense and
usually uninteresting, and computing it is where the time goes. Leave it at `0`
for a twenty-feature toy.

**The path selects nothing.** On its own it is a diagnostic: you look at the
beta-path plot and see the order in which features enter and whether any
coefficient changes sign. It becomes a selection rule only in combination with
`qiime classo predict` on a held-out table, which adds a test-error-versus-λ
curve to the same tab. That combination is covered in
[Predict & Summarize](../06_predict_and_summarize.md), including the trap that
the curve is degenerate if you predict on the training samples.

## Cross-validation — `--p-cv`

```bash
qiime classo regress \
    --i-features synthetic-xclr.qza \
    --i-c synthetic-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-cv \
    --p-cv-subsets 5 \
    --p-cv-nlam 100 \
    --p-cv-lamin 0.001 \
    --p-cv-seed 1 \
    --p-cv-one-se \
    --p-no-path \
    --p-no-stabsel \
    --p-no-lamfixed \
    --o-result ms-cv.qza
```

| Parameter | Default | What it controls |
|---|---|---|
| `--p-cv-subsets` | `5` | number of CV folds |
| `--p-cv-nlam` | `100` | grid points on the CV λ grid |
| `--p-cv-lamin` | `1e-3` | smallest λ as a fraction of λ_max |
| `--p-cv-logscale` | on | log-spaced grid; `--p-no-cv-logscale` for linear |
| `--p-cv-seed` | `1` | fold assignment |
| `--p-cv-one-se` | on | use the one-standard-error λ instead of the error minimum |
| `--p-cv-numerical-method` | `not specified` | solver |

Two λ values come out of every CV run and **both are always reported** in the
Cross-Validation tab: `lambda_min`, the grid point with the lowest mean error,
and `lambda_1SE`, the largest λ whose mean error is still within one standard
error of that minimum. `--p-cv-one-se` decides only which of the two the
displayed refit and the exported `CV-beta.csv` are computed at.

Prefer the default. The error curve near its minimum is usually flat and noisy,
so `lambda_min` is chasing sampling noise; the one-standard-error rule buys a
noticeably sparser model at a cost in mean error that is, by construction, within
the noise. Switch to `--p-no-cv-one-se` when you are optimizing prediction and
have enough samples that the curve has a genuine minimum — which is what the
worked Atacama commands in [Log-Contrast Regression](../03_regression/01_logcontrast.md)
do.

Note that the CV grid and the path grid are configured separately and by default
do not agree: 100 points down to `1e-3` for CV against 40 points down to `1e-2`
for the path. If you want the CV-selected λ to land exactly on a point you can
read off the beta-path plot, you have to say so explicitly.

`--p-cv-seed` matters more than it looks on small problems. With a hundred-odd
samples split into five folds, each fold is small enough that a different fold
assignment perturbs the error curve, and the selected λ can land on a different
grid point — which changes the support you report. Pin the seed, and report it.

### `--p-cv-nlam` and the deprecated `--p-cv--nlam`

The number of CV grid points is spelled **`--p-cv-nlam`**. That is the parameter
to use.

You will also find `--p-cv--nlam`, with two dashes, in older scripts and in
published versions of these tutorials. It is not a typo: the parameter was
originally registered as `cv__nlam` with a double underscore, which QIIME 2
renders on the command line as `--p-cv--nlam`. The old spelling is still
accepted so that existing pipelines and any command already printed in a paper
keep working, but it emits a `DeprecationWarning` and will be removed.

The resolution rule, if both are present:

* Both given with **different** values → `ValueError`, refusing to guess.
* One given → that one is used; the deprecated spelling additionally warns.

```{important}
There is one hole in that rule. The check treats `cv_nlam == 100` as "not set",
because 100 is its default and the framework cannot distinguish an explicit
`--p-cv-nlam 100` from an omitted one. So
`--p-cv-nlam 100 --p-cv--nlam 50` does **not** raise — it silently uses 50.
Pass only one of the two spellings, and make it `--p-cv-nlam`.
```

## Stability selection — `--p-stabsel`

```bash
qiime classo regress \
    --i-features synthetic-xclr.qza \
    --i-c synthetic-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-stabsel \
    --p-stabsel-method first \
    --p-stabsel-b 50 \
    --p-stabsel-q 10 \
    --p-stabsel-percent-ns 0.5 \
    --p-stabsel-threshold 0.7 \
    --p-stabsel-seed 1 \
    --p-no-path \
    --p-no-cv \
    --p-no-lamfixed \
    --o-result ms-stabsel.qza
```

Stability selection refits the model on `--p-stabsel-b` random subsamples, each
containing a `--p-stabsel-percent-ns` fraction of the samples, and records how
often each feature is picked. A feature is called selected when that frequency
clears `--p-stabsel-threshold`.

What "picked" means on a single subsample depends on `--p-stabsel-method`:

| `--p-stabsel-method` | Per-subsample rule | Parameters that apply |
|---|---|---|
| `first` (default) | run the path and take the first `q` variables to enter | `--p-stabsel-q` |
| `lam` | fit at one fixed λ and take its support | `--p-stabsel-lam`, `--p-stabsel-true-lam` |
| `max` | run the path down to `lamin` and take the `q` largest coefficients | `--p-stabsel-q`, `--p-stabsel-lamin` |

The parameters that do not apply to your chosen method are simply ignored, with
no diagnostic. Setting `--p-stabsel-lamin` while leaving the method at `first`
does nothing at all, and that is the single most common way to be confused by
this block.

| Parameter | Default | Notes |
|---|---|---|
| `--p-stabsel-b` | `50` | more subsamples, smoother probabilities, linearly more time |
| `--p-stabsel-q` | `10` | variables per subsample; the main sparsity dial for `first`/`max` |
| `--p-stabsel-percent-ns` | `0.5` | subsample size as a fraction of n |
| `--p-stabsel-threshold` | `0.7` | selection frequency required; lower it to widen the support |
| `--p-stabsel-threshold-label` | `0.4` | stored in the artifact but **unused by the QIIME 2 visualization** — no bar labels are drawn; it affects only c-lasso's own matplotlib output |
| `--p-stabsel-lam` | `-1.0` | method `lam` only; any value ≤ 0 means "use the theoretical λ" |
| `--p-stabsel-true-lam` | on | method `lam` only; on = the number you give is λ itself, `--p-no-stabsel-true-lam` = it is λ/λ_max in [0, 1] |
| `--p-stabsel-lamin` | `1e-2` | method `max` only |
| `--p-stabsel-seed` | *unset* | see below |

```{important}
`--p-stabsel-seed` has **no default**. Unlike `--p-cv-seed`, which defaults to
`1`, leaving the stability-selection seed out means the subsamples are drawn
pseudo-randomly and the selection probabilities — and therefore possibly the
selected support — differ between two otherwise identical runs. Always set it.
```

The threshold and `q` are the two dials worth turning. `--p-stabsel-threshold 0.7`
is conservative and is the right starting point; the worked Atacama commands
elsewhere in this tier relax it to `0.5`, which admits any feature picked in more
than half the subsamples and so gives a wider support. Report whichever you
used — a stability-selection support is meaningless without its threshold.

## Fixed lambda — `--p-lamfixed`

```bash
qiime classo regress \
    --i-features synthetic-xclr.qza \
    --i-c synthetic-c.qza \
    --m-y-file randomy.tsv \
    --m-y-column col \
    --p-lamfixed \
    --p-lamfixed-lam -1.0 \
    --p-lamfixed-true-lam \
    --p-no-path \
    --p-no-cv \
    --p-no-stabsel \
    --o-result ms-lamfixed.qza
```

| Parameter | Default | Notes |
|---|---|---|
| `--p-lamfixed-lam` | `-1.0` | any value ≤ 0 means "use the theoretical λ" |
| `--p-lamfixed-true-lam` | on | on = the number given is λ; `--p-no-lamfixed-true-lam` = it is λ/λ_max |
| `--p-lamfixed-numerical-method` | `not specified` | solver |

The theoretical λ is a distribution-derived penalty level that requires no
resampling at all — no folds, no subsamples, one solve. It is the cheapest
possible answer and a good sanity anchor: if CV and stability selection both
land somewhere wildly different from it, look at your data before believing
either.

This block reports **two** coefficient vectors: the raw penalized β at that λ,
and a refit β obtained by re-estimating on the selected support without the
penalty. The refit removes the shrinkage bias that makes penalized coefficients
systematically too small in magnitude, and it is what you should quote as an
effect size. The raw β is what tells you how strongly the penalty is biting.

## Formulation parameters that change what you are selecting

These are not model-selection parameters, but they change the model that
selection is applied to — and therefore change λ_max, and therefore change the
meaning of every λ/λ_max fraction above. Two runs with different formulations
are not comparable on the λ scale unless you work in true λ.

### `--p-intercept`

On by default. c-lasso prepends an unpenalized `intercept` entry to the
coefficient vector, so β has one more element than you have features. Three
consequences you will actually notice: the beta bar plots skip the first entry,
`qiime classo predict` fills that column with ones, and `--p-maxplot` needs to
be at least *features + 1* to avoid truncation — see
[Predict & Summarize](../06_predict_and_summarize.md).

Use `--p-no-intercept` only if you have centered the response yourself and want
the fit forced through the origin.

### `--p-do-yshift`

Off by default, and available on `regress` only — `classify` has no such
parameter. It centers the response by subtracting its mean before solving.

With `--p-intercept` on, this is close to a no-op: the intercept already absorbs
any location shift in y. It earns its keep in the other combination, with
`--p-no-intercept`, where an uncentered response forces the constrained
coefficients to carry the mean and the penalty then distorts the whole fit.

### `--p-rho`

Defaults to `1.345` on `regress` and `0.0` on `classify`. It has **no effect
unless `--p-huber` is also set** — it is the transition point of the Huber loss,
in residual units, below which the loss is quadratic and above which it is
linear. Small ρ means aggressive down-weighting of large residuals; large ρ
approaches ordinary least squares. The regression default of 1.345 is the
classical value giving 95% efficiency under Gaussian errors.

On `classify` the value is wired to the separate `rho_classification` setting of
the underlying solver, which is why its default differs. Robust classification
is the topic of
[Advanced: Concomitant Formulation](01_concomitant_formulation.md), which also
explains why `--p-concomitant` exists for `regress` and not for `classify`.

## Running all four at once

This is the command the rest of the tier uses, and the one to start from:

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
    --p-stabsel-threshold 0.7 \
    --p-lamfixed \
    --o-result ms-all.qza \
    --verbose
```

```bash
qiime classo summarize \
    --i-problem ms-all.qza \
    --p-maxplot 25 \
    --o-visualization ms-all.qzv
```

Each block adds a tab; each tab reports its own running time, so the `.qzv`
itself tells you which of the four dominated the cost on your data.

```{note}
No timings, λ values, selected supports or error curves are quoted anywhere in
this chapter. None of these commands has been re-run for this text, and the
synthetic problem is unseeded, so a verified run would not reproduce for you
regardless. Compare the four tabs *within your own run* against the planted
support printed by `generate-data`.
```

## Numerical methods

All four blocks take a `--p-*-numerical-method` parameter:
`--p-path-numerical-method`, `--p-cv-numerical-method`,
`--p-stabsel-numerical-method`, `--p-lamfixed-numerical-method`. The recognised
values are `Path-Alg`, `P-PDS`, `PF-PDS` and `DR`.

Leave them alone. The default is the literal string `not specified`, which the
solver reads as "choose for me" and resolves against the formulation you asked
for — and the formulation is a better guide to the right algorithm than intuition
is.

```{important}
Because the default is itself a non-canonical string, an unrecognised value here
is silently accepted rather than rejected: a misspelled `--p-cv-numerical-method
Path-Algo` does not raise, it just gives you the automatic choice. Neither plugin
declares `Choices()` on any string parameter. See
[Troubleshooting](../../90_reference/04_troubleshooting.md).
```

## Gotchas

* **`--p-stabsel-lamin` with the default method does nothing.** It applies only
  to `--p-stabsel-method max`. The same holds for `--p-stabsel-lam` and
  `--p-stabsel-true-lam`, which apply only to method `lam`.
* **`--p-stabsel-threshold-label` is not a selection threshold.** It changes
  nothing in the `.qzv`. Changing it will not change which features are
  selected, however much it looks like it should.
* **`--p-stabsel-seed` is unset by default** — stability selection is not
  reproducible until you set it.
* **A misspelled `--p-stabsel-method` fails at runtime, not at parse time.** No
  `Choices()` is declared.
* **`--p-no-cv-logscale` used to crash** with `NameError: name 'xGraph' is not
  defined`, from a typo on the linear-scale branch of the CV plot. It is fixed in
  current q2-classo; if you hit it, your install predates the fix.
* **λ/λ_max fractions are not portable across formulations.** Turning
  `--p-concomitant` or `--p-huber` on or off moves λ_max, so the same
  `--p-lamfixed-lam 0.1` with `--p-no-lamfixed-true-lam` is a different penalty.

## Next

* [Predict & Summarize](../06_predict_and_summarize.md) — how each of these four
  blocks is rendered, and what `predict` adds to them.
* [Interpretation](../07_interpretation.md) — turning a selected support into a
  statement about the data.
* [q2-classo Parameter Reference](../../90_reference/03_classo_parameters.md) —
  the complete parameter list with types and defaults.
