# q2-classo Parameter Reference

Every parameter of every registered `q2-classo` action, with its CLI spelling,
its type as registered with QIIME 2, its default as declared in the function
signature, and the chapter where the tutorial exercises it.

`regress` and `classify` carry most of the surface area — 36 and 34 registered
parameters respectively — and almost none of them describe the model itself.
They describe the four **model-selection procedures** c-lasso can run on the same
fitted path: PATH, CV, StabSel and LAMfixed. Each procedure has its own
independent prefix, its own on/off switch, and its own numerical method. Once you
see that structure the parameter list stops being intimidating: you are reading
five small tables, not one enormous one.

Use [Troubleshooting](04_troubleshooting.md) when a flag misbehaves, and
[Command Coverage Matrix](01_command_coverage.md) for which actions the tutorial
covers.

```{important}
These tables are maintained **by hand** today, which means they can drift away
from the plugin. The intended end state is to capture
`qiime classo <action> --help` into `docs/_data/help/classo-<action>.txt` at
build time, render it with `{literalinclude}`, and add a CI check that asserts
the documented parameter set equals the set registered in
`q2_classo/plugin_setup.py` and `q2_classo/_dict.py`. **That generation step is
not wired up yet.** Until it is, treat `--help` on your own install as the final
authority and open an issue if this page disagrees with it.
```

```{note}
Every name, type and default below was read directly from
`q2_classo/plugin_setup.py`, `q2_classo/_dict.py` and the function signatures in
`q2_classo/_func.py`. They have **not** yet been checked against captured
`--help` output from a QIIME 2 2026.7 build, because that environment does not
exist yet. Nothing on this page is a recorded command output.
```

## How to read the tables

QIIME 2 derives every CLI flag mechanically: a parameter `foo_bar` becomes
`--p-foo-bar`, an input `--i-foo-bar`, an output `--o-foo-bar`, a `Metadata`
parameter `--m-foo-bar-file`, and a `MetadataColumn` parameter the *pair*
`--m-foo-bar-file` plus `--m-foo-bar-column`. Underscores are not collapsed,
which is exactly why the deprecated `cv__nlam` surfaces as `--p-cv--nlam` with
two dashes.

Conventions in the **Notes** column:

- **deprecated** — still accepted, emits a `DeprecationWarning`, will be removed.
- **no `Choices()`** — the plugin does not constrain the accepted strings; see
  the last section.

The **Demonstrated in** column names the chapter that exercises the parameter,
or the closest chapter that discusses it. `not demonstrated` means the parameter
is never passed in a runnable command in this tutorial — it may still be
discussed in prose or listed in a chapter's own parameter table; it is not a
statement about whether the parameter works.

---

## `generate-data`

Builds a synthetic log-contrast problem with a known ground-truth support. Useful
for checking that your invocation works before pointing it at real data.

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `taxa` | `--i-taxa` | `FeatureData[Taxonomy]` | `None` | [Generate Data](../03_lowdim_classo/01_generate_data.md) | If given, feature labels are drawn from the taxonomy and a tree matrix is built. |
| `n` | `--p-n` | `Int` | `100` | [Generate Data](../03_lowdim_classo/01_generate_data.md) | Number of samples. |
| `d` | `--p-d` | `Int` | `80` | [Generate Data](../03_lowdim_classo/01_generate_data.md) | Number of features. |
| `d_nonzero` | `--p-d-nonzero` | `Int` | `5` | [Generate Data](../03_lowdim_classo/01_generate_data.md) | Size of the true support. |
| `classification` | `--p-classification` | `Bool` | `False` | [Generate Data](../03_lowdim_classo/01_generate_data.md) | When `True`, the response is binary. |
| `x` | `--o-x` | `FeatureTable[Design]` | required | [Generate Data](../03_lowdim_classo/01_generate_data.md) | |
| `c` | `--o-c` | `ConstraintMatrix` | required | [Generate Data](../03_lowdim_classo/01_generate_data.md) | Zero-sum constraint. |

```{warning}
**Known bug: `generate-data` writes `randomy.tsv` into your current working
directory.** The generated response vector is not returned as an artifact; it is
dumped as a side effect to a fixed relative path. Run this action from a
directory you do not mind writing into, and expect the file to be overwritten on
every call. The response you then pass to `regress` comes from that file.
```

---

## `transform-features`

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `features` | `--i-features` | `FeatureTable[Composition \| Frequency \| Design]` | required | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | |
| `transformation` | `--p-transformation` | `Str` | `clr` | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | **`clr` is the only accepted value.** Anything else raises `ValueError`. **No `Choices()`**. |
| `coef` | `--p-coef` | `Float` | `0.5` | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Pseudocount substituted for non-positive entries before the log. |
| `x` | `--o-x` | `FeatureTable[Design]` | required | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | |

Note that this is a *different* CLR implementation from
[`qiime gglasso transform-features`](02_gglasso_parameters.md): it takes a
`coef` rather than a `pseudo_count`, has no `mclr` option, has no metadata
handling, and preserves the orientation it is given, centring along rows — it
assumes samples in rows, which is what the QIIME 2 `FeatureTable` view supplies
and what `regress` expects. It does not transpose, so a feature-major input is
silently CLR-ed along the wrong axis. The two are not interchangeable.

---

## `add-taxa`

Replaces the feature matrix `log(X)` with `log(X)A`, where `A` encodes the
taxonomic tree — the change of basis that turns a plain log-contrast model into
`trac`.

| Name | CLI flag | Type | Notes |
|---|---|---|---|
| `features` | `--i-features` | `FeatureTable[Design \| Frequency]` | |
| `weights` | `--i-weights` | `Weights` | Optional; defaults to the all-ones vector. |
| `taxa` | `--i-taxa` | `FeatureData[Taxonomy]` | Converted to a tree internally. |
| `x` | `--o-x` | `FeatureTable[Design]` | Columns are now internal nodes plus leaves. |
| `aweights` | `--o-aweights` | `Weights` | Each node's weight divided by its number of leaves. |

**This action has no parameters.** The leaf-count rescaling of the weights is
what makes coarse taxonomic ranks comparable with fine ones, and it is not
optional.

---

## `add-covariates`

Appends metadata columns to the feature matrix and extends the constraint matrix
and weight vector to match. This is how a covariate is admitted to the model
*without* being subject to the zero-sum constraint.

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `features` | `--i-features` | `FeatureTable[Design \| Frequency]` | `None` (optional) | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Registered as optional, but required in practice — see the note below. |
| `c` | `--i-c` | `ConstraintMatrix` | `None` (optional) | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Genuinely optional; falls back to an all-ones constraint row. |
| `weights` | `--i-weights` | `Weights` | `None` | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Defaults to all ones, length = number of features. |
| `covariates` | `--m-covariates-file` | `Metadata` | **required** | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | The metadata file the columns come from. |
| `to_add` | `--p-to-add` | `List[Str]` | **required** | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Column names to append; repeat the flag once per column. |
| `rescale` | `--p-rescale` | `List[Bool]` | `None` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | One value per entry of `to_add`, or a length mismatch raises `ValueError`. Defaults to all `False`. |
| `w_to_add` | `--p-w-to-add` | `List[Float]` | `None` | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Penalty weight per added covariate; same length rule. Defaults to all `1.0`. |
| `new_features` | `--o-new-features` | `FeatureTable[Design]` | required | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | |
| `new_c` | `--o-new-c` | `ConstraintMatrix` | required | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | Zeros in the new columns, i.e. the added covariates are exempt from the constraint. |
| `new_w` | `--o-new-w` | `Weights` | required | [Data Preparation](../03_lowdim_classo/02_data_preparation.md) | |

```{note}
`--i-features` and `--i-c` are both registered as **optional** inputs (`--help`
prints `[optional]` for both), but only `--i-c` is genuinely optional: omitting
it falls back to an all-ones constraint row. Omitting `--i-features` fails with
`AttributeError: 'NoneType' object has no attribute 'columns'` rather than a CLI
usage error, so treat it as required. Only `--m-covariates-file` and
`--p-to-add` are marked `[required]`.
```

Two behaviours worth knowing before you build a design matrix:

**Categorical columns are one-hot expanded**, and the generated column labels
have the form `<name> = <value>` — with spaces around the equals sign. Those
labels are what appear in the `summarize` coefficient plots, so a categorical
covariate contributes several rows to the output, not one.

**`rescale` only applies to numeric columns**, where it centres the column and
divides by its norm. Since the log-contrast features are CLR values of order
one, an unrescaled covariate measured in metres or in cells per gram will
dominate or vanish relative to them under a shared penalty. Setting the
corresponding `w_to_add` entry is the other lever on the same problem.

---

## `regress`

Constrained sparse log-contrast **regression**. 36 registered parameters: 35
current plus one deprecated alias.

### Inputs and outputs

| Name | CLI flag | Type | Notes |
|---|---|---|---|
| `features` | `--i-features` | `FeatureTable[Design \| Frequency]` | Samples in rows. |
| `c` | `--i-c` | `ConstraintMatrix` | Defaults to the zero-sum constraint when omitted. |
| `weights` | `--i-weights` | `Weights` | Per-feature penalty weights. Shorter than the feature count is padded with ones; longer is truncated. |
| `result` | `--o-result` | `CLASSOProblem` | Zarr store with every model selection that was run. |

### Response and formulation

These six choose the loss and the constraint structure — the actual statistical
model. Everything after this section only affects how lambda is chosen.

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `y` | `--m-y-file` / `--m-y-column` | `MetadataColumn[Numeric]` | required | [Log-contrast Regression](../03_lowdim_classo/03_regression/01_logcontrast.md) | Samples missing from `y` are dropped after an inner join and reported. |
| `do_yshift` | `--p-do-yshift` | `Bool` | `False` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Centres `y` before fitting. Not present on `classify`. |
| `concomitant` | `--p-concomitant` | `Bool` | `True` | [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md) | Joint M-estimation of the noise level sigma. **Not present on `classify`.** |
| `huber` | `--p-huber` | `Bool` | `False` | [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md) | Robust loss; combine with `concomitant` for the Huber-concomitant formulation. |
| `rho` | `--p-rho` | `Float` | `1.345` | [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md) | Huber transition point. Only meaningful when `huber` is `True`. Note `classify` defaults this to `0.0`. |
| `intercept` | `--p-intercept` | `Bool` | `True` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Adds an unpenalised intercept, which appears in the output as a coefficient labelled `intercept`. |

The `concomitant` default of `True` is the one to be deliberate about. It
estimates sigma jointly with beta, which makes the selected lambda scale-free
with respect to the noise level — desirable, but it changes which numerical
methods are available and therefore what `*_numerical_method` can legally be.

### PATH parameters

Computes the full regularization path. On by default.

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `path` | `--p-path` | `Bool` | `True` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Registered default is `True` despite the parameter description saying `False`. |
| `path_numerical_method` | `--p-path-numerical-method` | `Str` | `not specified` | not demonstrated | `Path-Alg`, `P-PDS`, `PF-PDS` or `DR`. The default string means "let c-lasso choose". **No `Choices()`** — and because the default is itself an unrecognised literal, a typo is silently accepted. |
| `path_n_active` | `--p-path-n-active` | `Int` | `0` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Stop once this many variables are active. `0` means no early stop. |
| `path_nlam_log` | `--p-path-nlam-log` | `Int` | `40` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Number of lambdas on the log-spaced path. |
| `path_lamin_log` | `--p-path-lamin-log` | `Float` | `1e-2` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Smallest lambda as a fraction of lambda_max. |

### CV parameters

K-fold cross-validation over the path. On by default.

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `cv` | `--p-cv` | `Bool` | `True` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | |
| `cv_numerical_method` | `--p-cv-numerical-method` | `Str` | `not specified` | not demonstrated | As `path_numerical_method`. **No `Choices()`**. |
| `cv_seed` | `--p-cv-seed` | `Int` | `1` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Fold assignment seed. Fixed by default, so repeated runs agree — change it to check fold stability. |
| `cv_one_se` | `--p-cv-one-se` | `Bool` | `True` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Select lambda by the one-standard-error rule rather than the CV minimum. |
| `cv_subsets` | `--p-cv-subsets` | `Int` | `5` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Number of folds. |
| `cv_nlam` | `--p-cv-nlam` | `Int` | `100` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Lambdas on the CV path. **Current spelling.** |
| `cv__nlam` | `--p-cv--nlam` | `Int` | `None` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | **Deprecated** alias of `cv_nlam` — note the double underscore and the double dash. Still works, emits a `DeprecationWarning`. Passing both raises only when `cv_nlam` was changed from its default of 100 and the two values differ; `--p-cv-nlam 100 --p-cv--nlam 50` does not raise — it silently uses 50. |
| `cv_lamin` | `--p-cv-lamin` | `Float` | `1e-3` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Smallest lambda on the CV path. |
| `cv_logscale` | `--p-cv-logscale` | `Bool` | `True` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | Log-spaced CV path. |

```{note}
`--p-cv-nlam` and `--p-cv-lamin` are a separate grid from `--p-path-nlam-log`
and `--p-path-lamin-log`. Changing the PATH grid does not change what CV
searches over, which is a common source of confusion when the CV-selected lambda
does not appear on the plotted path.
```

```{warning}
**Open docs issue.** The only runnable commands that pass this grid's size use
the *deprecated* spelling: [Cross-Validation](../04_highdim_atacama/05_classo_cv.md)
invokes `--p-cv--nlam`, not `--p-cv-nlam`. Rule 5 of
[Command Coverage Matrix](01_command_coverage.md) says `--p-cv--nlam` must not
appear in a runnable command outside the quarantined chapters, so those two
commands should be migrated to `--p-cv-nlam`.
```

### StabSel parameters

Stability selection: refit on many subsamples and keep the features selected
often enough. On by default, and generally the most trustworthy of the four
procedures for microbiome data.

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `stabsel` | `--p-stabsel` | `Bool` | `True` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | |
| `stabsel_numerical_method` | `--p-stabsel-numerical-method` | `Str` | `not specified` | not demonstrated | As above. **No `Choices()`**. |
| `stabsel_seed` | `--p-stabsel-seed` | `Int` | `None` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Unlike `cv_seed`, this defaults to unset, so subsampling is **not** reproducible unless you set it. Set it for anything you intend to report. |
| `stabsel_method` | `--p-stabsel-method` | `Str` | `first` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | `first`, `lam` or `max`. Decides what is recorded per subsample. **No `Choices()`**. |
| `stabsel_lam` | `--p-stabsel-lam` | `Float` | `-1.0` | not demonstrated | Only used when `stabsel_method` is `lam`. A negative value means "use the theoretical lambda". |
| `stabsel_true_lam` | `--p-stabsel-true-lam` | `Bool` | `True` | not demonstrated | Only used when `stabsel_method` is `lam`. `True` = the value given is a real lambda; `False` = it is lambda/lambda_max in [0, 1]. |
| `stabsel_b` | `--p-stabsel-b` | `Int` | `50` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Number of subsamples. The main runtime knob of the whole action. |
| `stabsel_q` | `--p-stabsel-q` | `Int` | `10` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Variables selected per subsample. |
| `stabsel_percent_ns` | `--p-stabsel-percent-ns` | `Float` | `0.5` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Subsample size as a fraction of n. |
| `stabsel_lamin` | `--p-stabsel-lamin` | `Float` | `1e-2` | not demonstrated | Only used when `stabsel_method` is `max`. |
| `stabsel_threshold` | `--p-stabsel-threshold` | `Float` | `0.7` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Selection frequency above which a feature is reported as selected. |
| `stabsel_threshold_label` | `--p-stabsel-threshold-label` | `Float` | `0.4` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Frequency above which a feature is *labelled* in c-lasso's own matplotlib plot. Recorded in the artifact; not used by the QIIME 2 visualization. |

`stabsel_b`, `stabsel_q` and `stabsel_threshold` interact. Raising `q` makes each
subsample select more variables, which raises everyone's selection frequency, so
a `threshold` that was strict at `q=10` becomes permissive at `q=30`. Change one
at a time, and report all three alongside any selected feature set.

### LAMfixed parameters

A single fit at one fixed lambda. On by default.

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `lamfixed` | `--p-lamfixed` | `Bool` | `True` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | |
| `lamfixed_numerical_method` | `--p-lamfixed-numerical-method` | `Str` | `not specified` | not demonstrated | As above. **No `Choices()`**. |
| `lamfixed_lam` | `--p-lamfixed-lam` | `Float` | `-1.0` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | Negative means "use the theoretical lambda once it is computed". |
| `lamfixed_true_lam` | `--p-lamfixed-true-lam` | `Bool` | `True` | [Model Selection](../03_lowdim_classo/05_advanced/02_model_selection.md) | `True` = the value is a real lambda; `False` = it is lambda/lambda_max in [0, 1]. With `True` and `lam = -1`, the value becomes `n * theoretical_lam`. |

```{tip}
All four procedures default to on, so a bare `qiime classo regress` runs PATH,
CV, StabSel and LAMfixed in one call — StabSel refits the model `stabsel_b` times
(50 by default), so it usually dominates the cost. Turn off what you are not
going to read: `--p-cv False`, `--p-stabsel False`, `--p-lamfixed False`.
```

```{warning}
**`--help` misreports several of these defaults.** QIIME 2 renders the
registered description text verbatim, so `--help` contradicts its own
`[default: ...]` marker wherever the description in `q2_classo/_dict.py` has
drifted from the function signature. The marker is right and the prose is wrong
in every case below; the defaults in the tables above follow the signatures.

| Parameter | `--help` prose says | Signature default |
|---|---|---|
| `path` | `Default Value = False` | `True` |
| `cv` | `Default Value = False` | `True` |
| `lamfixed` | `Default Value = False` | `True` |
| `cv_seed` | `Default value : None` | `1` |
| `path_numerical_method` | `Default value : 'choose'` | `not specified` |
| `cv_numerical_method` | `Default value : 'choose'` | `not specified` |
| `stabsel_numerical_method` | `Default value : 'choose'` | `not specified` |
| `lamfixed_numerical_method` | `Default value : 'choose'` | `not specified` |
| `path_n_active` | `Dafault value : False` | `0` |
```

---

## `classify`

Constrained sparse **classification**. 34 registered parameters: 33 current plus
the same deprecated `cv__nlam` alias.

The PATH, CV, StabSel and LAMfixed blocks are **identical** to `regress` — same
names, same types, same defaults — so they are not repeated here. Read the four
tables above and substitute `classify` for `regress`.

### Inputs and outputs

Identical to `regress`: `--i-features`, `--i-c`, `--i-weights`, `--o-result`.

### What differs from `regress`

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `y` | `--m-y-file` / `--m-y-column` | `MetadataColumn[Categorical]` | required | [Log-contrast Classification](../03_lowdim_classo/04_classification/01_logcontrast.md) | **Categorical**, not numeric. Must be binary; a non-binary column is rejected. |
| `huber` | `--p-huber` | `Bool` | `False` | [Log-contrast Classification](../03_lowdim_classo/04_classification/01_logcontrast.md) | Huber hinge loss. |
| `rho` | `--p-rho` | `Float` | `0.0` | [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md) | **Different default from `regress`** (`1.345`). Also passed on a `classify` call in [Gut-to-Soil Regression](../05_metagenomics/01_gut_to_soil/03_regression.md). |
| `intercept` | `--p-intercept` | `Bool` | `True` | [Cross-Validation](../04_highdim_atacama/05_classo_cv.md) | |
| `do_yshift` | — | — | — | — | **Does not exist on `classify`.** Centring a categorical response is meaningless. |
| `concomitant` | — | — | — | — | **Does not exist on `classify`.** See the warning below. |

```{warning}
**`qiime classo classify --p-concomitant` does not exist.** The concomitant
formulation is unavailable for classification: the parameter is not registered,
and the solver forces `formulation.concomitant = False` for classification
problems regardless. Passing the flag is a CLI error, not a silently ignored
option.

If you wanted the robustness that motivated `concomitant`, use the Huber hinge
loss instead — `--p-huber True`, tuned via `--p-rho`. Note that `rho` defaults to
`0.0` here rather than `1.345` — it is wired to c-lasso's
`formulation.rho_classification`, a different field from `regress`'s
`formulation.rho`, and `0.0` is a legal value for it (c-lasso requires only that
it be strictly less than 1). The `--help` text for `classify` reports `1.345`
(`_dict.py:286`) and is wrong; the registered default is `0.0`, which is what
every worked `classify` command in this book runs with. See
[Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md).
```

---

## `predict`

| Name | CLI flag | Type | Notes |
|---|---|---|---|
| `features` | `--i-features` | `FeatureTable[Design \| Frequency]` | Columns are matched to the fitted model by label; an `intercept` column is synthesised. |
| `problem` | `--i-problem` | `CLASSOProblem` | Output of `regress` or `classify`. |
| `predictions` | `--o-predictions` | `CLASSOProblem` | One prediction set per model selection that was computed. |

**This action has no parameters.** It emits a prediction for every model
selection present in `problem`, so if you turned off CV and StabSel at fit time
you get correspondingly fewer prediction sets here.

---

## `summarize` (visualizer)

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `problem` | `--i-problem` | `CLASSOProblem` | required | [Predict and Summarize](../03_lowdim_classo/06_predict_and_summarize.md) | |
| `taxa` | `--i-taxa` | `FeatureData[Taxonomy]` | `None` | [Predict and Summarize](../03_lowdim_classo/06_predict_and_summarize.md) | Used to label coefficients taxonomically. |
| `predictions` | `--i-predictions` | `CLASSOProblem` | `None` | [Predict and Summarize](../03_lowdim_classo/06_predict_and_summarize.md) | Output of `predict`. Omit it and the prediction panes are simply absent. |
| `maxplot` | `--p-maxplot` | `Int` | `200` | [Predict and Summarize](../03_lowdim_classo/06_predict_and_summarize.md) | Maximum number of coefficients drawn in a StabSel profile or beta bar plot. Raise it on wide problems or the plot silently truncates. |

---

## String parameters have no validation

`q2_classo/plugin_setup.py` imports `Choices` but never applies it. QIIME 2
therefore accepts **any** string for the parameters below, and the check — if
there is one — happens inside the function:

| Flag | Accepted values | Behaviour on a typo |
|---|---|---|
| `--p-transformation` | `clr` only | `ValueError: Unknown transformation name, use clr and not '...'` |
| `--p-path-numerical-method` | `Path-Alg`, `P-PDS`, `PF-PDS`, `DR` | **Silently accepted** — the default is itself the unrecognised literal `not specified`, which means "choose automatically" |
| `--p-cv-numerical-method` | same set | Silently accepted |
| `--p-stabsel-numerical-method` | same set | Silently accepted |
| `--p-lamfixed-numerical-method` | same set | Silently accepted |
| `--p-stabsel-method` | `first`, `lam`, `max` | Reaches c-lasso; behaviour depends on the solver |

The `*_numerical_method` family is the dangerous one. Because the sentinel
default is not a valid method name, there is no way for the code to distinguish
"user asked for automatic selection" from "user misspelled `Path-Alg`" — both are
just strings that are not in the recognised set, and both fall through to
automatic selection. You will get an answer, and it will not be the method you
asked for.

---

## See also

- [q2-gglasso Parameter Reference](02_gglasso_parameters.md)
- [Troubleshooting & Known Gotchas](04_troubleshooting.md)
- [Command Coverage Matrix](01_command_coverage.md)
