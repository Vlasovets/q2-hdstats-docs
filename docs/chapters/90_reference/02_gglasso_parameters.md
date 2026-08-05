# q2-gglasso Parameter Reference

Every parameter of every registered `q2-gglasso` action, with its CLI spelling,
its type as registered with QIIME 2, its default as declared in the function
signature, and the chapter where the tutorial exercises it.

Use this page when you already know *what* you want to do and need the exact
flag. Use [Troubleshooting](04_troubleshooting.md) when a flag does something
other than what you expected, and
[Command Coverage Matrix](01_command_coverage.md) to see which actions the
tutorial covers at all.

```{important}
These tables are maintained **by hand** today, which means they can drift away
from the plugin. The intended end state is to capture
`qiime gglasso <action> --help` into `docs/_data/help/gglasso-<action>.txt` at
build time, render it with `{literalinclude}`, and add a CI check that asserts
the documented parameter set equals the set registered in
`q2_gglasso/plugin_setup.py` and `q2_gglasso/_dict.py`. **That generation step is
not wired up yet.** Until it is, treat `--help` on your own install as the final
authority and open an issue if this page disagrees with it.
```

```{note}
Every name, type and default below was read directly from
`q2_gglasso/plugin_setup.py`, `q2_gglasso/_dict.py` and the function signatures in
`q2_gglasso/_func.py`. They have **not** yet been checked against captured
`--help` output from a QIIME 2 2026.7 build, because that environment does not
exist yet. Nothing on this page is a recorded command output.
```

## How to read the tables

QIIME 2 derives every CLI flag mechanically from the registered name: a
parameter `foo_bar` becomes `--p-foo-bar`, an input becomes `--i-foo-bar`, an
output `--o-foo-bar`, and a `Metadata` parameter becomes `--m-foo-bar-file`.
Underscores in the registered name are *not* collapsed, which is why the
deprecated `cv__nlam` in the sibling plugin surfaces as `--p-cv--nlam`.

`List[...]` parameters are repeatable options on the command line: pass the flag
once per value (`--p-lambda1-path 1.0 --p-lambda1-path 0.1 ...`). Several
`List[...]` parameters here — the `lambda*_min`/`lambda*_max`/`mu1_*` bounds —
are registered as lists but consumed as scalars, so in practice you pass them
exactly once. (`n_samples` is the exception: pass it once per instance.)

Two conventions used in the **Notes** column:

- **always raises** — the parameter is registered and accepted by the CLI, but
  every call that sets it fails. It exists so the limitation is explicit rather
  than silently ignored.
- **no `Choices()`** — the plugin does not constrain the accepted strings, so a
  typo passes argument parsing and fails (or is silently accepted) inside the
  function. See the last section of this page.

The **Demonstrated in** column names the chapter that exercises the parameter,
or the closest chapter that discusses it. `not demonstrated` means the parameter
is never passed in a runnable command in this tutorial (it may still be
discussed in prose); it is not a statement about whether the parameter works.

---

## `transform-features`

Turns a count table into a CLR- or mCLR-transformed table suitable for
covariance estimation. Optionally joins scaled sample metadata onto the table as
extra rows, which is how continuous covariates enter the network.

### Inputs and outputs

| Name | CLI flag | Type | Required | Notes |
|---|---|---|---|---|
| `table` | `--i-table` | `FeatureTable[Composition \| Frequency \| Design]` | yes | |
| `taxonomy` | `--i-taxonomy` | `FeatureData[Taxonomy]` | yes | **Required but never read.** The function body ignores it entirely. |
| `transformed_table` | `--o-transformed-table` | `FeatureTable[Frequency]` | yes | Output is stored feature-major (p x n). |

### Parameters

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `sample_metadata` | `--m-sample-metadata-file` | `Metadata` | `None` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | Only read when `add_metadata` is `True`. |
| `transformation` | `--p-transformation` | `Str` | `clr` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | `clr` or `mclr`. **No `Choices()`** — a typo is a runtime `ValueError`. |
| `pseudo_count` | `--p-pseudo-count` | `Int` | `1` | [The 300-ASV Dataset](../04_highdim_atacama/01_data.md) | Zero replacement before the log. Only meaningful for `clr`; `mclr` handles zeros itself. |
| `scale_metadata` | `--p-scale-metadata` | `Bool` | `True` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | Standardises metadata columns before joining. |
| `add_metadata` | `--p-add-metadata` | `Bool` | `False` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | Appends metadata columns as additional variables of the network. Only **numeric** columns are used — `filter_columns(column_type='numeric')` silently drops categorical ones, so a categorical covariate never becomes a node. Missing values in the retained columns are filled with `0` behind a `Missing values are imputed with 0!` warning. |
| `keep_original_id` | `--p-keep-original-id` | `Bool` | `True` | [The 300-ASV Dataset](../04_highdim_atacama/01_data.md) | Rows are always reordered by **ascending** total abundance, at either setting, with ties broken on the feature ID so the order depends only on the table's contents. With `False` the index is additionally replaced by `ASV-1 … ASV-p`, so `ASV-1` is the *least* abundant feature. **`ASV-k` is a position, not an identifier** — it is only meaningful inside the artifact that defines it and must never be used as a join key across artifacts. Ties are common (209 of the 300 Atacama features share a total-abundance value), so before ties were broken deterministically the same table could yield different `ASV-k` assignments. Prefer the default. |

Two of these deserve more than a table row. `add_metadata` is the switch that
decides *what kind of question the network answers*: with it off you get taxon-taxon
conditional dependencies only; with it on, environmental covariates become nodes
and a taxon-covariate edge means the association survives conditioning on every
other taxon. `scale_metadata` matters because the CLR-transformed taxa and a raw
covariate such as elevation live on wildly different scales, and the graphical
lasso applies one global `lambda1` to all of them — leaving a covariate unscaled
effectively changes its penalty relative to everything else.

```{warning}
`--i-taxonomy` is required by the registration but unused by the implementation.
Pass any valid `FeatureData[Taxonomy]` artifact; its contents cannot change the
result. This is tracked in [Troubleshooting](04_troubleshooting.md).
```

---

## `build-groups`

Builds the bookkeeping array that a non-conforming Group Graphical Lasso needs
when the K instances do not share the same feature set.

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `tables` | `--i-tables` | `List[FeatureTable[Frequency]]` | required | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | Repeat the flag once per instance. |
| `check_groups` | `--p-check-groups` | `Bool` | `True` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | Validates the constructed overlap groups. |
| `group_array` | `--o-group-array` | `TensorData` | required | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | `(2, L, K)` index array: L overlap groups across K instances. |

```{warning}
**Known gap: `build-groups` does not chain into `solve-problem`.** This action
emits a `TensorData` *artifact*, while `solve-problem` accepts `group_array` as a
`List[Int]` *parameter*. The QIIME 2 type system cannot connect the two, so you
have to export the artifact and pass the indices explicitly with
`--p-group-array`. The workaround is written out in
[Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md).
```

---

## `calculate-covariance`

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `table` | `--i-table` | `FeatureTable[Frequency]` | required | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | Expects the p x n output of `transform-features`. |
| `method` | `--p-method` | `Str` | `scaled` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | `scaled` divides by the square root of the diagonal (a correlation matrix); `unscaled` returns the raw covariance. **No `Choices()`**. |
| `bias` | `--p-bias` | `Bool` | `True` | [The 300-ASV Dataset](../04_highdim_atacama/01_data.md) | `True` normalises by N, `False` by N-1. |
| `covariance_matrix` | `--o-covariance-matrix` | `PairwiseFeatureData` | required | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) | p x p, rounded to 10 decimals. |

`method` is not cosmetic. A single `lambda1` penalises every entry of the
precision matrix equally, so on an unscaled covariance the taxa with the largest
variance absorb the penalty budget. `scaled` is the default for that reason, and
it is what the tutorial uses throughout.

---

## `solve-problem`

The core solver. Twenty parameters, but they fall into five groups: how big the
problem is, the sparsity grid, the low-rank grid, the multiple-instance
settings, and model selection.

### Inputs and outputs

| Name | CLI flag | Type | Notes |
|---|---|---|---|
| `covariance_matrix` | `--i-covariance-matrix` | `PairwiseFeatureData` | p x p, positive semi-definite. |
| `solution` | `--o-solution` | `GGLassoProblem` | Zarr store holding the solution and the hyperparameters used. |

### Parameters

| Parameter | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `n_samples` | `--p-n-samples` | `List[Int]` | **required** | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) | Number of samples behind the covariance estimate; one value per instance. The only parameter with no default. |
| `lambda1_min` | `--p-lambda1-min` | `List[Float]` | `None` | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) | Lower bound of the sparsity grid. If unset while `lambda1_max` is set, falls back to `1e-3`. |
| `lambda1_max` | `--p-lambda1-max` | `List[Float]` | `None` | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) | Upper bound of the sparsity grid. If unset while `lambda1_min` is set, falls back to `1`. |
| `n_lambda1` | `--p-n-lambda1` | `Int` | `1` | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) | Grid points between min and max. Leaving it at `1` only limits the `lambda1` grid to one point; on its own it does not turn model selection off — see "When model selection actually runs" below. |
| `lambda1_path` | `--p-lambda1-path` | `List[Float]` | `None` | [Lambda Paths](../02_lowdim_gglasso/05_lambda_paths.md) | Explicit grid. Overrides `lambda1_min`/`lambda1_max`/`n_lambda1` *and* `path_scale`; values are used verbatim. |
| `lambda2_min` | `--p-lambda2-min` | `List[Float]` | `None` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | Across-instance penalty, MGL only. Same `1e-3` fallback. |
| `lambda2_max` | `--p-lambda2-max` | `List[Float]` | `None` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | Same `1` fallback. |
| `n_lambda2` | `--p-n-lambda2` | `Int` | `1` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | |
| `mu1_min` | `--p-mu1-min` | `List[Float]` | `None` | [Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) | Low-rank penalty. Only used when `latent` is `True`. |
| `mu1_max` | `--p-mu1-max` | `List[Float]` | `None` | [Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) | |
| `n_mu1` | `--p-n-mu1` | `Int` | `1` | [Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md) | |
| `mu1_path` | `--p-mu1-path` | `List[Float]` | `None` | [Lambda Paths](../02_lowdim_gglasso/05_lambda_paths.md) | Explicit low-rank grid; overrides `mu1_min`/`mu1_max`/`n_mu1` and `path_scale`. Ignored unless `latent` is `True`. |
| `path_scale` | `--p-path-scale` | `Str` | `log` | [Lambda Paths](../02_lowdim_gglasso/05_lambda_paths.md) | `log` or `linear` spacing when a grid is built from min/max/count. Ignored when an explicit path is given. **No `Choices()`** — an unrecognised value raises `ValueError`, but only once a grid is actually built from bounds: with every `lambda1`/`lambda2`/`mu1` bound unset the misspelling is never reached. |
| `latent` | `--p-latent` | `Bool` | `None` | [Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) | Switches on the sparse + low-rank decomposition. Required for `pca`. |
| `rank` | `--p-rank` | `Int` | `None` | [Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md) | **Always raises.** See the box below. |
| `weights` | `--p-weights` | `List[Str]` | `None` | [Adaptive Graphical Lasso](../02_lowdim_gglasso/04_adaptive_glasso.md) | `entry weight` pairs, applied element-wise to `lambda1`. Later entries override earlier ones, so order matters. SGL only. |
| `non_conforming` | `--p-non-conforming` | `Bool` | `None` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | Instances with differing feature sets; requires `group_array`. |
| `group_array` | `--p-group-array` | `List[Int]` | `None` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | A **parameter**, not an input — see the `build-groups` gap above. |
| `reg` | `--p-reg` | `Str` | `GGL` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) | `GGL` (group) or `FGL` (fused). MGL only. **No `Choices()`**. |
| `gamma` | `--p-gamma` | `Float` | `0.01` | [Model Selection](../04_highdim_atacama/02_model_selection.md) | eBIC parameter in [0, 1]. Larger values push model selection toward sparser solutions. |

```{warning}
**`--p-rank` always raises** — `ValueError` if `--p-latent` is not set,
`NotImplementedError` otherwise. No released GGLasso
(up to and including 0.3.0) can fix the rank of the low-rank component; it
exposes only the continuous `mu1` penalty and reports the achieved rank as an
*output*. The parameter is registered and guarded so that setting it fails
loudly instead of being ignored.

Size the low-rank block through `mu1` instead: **a larger `mu1` gives a smaller
rank.** Scout a small `mu1` grid and read the achieved rank out of the solution —
[Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md) walks through
this. Setting `rank` without `latent=True` raises `ValueError` before the
version check is even reached.
```

### When the defaults kick in

If you leave a grid completely unset, the solver substitutes a built-in path and
emits a warning rather than failing:

| Grid | Substituted when both bounds are unset | Warning |
|---|---|---|
| `lambda1` | `np.logspace(0, -4, 15)` | `Default values for lambda1 have been used.` |
| `lambda2` | `np.logspace(-1, -4, 5)` | `Default values for lambda2 have been used.` |
| `mu1` (only if `latent=True`) | `np.logspace(2, -1, 10)` | `Default values for mu1 have been used.` |

Setting only one of the two bounds does *not* trigger this. A missing lower
bound becomes `1e-3` and a missing upper bound becomes `1`, silently, which is
easy to mistake for a grid you specified yourself.

### When model selection actually runs

Model selection runs whenever **at least one grid ends up with more than one
value** — and the grids the solver sees are the ones *after* the defaults in the
previous section have been substituted. That substitution is what makes single
fits harder to get than the flags suggest:

- **Non-latent (SGL/MGL).** A run counts as a single fit only if `lambda1`
  collapses to one value **and** at least one of `--p-lambda2-min`/
  `--p-lambda2-max` is set (with `--p-n-lambda2 1`). With `lambda2` left unset
  the solver substitutes the 5-point default grid even for SGL, so
  `model_selection` is on and `modelselect_stats` is written — pinning
  `--p-lambda1-min`/`--p-lambda1-max` to the same value, or passing a
  one-element `--p-lambda1-path`, is not enough on its own.
- **Latent (`--p-latent True`).** `lambda1`, `lambda2` *and* `mu1` must all
  collapse to a single value. Supplying a `lambda1` range while forgetting `mu1`
  — or the reverse — produces a grid search when you expected one fit.

To check what happened, look for a `modelselect_stats` group in the solution: a
genuine single fit has none, and `summarize` omits the statistics accordingly.
Because of the `lambda2` substitution above, most non-latent runs *will* have
one.

```{note}
This describes current solver behaviour, not a design decision. Applying the
`lambda2` default only to MGL problems would be the cleaner fix, but that is a
behaviour change in `q2_gglasso/utils.py` and needs a maintainer decision; the
same inaccurate "n_lambda1 == 1 means a single fit" rule is also repeated as a
code comment in `q2_gglasso/_summarize/_visualizer.py`.
```

---

## `pca` (visualizer)

Projects samples onto the principal components of the **low-rank** component of
an SLR solution.

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `table` | `--i-table` | `FeatureTable[Frequency]` | required | [PCA of the Latent Space](../02_lowdim_gglasso/07_pca.md) | The transformed table, p x n. |
| `solution` | `--i-solution` | `GGLassoProblem` | required | [PCA of the Latent Space](../02_lowdim_gglasso/07_pca.md) | **Must have been solved with `--p-latent True`.** |
| `sample_metadata` | `--m-sample-metadata-file` | `Metadata` | `None` | [PCA of the Latent Space](../02_lowdim_gglasso/07_pca.md) | **Optional in the signature, required in practice** — the visualizer dereferences it unconditionally. |
| `n_components` | `--p-n-components` | `Int` | `3` | [Latent PCA](../04_highdim_atacama/04_latent_pca.md) | Number of components to plot. |
| `color_by` | `--p-color-by` | `Str` | `None` | [PCA of the Latent Space](../02_lowdim_gglasso/07_pca.md) | Name of a metadata column. Must exist in the file passed above. |

```{warning}
Two failure modes, both easy to hit:

- Running `pca` on a sparse-only (SGL) solution fails, because the visualizer
  reads `solution/lowrank_`, which only exists when `--p-latent True` was used.
- Omitting `--m-sample-metadata-file` crashes with an `AttributeError` despite
  the parameter being optional in the signature. Treat it as required.
```

---

## `summarize` (visualizer)

| Name | CLI flag | Type | Default | Demonstrated in | Notes |
|---|---|---|---|---|---|
| `solution` | `--i-solution` | `GGLassoProblem` | required | [Summarize](../02_lowdim_gglasso/08_summarize.md) | Works for sparse-only and latent solutions alike. |
| `width` | `--p-width` | `Int` | `1500` | [Summarize](../02_lowdim_gglasso/08_summarize.md) | Plot width in pixels. |
| `height` | `--p-height` | `Int` | `1500` | [Summarize](../02_lowdim_gglasso/08_summarize.md) | Plot height in pixels. |
| `label_size` | `--p-label-size` | `Str` | `5pt` | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) | A bokeh font-size string, e.g. `25pt`. Raise it for small networks; the default is sized for hundreds of nodes. |
| `n_cov` | `--p-n-cov` | `Int` | `None` | [Summarize](../02_lowdim_gglasso/08_summarize.md) | Number of *trailing* variables that are covariates rather than taxa. When set, the heatmaps cluster the taxon block and the covariate block separately instead of mixing them. Set this if you used `--p-add-metadata True`. |

---

## String parameters have no validation

Neither plugin declares `Choices()` on any string parameter. QIIME 2 therefore
accepts **any** string, and the check happens inside the function — if it
happens at all. In `q2-gglasso` this affects:

| Flag | Accepted values | Behaviour on a typo |
|---|---|---|
| `--p-transformation` | `clr`, `mclr` | `ValueError` inside `transform-features` |
| `--p-method` | `scaled`, `unscaled` | `ValueError` inside `calculate-covariance` |
| `--p-path-scale` | `log`, `linear` | `ValueError` inside grid construction — but only if at least one of the `lambda1`/`lambda2`/`mu1` bounds is set; with every grid unbounded (or driven entirely by an explicit `--p-lambda1-path`/`--p-mu1-path`) the misspelling is silently ignored |
| `--p-reg` | `GGL`, `FGL` | Reaches the solver; behaviour depends on the GGLasso build |
| `--p-label-size` | any CSS size string | No validation; a bad value degrades the plot silently |

The practical consequence is that a typo costs you a whole run rather than being
caught at argument-parsing time. On a 13-ASV problem that is a minor annoyance;
on the 300-ASV problem in Tier 2 it is the entire solve. Check the spelling
before you submit the job.

---

## See also

- [q2-classo Parameter Reference](03_classo_parameters.md)
- [Troubleshooting & Known Gotchas](04_troubleshooting.md)
- [Command Coverage Matrix](01_command_coverage.md)
