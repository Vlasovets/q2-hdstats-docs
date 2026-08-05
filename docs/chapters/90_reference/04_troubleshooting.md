# Troubleshooting & Known Gotchas

Every trap in one place. Each entry gives the symptom you will actually see, the
cause, and the workaround.

---

## Installation

### `conda env create` fails on `deblur` → `sortmerna`

```
package deblur-1.1.1 requires sortmerna 2.0, but none of the providers
can be installed
```

**Cause.** A defect in the upstream QIIME 2 2026.7 `linux-64` distribution file:
it pins `zlib=1.3.2`, while every `sortmerna` 2.0 build in bioconda requires
`zlib >=1.2.11,<1.3.0a0`. The two cannot be satisfied together.

**Workaround.** `deblur` is not used anywhere in this tutorial. Download the
environment file, delete the `deblur`, `q2-deblur` and `sortmerna` lines, and
create the environment from your edited copy.

### `qiime gglasso pca` fails with `TemplateNotFound`

**Cause.** A non-editable install built before the fix did not ship
`_pca/assets/index.html` — it was missing from `package_data`.

**Workaround.** Reinstall from a current checkout, or use `pip install -e .`.

### pip silently breaks a working environment

**Cause.** pip cannot see conda's pins and will happily install a wheel over
`numpy=2.4.2`.

**Workaround.** Always `pip install --no-deps` inside a QIIME 2 environment.
Check with `conda list numpy pandas scipy` afterwards.

---

## q2-gglasso

### `--p-rank` always fails

```
NotImplementedError: Explicit rank selection requires a GGLasso version that
exposes an explicit low-rank rank parameter...
```

**Cause.** No released GGLasso (up to and including 0.3.0) can fix the rank of
the low-rank component. It exposes only the continuous `mu1` penalty and reports
the achieved rank as an *output*. The parameter is registered but guarded so it
raises rather than being silently ignored.

**Workaround.** Size the low-rank block through `mu1`: a **larger `mu1` gives a
smaller rank**. Scout a small `mu1` grid and read the achieved rank out of the
solution — see [Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md).
Explicit rank selection becomes available when GGLasso PR #50
(`fix_latent_rank`) is merged and released.

### `qiime gglasso pca` crashes on a solution that has no low-rank part

**Cause.** `pca` reads `solution/lowrank_`, which only exists when the problem
was solved with `--p-latent True`.

**Workaround.** Run `solve-problem` with `--p-latent` first. A sparse-only (SGL)
solution cannot be used with `pca`.

### `qiime gglasso pca` crashes with `AttributeError` on metadata

**Cause.** `--m-sample-metadata-file` is optional in the signature but is
dereferenced unconditionally.

**Workaround.** Always pass `--m-sample-metadata-file`. Treat it as required.

### `transform-features` demands a taxonomy I don't have

**Cause.** `--i-taxonomy` is a required input but the function body never uses
it. This is a registration wart, not a real dependency.

**Workaround.** Pass any valid `FeatureData[Taxonomy]` artifact; its contents do
not affect the result.

### `build-groups` output cannot be fed to `solve-problem`

**Cause.** `build-groups` emits a `TensorData` **artifact**, but `solve-problem`
takes `group_array` as a `List[Int]` **parameter**. They do not chain through the
QIIME 2 type system.

**Workaround.** Export the artifact and pass the group index explicitly as
`--p-group-array`. See
[Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md).

### The `.qzv` opens as a blank white page

**Cause.** A bokeh version mismatch: the visualizer emits a bokeh 3 document
while the page loads a bokeh 2.4.3 runtime. There is no Python traceback — the
failure is entirely client-side.

**Workaround.** Fixed in current q2-gglasso, which injects version-matched bokeh
resources instead of hardcoding a CDN URL. If you see it, you are on an older
build. Confirm by unzipping the `.qzv` and grepping `index.html` for `bokeh-2.`.

### Model selection ran even though I only wanted a single fit

**Cause.** The solver performs model selection whenever at least one grid has
more than one value, so a single fit requires `lambda1` *and* `lambda2` — plus
`mu1` for a **latent** problem — to each resolve to exactly one value. Leaving
`--p-lambda2-min` / `--p-lambda2-max` unset does not count as one value: an unset
pair of bounds expands to the 5-point default `np.logspace(-1, -4, 5)` and warns
`Default values for lambda2 have been used.`. That alone is enough to switch
model selection on, which is why a lone `--p-lambda1-min 0.5 --p-lambda1-max 0.5`
still gives you a path. (`lambda2` is not used by the single-instance solver
otherwise; it only decides this branch.)

**Workaround.** Pin `λ₂` explicitly — pass the same value to
`--p-lambda2-min` and `--p-lambda2-max` — and do the same for `--p-mu1-min` /
`--p-mu1-max` on a latent problem. Check the result by looking for a
`modelselect_stats` group in the solution: a single fit has none, and `summarize`
reduces its statistics tab accordingly.

---

## q2-classo

### `--p-cv--nlam` has two dashes — is that a typo?

**No.** The parameter really was named `cv__nlam` with a double underscore, which
QIIME 2 renders as `--p-cv--nlam`.

It is now spelled **`cv_nlam`** (`--p-cv-nlam`). The old spelling still works but
emits a `DeprecationWarning`. Passing both with *different* values raises — but
only if `cv_nlam` was changed from its registered default of `100`. QIIME 2 fills
that default in whether or not you typed it, so `100` is indistinguishable from
"not given" and is treated as unset: `--p-cv-nlam 100 --p-cv--nlam 50` does *not*
raise, it silently uses 50.

### `qiime classo classify --p-concomitant True` is rejected

**Cause.** The concomitant formulation is not available for classification.
`classify` has no `concomitant` parameter, and the underlying solver forces
`formulation.concomitant = False` for classification problems.

**Workaround.** Use the Huber hinge loss: `--p-huber True`, tuned via `--p-rho`.
See [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md).

### `qiime classo --help` lists two actions called "regress"

**Cause.** `classify` was registered with `name="regress"`.

**Workaround.** Fixed in current q2-classo. If you see it, your install predates
the fix — the action still works, it is only mislabelled.

### Every plot pane in the `summarize` `.qzv` is blank

**Cause.** The plot-writing calls were commented out while the templates still
referenced 16 `<iframe>` files. No error was raised; the panes just rendered
empty.

**Workaround.** Fixed in current q2-classo. A regression test now asserts that
every iframe the templates reference has something that writes it.

### `NameError: name 'xGraph' is not defined`

**Cause.** A typo (`xGrpah`) on the non-log-scale branch, reachable via
`--p-cv-logscale False`.

**Workaround.** Fixed in current q2-classo.

---

## Both plugins

### A misspelled enum value fails at runtime, not on the command line

```
ValueError: Unknown transformation name, use clr and not 'clrr'
```

**Cause.** Neither plugin declares `Choices()` on its string parameters, so QIIME 2
accepts any string and the check happens inside the function. Affects
`--p-transformation` (`clr`/`mclr`), `--p-method` (`scaled`/`unscaled`),
`--p-reg` (`GGL`/`FGL`), `--p-path-scale` (`log`/`linear`) and the various
`--p-*-numerical-method` parameters.

**Workaround.** Check spelling against the
[parameter reference](02_gglasso_parameters.md). Note that the
`*_numerical_method` parameters default to the literal string
`"not specified"`, so a typo there is silently accepted.

### The CLR transform produces obviously wrong values

**Cause (historical).** Zero imputation used a chained pandas assignment that
becomes a silent no-op under Copy-on-Write, so no pseudo count was ever added.
The failure surfaced far from its origin, as
`AssertionError: Add pseudo count before using clr`.

**Workaround.** Fixed in current q2-gglasso.
