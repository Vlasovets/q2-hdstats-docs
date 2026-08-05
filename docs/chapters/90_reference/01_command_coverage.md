# Command Coverage Matrix

"Comprehensive coverage of every command" is a claim, and claims about
documentation rot quietly. This page exists to make the claim **checkable**: it
maps every registered action of both plugins, broken down by parameter group, to
the one chapter that owns it and to the chapters that reuse it.

Read it two ways. Forwards, it is a table of contents organised by command rather
than by narrative — if you want to know where `--p-stabsel-b` is explained, look
it up here. Backwards, it is the input to a CI check: if a new parameter is
registered in a plugin and nothing here claims it, the build should fail.

There are **6 q2-gglasso actions** and **8 q2-classo actions**. Every one appears
below, together with every **parameter** (`--p-`) and **input** (`--i-`, `--m-`)
group listed in [q2-gglasso Parameter Reference](02_gglasso_parameters.md) and
[q2-classo Parameter Reference](03_classo_parameters.md).

```{note}
**Outputs are deliberately out of scope here.** Each action has a fixed output
signature that does not vary by chapter, so tracking `--o-` flags row by row adds
noise without adding coverage. The one exception below (`--o-group-array`) is
listed because the artifact it produces cannot be fed back into `solve-problem`
through the type system, which is a fact about the workflow rather than about the
output itself. Full output signatures live in the two reference pages, and CI
check 2 below is scoped to parameters accordingly.
```

## How to read the matrix

Each row is an **action x parameter group** — not an action and not a single
parameter. Grouping is what makes the table usable: `--p-stabsel-b`,
`--p-stabsel-q` and `--p-stabsel-threshold` interact so tightly that documenting
them apart would be worse than useless, so they share a row and a primary
chapter.

Two columns carry the traceability:

- **Primary** — the single chapter that *introduces* the group, explains why the
  parameters matter, and is responsible for keeping the explanation correct.
  Exactly one chapter per row. If you are fixing an error about a parameter, fix
  it here first.
- **Also in** — chapters that use the group again at a different scale or on a
  different dataset without re-explaining it. These should link back to the
  primary rather than duplicate it.

The book's structure intends tier 1 to be the reference tier: every action gets
its canonical demonstration on the 13-ASV Atacama table, and the later tiers
introduce new *values* and new *questions*, not new commands. Where a row's
primary chapter is not in tier 1, that is a deliberate exception or an outstanding
gap — {ref}`coverage-debt` lists them all.

### Chapter keys

The matrix uses short keys so the tables stay narrow.

| Key | Chapter |
|---|---|
| `DL` | [Download the Tutorial Data](../00_getting_started/03_download_data.md) |
| `INST-GG` | [Installing q2-gglasso](../01_installation/02_q2_gglasso.md) |
| `INST-CL` | [Installing q2-classo](../01_installation/03_q2_classo.md) |
| `VERIFY` | [Verifying Your Installation](../01_installation/04_verify.md) |
| `G-PREP` | [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md) |
| `G-SGL` | [Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md) |
| `G-SLR` | [Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) |
| `G-ADAPT` | [Adaptive Graphical Lasso](../02_lowdim_gglasso/04_adaptive_glasso.md) |
| `G-PATH` | [Regularization Paths & Model Selection](../02_lowdim_gglasso/05_lambda_paths.md) |
| `G-MGL` | [Multiple Graphical Lasso](../02_lowdim_gglasso/06_multiple_graphical_lasso.md) |
| `G-PCA` | [Latent-Component PCA](../02_lowdim_gglasso/07_pca.md) |
| `G-SUM` | [Summarizing a Solution](../02_lowdim_gglasso/08_summarize.md) |
| `G-INT` | [Interpretation (gglasso)](../02_lowdim_gglasso/09_interpretation.md) |
| `C-GEN` | [Synthetic Data with Known Truth](../03_lowdim_classo/01_generate_data.md) |
| `C-PREP` | [Data Preparation (classo)](../03_lowdim_classo/02_data_preparation.md) |
| `C-REG` | [Log-Contrast Regression](../03_lowdim_classo/03_regression/01_logcontrast.md) |
| `C-RTRAC` | [Regression with trac](../03_lowdim_classo/03_regression/02_trac.md) |
| `C-CLF` | [Log-Contrast Classification](../03_lowdim_classo/04_classification/01_logcontrast.md) |
| `C-CTRAC` | [Classification with trac](../03_lowdim_classo/04_classification/02_trac.md) |
| `C-CONC` | [Concomitant Formulation](../03_lowdim_classo/05_advanced/01_concomitant_formulation.md) |
| `C-MSEL` | [Choosing a Model](../03_lowdim_classo/05_advanced/02_model_selection.md) |
| `C-PRED` | [Predict & Summarize](../03_lowdim_classo/06_predict_and_summarize.md) |
| `C-INT` | [Interpretation (classo)](../03_lowdim_classo/07_interpretation.md) |
| `H-DATA` | [The 300-ASV Dataset](../04_highdim_atacama/01_data.md) |
| `H-LAM` | [Selecting lambda](../04_highdim_atacama/02_model_selection.md) |
| `H-RANK` | [Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md) |
| `H-PCA` | [Latent Components & Covariates](../04_highdim_atacama/04_latent_pca.md) |
| `H-CV` | [Log-Contrast Models at Scale](../04_highdim_atacama/05_classo_cv.md) |
| `H-INT` | [Interpretation (tier 2)](../04_highdim_atacama/06_interpretation.md) |
| `M-DATA` | [Gut-to-Soil: The Dataset](../05_metagenomics/01_gut_to_soil/01_data.md) |
| `M-NET` | [Gut-to-Soil: Network Inference](../05_metagenomics/01_gut_to_soil/02_network.md) |
| `M-REG` | [Gut-to-Soil: Log-Contrast Regression](../05_metagenomics/01_gut_to_soil/03_regression.md) |
| `R-GG` | [q2-gglasso Parameter Reference](02_gglasso_parameters.md) |
| `R-CL` | [q2-classo Parameter Reference](03_classo_parameters.md) |
| `R-TS` | [Troubleshooting & Known Gotchas](04_troubleshooting.md) |

---

## q2-gglasso

### `transform-features`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Input table | `--i-table` | `G-PREP` | `G-ADAPT`, `G-MGL`, `H-DATA`, `M-DATA`, `VERIFY` |
| Choice of transform | `--p-transformation` | `G-PREP` | `G-ADAPT`, `G-MGL`, `H-DATA`, `M-DATA`, `VERIFY` |
| Zero handling | `--p-pseudo-count` | `H-DATA` | — |
| Metadata as network nodes | `--m-sample-metadata-file`, `--p-add-metadata`, `--p-scale-metadata` | `G-PREP` | `G-ADAPT`, `G-MGL`, `H-DATA` |
| Feature relabelling | `--p-keep-original-id` | `H-DATA` | `G-INT` |
| The unused required input | `--i-taxonomy` | `G-PREP` | `DL`, `R-TS` |

`--p-transformation` and the metadata switches are the two decisions that change
what the network *means*, which is why both are owned by `G-PREP` rather than
being scattered. `--i-taxonomy` gets a row of its own because it is a trap, not a
parameter: it is required by the registration and never read by the function, so
readers hit it before they hit anything statistical.

### `build-groups`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Instance tables and validation | `--i-tables`, `--p-check-groups`, `--o-group-array` | `G-MGL` | — |
| The `TensorData` -> `List[Int]` gap | (export workaround) | `G-MGL` | `R-GG`, `R-TS` |

```{warning}
`build-groups` emits a `TensorData` **artifact** while `solve-problem` accepts
`group_array` as a `List[Int]` **parameter**. They do not chain through the QIIME 2
type system, so the only way to connect them is to export the artifact and pass
the indices by hand with `--p-group-array`. This is a **known gap**, not a
documentation shortcut, and `G-MGL` is the one chapter responsible for spelling
out the workaround.
```

### `calculate-covariance`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Input table (the p x n transform output) | `--i-table` | `G-PREP` | `G-ADAPT`, `G-MGL`, `H-DATA`, `M-NET`, `VERIFY` |
| Scaling | `--p-method` | `G-PREP` | `G-ADAPT`, `G-MGL`, `H-DATA`, `M-NET`, `VERIFY` |
| Normalisation denominator | `--p-bias` | `H-DATA` | — |

### `solve-problem`

The largest surface in either plugin, and the reason this page is organised by
group. Its twenty parameters divide into eight concerns, and no chapter tries to
cover more than two of them at once.

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Input covariance | `--i-covariance-matrix` | `G-SGL` | `G-SLR`, `G-ADAPT`, `G-PATH`, `G-MGL`, `H-LAM`, `H-RANK`, `M-NET` |
| Problem size | `--p-n-samples` | `G-SGL` | `G-SLR`, `G-ADAPT`, `G-PATH`, `G-MGL`, `H-LAM`, `H-RANK`, `M-NET` |
| Sparsity grid | `--p-lambda1-min`, `--p-lambda1-max`, `--p-n-lambda1` | `G-SGL` | `G-SLR`, `G-ADAPT`, `G-PATH`, `G-MGL`, `H-LAM`, `H-RANK`, `M-NET` |
| Explicit grids and spacing | `--p-lambda1-path`, `--p-mu1-path`, `--p-path-scale` | `G-PATH` | `G-MGL`, `H-LAM`, `H-RANK` |
| Low-rank block | `--p-latent`, `--p-mu1-min`, `--p-mu1-max`, `--p-n-mu1` | `G-SLR` | `G-PATH`, `G-PCA`, `H-LAM`, `H-RANK`, `H-PCA` |
| Explicit rank (always raises) | `--p-rank` | `H-RANK` | `G-PATH`, `G-PCA`, `R-TS` |
| Adaptive penalty weights | `--p-weights` | `G-ADAPT` | — |
| Multiple instances | `--p-reg`, `--p-lambda2-min`, `--p-lambda2-max`, `--p-n-lambda2`, `--p-non-conforming`, `--p-group-array` | `G-MGL` | `G-PATH`, `VERIFY` |
| Model-selection criterion | `--p-gamma` | `G-PATH` | `G-SGL`, `G-SLR`, `G-ADAPT`, `G-MGL`, `H-LAM`, `H-RANK` |

```{important}
Two behaviours cut across the grid rows above and must be repeated wherever a
grid is set, because getting them wrong produces a plausible-looking result
rather than an error.

**Defaults appear silently.** Leaving a grid entirely unset substitutes a
built-in path and emits a warning; setting only one bound substitutes the other
one with no warning at all.

**Model selection runs only if at least one grid has more than one value** — and
for a latent problem, `lambda1`, `lambda2` *and* `mu1` must all be singletons
before the run counts as a single fit. `G-PATH` owns the full explanation;
`H-LAM` and `H-RANK` reuse it.
```

`--p-rank` is the one row whose primary chapter exists mainly to talk you out of
the parameter. It **always raises** on every released GGLasso up to and including
0.3.0 — `ValueError` if `--p-latent` is not set, `NotImplementedError`
otherwise — because no release can fix the rank of the low-rank component.
`H-RANK` therefore documents the alternative — steer the rank through `mu1`,
where a **larger `mu1` gives a smaller rank** — and reads the achieved rank back
out of the solution.

### `pca` (visualizer)

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Inputs and the required metadata file | `--i-table`, `--i-solution`, `--m-sample-metadata-file` | `G-PCA` | `H-PCA`, `R-TS` |
| Projection and colouring | `--p-n-components`, `--p-color-by` | `G-PCA` | `H-PCA` |

```{warning}
`pca` has two prerequisites that the signature does not state.

The solution must have been produced with `--p-latent True` — the visualizer
reads `solution/lowrank_`, which a sparse-only SGL solution does not have.

And `--m-sample-metadata-file` is **optional in the signature but required in
practice**: it is dereferenced unconditionally, so omitting it crashes with an
`AttributeError`. `G-PCA` states both before its first command.
```

### `summarize` (visualizer)

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Input solution | `--i-solution` | `G-SUM` | `G-SGL`, `G-SLR`, `G-ADAPT`, `G-PATH`, `H-RANK` |
| Label sizing | `--p-label-size` | `G-SUM` | `G-SGL`, `G-SLR`, `G-ADAPT`, `G-PATH`, `H-RANK` |
| Canvas size | `--p-width`, `--p-height` | `G-SUM` | `H-RANK` |
| Covariate block separation | `--p-n-cov` | `G-SUM` | `G-ADAPT`, `G-INT`, `H-DATA` |

`--p-n-cov` is paired with `--p-add-metadata`: it tells the heatmaps how many
*trailing* variables are covariates rather than taxa, so the two blocks cluster
separately. Anyone who turned metadata into nodes in `G-PREP` needs it here.

---

## q2-classo

### `generate-data`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Problem shape | `--p-n`, `--p-d`, `--p-d-nonzero` | `C-GEN` | — |
| Response type | `--p-classification` | `C-GEN` | — |
| Taxonomy-derived labels and tree | `--i-taxa` | `C-GEN` | — |
| The `randomy.tsv` side effect | (no flag) | `C-GEN` | `R-CL` |

```{warning}
`generate-data` writes `randomy.tsv` into the **current working directory** — the
generated response is not returned as an artifact. It is overwritten on every
call. `C-GEN` is the only chapter that runs this action, and it says where to run
it from.
```

### `transform-features`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Input features | `--i-features` | `C-PREP` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `M-REG` |
| CLR transform and pseudocount | `--p-transformation`, `--p-coef` | `C-PREP` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `M-REG` |

This is a **different** implementation from `qiime gglasso transform-features`:
`coef` rather than `pseudo_count`, no `mclr`, no metadata handling, and a
sample-major output because `regress` wants samples in rows. `C-PREP` says so
explicitly, because the shared action name is the single most common source of
confusion between the two plugins.

### `add-taxa`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Tree change of basis | `--i-features`, `--i-weights`, `--i-taxa` (no parameters) | `C-RTRAC` | `C-PREP`, `C-CTRAC`, `C-PRED` |

### `add-covariates`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Inputs, column selection and one-hot expansion | `--i-features`, `--i-c`, `--i-weights`, `--m-covariates-file`, `--p-to-add` | `C-PREP` | `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `H-CV` |
| Per-covariate penalty weight | `--p-w-to-add` | `C-PREP` | `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `H-CV` |
| Rescaling numeric covariates | `--p-rescale` | `H-CV` | — |

Categorical columns are expanded to one-hot indicators labelled
`<name> = <value>`, spaces included, and those labels are what appear in the
`summarize` coefficient plots — so one categorical covariate contributes several
rows to the output. `C-PREP` owns that fact.

### `regress`

The four model-selection procedures — PATH, CV, StabSel, LAMfixed — are all **on
by default**, each has its own prefix, and each has its own numerical method. That
structure is why the rows below look repetitive: they are genuinely four parallel
blocks over the same fitted path.

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Inputs | `--i-features`, `--i-c`, `--i-weights` | `C-REG` | `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-CONC`, `C-MSEL`, `H-CV`, `M-REG` |
| Numeric response | `--m-y-file`, `--m-y-column` | `C-REG` | `C-RTRAC`, `C-CONC`, `C-MSEL`, `H-CV`, `M-REG` |
| Response shift | `--p-do-yshift` | `C-MSEL` | `H-CV` |
| Intercept | `--p-intercept` | `C-MSEL` | `C-PRED`, `H-CV` |
| Loss and noise model | `--p-concomitant`, `--p-huber`, `--p-rho` | `C-CONC` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-PRED`, `H-CV` |
| PATH | `--p-path`, `--p-path-nlam-log`, `--p-path-lamin-log`, `--p-path-n-active`, `--p-path-numerical-method` | `C-MSEL` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-PRED`, `H-CV` |
| CV | `--p-cv`, `--p-cv-subsets`, `--p-cv-nlam`, `--p-cv-lamin`, `--p-cv-logscale`, `--p-cv-one-se`, `--p-cv-seed`, `--p-cv-numerical-method` | `C-MSEL` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-PRED`, `H-CV` |
| Deprecated CV alias | `--p-cv--nlam` | `C-MSEL` | `R-CL`, `R-TS` |
| StabSel | `--p-stabsel`, `--p-stabsel-b`, `--p-stabsel-q`, `--p-stabsel-threshold`, `--p-stabsel-threshold-label`, `--p-stabsel-seed`, `--p-stabsel-method`, `--p-stabsel-lam`, `--p-stabsel-true-lam`, `--p-stabsel-lamin`, `--p-stabsel-percent-ns`, `--p-stabsel-numerical-method` | `C-MSEL` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-PRED`, `H-CV` |
| LAMfixed | `--p-lamfixed`, `--p-lamfixed-lam`, `--p-lamfixed-true-lam`, `--p-lamfixed-numerical-method` | `C-MSEL` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `H-CV` |

```{note}
`--p-cv--nlam` — two dashes — is not a typo in this book. The parameter was
originally registered as `cv__nlam` with a double underscore, which QIIME 2
renders literally. **`--p-cv-nlam` is the current spelling**; the old one still
works and emits a `DeprecationWarning`. New commands should use `--p-cv-nlam`,
and the deprecated form should appear only where it is being explained.
```

### `classify`

`classify` shares the PATH, CV, StabSel and LAMfixed blocks with `regress` — same
names, same defaults — so those rows are not repeated. Only the differences are
owned separately.

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Categorical response | `--m-y-file`, `--m-y-column` | `C-CLF` | `C-GEN`, `C-CTRAC`, `C-CONC` |
| Hinge loss and its transition point | `--p-huber`, `--p-rho` | `C-CLF` | `C-CTRAC`, `C-CONC` |
| Intercept | `--p-intercept` | `C-CLF` | `C-MSEL` |
| What `classify` does **not** have | (`--p-concomitant`, `--p-do-yshift`) | `C-CONC` | `C-CLF`, `R-CL`, `R-TS` |
| Selection procedures | as `regress` | `C-MSEL` | `C-CLF`, `C-CTRAC` |

```{warning}
**`qiime classo classify --p-concomitant` does not exist.** The parameter is not
registered on `classify`, and the solver forces the concomitant formulation off
for classification problems regardless. Passing the flag is a command-line error.
Use the Huber hinge loss instead — `--p-huber True` with an explicit `--p-rho`,
since `rho` defaults to `0.0` here rather than the `1.345` used by `regress`.
`C-CONC` owns this comparison.
```

### `predict`

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Inputs (no parameters) | `--i-features`, `--i-problem` | `C-PRED` | `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-MSEL` |

`predict` emits one prediction set per model selection present in the problem, so
switching CV or StabSel off at fit time silently reduces what you get here. That
coupling is `C-PRED`'s to explain.

### `summarize` (visualizer)

| Parameter group | Flags | Primary | Also in |
|---|---|---|---|
| Inputs | `--i-problem`, `--i-taxa`, `--i-predictions` | `C-PRED` | `C-GEN`, `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-MSEL` |
| Plot truncation | `--p-maxplot` | `C-PRED` | `C-GEN`, `C-MSEL` |

---

## Supporting QIIME 2 commands

The tutorial does not run in a vacuum. These commands are not part of either
plugin, but a reader who skips them cannot complete the chapters, so they get the
same treatment: one owning chapter each.

| Command | What the tutorial uses it for | Primary | Also in |
|---|---|---|---|
| `qiime sample-classifier split-table` | Train/test split before `regress` or `classify` (`--p-test-size`, `--p-random-state`, `--p-stratify`) | `C-PREP` | `C-REG`, `C-RTRAC`, `C-CLF`, `C-CTRAC`, `C-PRED` |
| `qiime feature-table filter-features` | Restricting a table to a shared feature set before building multiple graphical-lasso instances | `G-MGL` | — |
| `qiime feature-table filter-samples` | Splitting one table into the K per-group instances (`--p-where`) | `G-MGL` | — |
| `qiime feature-table summarize` | Reading off the sample count that `--p-n-samples` needs | `G-PATH` | — |
| `qiime metadata tabulate` | Inspecting a grouping variable before splitting on it | `G-MGL` | — |
| `qiime tools export` | The `build-groups` workaround, and reading the achieved rank out of a solution | `G-MGL` | `H-RANK` |
| `qiime tools view` | Opening a `.qzv` without a browser round-trip | `G-SUM` | `G-PCA` |
| `qiime tools peek` | Confirming a downloaded artifact's type and UUID | `DL` | `VERIFY` |
| `qiime dev refresh-cache` | Making a freshly installed plugin visible to the CLI | `INST-GG` | `INST-CL`, `VERIFY` |
| `qiime gglasso --help`, `qiime classo --help` | The authoritative parameter list on *your* install | `VERIFY` | `INST-GG`, `INST-CL` |

```{tip}
`qiime feature-table summarize` earns its place because `--p-n-samples` is the
only `solve-problem` parameter with no default, and the value you give is passed
straight through as the sample size `N` of the underlying problem — the same `N`
that the model-selection criterion is computed against. Read it off the table
rather than from memory.
```

(coverage-debt)=
## Coverage debt

Rows whose primary chapter is not in tier 1, i.e. where the "tier 1 is the
reference tier" rule is currently broken. These are tracked deliberately; each is
either a justified exception or work outstanding.

| Group | Current primary | Why, or what is missing |
|---|---|---|
| `--p-pseudo-count` | `H-DATA` | Not exercised in tier 1. Zero handling only becomes visible on a sparse table, but tier 1 should still name it. **Outstanding.** |
| `--p-keep-original-id` | `H-DATA` | The tier 1 table has only 13 features, so relabelling MD5 hashes to `ASV-n` buys little; relabelling 300 features is where it earns its place. Justified. |
| `--p-bias` | `H-DATA` | `N` versus `N-1` changes nothing structural about the network. It is discussed once, in the chapter where the covariance estimate itself is under scrutiny. Justified. |
| `--p-rescale` | `H-CV` | Not exercised in tier 1, although the tier 1 covariates (`elevation`, `ph`) are exactly the case that needs it. **Outstanding.** |
| `--p-rank` | `H-RANK` | The parameter always raises, so its owning chapter is the one about choosing a rank the working way. Justified. |

```{note}
Most of this book has not yet been re-run against QIIME 2 2026.7, so the "Also in"
columns record where a flag is *written*, not where it has been *observed to
work*. The distinction disappears once the CI checks below run.
```

## How this page gets checked

The two parameter reference pages are meant to be **generated, not written**:
`qiime gglasso <action> --help` and `qiime classo <action> --help` captured into
`docs/_data/help/<plugin>-<action>.txt` at build time and rendered with
`{literalinclude}`. Once the flag list on those pages comes from the plugin rather
than from a human, this matrix becomes machine-checkable, because both sides of
every comparison are then mechanical.

The checks worth wiring, in rough order of value:

1. **No invented flags.** Every `--p-`, `--i-`, `--o-` and `--m-` token appearing
   in a fenced `bash` block anywhere under `docs/chapters/` must appear in one of
   the captured help files. This is the check that catches a documented parameter
   that does not exist.
2. **No missing parameters.** Every `--p-`, `--i-` and `--m-` flag in a captured
   help file must appear in at least one row of this matrix. This is the check
   that turns "comprehensive coverage" from a claim into a build failure. Scoped
   to parameters and inputs — outputs are covered by the reference pages, per the
   note at the top of this page.
3. **Primaries resolve.** Every chapter key used in a **Primary** cell must
   resolve to a file listed in `docs/_toc.yml`, and each row must name exactly
   one.
4. **Reference pages agree with the matrix.** The chapter named in the
   *Demonstrated in* column of `R-GG` and `R-CL` must appear in this matrix as
   either the primary or an "Also in" chapter for that parameter's group.
5. **Deprecated spellings stay quarantined.** `--p-cv--nlam` must not appear in a
   runnable command outside `C-MSEL`, `R-CL` and `R-TS`, where it is being
   explained rather than recommended.
6. **Links resolve.** Every relative link in every chapter points at a file that
   exists — cheap, and the failure mode readers notice first.

```{important}
**None of this is wired up yet.** The `--help` capture step is not wired into the
build, and the parameter reference pages are maintained by hand in the meantime.
Until that changes, treat
`qiime <plugin> <action> --help` on your own installation as the final authority,
this matrix as the intent, and any disagreement between them as a bug worth
filing.
```

## See also

- [q2-gglasso Parameter Reference](02_gglasso_parameters.md) — every flag, type
  and default for the 6 gglasso actions.
- [q2-classo Parameter Reference](03_classo_parameters.md) — the same for the 8
  classo actions, organised by model-selection procedure.
- [Troubleshooting & Known Gotchas](04_troubleshooting.md) — the traps referenced
  throughout this page, with the symptom you will actually see.
