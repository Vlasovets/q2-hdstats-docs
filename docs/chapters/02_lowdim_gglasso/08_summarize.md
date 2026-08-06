# Summarizing a Solution

Chapters [02](02_sgl.md), [03](03_slr.md) and [05](05_lambda_paths.md) end with a
`qiime gglasso summarize` call that sets only `--p-label-size`; chapter
[04](04_adaptive_glasso.md) adds `--p-n-cov`. Chapters 01, 06 and 07 do not call
it at all. This chapter finishes the action: what the other three parameters do,
and — more usefully — what each of the four tabs in the `.qzv` actually contains,
since two of them change shape depending on how the solution was fitted.

`summarize` accepts any **single-instance** `GGLassoProblem` artifact, sparse or
latent. It is the default way to look at a solution before you export anything.

```{note}
Multi-group solutions — MGL and the non-conforming variant — are not handled. The
visualizer builds its tables by wrapping the stored covariance and precision
arrays in a `DataFrame` with no dimension check, which a $K \times p \times p$
stack cannot satisfy. This is academic in practice: as
[Multiple Graphical Lasso](06_multiple_graphical_lasso.md) sets out in its gap 3,
a stack of $K$ covariance matrices cannot be supplied through the QIIME 2
interface in the first place, so no such solution artifact can be produced to
pass here. **Pending verification against QIIME 2 2026.7.**
```

## Running the visualizer

```bash
qiime gglasso summarize \
    --i-solution data/atacama-solution-slr.qza \
    --p-width 900 \
    --p-height 900 \
    --p-label-size 12pt \
    --o-visualization data/slr-summary-annotated.qzv
```

**Explanation:**

- `--i-solution`: the solution artifact from `solve-problem`.
- `--p-width` / `--p-height`: heatmap size in pixels, `1500` each by default.
- `--p-label-size`: tick-label font, `"5pt"` by default.
- `--p-n-cov`: number of trailing covariate columns; unset by default. See below
  — this is the one that can quietly corrupt the figure.
- `--o-visualization`: open with `qiime tools view` or at
  [QIIME 2 View](https://view.qiime2.org/).

```{note}
No screenshots or statistics from this `.qzv` are reproduced in this chapter. No
QIIME 2 2026.7 environment exists yet, so nothing has been regenerated; the
rendered output is pending verification against QIIME 2 2026.7.
```

## Sizing: `--p-width`, `--p-height`, `--p-label-size`

The heatmap is drawn at a *fixed pixel size* regardless of how many features it
contains, so the three sizing parameters interact with `p` rather than being
cosmetic. At the default 1500 × 1500 the 13-ASV Atacama matrix of this tier
gives enormous cells and a lot of whitespace, which is why the earlier chapters
pair it with a large `--p-label-size 25pt`; the same 1500 px spread over the
300-ASV table of [Tier 2](../04_highdim_atacama/00_index.md) gives cells a few
pixels wide, where the default `"5pt"` labels are the only ones that fit.

The rule of thumb is to keep width and height equal — the matrix is square and a
non-square canvas distorts the diagonal — and to scale the label size so that
`p × label height` stays under the canvas height. If labels still collide, the
plot is zoomable and the hover tooltip reports `taxa_x`, `taxa_y` and the cell
value, so an unreadable axis is an inconvenience rather than a loss of
information.

`--p-label-size` is a string handed straight to bokeh and is not validated by the
plugin, so give it a proper CSS-style font size with a unit (`"5pt"`, `"12pt"`,
`"25pt"`). A bare number is rejected by bokeh, not by QIIME 2, so the error
arrives late and mentions a property name rather than your flag.

## `--p-n-cov` must match the appended covariates

This is the parameter to get right.

When you run `transform-features --p-add-metadata True`, the numeric metadata
columns are standardized and joined onto the transformed table, and they end up
as the **last** rows and columns of the covariance matrix that
`calculate-covariance` produces. The matrix is therefore blocked: an ASV × ASV
block, an ASV × covariate block, and a small covariate × covariate block.

`summarize` reorders rows and columns by hierarchical clustering (average
linkage, Euclidean distance) so that block structure is visible. Reordering is
not optional — there is no flag to switch it off. `--p-n-cov` tells it how many
trailing entries are covariates: with the parameter set, only the ASV block is
clustered and reordered, and the covariate rows and columns are pinned at the
edge of the matrix in their original order. This is what you want — covariates
shuffled into the middle of a taxon dendrogram are unreadable.

```{note}
If `--p-n-cov` does not equal the number of columns actually appended, the split
happens in the wrong place. Some ASVs get pinned as if they were covariates,
some covariates get clustered as if they were ASVs, and the axis labels stop
corresponding to the cells they sit next to. Nothing errors — you get a
plausible-looking heatmap that is wrong.
```

Count the covariates the way `transform-features` does: only **numeric**
metadata columns are appended, because the metadata is filtered by column type
first. A categorical column in your metadata file is not appended and must not be
counted. That is why [Adaptive Graphical Lasso](04_adaptive_glasso.md), whose
metadata file `selected-atacama-sample-metadata.tsv` carries exactly four numeric
covariates — `ph`, `average-soil-relative-humidity`, `elevation` and
`average-soil-temperature` — passes `--p-n-cov 4`. The file that does carry a
categorical `vegetation` column is `atacama-selected-covariates-veg.tsv`, but that
is the q2-classo covariate file, not the one chapter 04 hands to
`transform-features`.

If you did not use `--p-add-metadata`, **leave `--p-n-cov` unset** rather than
passing `0`. Unset means "one block, cluster everything", which is correct for a
pure ASV covariance. Zero is a different code path: it is not the same as unset,
and because `-0 == 0`, the ASV slice `iloc[:-0, :-0]` comes out empty while the
covariate slice `iloc[-0:, -0:]` swallows the whole matrix, so the clustering step
is handed an empty frame.

## The four tabs

### Sample covariance

The input matrix `S` — exactly what `calculate-covariance` produced — after
clustering. Use it as the "before" picture: whatever structure the graphical
lasso later attributes to direct edges should be visible here as marginal
correlation.

The clustering order is computed from this matrix and then reused verbatim for
the precision heatmap, so the two tabs are cell-for-cell comparable. That is the
whole reason to look at them side by side.

The colour scale is **fixed to `[-1, 1]`**, red for positive and blue for
negative, and is not rescaled to the data. With
`calculate-covariance --p-method scaled` (the default, i.e. a correlation
matrix) that is exactly right. With `--p-method unscaled` the entries are
unbounded, everything beyond ±1 saturates to solid red or blue, and the tab
becomes uninformative. If you need to view an unscaled covariance, expect that
limitation.

### Estimated inverse covariance

The fitted precision matrix `Θ`, which is the actual result: an off-diagonal
entry that is nonzero is an edge, and a zero is a conditional independence.

Note the sign convention. The tab is titled "Estimated inverse covariance", but
the matrix drawn is **negated** — the plot title inside the tab says
"Estimated (negative) inverse covariance". Precision entries carry the opposite
sign to the corresponding partial correlation, so negating restores the intuitive
reading: **red means a positive partial association**, consistent with the sample
covariance tab next to it. Comparing the two tabs directly is the point of the
convention; without the flip, every association would appear to change sign
between them.

### Low-rank

The low-rank component `L`. This tab exists **only for latent solutions** — those
fitted with `--p-latent True`, as in [Sparse + Low-Rank](03_slr.md). For a
sparse-only SGL solution there is no `lowrank_` group and the tab is simply
absent; the visualizer logs `NO low-rank solution has been found.`, which you
will see if you pass `--verbose`.

Read it as the counterpart to the sparse tab: dense, smooth, no isolated cells.
Bands of taxa moving together here are the global structure that the sparse
component was relieved of. To turn that into a per-sample statement, project it
with [`pca`](07_pca.md).

```{note}
The tab is built inside a `try` block whose `except` catches everything and emits
the same `NO low-rank solution has been found.` message. Any failure while
rendering the low-rank heatmap — not just a genuinely missing component — is
therefore reported as an absent low-rank solution, and the tab disappears from a
latent solution with no traceback. Reviewing the current sources, the clustered
low-rank branch calls the reordering helper with an argument name that helper
does not accept, and it derives a fresh clustering order from the *already
reordered* covariance matrix rather than reusing the one applied to the other
tabs. Both would be swallowed in exactly this way. Whether the tab renders, and
in what order, is **pending verification against QIIME 2 2026.7**. If a latent
solution shows only three tabs, this is the first thing to check.
```

### Statistics

Tables rather than heatmaps, and the tab whose contents depend most on how you
fitted the model.

If the solution came from a **model-selection run** — a grid with more than one
value — the tab shows one row per grid point with `sparsity`, `lambda`, `mu` and
`rank`, followed by a summary table giving the selected `best lambda`,
`best mu`, and the percentage of positive edges. That per-grid-point table is
where you read the achieved rank needed by
[`pca --p-n-components`](07_pca.md), and where you see how sparsity responds to
`λ₁` before committing to a value; the fuller treatment is in
[Regularization Paths & Model Selection](05_lambda_paths.md).

If the solution came from a **single fit**, none of that exists. The solver
writes no `modelselect_stats` group, so the grid table and the best-parameter
table are both omitted and the tab is reduced to the positive-edge percentage
alone.

```{important}
A single fit is harder to get than it looks. Model selection runs unless *every*
grid collapses to exactly one value: `lambda1` and `lambda2` for a sparse
problem, plus `mu1` for a **latent** one. Crucially, leaving
`--p-lambda2-min` / `--p-lambda2-max` unset does **not** count as one value — an
unset pair of bounds expands to the 5-point default `np.logspace(-1, -4, 5)` and
warns `Default values for lambda2 have been used.`, which is more than one value
and so forces model selection on. Since none of the commands in these chapters
pin `λ₂`, all of them run model selection and all of them get the full statistics
tab.

The realistic accident is therefore the opposite of what you might expect: you
ask for a single fit and get a path. To actually obtain one, pin `λ₂` explicitly
(and `μ₁` too, if latent) to a single value each. Note that `λ₂` is not used by
the single-instance solver at all — it only decides this branch. A reduced
statistics tab is the symptom of a genuine single fit. See
[Troubleshooting](../90_reference/04_troubleshooting.md).
```

The table to the right of the grid and best-parameter tables — not below them —
lists the precision-matrix entries pair by pair, with the diagonal removed and
the taxon labels substituted in. It is sortable, and it is the quickest way to
get a ranked edge list out of the visualization without exporting the artifact.

```{important}
Treat the taxon names in that pair table with suspicion until they are checked.
Reading the current sources, the table is built from the precision matrix as
stored — in the original feature order — while the labels substituted into it are
the ones renumbered by the clustering applied to the heatmaps. Wherever
clustering changes the order, which is almost always, the two would not line up,
and the values would carry the wrong pair of names. The magnitudes themselves are
unaffected. **Pending verification against QIIME 2 2026.7**; until then, confirm
any edge you intend to report against the heatmap hover tooltip, which is built
from the reordered matrix and its matching labels.
```

The complete parameter list for all six actions is in the
[q2-gglasso parameter reference](../90_reference/02_gglasso_parameters.md).
