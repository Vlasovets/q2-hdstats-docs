# Latent-Component PCA

`qiime gglasso pca` is not a PCA of the feature table. It projects your samples
onto the eigenvectors of the **low-rank component** `L` of a Sparse + Low-Rank
solution. The sparse component `Θ` answers "which taxa are directly linked?";
`L` absorbs the dense, global structure that no sparse graph can represent —
environmental gradients, sequencing depth, batch. This visualizer is how you
look at that structure and ask what it corresponds to.

That framing explains both of the traps below: without a low-rank block there is
nothing to project onto, and without metadata there is nothing to interpret the
projection against.

## Two prerequisites, before you run anything

```{important}
**1. The solution must come from `--p-latent True`.** `pca` reads
`solution/lowrank_` out of the solution artifact, and that group only exists for
a latent problem. An SGL solution from [Single Graphical Lasso](02_sgl.md) will
fail here — not with a helpful message about missing latent variables, but with
a lookup error on the solution's internals.

**2. `--m-sample-metadata-file` is effectively required.** It is registered as an
optional `Metadata` parameter, but the function body dereferences it
unconditionally — the first thing it does with the metadata is filter it to
numeric columns — so omitting it raises an `AttributeError` on `NoneType`. Treat
it as mandatory.
```

Both are listed in [Troubleshooting](../90_reference/04_troubleshooting.md); they
are repeated here because they are the two ways almost everyone's first `pca`
call fails.

## Running the visualizer

This continues from the latent solution fitted in
[Sparse + Low-Rank](03_slr.md):

```bash
qiime gglasso pca \
    --i-table data/atacama-table-mclr.qza \
    --i-solution data/atacama-solution-slr.qza \
    --m-sample-metadata-file data/selected-atacama-sample-metadata.tsv \
    --p-n-components 2 \
    --p-color-by ph \
    --o-visualization data/slr-pca.qzv
```

**Explanation:**

- `--i-table`: the **transformed** table — the same artifact you passed to
  `calculate-covariance` in [Data Preparation](01_data_preparation.md). The
  projection is a matrix product between the table and the eigenvectors of `L`,
  so the features must be the same set, in the same order, as the ones behind
  the covariance matrix. Passing the raw counts, or a table filtered
  differently, produces either a shape error or a silently meaningless plot.
- `--i-solution`: a latent solution, per the prerequisite above.
- `--m-sample-metadata-file`: sample metadata, used both to colour the points and
  to populate the axis-selection dropdowns.
- `--p-n-components`: how many principal components enter the pair-plot grid
  (default `3`). See the constraint below — this one is not free, and the
  command above will only run if the fitted low-rank block has rank 3 or more.
- `--p-color-by`: the metadata column used for the colour scale.
- `--o-visualization`: the `.qzv`, viewable with `qiime tools view` or at
  [QIIME 2 View](https://view.qiime2.org/).

```{note}
The rendered figure is not reproduced here, and neither is the rank this
particular fit achieves. No QIIME 2 2026.7 environment exists yet, so nothing has
been re-run; the achieved rank, the variance shares on the axis labels and the
appearance of the tabs are all pending verification against QIIME 2 2026.7. Read
the rank off your own solution before fixing `--p-n-components`.
```

## `--p-n-components` must not exceed the achieved rank

The constraint is `n_components <= rank(L)` — asking for exactly as many
components as the achieved rank is legitimate. `pca` validates it up front and
raises

```
n_components (N) exceeds the rank of the low-rank component (R).
Pass --p-n-components R or lower.
```

The pair-plot carries the same check as a defensive assertion for direct callers,
worded `n_components (N) is greater than the rank of the low-rank component (R)`.

With the default of `3` you therefore need a low-rank block of rank 3 or more.
That is a real constraint, not a formality: the rank is set indirectly by `μ₁`,
and a **larger `μ₁` gives a smaller rank**. If you tuned `μ₁` for a compact,
interpretable latent block you may well end up with rank 2, at which point the
largest admissible value is `--p-n-components 2` — enough for exactly one pair in
the grid. Rank 1 is unusable outright: the single-plot tab defaults to PC1
against PC2 and needs at least two components to exist.

You cannot set the rank directly. `--p-rank` raises `NotImplementedError` on
every released GGLasso; see
[Troubleshooting](../90_reference/04_troubleshooting.md) and, for the practical
`μ₁` scouting procedure, [Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md).

To read the rank you actually got, open the solution with
[`summarize`](08_summarize.md) and look at the `rank` column of the statistics
tab. That column comes from the model-selection path, so it is only present if
the solution was fitted over a grid rather than as a single fit — see
[Regularization Paths & Model Selection](05_lambda_paths.md).

```{important}
Do not size `--p-n-components` exactly to the number you read off that tab. Two
things stand between it and the components you actually get. First, the tab
prints **one rank per grid point**, so you have to match the row to the selected
`λ₁`/`μ₁` from the best-parameter table yourself; the best-parameter table
reports the chosen hyperparameters, not the rank. Second, the limit `pca`
enforces is `np.linalg.matrix_rank` recomputed on the *stored* `L`, whose default
tolerance is far finer than the hard `1e-9` eigenvalue cut the projection itself
applies — so `matrix_rank` can count eigenvalues that the projection then
discards, and the guard can pass while fewer than `n_components` components
exist. Leave yourself a margin below the reported rank rather than sizing to it
exactly.
```

## `--p-color-by` sees numeric columns only

Before anything is plotted, the metadata is filtered to numeric columns. A
categorical column — a yes/no vegetation flag, say, as carried by tier 2's
`sample-metadata.tsv` — is dropped, and naming it in `--p-color-by` fails with a
lookup error on a column that is no longer there. The same filter governs the
axis dropdowns, so categorical variables cannot be plotted against a component
either. If you need to see a grouping variable, encode it numerically in the
metadata file first. Note that the tier-1 file used above,
`selected-atacama-sample-metadata.tsv`, holds four numeric columns and no
categorical one, so on this data the filter is a no-op and nothing is dropped.

The filtered metadata is then re-indexed onto the sample IDs of the table. Any
sample ID present in the table but absent from the metadata file, or spelled
differently in the two, becomes a missing value rather than an error. Missing
values in the colouring column are load-bearing: the pair-plot drops those
samples before projecting, so it can show fewer points than the single plot, and
the single plot derives its colour-scale bounds from the raw minimum and maximum
of the column, which missing entries can poison. Check that your colouring
column is complete for exactly the samples in the table.

```{tip}
`pca` builds a complete pair-plot grid for **every** numeric metadata column and
then displays only the one named by `--p-color-by`. Wide metadata therefore costs
real time for output you never see. Trim the metadata file to the columns you
care about before running it on anything large.
```

### The `seq-depth` fallback

If you omit `--p-color-by`, the visualizer computes a column itself: per-sample
sequencing depth from the input table, rescaled to `[0, 1]` and added as
`seq-depth`. Conceptually this is the single most useful first plot you can
make here — if the leading component tracks `seq-depth`, the low-rank block is
modelling a technical gradient, and any biological reading of it is unsafe.

```{note}
The orientation is not guessed. `pca` first lines the table up against the
low-rank component — whichever axis has $p$ entries becomes the feature axis,
and the table is transposed if needed — and raises if neither axis matches. The
helper then sums over features unconditionally, so the derived column is always
indexed by sample and joins cleanly to the sample-indexed metadata at any $p/N$
ratio, including the 13 ASVs of this tier. `--p-color-by` remains the better
choice whenever you have a variable in mind; the fallback is the default because
it is the right *first* plot, not because it is a wide-table-only convenience.
```

## The two tabs

**Single plot.** One large scatter of the samples, PC1 against PC2 by default,
with a colour bar for the selected variable and hover tooltips carrying the
coordinates. The two `X-Axis` / `Y-Axis` dropdowns to the right of the plot
re-map the axes client-side to any other component *or any numeric metadata
column*. That second option is the point of the tab: plotting PC1 directly
against `ph` or `elevation` is the quickest way to see whether a latent axis is
a measured gradient in disguise.

Axis labels carry the share of variance each component explains, computed from
the eigenvalues of `L` — that is, the proportion of the *latent* structure, not
of the total variance in the table. Only eigenvalues above a numerical tolerance
are kept, so the number of components on offer is the numerical rank of `L`.
Note also that the projection is whitened by those eigenvalues, so scores on
different components are on comparable scales and a wide spread on a late
component does not mean it dominates.

**Pair-plot.** A triangular grid of every pair among the first
`--p-n-components` components, all coloured by `--p-color-by`, with the colour
bar in the corner cell; the diagonal and one triangle stay empty. Panels are
small and fixed-size, so the grid is for spotting which pair separates your
samples, after which you go back to the single plot to look at that pair
properly.

The two tabs use different colour palettes for the same variable — a blue ramp
in the single plot, viridis in the pair-plot. Compare positions across tabs, not
shades.

## Reading latent components as putative drivers

The useful question is not "what is PC1?" but "does PC1 line up with anything I
measured?".

If a component correlates visibly with a measured covariate, the low-rank block
is standing in for that covariate. The follow-up is to stop treating it as
latent: append the covariate to the table with
`transform-features --p-add-metadata True` and re-fit, as in
[Adaptive Graphical Lasso](04_adaptive_glasso.md), so the variable is modelled
explicitly and the remaining sparse edges are conditional on it. Comparing the
edge set before and after is what turns "there is confounding" into "these
specific edges were environment-mediated" — the comparison set out in
[Interpretation](09_interpretation.md).

If a component lines up with nothing you measured, you have a candidate
unmeasured driver: an unrecorded gradient, a collection batch, an extraction
run. That is a hypothesis to check against your study design, not a result.

```{important}
Do not over-read individual axes. The eigenvectors of `L` are determined only up
to sign, so the direction of an axis carries no meaning, and components with
close eigenvalues are only defined up to rotation within their shared subspace —
which of two near-tied components a sample loads on is not interpretable.
Interpret the subspace spanned by the leading components, and treat a
component-by-component story as unsupported unless the eigenvalues are clearly
separated.
```

A worked, higher-dimensional version of this analysis — several ranks, the same
metadata — is in [Latent PCA on the 300-ASV table](../04_highdim_atacama/04_latent_pca.md).
The full parameter list is in the
[q2-gglasso parameter reference](../90_reference/02_gglasso_parameters.md).
