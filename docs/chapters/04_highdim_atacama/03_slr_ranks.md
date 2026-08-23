# Choosing the Latent Rank

[Selecting lambda](02_model_selection.md) fixed the sparsity penalty at
$\lambda_1 = 0.8$. The second penalty, $\mu_1$, controls the size of the latent
block in the sparse + low-rank decomposition

$$
\hat{\Theta} = \hat{\Theta}_S - \hat{L},
$$

where the rank of $\hat{L}$ counts the unobserved factors the model is allowed to
invoke {cite}`chandrasekaran2010latent,kurtz2019disentangling` (see the
[appendix](../99_appendix/01_math.md)).

The map from $\mu_1$ to achieved rank is specific to the dataset and to the
$\lambda$ it was fitted at, so the values quoted below hold for this table at
$\lambda_1 = 0.8$ and nowhere else. The method transfers to your own data; the
numbers do not transfer.

## You cannot set the rank

````{important}
`--p-rank` raises on every released GGLasso (up to and including
0.3.0):

```
NotImplementedError: Explicit rank selection requires a GGLasso version that
exposes an explicit low-rank rank parameter...
```

The parameter is registered but guarded, so it fails loudly rather than being
ignored. Released GGLasso exposes only the continuous $\mu_1$ penalty and
reports the achieved rank as an output. (Passing `--p-rank` together with
`--p-no-latent` raises a `ValueError` first, since a rank is meaningless without
a low-rank block.)

Explicit rank selection becomes available when GGLasso PR #50
(`fix_latent_rank`) is merged and released. Until then, use the scouting
procedure below.
[Troubleshooting](../90_reference/04_troubleshooting.md) records the same
behaviour.
````

The control you have instead is the nuclear-norm penalty $\mu_1$, together with
one reliable qualitative fact: a larger $\mu_1$ gives a smaller rank. The nuclear
norm is the $\ell_1$ norm of the eigenvalues, so penalizing it harder drives more
eigenvalues of $\hat{L}$ to exactly zero, and the rank counts the ones that
survive.

The relationship is monotone in intent but not calibrated: no formula takes a
target rank to a $\mu_1$. Fit, and read the rank back.

```{important}
The $\mu_1 \rightarrow \text{rank}$ map is specific to this covariance matrix at
this $\lambda_1$. The sparse and low-rank blocks compete to explain the same
structure, so changing $\lambda_1$ reallocates covariance between them and moves
the whole map. Scout $\mu_1$ at the $\lambda_1$ you intend to report, and rescout
if you change it. The superseded $\lambda = 0.95$ / $\mu = 10.5$ bundle described
in the [overview](00_index.md) therefore cannot be compared with the numbers
here.
```

## Procedure A: one model-selection pass over a $\mu_1$ grid

Build the map in a single run: pin $\lambda_1$ to one value and sweep $\mu_1$
over a grid. More than one grid then holds more than one value, so the solver
runs model selection and records per-grid-point statistics.

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix atacama-top-300-correlation.qza \
    --p-n-samples 54 \
    --p-latent \
    --p-lambda1-min 0.8 --p-lambda1-max 0.8 --p-n-lambda1 1 \
    --p-mu1-min 5 --p-mu1-max 20 --p-n-mu1 7 \
    --p-path-scale linear \
    --p-gamma 0.3 \
    --o-solution atacama-top-300-slr-lambda0.8-mu-scout.qza \
    --verbose
```

**Explanation:**

- `--p-lambda1-min 0.8 --p-lambda1-max 0.8 --p-n-lambda1 1` pins the sparsity
  penalty to the selected value. Pass both bounds — the trap below covers what
  the solver does when you leave them unset.
- `--p-mu1-min 5 --p-mu1-max 20 --p-n-mu1 7` with `--p-path-scale linear` gives
  the grid `5, 7.5, 10, 12.5, 15, 17.5, 20`. A linear grid is the right shape for
  a first scout because you do not yet know the order of magnitude at which the
  rank starts collapsing; switch to `--p-path-scale log` once you do.
- `--p-gamma 0.3` is the eBIC penalty. It matters here because model selection
  picks a *best* grid point by eBIC — see the caveat below.
- `--p-path-scale` also applies to the $\lambda_1$ grid, but that grid has one
  point, so the choice is inert for it.

Then read the map out of the visualization:

```bash
qiime gglasso summarize \
    --i-solution atacama-top-300-slr-lambda0.8-mu-scout.qza \
    --p-width 1500 --p-height 1500 \
    --p-label-size 4pt \
    --o-visualization atacama-top-300-slr-mu-scout-summary.qzv
```

The **Statistics** tab of the resulting `.qzv` carries one row per grid point with
`sparsity`, `lambda`, `mu` and `rank`. Those four columns are the
$\mu_1 \rightarrow \text{rank}$ map, and the same rows give the sparsity of the
corresponding sparse block.
[Summarizing a Solution](../02_lowdim_gglasso/08_summarize.md) describes the tab
in full.

```{note}
The achieved ranks the Statistics tab reports for the Atacama data are tabulated
under [Building the map](#building-the-map) below, re-run under QIIME 2 2026.7.
```

```{important}
A model-selection run returns one solution — the eBIC-best grid point — not one
per grid point. The statistics table is the map. The `solution/` group in the
artifact belongs to the selected $\mu_1$ only. When you need a fitted artifact
at each $\mu_1$ — to run `pca` on, or to count edges — use Procedure B.

The eBIC-selected $\mu_1$ is not automatically the $\mu_1$ to report. eBIC scores
fit against complexity and has no notion of whether a rank-8 latent block is
interpretable. Treat it as one input to the argument below rather than as the
decision.
```

## Procedure B: one single fit per $\mu_1$

When you need the individual solutions, fit them one at a time with every grid
collapsed to a single value:

```bash
for MU in 5 7.5 10 12.5 15 17.5 20; do
  qiime gglasso solve-problem \
      --i-covariance-matrix atacama-top-300-correlation.qza \
      --p-n-samples 54 \
      --p-latent \
      --p-lambda1-min 0.8 --p-lambda1-max 0.8 --p-n-lambda1 1 \
      --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
      --p-mu1-min ${MU} --p-mu1-max ${MU} --p-n-mu1 1 \
      --p-gamma 0.3 \
      --o-solution atacama-top-300-slr-lambda0.8-mu${MU}.qza
done
```

`lambda2` is inert for this single-graph problem, but leave the `--p-lambda2-*`
triple unset and the solver expands it to a five-point default path, which flips
the run into model selection — see the trap on singleton grids below. Any
positive value collapses it; `0.1` is used here because the default
`--p-path-scale log` takes a $\log_{10}$ of the bounds, so zero is not
available.

A single fit writes no `modelselect_stats` group, so `summarize` shows a reduced
Statistics tab with no `rank` column. Read the rank from the low-rank component
directly instead:

```bash
qiime tools export \
    --input-path atacama-top-300-slr-lambda0.8-mu15.qza \
    --output-path slr-mu15-export
```

The export contains a single file, `problem.zip`, which is a Zarr store:

```python
import numpy as np
import zarr

store = zarr.ZipStore("slr-mu15-export/problem.zip", mode="r")
root = zarr.open(store=store)

L = np.asarray(root["solution/lowrank_"])
Theta = np.asarray(root["solution/precision_"])

rank = np.linalg.matrix_rank(L)
n_edges = (np.count_nonzero(Theta) - Theta.shape[0]) // 2

print("achieved rank:", rank)
print("off-diagonal edges:", n_edges)
```

`np.linalg.matrix_rank` is what the plugin's own `pca` visualizer uses to decide
how many components exist, so this number governs `--p-n-components` in
[Latent Components & Covariates](04_latent_pca.md). The edge count follows the
same convention as the eBIC in the [appendix](../99_appendix/01_math.md): the
off-diagonal nonzeros of the precision matrix, halved for symmetry.

```{note}
`matrix_rank` thresholds singular values numerically, so the "achieved rank" is a
numerical statement about $\hat{L}$, not an exact algebraic one. Eigenvalues that
the solver drove very close to — but not exactly — zero can flip the reported rank
by one. If a rank sits on a boundary you care about, inspect the eigenvalue
spectrum of `L` with `np.linalg.eigvalsh` rather than trusting the integer.
```

## Two traps while scouting

```{important}
**Leaving a grid unset does not switch it off.** Omit `--p-mu1-min` and
`--p-mu1-max` on a latent problem and the solver substitutes
`np.logspace(2, -1, 10)` — a ten-point $\mu_1$ path spanning three orders of
magnitude — and emits a warning. The same applies to $\lambda_1$, which defaults
to `np.logspace(0, -4, 15)`. You get a long model-selection run you did not ask
for, at parameter values you did not choose. Pass the grids explicitly, and run
the scout with `--verbose` so you see the warning if you leave one unset.
```

```{note}
**A single fit needs every grid to be a singleton.** For a latent problem,
`lambda1`, `lambda2` and `mu1` must all collapse to one value before the solver
treats the run as a single fit. Pinning $\lambda_1$ and $\mu_1$ while leaving
`lambda2` to its default is enough to trigger a model-selection run.
Conversely, forgetting `--p-n-mu1` on what you intended to be a sweep gives you a
single fit at `mu1_min`. The symptom either way is the Statistics tab: a single
fit has no per-grid-point table.
```

## Building the map

Record the scout as a table. These are measured values, read out of the fitted
solutions at the Gate-C1-selected $\lambda = 0.8$:

```{csv-table} Achieved rank as a function of $\mu_1$ (at $\lambda = 0.8$)
:file: ../../_data/atacama-mu-rank-map.tsv
:delim: tab
:header-rows: 1
:widths: 20, 25, 25, 30
```

Read it the way the penalty works: a larger $\mu_1$ shrinks the low-rank block
harder, so it buys you a smaller rank, and a smaller rank leaves more sparse
edges behind.

The Atacama analysis uses three of these points and reports the fit at rank 2:

| Target rank | $\mu_1$ | Role |
|---|---|---|
| **2** | **15** | the canonical model, at rank 2 |
| 5 | 10 | comparison, and a usable input for `pca` |
| 10 | 7.5 | comparison |

```{note}
`analysis/slurm/03_mu_rank_map.sh`, in the repository accompanying this book,
generates this table: it fits each $\mu_1$ and writes
`analysis/results/tables/mu-rank-map.tsv` from the solution artifacts.
`docs/_data/atacama-mu-rank-map.tsv` is a copy of that file with the same four
columns, and the two are kept in step by hand. The script also asserts that each
fit really was a single fit — with `--p-lambda2-*` pinned, none of the three
carries a `modelselect_stats` group. $\mu_1$ *decreasing* from 15 to 7.5 makes
the rank *increase* from 2 to 10, the expected sign but a fast response over a
narrow range: verify rather than interpolate.
```

Every dimension you give the latent block removes covariance that the sparse
block would otherwise have explained with edges, so the edge count falls as the
rank rises. You are choosing where to draw the line between "an edge" and "a
shared latent driver", and you cannot see that trade-off from the rank alone.

## The argument for a low rank

Once you have the map, the question is which point on it to report.

**The sample budget.** A rank-$r$ latent block over $p = 300$ features introduces
on the order of $rp$ additional free parameters, against $n = 54$ samples. At
rank 2 that is a small extension of the model; at rank 10 the latent block alone
has more parameters than the data has observations several times over, and its
individual directions are not estimable in any meaningful sense — only the
subspace, and only loosely.

**Identifiability.** The sparse-plus-low-rank decomposition is only identifiable
when the sparse part is not itself low-rank-like and the low-rank part is not
sparse-like {cite}`chandrasekaran2010latent,candes2011robust`. That condition
degrades as the rank grows: a high-rank $\hat{L}$ can start absorbing structure
that is genuinely a set of direct edges, and you have no way to tell from the fit
which happened. A low rank keeps the two blocks doing recognisably different
jobs.

**Interpretability.** A latent dimension is only useful if you can eventually say
what it might be. In this study the plausible global drivers are few — the
aridity/elevation transect, soil moisture and temperature, sampling depth — and
they are strongly inter-correlated, so they do not span many independent
directions. A rank commensurate with the number of plausible drivers is a model
you can argue about; a rank of 10 is a model in which no individual axis can be
defended.

**Rotation ambiguity.** Components with close eigenvalues are defined only up to
rotation within their shared subspace, so a "PC7 corresponds to pH" story is
unsupportable unless the eigenvalues are clearly separated. Higher ranks produce
more near-ties, and therefore more axes about which nothing can be said.

**The external check.** Validate the rank against something outside the
graphical model: whether the latent subspace is the same structure the
log-contrast regressions are using. The [appendix](../99_appendix/01_math.md)
defines the two quantities for this — the strongest correlation $m_t$ between any
robust principal component and outcome $t$, and the fraction $q_t$ of the
log-contrast coefficient vector lying in the latent subspace. If a rank-2
subspace already accounts for the tasks the regressions can predict, adding
dimensions adds detail nothing downstream uses. The rank-2 choice for this
dataset rests on that comparison. The correlation itself is **pending recompute**,
and the reading of it belongs to [Interpretation](06_interpretation.md).

**The counter-argument.** A rank that is too *small* leaves confounding in the
sparse block, where it shows up as implausibly dense neighbourhoods and blocks of
taxa all connected to each other. Compare the sparse edge sets across your
scouted ranks. The rank to report is the smallest one at which two things have
stopped happening: the edge set stops changing materially as you add a dimension,
and the additional components stop tracking any measured covariate. If the edge
set is still moving at rank 2, the low-rank block is too small, whatever the
parsimony argument says.

```{note}
A low rank costs you something in the *visualizer*, not in the model.
`qiime gglasso pca` requires `--p-n-components` to be no greater than the
achieved rank, so a rank-2 solution admits at most `--p-n-components 2` — one
pair-plot panel, showing the same PC1-against-PC2 scatter the single-plot tab
already draws. That is a reason to keep the $\mu_1 = 10$ and $\mu_1 = 7.5$ fits
around for inspection. Report only a rank you can defend.
See [Latent Components & Covariates](04_latent_pca.md).
```

## Fitting the selected model

Fit the canonical model at the selected $\mu_1$:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix atacama-top-300-correlation.qza \
    --p-n-samples 54 \
    --p-latent \
    --p-lambda1-min 0.8 --p-lambda1-max 0.8 --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min 15 --p-mu1-max 15 --p-n-mu1 1 \
    --p-gamma 0.3 \
    --o-solution atacama-top-300-slr-lambda0.8-rank2.qza
```

Then project its samples onto the latent axes and ask what they correspond to:
[Latent Components & Covariates](04_latent_pca.md).
