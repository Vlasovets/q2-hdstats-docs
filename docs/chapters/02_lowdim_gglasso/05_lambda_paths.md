# Regularization Paths & Model Selection

Every chapter so far has handed `solve-problem` a range of $\lambda_1$ values and
accepted whatever came back. This chapter opens that box. Two separate decisions
are being made on your behalf:

1. **Which values are tried** — the *regularization path*, built from
   `lambda1_min` / `lambda1_max` / `n_lambda1` and the spacing rule
   `path_scale`, or supplied verbatim through `lambda1_path`.
2. **Which value is kept** — the extended BIC (eBIC), whose behaviour is
   controlled by a single parameter, `gamma`.

Getting either one wrong is quiet rather than loud: you still get a solution
artifact, it is just not the one you thought you asked for. This chapter is the
prerequisite for [Selecting lambda](../04_highdim_atacama/02_model_selection.md),
where the same machinery runs on 300 features and the choice actually changes
the biology.

## How the grid is built

For each penalty, `solve-problem` builds a grid from three numbers and a spacing
rule. With `--p-path-scale log` (the default) the grid is log-spaced between the
bounds; with `--p-path-scale linear` it is evenly spaced. The same
`path_scale` applies to $\lambda_1$, $\lambda_2$ and $\mu_1$ — there is one
spacing rule per invocation, not one per penalty.

Three behaviours are worth committing to memory, because none of them is what a
casual reading of the flags suggests:

**Omitting one bound does not omit the grid.** If you give `lambda1_min` but not
`lambda1_max`, the missing upper bound silently becomes `1`; a missing lower
bound becomes `1e-3`. You get a full grid built against a default endpoint you
never chose.

**Omitting both bounds gives you a default path, not a single fit.** When both
ends of a penalty range are unset, the plugin substitutes a built-in path and
emits a warning:

| Penalty | Default path when both bounds are unset | Values |
|---|---|---|
| $\lambda_1$ | `np.logspace(0, -4, 15)` | 15, from `1` down to `1e-4` |
| $\lambda_2$ | `np.logspace(-1, -4, 5)` | 5, from `0.1` down to `1e-4` |
| $\mu_1$ (only when `--p-latent`) | `np.logspace(2, -1, 10)` | 10, from `100` down to `0.1` |

Each substitution prints a warning of the form *"Default values for lambda1 have
been used."* — run with `--verbose` or you will not see it.

**An explicit path wins over everything.** `--p-lambda1-path` and
`--p-mu1-path` take a list of values that is used exactly as given. They
override `*_min`, `*_max`, `n_*` **and** `path_scale`. There is deliberately no
`lambda2_path`; if you need a hand-built $\lambda_2$ grid you have to express it
through `lambda2_min` / `lambda2_max` / `n_lambda2`.

```{note}
Neither plugin declares `Choices()` on its string parameters, so
`--p-path-scale logg` is accepted by the command line and fails inside the
function with `ValueError: Unknown scale 'logg', use 'log' or 'linear'.`
See [Troubleshooting](../90_reference/04_troubleshooting.md).
```

## Single fit versus model selection — the exact rule

`solve-problem` runs a model-selection search **whenever any relevant grid holds
more than one value**. The rule is evaluated on the *final* grids, after the
defaults above have been substituted, and the set of grids consulted depends on
`latent`:

- **Without `--p-latent`**: a single fit requires that the $\lambda_1$ grid *and*
  the $\lambda_2$ grid each hold exactly one value. The $\mu_1$ grid is not
  consulted.
- **With `--p-latent`**: a single fit requires that $\lambda_1$, $\lambda_2$
  **and** $\mu_1$ each hold exactly one value. All three.

Combine this with the default-path rule above and one consequence follows that
catches almost everyone: **leaving $\lambda_2$ unset does not make it a
singleton.** An unset $\lambda_2$ range is replaced by the five-value default
path, the grid therefore holds more than one value, and model selection runs no
matter what you did to $\lambda_1$. On a single-instance problem $\lambda_2$ has
no effect on the answer, but it still decides which branch you are in, so you
have to pin it explicitly to get one fit.

The latent case is where people get caught. Pinning $\lambda_1$ to one value and
leaving $\mu_1$ unset does *not* give you one model at that $\lambda_1$ — the
$\mu_1$ default path kicks in with ten values, model selection runs, and eBIC is
free to move away from the $\lambda_1$ you thought you had fixed. To fit exactly
one latent model you must pin all three:

```bash
# exactly one model: one lambda1, one lambda2, one mu1
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 \
    --p-latent True \
    --p-lambda1-min 0.1 --p-lambda1-max 0.1 --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min 0.5 --p-mu1-max 0.5 --p-n-mu1 1 \
    --o-solution data/atacama-solution-single-fit.qza \
    --verbose
```

```{note}
`--p-n-samples 50` matches the value used by the other Tier 1 chapters. It must
equal the number of samples behind the covariance matrix, and for the 13-ASV
subset that is settled: the table is **13 features × 50 samples**, read directly
off the artifact (see
[Downloading the Data](../00_getting_started/03_download_data.md)). On your own
table, read the count off your `qiime feature-table summarize` output rather than
copying a published number.
```

To tell after the fact which branch ran, look for a `modelselect_stats` group in
the solution: a single fit has none, and
[`summarize`](08_summarize.md) omits the statistics accordingly.

## The same problem, three ways

The three runs below differ only in how the $\lambda_1$ path is constructed.
Everything else — the covariance matrix, the sample count, `gamma`, the absence
of latent variables — is held fixed, so any difference in the selected model is
attributable to the path alone.

**A. Log-spaced (the default).** Fifteen values between `0.001` and `1`, each
about 1.64× the previous one:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 \
    --p-latent False \
    --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 15 \
    --p-path-scale log \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-path-log.qza \
    --verbose
```

**B. Linearly spaced.** Same endpoints, same number of points, uniform steps of
about `0.071`:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 \
    --p-latent False \
    --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 15 \
    --p-path-scale linear \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-path-linear.qza \
    --verbose
```

**C. An explicit path.** Ten values chosen by hand; `path_scale` is ignored:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 \
    --p-latent False \
    --p-lambda1-path 1.0 0.5 0.25 0.1 0.05 0.025 0.01 0.005 0.0025 0.001 \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-path-explicit.qza \
    --verbose
```

These are genuinely different questions, not cosmetic variants. Nine of the
fifteen log-spaced values in run A lie below the *second* linear value in run B
(`≈0.072`), so run A spends most of its budget resolving weak penalties and dense
networks, while run B spends most of its budget on the sparse end. If the eBIC
optimum sits near $\lambda_1 = 0.01$, run B cannot find it — its nearest grid
point is an order of magnitude away. Conversely, if the optimum sits near
$\lambda_1 = 0.6$, run A has only two candidates up there.

Run C exists for the case where you already know roughly where the optimum is
and want a reproducible, human-readable grid in the provenance record rather
than a triple of bounds that a reader has to re-derive.

```{note}
These three runs have not yet been recomputed against QIIME 2 2026.7, so the
selected $\lambda_1$, the eBIC value and the resulting edge counts are not
quoted here.
```

## Comparing the three selections

Visualize each solution and compare them on three axes:

```bash
qiime gglasso summarize \
    --i-solution data/atacama-solution-path-log.qza \
    --p-label-size 25pt \
    --o-visualization data/path-log-summary.qzv

qiime gglasso summarize \
    --i-solution data/atacama-solution-path-linear.qza \
    --p-label-size 25pt \
    --o-visualization data/path-linear-summary.qzv

qiime gglasso summarize \
    --i-solution data/atacama-solution-path-explicit.qza \
    --p-label-size 25pt \
    --o-visualization data/path-explicit-summary.qzv
```

- **The selected $\lambda_1$.** If the three runs land on values that are close
  in absolute terms, the optimum is well identified and the path only affects
  precision. If they land in different regimes, the eBIC surface is flat and you
  are choosing a model by grid design, not by evidence.
- **The sparsity level.** Two paths can select different $\lambda_1$ values and
  still produce nearly the same network. Edge count is the quantity you actually
  care about; $\lambda_1$ is only the dial that produces it.
- **Where the optimum sits on the path.** This is the one people skip, and it is
  the one that invalidates results.

```{important}
If the selected $\lambda_1$ is the smallest or largest value on the path,
`solve-problem` emits two warnings — a directional one naming the fix, followed
by a general one:

> `lambda is on the edge of the interval, try SMALLER lambda1`
>
> `lambda is on the edge of the interval, the solution might have not reached global minimum!`

(`try BIGGER lambda1` when the selection sits at the top of the range.) The
optimum is outside your search range. Widen the bounds in the direction named
and re-run. A selection at the boundary is not a selection — it is the grid
running out.

Note what is *not* checked: the boundary test covers $\lambda_1$ and, for
multiple-instance problems, $\lambda_2$. There is no equivalent check on
$\mu_1$, so a latent solution can silently select the smallest or largest
$\mu_1$ on your path with no warning at all. Read the selected $\mu_1$ out of
[`summarize`](08_summarize.md) yourself and confirm it is interior.
```

## `gamma`: what eBIC is actually optimizing

Model selection scores each fitted model with the extended BIC of Foygel and
Drton {cite}`foygel2010extended`. Ordinary BIC penalizes the number of edges; the *extended* version adds a
second penalty term scaled by $\gamma \in [0, 1]$ that accounts for the size of
the model space. Larger $\gamma$ means a heavier penalty on extra edges, hence a
sparser selection.

Three values appear in this tutorial, and the spread between them is not a
rounding difference:

| $\gamma$ | Where it is used | Character |
|---|---|---|
| `0.01` | the plugin default, and every Tier 1 chapter | almost plain BIC; the most permissive of the three |
| `0.3` | the 300-ASV analysis in [Tier 2](../04_highdim_atacama/02_model_selection.md) | a deliberate middle ground for a high-dimensional problem |
| `0.5` | the conventional choice in the eBIC literature | conservative; the usual default outside this tutorial |

Re-running the log path under all three makes the sensitivity visible:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 --p-latent False \
    --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 15 \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-gamma-001.qza

qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 --p-latent False \
    --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 15 \
    --p-gamma 0.3 \
    --o-solution data/atacama-solution-gamma-03.qza

qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 --p-latent False \
    --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 15 \
    --p-gamma 0.5 \
    --o-solution data/atacama-solution-gamma-05.qza
```

```{figure} ../../images/png/generated/toy-lambda-path-gamma.png
:name: fig-toy-gamma-path
:width: 100%

eBIC against $\lambda_1$ on the 13-ASV table, scored at five values of $\gamma$.
Open circles mark each criterion's choice. The curves separate vertically — a
larger $\gamma$ really does charge more for the same model — but they are
minimised at almost the same place, so the five criteria between them pick only
**two** distinct $\lambda_1$ values.

Generated by `analysis/slurm/30_tier1_figures.sh` from a single 30-point path,
not from the three commands above: GGLasso scores every candidate model under a
fixed internal set of $\gamma$, so one `solve-problem` run already carries all
five curves in `modelselect_stats/BIC`. The three-command form is the explicit
way to get one *selected solution* per $\gamma$; run at 15 grid points it lands
on neighbouring $\lambda_1$ values rather than these exact ones.
```

**On this dataset $\gamma$ barely matters, and that is the lesson.** The five
criteria collapse onto two adjacent grid points, and both are networks with
essentially nothing in them:

| $\gamma$ | selected $\lambda_1$ | edges (of 78 possible) |
|---|---|---|
| 0.01, 0.1 | 0.4894 | 1 |
| 0.3, 0.5, 0.7 | 0.6210 | 0 |

With 13 features and 50 samples the problem is not high-dimensional, eBIC's
model-space penalty has little to bite on, and the criterion is decisive well
before $\gamma$ gets a say. Turning $\gamma$ up from 0.01 to 0.7 costs you the
single surviving edge. This is why the Tier 1 chapters mostly show fits at a
*fixed* $\lambda_1$: on a table this small, letting eBIC choose gives you an
empty graph and nothing to talk about.

Contrast [Tier 2](../04_highdim_atacama/02_model_selection.md), where $p = 300$
and $n = 54$. There the same parameter flips the answer between 1405 edges, 216
edges and the empty graph across the window $\gamma \in [0.29, 0.32]$ — three
qualitatively different scientific conclusions, a hundredth apart.

Two practical consequences. First, **`gamma` does nothing at all in a single
fit** — it is only consumed by the model-selection routine, so pinning every
grid to one value makes the flag inert. Second, **`gamma` is a modelling choice
that must be reported**, in the same way a significance threshold is. A network
reported without its $\gamma$ is not reproducible, and as the two cases above
show, you cannot predict from the value alone whether it mattered — you have to
look.

```{tip}
Do not tune `gamma` until you are satisfied with the path. A `gamma` sweep on a
grid whose optimum sits at the boundary tells you about the grid, not about the
data.
```

## Latent paths: the $\mu_1$ dimension

With `--p-latent` there is a second path to design. `--p-mu1-min`,
`--p-mu1-max` and `--p-n-mu1` build it from bounds under the shared
`path_scale`; `--p-mu1-path` supplies it verbatim and overrides all three. The
search is over the *product* of the $\lambda_1$ and $\mu_1$ grids, so a 15 × 10
specification is 150 fits — cheap on 13 features, much less so on 300.

```bash
# scout a coarse mu1 path at a fixed lambda1
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-table-corr.qza \
    --p-n-samples 50 \
    --p-latent True \
    --p-lambda1-min 0.1 --p-lambda1-max 0.1 --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-path 10.0 5.0 2.0 1.0 0.5 0.2 \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-mu1-scout.qza \
    --verbose
```

$\mu_1$ controls the rank of the low-rank component: **a larger $\mu_1$ yields a
smaller rank**. This is the only handle available — `--p-rank` is registered but
never works: without `--p-latent` it raises `ValueError` ("the `rank` parameter
is only meaningful for the sparse + low-rank model"), and with `--p-latent` it
raises `NotImplementedError` on every released GGLasso up to and including
0.3.0, because no release can fix the rank directly. Scouting a $\mu_1$ path and
reading the achieved rank out of each solution is the supported workflow; see
[Sparse + Low-Rank](03_slr.md) and, at scale,
[Choosing the Latent Rank](../04_highdim_atacama/03_slr_ranks.md).

## A practical recipe

1. Start wide and log-spaced: `--p-lambda1-min 0.001 --p-lambda1-max 1
   --p-n-lambda1 15`. Log spacing is the right default because the interesting
   behaviour of an L1 penalty is multiplicative.
2. Run with `--verbose` and read the warnings. A "default values have been used"
   warning means you did not specify the grid you thought you did.
3. Check that the selection is interior, not at a boundary. Widen and re-run if
   it is not.
4. Refine with a narrower path — or an explicit `--p-lambda1-path` — around the
   optimum.
5. Only then vary `gamma`, and report the value you settled on.
6. For latent models, pin $\lambda_1$ and $\lambda_2$ while scouting $\mu_1$, and
   remember that leaving any of the three unpinned re-enables the full search.

The full parameter list, with types and defaults, is in the
[q2-gglasso Parameter Reference](../90_reference/02_gglasso_parameters.md).
