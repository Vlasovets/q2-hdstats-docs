# Overview

Tier 2 introduces **no new action**. Every command on the following pages has
already been demonstrated in
[Tier 1](../02_lowdim_gglasso/00_index.md) on a 13-ASV table you could read off
the screen. What changes here is the size of the problem — **300 ASVs across 54
samples** of the Atacama soil microbiome — and, with it, the fact that the
parameter values stop being arbitrary. At $p = 13$ you can pick $\lambda_1$ by
eye. At $p = 300$, with $p \gg n$, the penalty *is* the model: it decides how
many of the 44,850 possible edges you are willing to believe.

So this tier is about three things the small example cannot teach you:

1. **Choosing $\lambda_1$ with a criterion** rather than a guess
   ([Selecting lambda](02_model_selection.md)).
2. **Sizing the latent block** when you cannot set its rank directly
   ([Choosing the Latent Rank](03_slr_ranks.md)), and reading the resulting
   components against measured variables
   ([Latent Components & Covariates](04_latent_pca.md)).
3. **Asking whether the two plugins agree** — whether the structure the network
   attributes to unobserved drivers is the same structure the log-contrast
   models are exploiting ([Interpretation](06_interpretation.md)).

## The canonical parameterization

There is exactly **one** parameterization for this tier. Every page uses it, and
any figure or artifact that does not carry these values belongs to an older
bundle (see below).

| Setting | Value | Where it is argued |
|---|---|---|
| `--p-n-samples` | `54` | fixed by the table |
| `--p-lambda1-min` / `--p-lambda1-max` | `0.8` | [Selecting lambda](02_model_selection.md) |
| `--p-gamma` | `0.3` | [Selecting lambda](02_model_selection.md) |
| `--p-path-scale` | `linear` | [Selecting lambda](02_model_selection.md) |
| `--p-latent` | on, for the SLR models | [Choosing the Latent Rank](03_slr_ranks.md) |
| `--p-mu1-min` / `--p-mu1-max` | `15`, giving **rank 2** | [Choosing the Latent Rank](03_slr_ranks.md) |

Three of these deserve a sentence of *why*, because they are the ones people
change first and regret later.

**$\lambda_1 = 0.8$** is not a round number someone liked; it is the minimiser of
the extended BIC along a linear path. Move it down and the graph densifies fast —
the edge count is a steep function of $\lambda_1$ in this regime, so a value that
looks "only slightly less conservative" can multiply the edge set several-fold.

**$\gamma = 0.3$** is the eBIC's extra edge penalty. It is *not* a nuisance
constant: $\gamma$ and $\lambda_1$ are selected jointly, in the sense that
changing $\gamma$ changes which $\lambda_1$ wins. Foygel and Drton's conventional
$\gamma = 0.5$ {cite}`foygel2010extended` and the tier-1 default of `0.01` both
select a different network on this table. If you report a network, report the
$\gamma$ that produced it.

**$\mu_1 = 15$** is the nuclear-norm penalty on the low-rank block, and it is the
*only* way to influence the rank — `--p-rank` raises `NotImplementedError` on
every released GGLasso. Larger $\mu_1$ means a smaller rank. The value 15 was
found by scouting, at $\lambda_1 = 0.8$; it does not transfer to another
$\lambda_1$ or another table.

```{important}
Read the `--p-mu1-*` value as a *penalty*, never as a rank. The rank is an
output. The procedure for turning a target rank into a $\mu_1$ — and for reading
the rank you actually got — is [Choosing the Latent Rank](03_slr_ranks.md).
```

Two further SLR fits, at $\mu_1 = 10$ and $\mu_1 = 7.5$, appear on pages
[03](03_slr_ranks.md) and [04](04_latent_pca.md). They are **not** alternative
canonical models: they exist so that the rank-2 choice can be argued against
something, and — because `qiime gglasso pca` caps `--p-n-components` at the
achieved rank — so that there is a solution with room for more than the single
pair-plot panel a rank-2 fit allows.

## One parameterization, not two: closing the $\lambda = 0.95$ question

An earlier **exploratory bundle** of this analysis circulated with
$\lambda = 0.95$ and $\mu = 10.5$. Drafts, notebooks and slide figures carrying
those values exist. **They are superseded.** Do not mix them with anything on
these pages.

The reason is not that 0.95 is wrong in isolation — it is that the two bundles
answer different questions:

- The 0.95 bundle predates the linear-path eBIC sweep and the choice of
  $\gamma = 0.3$. Its $\lambda$ was a hand-picked point on the conservative end
  of the path, not a selected one.
- Its $\mu = 10.5$ was scouted **against $\lambda = 0.95$**. Because the sparse
  and low-rank blocks compete for the same covariance structure, the
  $\mu \rightarrow \text{rank}$ map is $\lambda$-specific. A rank obtained at
  $\lambda = 0.95$ is not the rank you get from the same $\mu$ at
  $\lambda = 0.8$, and the two are not comparable.

Concretely: if you find an artifact, figure or table that mentions
$\lambda = 0.95$ or $\mu = 10.5$, it is from the exploratory bundle, and its edge
counts, ranks and eBIC values cannot be quoted alongside the ones here. The
canonical bundle is $\lambda_1 = 0.8$, $\gamma = 0.3$, $\mu_1 = 15$.

## Most numbers in this tier are still pending recompute

```{warning}
**Recomputed and confirmed.** The single-graphical-lasso $\lambda$ path on
[Selecting lambda](02_model_selection.md) has been re-run under QIIME 2 2026.7:
the selected $\lambda = 0.8$, its 216 edges and the $\gamma$ sensitivity are
confirmed, and the eBIC table on that page is **generated** into
`docs/_data/atacama-lambda-path.tsv` by the recompute script rather than typed in.
The achieved ranks and the $\mu_1 \rightarrow \text{rank}$ map
([Choosing the Latent Rank](03_slr_ranks.md)) have likewise been re-run at the
selected $\lambda = 0.8$: $\mu_1 = 15/10/7.5$ give rank $2/5/10$, and that table
is generated into `docs/_data/atacama-mu-rank-map.tsv`.

**Still pending.** Nothing else has been re-run. The SLR rank-comparison edge
counts and the $m_t$-versus-$q_t$ Spearman
correlation on [Selecting lambda](02_model_selection.md), the
component-to-covariate correlations
([Latent Components & Covariates](04_latent_pca.md)), the cross-validated $R^2$
values and selected features ([Log-Contrast Models at Scale](05_classo_cv.md)),
and the SGL-versus-SLR edge-set comparison and named taxa
([Interpretation](06_interpretation.md)) are all **pending verification against
QIIME 2 2026.7**.

Treat the pending values as the *expected* result of the recompute, not as a
confirmed one. Where such a page would naturally show output, it says so rather
than printing a number.
```

Once the recompute has run, the numbers will not be typed into these pages by
hand. Each figure comes from **a pinned artifact bundle plus a recompute
script**: the bundle fixes the inputs (the count table, the taxonomy, the
metadata, and their checksums), the script fixes the commands and the parameter
values in the table above, and the page cites both. That is the only way a
reader can tell whether a number they disagree with came from a different input
or a different command.

```{warning}
The artifact bundle will be published as the Zenodo record
*q2-hdstats-tutorial-data* v1.0.0. **The DOI has not been minted yet.** Wherever
a download URL is required, these pages use the literal placeholder
`ZENODO_DOI_PENDING`, so every occurrence can be found with `grep` and replaced
in one pass when the record goes live.
```

## What each page does

| Page | What it adds | Actions used |
|---|---|---|
| [The 300-ASV Dataset](01_data.md) | provenance of the table, taxonomy and transformed artifacts | `transform-features`, `calculate-covariance` |
| [Selecting lambda](02_model_selection.md) | eBIC on a linear $\lambda$ path; the $\gamma$ sensitivity | `solve-problem` |
| [Choosing the Latent Rank](03_slr_ranks.md) | the $\mu_1 \rightarrow$ rank scouting procedure | `solve-problem`, `summarize` |
| [Latent Components & Covariates](04_latent_pca.md) | projecting samples on the latent axes and matching them to measured variables | `pca` |
| [Log-Contrast Models at Scale](05_classo_cv.md) | cross-validated log-contrast regression for each covariate | `regress`, `add-covariates` |
| [Interpretation](06_interpretation.md) | reading the sparse and low-rank parts together, and against the regression coefficients | — |

## Conventions on these pages

Commands assume you are working in the directory where the Tier 2 bundle was
unpacked, so artifacts are referred to by bare filename
(`atacama-top-300-correlation.qza`, not `data/atacama-top-300-correlation.qza`).
Adjust the paths if you keep them elsewhere; see
[Download the Tutorial Data](../00_getting_started/03_download_data.md) for the
layout the rest of the book assumes.

Everything here assumes both plugins are installed and registered — see
[Verifying Your Installation](../01_installation/04_verify.md) — and that you
have already worked through at least
[Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md),
[Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) and
[Latent-Component PCA](../02_lowdim_gglasso/07_pca.md), since none of them is
explained again here.

The sharp edges you are most likely to hit at this scale — `--p-rank` raising,
`pca` refusing a non-latent solution, a model-selection run that silently
collapsed to a single fit — are collected in
[Troubleshooting](../90_reference/04_troubleshooting.md), and the full parameter
lists are in the
[q2-gglasso](../90_reference/02_gglasso_parameters.md) and
[q2-classo](../90_reference/03_classo_parameters.md) references.
