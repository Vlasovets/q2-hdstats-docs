# Overview

The Atacama chapters introduce no new QIIME 2 action. Every command on the
following pages has already been demonstrated in
[Tier 1](../02_lowdim_gglasso/00_index.md) on a 13-ASV table you could read off
the screen. What changes is the size of the problem — 300 ASVs across 54
samples of the Atacama soil microbiome — and with it the need to argue for the
parameter values. At $p = 13$ you can pick $\lambda_1$ by eye. At $p = 300$,
with $p \gg n$, the penalty *is* the model: it decides how many of the 44,850
possible edges you are willing to believe.

Three problems come with the larger table:

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

One parameterization covers the Atacama analysis. Every page uses it, and any
figure or artifact carrying other values comes from an older bundle.

| Setting | Value | Where it is argued |
|---|---|---|
| `--p-n-samples` | `54` | fixed by the table |
| `--p-lambda1-min` / `--p-lambda1-max` | `0.8` | [Selecting lambda](02_model_selection.md) |
| `--p-gamma` | `0.3` | [Selecting lambda](02_model_selection.md) |
| `--p-path-scale` | `linear` | [Selecting lambda](02_model_selection.md) |
| `--p-latent` | on, for the SLR models | [Choosing the Latent Rank](03_slr_ranks.md) |
| `--p-mu1-min` / `--p-mu1-max` | `15`, giving rank 2 | [Choosing the Latent Rank](03_slr_ranks.md) |

Change these three last — they are the ones people change first and
regret later.

**$\lambda_1 = 0.8$** minimises the extended BIC along a linear path. Move it
down and the graph densifies fast: the edge count is a steep function of
$\lambda_1$ in this regime, so a value that looks "only slightly less
conservative" can multiply the edge set several-fold.

**$\gamma = 0.3$** is the eBIC's extra edge penalty. $\gamma$ and $\lambda_1$
are selected jointly, in the sense that changing $\gamma$ changes which
$\lambda_1$ wins. Foygel and Drton's conventional $\gamma = 0.5$
{cite}`foygel2010extended` and the tier-1 default of `0.01` both select a
different network on this table. If you report a network, report the $\gamma$
that produced it.

**$\mu_1 = 15$** is the nuclear-norm penalty on the low-rank block, and the only
way to influence the rank — `--p-rank` raises `NotImplementedError` on every
released GGLasso. Larger $\mu_1$ gives a smaller rank. The value 15 was found by
scouting at $\lambda_1 = 0.8$, and does not transfer to another $\lambda_1$ or
another table.

```{important}
Read the `--p-mu1-*` value as a penalty, never as a rank. The rank is an output.
[Choosing the Latent Rank](03_slr_ranks.md) gives the procedure for turning a
target rank into a $\mu_1$ and for reading back the rank you got.
```

Two further SLR fits, at $\mu_1 = 10$ and $\mu_1 = 7.5$, appear on
[03](03_slr_ranks.md) and [04](04_latent_pca.md). They give you something to
argue the rank-2 choice against, and — because `qiime gglasso pca` caps
`--p-n-components` at the achieved rank — they supply a solution with room for
more than the single pair-plot panel a rank-2 fit allows. Neither is an
alternative canonical model.

## The superseded $\lambda = 0.95$ bundle

An earlier exploratory bundle of this analysis circulated with
$\lambda = 0.95$ and $\mu = 10.5$, and drafts, notebooks and slide figures
carrying those values still exist. Do not mix them with anything on these pages.

Whatever the merits of 0.95 in isolation, the two bundles answer different
questions:

- The 0.95 bundle predates the linear-path eBIC sweep and the choice of
  $\gamma = 0.3$. Its $\lambda$ was hand-picked at the conservative end of the
  path rather than selected by a criterion.
- Its $\mu = 10.5$ was scouted against $\lambda = 0.95$. The sparse and low-rank
  blocks compete for the same covariance structure, so the
  $\mu \rightarrow \text{rank}$ map is $\lambda$-specific: a rank obtained at
  $\lambda = 0.95$ is not the rank the same $\mu$ gives at
  $\lambda = 0.8$, and the two are not comparable.

An artifact, figure or table mentioning $\lambda = 0.95$ or $\mu = 10.5$ comes
from the exploratory bundle, and its edge counts, ranks and eBIC values cannot
be quoted alongside the ones here. The canonical bundle is
$\lambda_1 = 0.8$, $\gamma = 0.3$, $\mu_1 = 15$.

## Recompute status

```{note}
**Recomputed and confirmed.** The single-graphical-lasso $\lambda$ path on
[Selecting lambda](02_model_selection.md) has been re-run under QIIME 2 2026.7:
the selected $\lambda = 0.8$, its 216 edges and the $\gamma$ sensitivity all
hold, and the recompute script generates the eBIC table on that page into
`docs/_data/atacama-lambda-path.tsv` instead of anyone typing it in.
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
([Interpretation](06_interpretation.md)) are all pending verification against
QIIME 2 2026.7.

Treat the pending values as the expected result of the recompute rather than a
confirmed one.
```

Once the recompute has run, the numbers will not be typed into these pages by
hand. Each figure comes from a pinned artifact bundle plus a recompute script:
the bundle fixes the inputs (the count table, the taxonomy, the metadata, and
their checksums), the script fixes the commands and the parameter values in the
table above, and the page cites both. Only then can a reader who disagrees with
a number tell whether it came from a different input or a different command.

```{note}
The artifact bundle will be published as the Zenodo record
*q2-hdstats-tutorial-data* v1.0.0. The DOI has not been minted yet. Wherever a
download URL is required, these pages carry the literal placeholder
`ZENODO_DOI_PENDING`; `grep` for it to find every occurrence and replace them in
one pass when the record goes live.
```

## The pages

| Page | What it adds | Actions used |
|---|---|---|
| [The 300-ASV Dataset](01_data.md) | provenance of the table, taxonomy and transformed artifacts | `transform-features`, `calculate-covariance` |
| [Selecting lambda](02_model_selection.md) | eBIC on a linear $\lambda$ path; the $\gamma$ sensitivity | `solve-problem` |
| [Choosing the Latent Rank](03_slr_ranks.md) | the $\mu_1 \rightarrow$ rank scouting procedure | `solve-problem`, `summarize` |
| [Latent Components & Covariates](04_latent_pca.md) | projecting samples on the latent axes and matching them to measured variables | `pca` |
| [Log-Contrast Models at Scale](05_classo_cv.md) | cross-validated log-contrast regression for each covariate | `regress`, `add-covariates` |
| [Interpretation](06_interpretation.md) | reading the sparse and low-rank parts together, and against the regression coefficients | — |

## Conventions

Commands assume you are working in the directory where the Tier 2 bundle was
unpacked, so artifacts appear as bare filenames
(`atacama-top-300-correlation.qza`, not `data/atacama-top-300-correlation.qza`).
Adjust the paths if you keep them elsewhere;
[Download the Tutorial Data](../00_getting_started/03_download_data.md) gives
the layout the rest of the book assumes.

These pages assume both plugins are installed and registered — see
[Verifying Your Installation](../01_installation/04_verify.md) — and that you
have worked through at least
[Single Graphical Lasso](../02_lowdim_gglasso/02_sgl.md),
[Sparse + Low-Rank](../02_lowdim_gglasso/03_slr.md) and
[Latent-Component PCA](../02_lowdim_gglasso/07_pca.md), none of which is
explained again here.

[Troubleshooting](../90_reference/04_troubleshooting.md) collects the failures
you are most likely to hit at this scale — `--p-rank` raising, `pca` refusing a
non-latent solution, a model-selection run that silently collapsed to a single
fit. The full parameter lists are in the
[q2-gglasso](../90_reference/02_gglasso_parameters.md) and
[q2-classo](../90_reference/03_classo_parameters.md) references.
