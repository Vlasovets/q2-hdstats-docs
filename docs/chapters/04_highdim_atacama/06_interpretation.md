# Interpretation

By this point Tier 2 has produced four things: a sparse network at
$\lambda_1 = 0.8$ ([Selecting lambda](02_model_selection.md)), a sparse + low-rank
decomposition at $\mu_1 = 15$ ([Choosing the Latent Rank](03_slr_ranks.md)), a
reading of what the latent axes correspond to
([Latent Components & Covariates](04_latent_pca.md)), and a set of
cross-validated log-contrast models, one per environmental outcome
([Log-Contrast Models at Scale](05_classo_cv.md)).

This page is about what those four things say **together**, and — at least as
importantly — what they do not say. It fits no models and introduces no
parameters.

```{note}
**Which numbers on this page are verified.** The recompute has run under QIIME 2
2026.7, but it did not cover everything, so the claims below are not uniformly
solid:

| quantity | status |
|---|---|
| $\lambda = 0.8$, 216 edges, eBIC 16130.0988 | reproduced through the CLI |
| $\mu_1 \rightarrow$ rank map, and the edge/node counts at each rank | reproduced |
| taxonomy of individual nodes | resolves for all 300 features since the bundle was rebuilt with real feature IDs |
| correlations between components and covariates | **not** re-run |
| $R^2$ values | **not** re-run — the recompute reports cross-validated *error*; converting it needs a held-out `predict` pass |
| named taxa in the narrative below | carried over from the reference analysis |

Read the unverified rows as "how to read this number once you have it", not as a
report of it.
```

## The two blocks are one model

It is tempting to treat $\hat{\Theta}_S$ as "the network" and $\hat{L}$ as
"the noise we removed". That is the wrong split. The model is

$$
\hat{\Theta} = \hat{\Theta}_S - \hat{L},
$$

and the two terms are estimated jointly, competing for the same covariance
{cite}`chandrasekaran2010latent,kurtz2019disentangling`. Neither is meaningful
without the other being stated.

Concretely, this changes what an edge *is*. In the SGL network of
[Selecting lambda](02_model_selection.md), a nonzero $\hat{\Theta}_{ij}$ means
taxa $i$ and $j$ are conditionally dependent **given the other 298 taxa in the
table**. In the SLR network, it means they are conditionally dependent given the
other 298 taxa **and given the $r$ latent directions**. Those are different
claims about the world, and an edge that survives the second is a stronger
statement than an edge that only satisfies the first.

Symmetrically, the loadings of $\hat{L}$ are not a nuisance parameter. A taxon
with a large loading on a latent component is one whose abundance is largely
explained by a global gradient shared across the community. It may have very few
edges precisely *because* the latent block already accounts for most of its
covariance. Low degree in $\hat{\Theta}_S$ plus high loading in $\hat{L}$ is a
recognisable and interpretable profile: a generalist responding to the
environment rather than to its neighbours. Reading only the sparse block would
record such a taxon as uninteresting.

## Reading a taxon's neighbourhood

The unit of interpretation is not the whole 300-node graph — at this size the
full picture is a hairball, and the `summarize` heatmap is for spotting block
structure, not for reading individual relationships. The unit is one taxon and
its neighbours.

Export the solution, convert the precision matrix to partial correlations, and
join the taxonomy:

```bash
qiime tools export \
    --input-path atacama-top-300-slr-lambda0.8-rank2.qza \
    --output-path slr-rank2-export

qiime tools export \
    --input-path atacama-top-300-clr.qza \
    --output-path clr-export

qiime tools export \
    --input-path atacama-taxonomy-silva138.qza \
    --output-path taxonomy-export
```

```python
import biom
import numpy as np
import pandas as pd
import zarr

store = zarr.ZipStore("slr-rank2-export/problem.zip", mode="r")
root = zarr.open(store=store)
Theta = np.asarray(root["solution/precision_"])

# partial correlations: -Theta_ij / sqrt(Theta_ii * Theta_jj)
d = np.sqrt(np.diag(Theta))
P = -Theta / np.outer(d, d)
np.fill_diagonal(P, 1.0)

# feature IDs, in the order used to build the covariance matrix:
# the solution artifact stores them in an ordered `labels/` group
ids = [str(root[f"labels/{i}"][()]) for i in range(Theta.shape[0])]
P = pd.DataFrame(P, index=ids, columns=ids)

# cross-check only: the transformed table should agree, in the same order
X = biom.load_table("clr-export/feature-table.biom").to_dataframe(dense=True)
if X.shape[1] != Theta.shape[0]:
    X = X.T
assert list(X.columns) == ids, "table and solution disagree on feature order"

tax = pd.read_csv("taxonomy-export/taxonomy.tsv", sep="\t", index_col=0)

degree = (P.abs() > 0).sum(axis=1) - 1
hub = degree.idxmax()
neighbours = P.loc[hub][P.loc[hub] != 0].drop(hub).sort_values()
print(tax.loc[hub, "Taxon"])
print(neighbours.to_frame("partial_r").join(tax["Taxon"]))
```

```{note}
`solution/precision_` holds the **sparse** block $\hat{\Theta}_S$; the low-rank
part is stored separately as `solution/lowrank_`. Reading edges off
`precision_` therefore gives you the sparse network conditional on the latent
directions, which is what you want here — but it is not the full precision
matrix, and the two should never be conflated in a figure caption.

The solution artifact does not store only matrices. Alongside the blocks it
carries an ordered `labels/` group — one entry per feature, in the order the
covariance matrix was built in — taken from the column names of the input
correlation matrix. That group is the authoritative ID list, which is why the
snippet reads IDs from it and uses the exported table only as a cross-check. If
the assertion fires, do not silently reorder: a mismatched order relabels the
whole network.
```

```{note}
**The taxonomy join needs matching key spaces, and the current bundle has them.**
`atacama-top-300-clr.qza` is built with `--p-keep-original-id`, so its features —
and therefore the `labels/` group of every solution derived from it — are the
same 32-character hexadecimal IDs that key
`atacama-taxonomy-silva138.qza`. `tax.loc[hub, "Taxon"]` works directly, and all
300 features resolve.

If you are working from an older copy whose features are `ASV-1` … `ASV-300`, the
join will not raise — it returns all-`NaN`, which is easy to miss. Rebuild with
`--p-keep-original-id` as [The 300-ASV Dataset](01_data.md) shows rather than
attempting to map the names back; the next box explains why mapping cannot work.
```

```{important}
**Do not recover the mapping from `top-300-asvs.tsv` by abundance rank.** It is
tempting — the file has `feature-id` / `total-abundance` / `abundance-rank`, the
relabelling helper sorts **ascending** by total abundance, so `ASV-1` is the
*least* abundant and the mapping looks like
$\texttt{ASV-}n \longleftrightarrow \texttt{abundance-rank} = 301 - n$.

**It does not work, and it fails silently.** Total abundance is not unique:
**209 of these 300 features share a total-abundance value with another feature**
(61 tie groups, the largest holding 13). Within a tie group the rank order is
arbitrary, so `abundance-rank` and the plugin's internal ordering are free to
disagree — and they do. Permuting the correlation matrix by this mapping fails to
reproduce the shipped one, off by 1.137. Only the **91** features with a unique
total abundance are placed correctly; the rest get a neighbour's taxonomy, and
nothing raises.

This is not a hypothetical drift. The helper originally ordered features with
`df.sort_index()`, whose default quicksort is **not stable**; a later change to a
stable sort moved **158 of the 300** features to different `ASV-k` labels. Every
published number was unaffected — the graphical-lasso objective is invariant
under permutation, so the λ path, the eBIC at every grid point and the 216 edges
are bit-identical — but every *feature identity* shifted.

The lesson generalises: with `--p-no-keep-original-id`, `ASV-k` is a position, not
an identifier. It is only meaningful within the single artifact that defines it,
and it is never a key you can join on across artifacts.
```

```{note}
The taxonomy file's columns are `Feature ID`, `Taxon` and `Consensus` — if you
ever need the classifier's confidence, the column is `Consensus`, not
`Confidence`.
```

Four rules for reading the result.

**Scale.** Work in partial correlations, not raw precision entries. Precision
entries depend on the scaling of the covariance and are not comparable between
rows; partial correlations are bounded in $[-1, 1]$ and comparable across the
matrix.

**Sign.** A negative partial correlation is *not* competition and a positive one
is *not* cooperation. It is a conditional association among log-ratio-transformed
abundances, and the compositional geometry alone can produce either sign
{cite}`gloor2017microbiome,aitchison1982statistical`. The honest reading is "these
two vary together after conditioning on everything else in the table", full stop.

**Degree.** A high-degree node is a candidate hub, but degree is a function of
$\lambda_1$ before it is a function of biology: at a smaller penalty every degree
goes up. Comparing degrees *within* one fitted network is legitimate; quoting a
degree as a property of the organism is not.

**The sub-composition.** Every edge is conditional on the top-300 table. Removing
or adding features changes the CLR reference and therefore the whole geometry —
the warning in [The 300-ASV Dataset](01_data.md) applies to every statement on
this page. "Taxon A and taxon B are conditionally associated" always carries the
silent suffix "within this sub-composition, at this $\lambda_1$ and this
$\gamma$".

## What the SGL-to-SLR difference buys you

The most informative comparison in this tier costs nothing extra: the edges
present in the SGL network but absent in the SLR network at the same
$\lambda_1$. Those are the associations that a small number of global directions
can explain — candidate environment-mediated or batch-mediated edges. The edges
present in both are the ones that survive conditioning on the latent subspace.

On this dataset the difference is small and highly structured. Adding two latent
dimensions removes 14 of the 216 edges and adds none, and the removed edges are
not scattered: they fall into three components — one of 6 nodes carrying 9 of them,
one of 4 nodes carrying 4, and a single isolated pair. The remaining 202 edges are
common to both models.

Two things make that difference interpretable rather than incidental.

**The removed edges were the weak ones.** Their partial correlations have median
$|r| = 0.021$ against $0.116$ for the edges that survive — and the strongest edge
the latent block removed, $|r| = 0.062$, is weaker than the *median* surviving
edge. The latent component is not competing with the strong structure in the
network; it is absorbing a haze at the bottom of the edge-weight distribution.

**They are clustered, not spread.** Fourteen edges over 12 nodes, nine of them
inside a single 6-node group that was nearly complete before the latent block was
added. That is the shape you expect when one unmeasured driver acts on a handful
of taxa at once and induces weak mutual correlation among all of them. A driver
acting on the whole community would have thinned edges everywhere instead.

Neither observation proves the driver is environmental — see
[What you cannot conclude](#what-you-cannot-conclude) — but together they say the
rank-2 block is doing something specific and local rather than shaving the
network uniformly, which is what makes those 12 nodes worth looking up in the
taxonomy.

Export the sparse-only solution alongside the one you already unpacked:

```bash
qiime tools export \
    --input-path atacama-top-300-sgl-linear-path.qza \
    --output-path sgl-export
```

```python
import numpy as np
import zarr

def edges(path):
    root = zarr.open(store=zarr.ZipStore(path, mode="r"))
    T = np.asarray(root["solution/precision_"])
    iu = np.triu_indices_from(T, k=1)
    return set(zip(*(idx[T[iu] != 0] for idx in iu)))

sgl = edges("sgl-export/problem.zip")
slr = edges("slr-rank2-export/problem.zip")

print("SGL only :", len(sgl - slr))
print("both     :", len(sgl & slr))
print("SLR only :", len(slr - sgl))
```

The third number is the one people forget to look at. Edges can also *appear*
when the latent block is added, because removing a dense confounding direction
can unmask a direct association that was previously cancelled out. A model where
the SLR network is a strict subset of the SGL network is a possible outcome, not
a guaranteed one, and an "SLR only" set of any size is worth inspecting
individually.

```{note}
The three counts, and the taxonomy of the edges in each set, are **pending
recompute**. The comparison figures that circulated with the superseded
$\lambda = 0.95$ bundle described in the [overview](00_index.md) cannot be used
for this — they were computed at a different penalty.
```

## The network and the log-contrast coefficients

[Log-Contrast Models at Scale](05_classo_cv.md) fits, for each environmental
outcome $t$, a sparse coefficient vector $\hat{\beta}^{(t)}$ subject to
$\mathbf{1}^\top \beta = 0$. That is a completely different estimator on the same
matrix, and the interesting question is whether the two agree.

They are not measuring the same thing, and the difference matters before any
comparison is attempted:

- An **edge** is a conditional dependence between two taxa. It involves no
  outcome variable at all.
- A **log-contrast coefficient** is a weight in a zero-sum contrast predicting an
  outcome. A taxon selected by the regression is selected *relative to the other
  taxa in the contrast* {cite}`aitchison1984log,lin2014variable,shi2016regression`;
  the coefficient of a single feature has no meaning in isolation.

With that stated, there are three comparisons worth making.

**1. Are the selected taxa neighbours?** Take the features with nonzero
$\hat{\beta}^{(t)}$ for one outcome and look them up in the partial-correlation
matrix. If they form a connected subgraph, the regression is picking up a
coherent module and the two methods are describing one structure from two
directions. If they are scattered isolated nodes, the regression is exploiting
marginal signal that the network sees as conditionally independent — which is
possible and not necessarily wrong, but it means the two results are separate
findings rather than mutual corroboration.

**2. Does the latent subspace explain what is predictable?** This is the
comparison the [appendix](../99_appendix/01_math.md) formalises, and it is the
external check the rank-2 choice rests on. For each task $t$ it defines

- $m_t$, the strongest correlation between any robust principal component and
  the outcome — computed with the procedure in
  [Latent Components & Covariates](04_latent_pca.md); and
- $q_t = \lVert U^\top \hat{\beta}^{(t)} \rVert_2^2 / \lVert \hat{\beta}^{(t)} \rVert_2^2$,
  the fraction of the coefficient vector lying in the latent subspace $U$.

If $m_t$ and $q_t$ are rank-correlated across tasks, then the outcomes the
regression can predict are exactly the outcomes aligned with the latent subspace,
and two dimensions are carrying the predictable structure. That is a statement
about the *rank*, which is why it belongs in the rank argument in
[Choosing the Latent Rank](03_slr_ranks.md) as well as here.

```{important}
The Spearman correlation between $m_t$ and $q_t$ that appears in the appendix and
in the earlier drafts of this tier is **pending recompute** and must not be
quoted until it has been recomputed at $\lambda_1 = 0.8$, $\gamma = 0.3$,
$\mu_1 = 15$. It is also a correlation over a small number of tasks with strongly
inter-correlated outcomes: report it with the number of tasks and a permutation
$p$-value, and do not treat it as an independent confirmation of anything.
```

**3. Watch for the same signal counted twice.** If a latent component tracks
elevation and a log-contrast model predicts elevation well, you have observed one
gradient through two instruments, not two independent lines of evidence. The
filtered-covariate analysis on the [q2-classo page](05_classo_cv.md) exists for
the same reason: when a covariate is nearly a proxy for the outcome, adding it as
a predictor inflates $R^2$ without adding microbial information. Apply the same
scepticism to a network-to-regression agreement.

## What you cannot conclude

```{important}
**Nothing here is causal.** Conditional dependence is not interaction; a
coefficient is not an effect. Every result in this tier is observational, from a
single sampling campaign, with $n = 54$.

**Nothing here is about absolute abundance.** The model lives in log-ratio
coordinates. An edge or a coefficient describes relative structure within the
top-300 sub-composition and can move when the composition changes for reasons
that have nothing to do with the taxa involved.

**Nothing here is penalty-free.** The edge set is a function of $\lambda_1$ and
$\gamma$; the rank is a function of $\mu_1$; the selected features are a function
of the cross-validation fold assignment and the one-standard-error rule. Report
all of them, or the result is not reproducible even in principle.

**The sample size limits what is checkable.** With 54 samples, 300 features and
strongly inter-correlated covariates, individual edges are not stably estimated.
Read modules and gradients, not single links. If a specific edge carries the
weight of an argument, it needs a stability check — refitting on subsamples and
recording selection frequency {cite}`meinshausen2010stability` — not a
correlation coefficient.
```

```{note}
Earlier drafts of this tier named a *Pseudarthrobacter* ASV as the leading
selected feature across the log-contrast models, a genus well described in
hyperarid Atacama soils
{cite}`finger2018pseudarthrobacter,horstmann2025subsurface,neilson2017significant`.
That reading is plausible, but the identity of the leading feature is an output
of a model that has not been re-run, and it is **pending verification against
QIIME 2 2026.7** like everything else on this page.
```

## A reporting checklist

If you take one thing from Tier 2 to your own data, take this. A network result
is reportable when it comes with:

1. the table it was estimated on — feature count, sample count, and the filtering
   rule that produced them;
2. the transform (`clr` or `mclr`) and the pseudo-count, if any;
3. `--p-method` for the covariance and whether it was scaled;
4. $\lambda_1$, $\gamma$, and how $\lambda_1$ was selected;
5. for an SLR model, $\mu_1$ **and** the achieved rank, since the rank is an
   output;
6. what the latent components correlated with, including the sequencing-depth
   check;
7. the software versions of both the plugin and GGLasso itself.

Points 5 and 7 are the ones most often missing, and they are the two that make a
sparse + low-rank result impossible to reproduce when they are.

## Where to go next

The method-level comparison of SGL, SLR and the adaptive
model is in the Tier 1 chapter
[Network Interpretation and Analysis](../02_lowdim_gglasso/09_interpretation.md);
the definitions behind $m_t$, $q_t$ and the eBIC are in
[Appendix: Mathematical Background](../99_appendix/01_math.md); and all the works
cited here are collected on the [References](../../references.md) page.
