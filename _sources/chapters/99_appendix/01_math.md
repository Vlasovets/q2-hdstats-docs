# Appendix: Mathematical Background

This appendix collects the models and selection criteria used throughout the
tutorial. It draws on the general log-contrast and perspective-estimation
framework of {cite}`combettes2021regression,combettes2020perspective`, the
tree-aggregation model of {cite}`bien2021tree`, and the latent (sparse +
low-rank) graphical model of {cite}`chandrasekaran2010latent,kurtz2019disentangling`.

## Notation

Let $X \in \mathbb{R}^{n \times p}$ be the transformed abundance matrix for
$n$ samples and $p$ features (ASVs). Compositional counts are mapped to
Euclidean space with the centered log-ratio (clr) or modified clr (mclr)
transform {cite}`yoon2019microbial`. For the tutorial's high-dimensional example
$n = 54$ and $p = 300$. We write $S$ for the empirical covariance (or
correlation) matrix of $X$ and $\hat{\Theta}$ for an estimated precision
(inverse covariance) matrix.

## Log-contrast regression

For a continuous outcome $y \in \mathbb{R}^n$, the sparse **log-contrast model**
{cite}`aitchison1984log,lin2014variable,shi2016regression` is

$$
\min_{\beta \in \mathbb{R}^p,\ \sigma}\ f(y - X\beta,\ \sigma) + \lambda \lVert \beta \rVert_1
\quad \text{s.t.} \quad \mathbf{1}_p^\top \beta = 0,
$$

where the zero-sum constraint $\mathbf{1}_p^\top \beta = 0$ makes the fit
invariant to the compositional scale. The loss $f$ selects the c-lasso
formulation:

- **R1 — standard (least squares):** $f = \lVert y - X\beta \rVert_2^2$, with
  $\sigma$ fixed.
- **R2 — Huber:** $f = h_\rho(y - X\beta)$, robust to outliers.
- **R3 — concomitant (least squares):**
  $\ \lVert y - X\beta \rVert_2^2 / (2\sigma) + (n/2)\,\sigma$, jointly estimating
  $\sigma$.
- **R4 — concomitant Huber:** the robust analogue, with the linear term
  $n\,\sigma$.

The concomitant formulations R3/R4 estimate the noise level $\sigma$ jointly
with $\beta$ and have **no tuning parameter for $\sigma$**: the coefficients
$n/2$ (R3) and $n$ (R4) are fixed constants that fall out of the *perspective
function* of the respective loss {cite}`combettes2020perspective,combettes2018perspective`,
not hyperparameters.

### trac: tree-aggregation

The **trac** model {cite}`bien2021tree` reparameterizes the log-contrast problem
on the taxonomic tree. With $A$ the binary ancestor matrix mapping the $p$ leaves
to their $t$ tree nodes, one solves for node coefficients $\gamma$ with
$\beta = A\gamma$ under a weighted $\ell_1$ penalty $\lVert w \odot \gamma
\rVert_1$, so that selection acts on internal nodes (aggregated clades) rather
than individual ASVs — "tree-aggregation of compositional data".

## Graphical lasso

Microbial association networks are estimated with the **graphical lasso**
{cite}`friedman2008sparse,dempster1972covariance`,

$$
\min_{\Theta \succ 0}\ -\log\det \Theta + \langle S, \Theta \rangle
+ \lambda \lVert \Theta \rVert_{1,\text{od}},
$$

where $\lVert \cdot \rVert_{1,\text{od}}$ is the off-diagonal $\ell_1$ penalty
and edges correspond to nonzero off-diagonal entries of $\hat{\Theta}$.

### Latent variables: sparse + low-rank (SLR)

Unobserved environmental/technical factors induce dense confounding. The
latent-variable graphical model {cite}`chandrasekaran2010latent,kurtz2019disentangling`
decomposes the precision as a **sparse** part $\hat{\Theta}_S$ minus a
**low-rank** part $\hat{L}$,

$$
\hat{\Theta} = \hat{\Theta}_S - \hat{L}, \qquad
\min_{\Theta_S,\ L}\ -\log\det(\Theta_S - L) + \langle S, \Theta_S - L\rangle
+ \lambda \lVert \Theta_S \rVert_{1,\text{od}} + \mu \lVert L \rVert_\star,
$$

with $\lVert \cdot \rVert_\star$ the nuclear norm. The sparse part carries the
direct microbial interactions; the rank of $\hat{L}$ counts the latent factors.
A larger $\mu$ shrinks the rank. The low-rank subspace can be recovered by robust
PCA {cite}`candes2011robust`.

## Model selection

### eBIC for the single graphical lasso

Along a path of $\lambda$ values we score each fitted $\hat{\Theta}$ with the
**extended BIC** {cite}`foygel2010extended`:

$$
\mathrm{eBIC}_\gamma(\hat{\Theta}) =
n\big\{\operatorname{tr}(S\hat{\Theta}) - \log\det \hat{\Theta}\big\}
+ |E| \log n + 4\gamma\,|E| \log p,
$$

where $|E| = \tfrac{1}{2}\big(\#\{\hat{\Theta}_{ij} \neq 0\} - p\big)$ is the
number of off-diagonal edges and $\gamma \in [0,1]$ tunes the extra
edge penalty. $\gamma = 0.5$ is the conventional choice; the tutorial uses
$\gamma = 0.3$, which selects $\lambda = 0.8$ for the 300-ASV network.

### Cross-validation for log-contrast regression

Log-contrast models are selected by $k$-fold cross-validation with the
one-standard-error rule, or alternatively by stability selection
{cite}`meinshausen2010stability` or a theoretically-derived fixed penalty
{cite}`shi2016regression`. The reported out-of-sample $R^2$ is averaged over the
folds.

## Comparing the latent subspace and the regression coefficients

To relate the graphical model's latent structure to the regression tasks, the
tutorial defines (with $\hat{L} = U\Lambda U^\top$ truncated to the top two
eigenvectors, giving loadings $\ell_j = u_j / \sqrt{\lambda_j}$):

- **Robust-PC / covariate association** $m_t = \max_j \lvert \operatorname{cor}(z_j, y^{(t)}) \rvert$,
  the strongest correlation of any robust principal component $z_j$ with task
  $t$'s outcome.
- **Rank-2 subspace energy** $q_t = \lVert U^\top \hat{\beta}^{(t)} \rVert_2^2 / \lVert \hat{\beta}^{(t)} \rVert_2^2$,
  the fraction of the log-contrast coefficient vector lying in the latent
  subspace, with a permutation $p$-value.

Across tasks, $m_t$ and $q_t$ are strongly rank-correlated (Spearman
$\rho = 0.90$), which is the quantitative basis for choosing **rank 2**.

All references are collected on the [References](../../references.md) page.
