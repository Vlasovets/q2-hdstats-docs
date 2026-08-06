# References

This page is a **guided bibliography** rather than a flat list. It supports the
three worked examples the book is built around: the low-dimensional toy example
that introduces the concepts, the 300-ASV Atacama example that puts them under
high-dimensional pressure, and the metagenomics tier that moves from networks to
log-contrast regression.

References are grouped by the **role they play in the tutorial** — the theory the
models rest on, the estimators we actually run, the data and pipelines the
examples consume, and the software that implements it all. The grouping is there
so you can decide what to read: most readers need one or two entries from a
section, not the whole section.

```{figure} _static/figs/pipeline_overview_placeholder.png
:name: fig-pipeline-overview
:width: 100%

**Placeholder.** This figure will trace the pipeline end to end: raw reads to
ASVs, ASVs to a log-ratio transformed table, and from there to both a sparse
network and a log-contrast regression on metagenomic features.
```

## How these references support the tutorial

- **Toy example (13 ASVs).** Start with {cite}`aitchison1982statistical` and
  {cite}`gloor2017microbiome` for why counts must be treated as compositions,
  then {cite}`friedman2008sparse` for the graphical lasso itself. Everything the
  toy example does is those two ideas composed; the dataset is small enough that
  you can check the estimator's output by eye.
- **300-ASV example (Atacama).** The methods section carries the weight here:
  {cite}`friedman2008sparse` for the estimator, {cite}`foygel2010extended` for
  the eBIC model selection that picks λ, {cite}`chandrasekaran2010latent` and
  {cite}`kurtz2019disentangling` for the latent/low-rank decomposition, and
  {cite}`meinshausen2010stability` for selection stability. For the data itself
  see {cite}`neilson2017significant`, and {cite}`callahan2016dada2` with
  {cite}`quast2012silva` for how the ASVs and their taxonomy were produced.
- **Metagenomics tier.** This is where compositional *regression* replaces
  network estimation: {cite}`aitchison1984log` for the log-contrast formulation,
  {cite}`lin2014variable` and {cite}`shi2016regression` for variable selection
  under the zero-sum constraint, and {cite}`combettes2021regression` with
  {cite}`mishra2022robust` for the proximal and robust formulations that
  `q2-classo` implements.

## Core compositional theory

Microbiome counts carry no absolute scale: sequencing depth is arbitrary, so only
*ratios* between features are meaningful. These references establish that, and
give the log-ratio and log-contrast machinery that lets ordinary multivariate
methods be applied to data living on the simplex. They underpin every
transformation in the book — `clr`, `mclr`, and the zero-sum constraint in the
regression tier.

- {cite}`aitchison1982statistical` — the foundational treatment of compositional
  data; the reason the toy example transforms counts before estimating anything.
- {cite}`aitchison1984log` — introduces log-contrast models, the direct ancestor
  of the constrained regression used in the metagenomics tier.
- {cite}`gloor2017microbiome` — the accessible argument for why compositionality
  is not optional in microbiome analysis; read this first if the others feel
  abstract.
- {cite}`lin2014variable`, {cite}`shi2016regression` — variable selection when
  covariates are compositional and coefficients must sum to zero.
- {cite}`combettes2021regression`, {cite}`mishra2022robust` — general
  log-contrast formulations and their robust variants, as implemented by
  `q2-classo`.

```{figure} _static/figs/compositional_regression_placeholder.png
:name: fig-compositional-regression
:width: 100%

**Placeholder.** This figure will depict the simplex and the log-contrast
constraint — why coefficients are required to sum to zero, and what that means
geometrically for the metagenomics tier.
```

```{bibliography}
:filter: keywords % "compositional"
:labelprefix: C
```

## High-dimensional methods used in this tutorial

With 300 features and 54 samples the sample covariance is singular, so every
estimate needs structure imposed on it. These references justify the modelling
choices made in both worked examples: sparsity in the inverse covariance,
a low-rank term for unobserved confounders, information criteria for choosing the
penalty, and stability-based control of what ends up selected. The optimisation
papers are here too — they matter when a fit fails to converge rather than when
it succeeds.

- {cite}`friedman2008sparse` — introduces the graphical lasso, the estimator
  behind both the toy and 300-ASV network examples.
- {cite}`danaher2014joint` — the joint graphical lasso, used when estimating
  several related networks at once (the multiple-graphical-lasso chapter).
- {cite}`foygel2010extended` — extended BIC for Gaussian graphical models; this
  is the γ-weighted criterion that selects λ = 0.8 in the tier-2 example.
- {cite}`chandrasekaran2010latent` — latent-variable graphical model selection,
  the sparse-plus-low-rank decomposition used in the SLR chapters.
- {cite}`kurtz2019disentangling` — applies that decomposition to microbiome data
  specifically, separating association from hidden environmental structure.
- {cite}`meinshausen2010stability` — stability selection, which we use to control
  which variables survive in the high-dimensional regime.
- {cite}`candes2011robust` — robust PCA, the background for the low-rank and
  latent-component views of the data.
- {cite}`kurtz2015sparse` — SPIEC-EASI; the compositionally-aware network
  inference this book's estimator is most directly compared against.

```{figure} _static/figs/toy_vs_hd_network_placeholder.png
:name: fig-toy-vs-hd
:width: 100%

**Placeholder.** This figure will place the 13-ASV toy network beside the 300-ASV
Atacama network at the same λ, showing what changes when the same estimator is
taken from an inspectable problem to a genuinely high-dimensional one.
```

```{bibliography}
:filter: keywords % "methods"
:labelprefix: M
```

## Microbiome data sources and pipelines

These provide the datasets the examples run on, and the processing context that
determines what an "ASV" actually is. Read at least the denoising and taxonomy
entries: the 300-ASV table did not arrive as a matrix, and choices made upstream
— denoising, reference database, taxonomy assignment — shape every downstream
network. The remaining entries are the wider tooling ecosystem, included because
the book positions its plugins relative to them.

- {cite}`neilson2017significant` — the Atacama soil study; **the source of the
  300-ASV dataset** used throughout tier 2.
- {cite}`meilander2025upcycling`, {cite}`caporaso2025guttosoil` — the
  gut-to-soil study and its supporting data; **the tier-3 metagenomics dataset**.
- {cite}`callahan2016dada2` — DADA2, which produces the amplicon sequence
  variants that are the features in every example here.
- {cite}`quast2012silva` — the SILVA reference database behind the taxonomy
  used to label network nodes.
- {cite}`bokulich2018optimizing` — how QIIME 2 assigns that taxonomy, and why
  classifier choice affects what a node is called.

```{bibliography}
:filter: keywords % "data"
:labelprefix: D
```

## Software and QIIME 2 plugins

These are the tools the tutorial actually invokes. They sit directly on top of
the previous two sections: the plugins wrap the estimators from *High-dimensional
methods* and apply them to data prepared according to *Core compositional
theory*.

- {cite}`Schaipp2021` — **GGLasso**, the solver underneath `q2-gglasso`; used in
  every network fit, toy and 300-ASV alike.
- {cite}`Simpson2021` — **c-lasso**, the solver underneath `q2-classo`; used in
  the metagenomics tier for constrained log-contrast regression.
- {cite}`bolyen2019reproducible` — QIIME 2 itself, the framework both plugins
  are registered against and the source of the artifact/provenance model the
  book relies on.
- {cite}`bokulich2018q2` — `q2-sample-classifier`, referenced when contrasting
  predictive classification with the interpretable models used here.
- {cite}`shaffer2023scnic` — SCNIC, an alternative compositional network tool,
  useful for orientation.
- {cite}`pedregosa2011scikit`, {cite}`scikit-learn` — scikit-learn, whose
  estimators appear in supporting analyses.

```{bibliography}
:filter: keywords % "software"
:labelprefix: S
```

## Everything else

Nothing should appear below this heading. It is a safety net: the sections above
select entries by keyword, so any reference added to `references.bib` **without**
a `keywords` field would silently vanish from this page rather than fail the
build. If entries show up here, tag them with one of `compositional`, `methods`,
`data` or `software` and they will move into place.

```{bibliography}
:filter: not (keywords % "compositional" or keywords % "methods" or keywords % "data" or keywords % "software")
:labelprefix: X
```
