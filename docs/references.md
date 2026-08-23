# References

The entries below are grouped by the role they play in the book: the theory the models
rest on, the estimators the chapters run, the data and pipelines the examples consume, and
the software that implements it all. Most readers need one or two entries from a section
rather than the whole section.

## Reading order for each worked example

- **The 13-ASV example.** Start with {cite}`aitchison1982statistical` and
  {cite}`gloor2017microbiome` for why counts must be treated as compositions, then
  {cite}`friedman2008sparse` for the graphical lasso itself. The example is those two
  ideas composed, on a dataset small enough to check the estimator's output by eye.
- **The 300-ASV Atacama example.** {cite}`friedman2008sparse` for the estimator,
  {cite}`foygel2010extended` for the eBIC model selection that picks λ,
  {cite}`chandrasekaran2010latent` and {cite}`kurtz2019disentangling` for the
  latent/low-rank decomposition, and {cite}`meinshausen2010stability` for selection
  stability. For the data see {cite}`neilson2017significant`, and
  {cite}`callahan2016dada2` with {cite}`quast2012silva` for how the ASVs and their
  taxonomy were produced.
- **The shotgun metagenomics example.** Compositional *regression* replaces network
  estimation here: {cite}`aitchison1984log` for the log-contrast formulation,
  {cite}`lin2014variable` and {cite}`shi2016regression` for variable selection under the
  zero-sum constraint, and {cite}`combettes2021regression` with {cite}`mishra2022robust`
  for the proximal and robust formulations `q2-classo` implements. The profiles come from
  mOTUs {cite}`ruscheweyh2022motus` and the study is {cite}`baitong2022map`.

## Core compositional theory

Microbiome counts carry no absolute scale: sequencing depth is arbitrary, so only *ratios*
between features are meaningful. These references establish that, and give the log-ratio
and log-contrast machinery that lets ordinary multivariate methods be applied to data
living on the simplex. They underpin every transformation in the book — `clr`, `mclr`, and
the zero-sum constraint in the regression chapters.

- {cite}`aitchison1982statistical` — the foundational treatment of compositional data, and
  the reason the 13-ASV example transforms counts before estimating anything.
- {cite}`aitchison1984log` — introduces log-contrast models, the direct ancestor of the
  constrained regression used in the shotgun chapter.
- {cite}`gloor2017microbiome` — the accessible argument for why compositionality is not
  optional in microbiome analysis. Read this first if the others feel abstract.
- {cite}`lin2014variable`, {cite}`shi2016regression` — variable selection when covariates
  are compositional and coefficients must sum to zero.
- {cite}`combettes2021regression`, {cite}`mishra2022robust` — general log-contrast
  formulations and their robust variants, as implemented by `q2-classo`.

```{bibliography}
:filter: keywords % "compositional"
:labelprefix: C
```

## High-dimensional methods used in this tutorial

With 300 features and 54 samples the sample covariance is singular, so every estimate
needs structure imposed on it. These references justify the modelling choices made in the
worked examples: sparsity in the inverse covariance, a low-rank term for unobserved
confounders, information criteria for choosing the penalty, and stability-based control of
what ends up selected. The optimisation papers matter when a fit fails to converge rather
than when it succeeds.

- {cite}`friedman2008sparse` — introduces the graphical lasso, the estimator behind both
  the 13-ASV and 300-ASV network examples.
- {cite}`danaher2014joint` — the joint graphical lasso, for estimating several related
  networks at once (the multiple-graphical-lasso chapter).
- {cite}`foygel2010extended` — extended BIC for Gaussian graphical models. This is the
  γ-weighted criterion that selects λ = 0.8 on the Atacama data.
- {cite}`chandrasekaran2010latent` — latent-variable graphical model selection, the
  sparse-plus-low-rank decomposition used in the SLR chapters.
- {cite}`kurtz2019disentangling` — applies that decomposition to microbiome data
  specifically, separating association from hidden environmental structure.
- {cite}`meinshausen2010stability` — stability selection, which controls which variables
  survive in the high-dimensional regime.
- {cite}`candes2011robust` — robust PCA, the background for the low-rank and
  latent-component views of the data.
- {cite}`kurtz2015sparse` — SPIEC-EASI, the compositionally-aware network inference this
  book's estimator is most directly compared against.

```{bibliography}
:filter: keywords % "methods"
:labelprefix: M
```

## Microbiome data sources and pipelines

These provide the datasets the examples run on, and the processing context that determines
what a feature actually is. Read at least the denoising and taxonomy entries: the 300-ASV
table did not arrive as a matrix, and choices made upstream — denoising, reference
database, taxonomy assignment — shape every downstream network.

- {cite}`neilson2017significant` — the Atacama soil study, and the source of the 300-ASV
  dataset.
- {cite}`baitong2022map` — the MAP preterm-infant study, and the source of the shotgun
  metagenomes profiled in
  [Shotgun Metagenomics](chapters/04_highdim_atacama/07_motus_shotgun.md).
- {cite}`callahan2016dada2` — DADA2, which produces the amplicon sequence variants that
  are the features in the 13-ASV and 300-ASV examples.
- {cite}`quast2012silva` — the SILVA reference database behind the taxonomy used to label
  network nodes.
- {cite}`bokulich2018optimizing` — how QIIME 2 assigns that taxonomy, and why classifier
  choice affects what a node is called.

```{bibliography}
:filter: keywords % "data"
:labelprefix: D
```

## Software and QIIME 2 plugins

These are the tools the chapters invoke. They sit on top of the previous two sections: the
plugins wrap the estimators from *High-dimensional methods* and apply them to data
prepared according to *Core compositional theory*.

- {cite}`Schaipp2021` — **GGLasso**, the solver underneath `q2-gglasso`, used in every
  network fit.
- {cite}`Simpson2021` — **c-lasso**, the solver underneath `q2-classo`, used for
  constrained log-contrast regression.
- {cite}`ruscheweyh2022motus` — **mOTUs**, the marker-gene profiler that produces the
  species-level table in the shotgun chapter.
- {cite}`bolyen2019reproducible` — QIIME 2 itself, the framework both plugins are
  registered against and the source of the artifact and provenance model the book relies
  on.
- {cite}`bokulich2018q2` — `q2-sample-classifier`, for contrasting predictive
  classification with the interpretable models used here.
- {cite}`shaffer2023scnic` — SCNIC, an alternative compositional network tool.
- {cite}`pedregosa2011scikit`, {cite}`scikit-learn` — scikit-learn, whose estimators
  appear in supporting analyses.

```{bibliography}
:filter: keywords % "software"
:labelprefix: S
```

## Everything else

Nothing should appear below this heading. It is a safety net: the sections above select
entries by keyword, so any reference added to `references.bib` **without** a `keywords`
field would silently vanish from this page rather than fail the build. If entries show up
here, tag them with one of `compositional`, `methods`, `data` or `software` and they will
move into place.

```{bibliography}
:filter: not (keywords % "compositional" or keywords % "methods" or keywords % "data" or keywords % "software")
:labelprefix: X
```
