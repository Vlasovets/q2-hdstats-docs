# Shotgun Metagenomics: mOTUs Profiles

Everything so far in this tier has been 16S amplicon data from soil. This page runs the
same two plugins on **shotgun metagenomes from human infants**, profiled to species
level with [mOTUs](https://github.com/motu-tool/mOTUs) {cite}`ruscheweyh2022motus`.

Nothing about the workflow changes. `transform-features`, `calculate-covariance`,
`solve-problem`, `add-taxa` and `regress` take the same arguments and mean the same
things. What changes is what a *feature* is — a marker-gene-defined species rather than
an exact sequence variant — and, more consequentially, how much data sits behind each
count. That difference turns out to decide which parameter values work, and it is the
reason this page exists rather than a sentence saying "shotgun tables also work".

## The study

| | |
|---|---|
| Qiita study | [13241](https://qiita.ucsd.edu/public/?study_id=13241) |
| ENA project | [PRJEB52147](https://www.ebi.ac.uk/ena/browser/view/PRJEB52147) |
| Publication | Bai-Tong et al., *Scientific Reports* 2022 {cite}`baitong2022map` |
| Design | preterm infants in two San Diego NICUs, 9 with a mother who has asthma and 9 without |

The **MAP study** (Microbiome, Atopic disease, and Prematurity) asked whether **maternal
asthma leaves a mark on the preterm infant gut** — in the microbiome, in the metabolome,
or both. Its answer was split: the authors found metabolomic differences between the two
groups, and reported that maternal asthma *"did not play an observable role in shaping
the infant gut microbiome during the study period"*.

That published null matters for this page, and we return to it at the end rather than
discovering it by accident.

```{note}
**Anonymous access uses a different URL than the one you will find first.**
`qiita.ucsd.edu/study/description/13241` redirects to the Qiita landing page when you
are not logged in, which looks like the study is private. It is not — the public view is
at `qiita.ucsd.edu/public/?study_id=13241`, and the download endpoints under
`qiita.ucsd.edu/public_download/` work without an account.
```

### What we took, and what we left

The study has **two sequencing preparations over the same samples**: 16S V4 amplicon on a
MiSeq, and whole-genome shotgun on a NovaSeq. Only the shotgun half can be profiled with
mOTUs. Pulling PRJEB52147 wholesale gets you 192 runs of which half are amplicon.

Of the shotgun half we kept only the **fecal** samples — the study also sequenced breast
milk and meconium, which are different niches and would confound a single covariance
estimate with sample type. That leaves **36 runs from 18 infants**, each contributing an
early (~11–16 d) and a late (~26–42 d) sample.

## How the table was made

The profiling is described here rather than asked of you: it needs 13 GiB of reads and
the 2.9 GiB mOTUs reference database. The stages that do it are
`analysis/slurm/32`–`35` in the book's repository.

1. **Select the runs.** Join the Qiita sample template to the ENA run report on the
   stripped `sample_alias`, keep `library_strategy == WGS` and `sample_type == feces`.
2. **Fetch and verify.** 72 FASTQ files, checked against the MD5s ENA publishes.
3. **Profile.** `qiime motus profile` per sample, then merge. mOTUs counts only reads
   landing on ten universal marker genes, so a 2 M-read run yields a few hundred assigned
   inserts — this is the single most important property of the resulting table.
4. **Drop `unassigned`.** mOTUs emits a literal `unassigned` row for reads it could not
   place. It is the complement of everything profiled, and including it in a log-ratio
   makes the "composition" a mixture of taxa and not-taxa.
5. **Filter.** Two samples fall below 100 assigned counts and are dropped, leaving 34.
   Features are then ranked **by total abundance** and the top 100 kept.

```{important}
**Rank by abundance, not by prevalence.** Both are one-line filters and they behave
completely differently here.

The unfiltered table is 84% zeros. Two features that are absent in the same 28 samples
have a large sample correlation driven entirely by shared zeros, and a *prevalence*
filter keeps such a pair as long as each clears the threshold. An *abundance* filter
removes them, because a feature whose entire signal is a handful of 1s carries almost no
counts.

This is not a stylistic preference. Under prevalence filtering the eBIC criterion on this
dataset selects the **empty graph** at every threshold tried. Under abundance filtering it
selects a network that beats a permutation null by three orders of magnitude. Same data,
same commands, different filter.
```

## The shape of the problem

```{csv-table} What you are modelling
:file: ../../../analysis/results/tables/motus-tutorial-shape.tsv
:delim: tab
:header-rows: 1
:widths: 45, 25
```

$p = 100$ against $n = 34$ is squarely the regime this tier is about: the sample
covariance has rank at most 33, so it is singular, and the penalty decides the model.

Two properties differ sharply from the Atacama table and both are worth carrying:

- **The table is sparse.** 85% of cells are zero, against a dense soil table. This is why
  `mclr` rather than `clr` — the modified transform is built for exactly this.
- **Samples are paired within infants.** Two samples per infant, 18 infants. Nothing in
  the network estimate cares, but the regression very much does, and
  [Reading the result honestly](#reading-the-result-honestly) is where that bill comes due.

## Getting the files

Three small artifacts, committed with the book:

```bash
BASE=https://raw.githubusercontent.com/Vlasovets/q2-hdstats-docs/main/docs/_data/motus
curl -L -O "${BASE}/motus-top100-table.qza"      # 20 KB, FeatureTable[Frequency]
curl -L -O "${BASE}/motus-top100-taxonomy.qza"   # 16 KB, FeatureData[Taxonomy]
curl -L -O "${BASE}/motus-outcomes.tsv"          # 1 KB, sample metadata
```

## Network inference with q2-gglasso

### Transform and build the covariance

```bash
qiime gglasso transform-features \
    --i-table motus-top100-table.qza \
    --i-taxonomy motus-top100-taxonomy.qza \
    --p-transformation mclr \
    --p-keep-original-id \
    --o-transformed-table mclr.qza

qiime gglasso calculate-covariance \
    --i-table mclr.qza \
    --p-method scaled \
    --o-covariance-matrix correlation.qza
```

```
Saved FeatureTable[Frequency] to: mclr.qza
Saved PairwiseFeatureData to: correlation.qza
```

`--i-taxonomy` is required by `transform-features` even though the function body never
reads it; the taxonomy earns its keep later, in `add-taxa`.

### Select $\lambda_1$

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix correlation.qza \
    --p-n-samples 34 \
    --p-no-latent \
    --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 \
    --p-gamma 0.15 \
    --o-solution motus-sgl-path.qza
```

```
Saved GGLassoProblem to: motus-sgl-path.qza
```

```{figure} ../../images/png/generated/motus-lambda-path.png
:name: fig-motus-lambda
:width: 100%

eBIC along the $\lambda_1$ path at $\gamma = 0.15$ (solid, left axis) with the edge count
(dashed, right axis). The minimum is interior — both neighbouring grid points score
worse — and sits at $\lambda_1 = 0.30$, giving **481 edges** among 4,950 possible pairs,
a density of 9.7%.
```

**$\gamma = 0.15$, not the $\gamma = 0.3$ used elsewhere in this tier.** That is a
reported modelling choice, not a default someone forgot to change. At $\gamma \geq 0.20$
this dataset selects the empty graph; the criterion's model-space penalty overwhelms the
likelihood at $n = 34$. Report $\gamma$ with any network, always, and especially here.

```bash
qiime gglasso summarize \
    --i-solution motus-sgl-path.qza \
    --p-label-size 6pt \
    --o-visualization motus-sgl-summary.qzv
```

📊 **[motus-sgl-summary.qzv](../../_static/qzv/motus-sgl-summary.qzv)** — the solution
path, sparsity and eBIC at every grid point. Open with
[view.qiime2.org](https://view.qiime2.org).

### Add a latent block

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix correlation.qza \
    --p-n-samples 34 \
    --p-latent \
    --p-lambda1-min 0.30 --p-lambda1-max 0.30 --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min 5 --p-mu1-max 5 --p-n-mu1 1 \
    --o-solution motus-slr-mu5.qza
```

```
Saved GGLassoProblem to: motus-slr-mu5.qza
```

All three penalties are pinned to single values. Leaving `--p-lambda2-*` unset
substitutes a five-point default path, which silently turns this single fit into a
model-selection run that can wander away from the $\lambda_1$ you just chose.

```{figure} ../../images/png/generated/motus-mu-rank.png
:name: fig-motus-mu-rank
:width: 100%

The rank is an **output**, not a setting — `--p-rank` is registered but raises on every
released GGLasso, so you reach a rank by tuning $\mu_1$ and reading back what you got.
$\mu_1 = 5$ gives rank 3 and leaves 217 of the 481 edges in the sparse component.
```

**Three latent dimensions absorb 55% of the network.** That is a much larger effect than
the Atacama example, where a rank-2 block removes 14 of 216 edges (6.5%). A dominant
latent gradient is what a maturing infant gut should look like — but do not carry
[Choosing the Latent Rank](03_slr_ranks.md)'s framing across, because there almost
nothing was explained away and here most of it is.

The eigenvalues of the latent block at $\mu_1 = 5$ are 2.428, 0.711 and 0.029. The third
is small enough that this is effectively rank 2 with a marginal third direction.

```bash
qiime gglasso pca \
    --i-solution motus-slr-mu5.qza \
    --i-table mclr.qza \
    --m-sample-metadata-file motus-outcomes.tsv \
    --p-n-components 3 \
    --o-visualization motus-latent-pca.qzv
```

📊 **[motus-latent-pca.qzv](../../_static/qzv/motus-latent-pca.qzv)** — samples projected
onto the latent components, coloured by metadata.

`--p-n-components` must not exceed the achieved rank, and `pca` fails without
`--m-sample-metadata-file` despite the signature marking it optional.

## Log-contrast regression with q2-classo

The outcome is **`host_age_days`**, postnatal age at sampling, 10–47 days.

It is chosen because it is the only variable that **varies within an infant**. Every
other numeric column — gestational age at birth, birth weight, maternal age — is constant
per infant, so 34 samples would carry no more information than 18. Age varies by 13–29
days within each infant, which is what makes the paired design informative rather than
merely repetitive.

```bash
qiime classo transform-features \
    --i-features motus-top100-table.qza \
    --o-x classo-x.qza

qiime classo add-taxa \
    --i-features classo-x.qza \
    --i-taxa motus-top100-taxonomy.qza \
    --o-x classo-x-trac.qza \
    --o-aweights classo-w-trac.qza
```

```
Saved FeatureTable[Design] to: classo-x.qza
Saved FeatureTable[Design] to: classo-x-trac.qza
Saved Weights to: classo-w-trac.qza
```

`add-taxa` rebuilds the design as $\log(X)A$, where $A$ aggregates features up the
taxonomy, so coefficients attach to **clades** instead of single species. The design
grows from **100 columns to 133** — 33 internal nodes, one per taxonomic group present.

This works here because the mOTUs lineage has a property trac needs and not every
taxonomy has: its leaf rank is `m__<mOTU_id>`, so all **100 leaves are unique** even
where two mOTUs share a species name. A taxonomy with duplicate leaves silently merges
distinct features into one node.

```bash
qiime classo regress \
    --i-features classo-x-trac.qza \
    --i-weights classo-w-trac.qza \
    --m-y-file motus-outcomes.tsv \
    --m-y-column host_age_days \
    --p-concomitant \
    --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
    --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed \
    --o-result motus-regress-age.qza
```

```
Saved CLASSOProblem to: motus-regress-age.qza
```

```{figure} ../../images/png/generated/motus-trac-coefficients.png
:name: fig-motus-trac
:width: 100%

The three clades the cross-validated refit keeps, of 134 candidates. Proteobacteria
rising and Actinobacteria falling with postnatal age is the textbook preterm gut
trajectory, which is a good sign: the method recovers something a neonatologist would
recognise.
```

**`--p-no-cv-one-se` is doing real work here.** The one-standard-error rule returns the
sparsest model within one standard error of the cross-validated minimum, and on this
dataset that is the intercept-only model — it selects **1** coefficient where the CV
minimum selects **9** on an identical CV curve. Neither is wrong; the rule is more
conservative than this sample size can afford.

```bash
qiime classo summarize \
    --i-problem motus-regress-age.qza \
    --p-maxplot 40 \
    --o-visualization motus-regress-summary.qzv
```

📊 **[motus-regress-summary.qzv](../../_static/qzv/motus-regress-summary.qzv)** — the
$\lambda$ path, the CV curve and the selected coefficients. `--p-maxplot` matters after
`add-taxa`, because the design is wider than the table you started from.

## Reading the result honestly

Four things about these numbers, none of which should be discovered by a reader later.

**The network is real, but it is dense.** 481 edges is 9.7% of all pairs, against 0.48%
for the 300-ASV Atacama network. Estimating 481 edges from 34 samples is
over-parameterised however good the criterion. What justifies believing there is
*structure* is not the count but a null: re-solving five times with every feature column
permuted independently across samples — destroying all between-feature dependence while
preserving each feature's exact marginal, zeros included — yields a mean of **0.2 edges**
and a maximum of 1.

**Cross-validation here leaks, and the direction is not what you expect.** Each infant
contributes two samples and `--p-cv-subsets` splits at random with no grouping option, so
one of a pair lands in training and the other in testing. Holding out whole infants
instead:

| outcome | | baseline | subject-holdout | sample-holdout | leakage |
|---|---|---|---|---|---|
| `host_age_days` | R² | 0.0 | **+0.164** | −0.041 | −0.205 |
| `diagnosis` | accuracy | 0.500 | **0.471** | 0.588 | +0.118 |

The signs are opposite, and the reason is structural. `diagnosis` is a property of the
*infant*, so a paired sample gives the label away and random folds flatter the classifier
by +0.118. `host_age_days` varies *within* an infant, so a paired sample reveals nothing
and the leaky estimate is merely noisier. "Grouped CV always inflates" is the intuition
most people carry, and it is wrong.

```{important}
q2-classo exposes `--p-cv-subsets` and `--p-cv-seed` but **no grouping parameter**. The
subject-held-out figures above were produced by driving the `classo` API directly
(`analysis/slurm/42_subject_holdout.sh`). If your design has repeated measures, the
plugin's cross-validation cannot express what you need, and you should say so in your
methods rather than quote its number.
```

**The classification reproduces the published null.** Held out by infant, predicting
maternal asthma from the microbiome scores 0.471 against a 0.500 majority baseline —
worse than guessing. With 9 versus 9 infants a 95% interval on any accuracy here spans
roughly ±0.23, so this is not evidence of *absence*; it is an absence of evidence,
consistent with what {cite}`baitong2022map` reported from a much fuller analysis. It is
included because a tutorial that only ever shows the method succeeding teaches the wrong
reflex.

**Eighteen infants is the binding constraint**, not the profiler and not the filter. It
is why $\gamma$ has to come down to 0.15, why the network is dense, and why the R² of
0.164 is a real result rather than a strong one.

## What you should have now

| artifact | what it is |
|---|---|
| `motus-sgl-path.qza` | the $\lambda_1$ path, eBIC at 20 grid points |
| `motus-slr-mu5.qza` | the selected model: $\lambda_1 = 0.30$, rank 3 |
| `motus-regress-age.qza` | the trac regression of postnatal age |
| three `.qzv` files | linked above |

The headline: **481 edges at $\lambda_1 = 0.30$**, reduced to **217** once three latent
factors are allowed, and **three clades** carrying postnatal age at a
subject-held-out R² of **+0.164**.

Compare [Interpretation](06_interpretation.md), which asks the same closing question of
the Atacama data — whether the structure the network attributes to unobserved drivers is
the structure the log-contrast model exploits.
