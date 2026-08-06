# The Dataset

This page is bookkeeping and honesty: where the table comes from, how to fetch
it, how to look at it before modelling it, and why its depth limits what the next
two pages can claim. No model is fitted here — that starts in
[Network Inference](02_network.md).

## Provenance

The data come from the **gut-to-soil axis** composting study
{cite}`meilander2025upcycling`, whose supporting data are deposited on Zenodo
{cite}`caporaso2025guttosoil`, and which tracks the microbial community of human
faecal material through composting and into soil. Sequencing is **16S rRNA V4 amplicon** with the
515F/806R primer pair following the Earth Microbiome Project protocol, and the
published tutorial processes the reads into amplicon sequence variants with DADA2
{cite}`callahan2016dada2` inside QIIME 2 {cite}`bolyen2019reproducible`.

What this book uses is the table from the study's own **QIIME 2 tutorial**, not
the full study. That tutorial publishes two feature tables; the one used here is
`asv-table-ms2.qza`, with **335 features across 99 samples**.

```{important}
`asv-table-ms2.qza` is a **10% subsample** of the published data. It was made
small so that a tutorial runs quickly on a laptop, not so that it supports
inference. Everything on the following pages inherits that limitation.
```

## Getting the files

[Download the Tutorial Data](../../00_getting_started/03_download_data.md) is the
canonical download page for the whole book. The tier 3 commands are repeated here
because they are the only ones in the book that work today — tiers 1 and 2 wait
on a Zenodo DOI, tier 3 is served from the upstream tutorial site.

```bash
mkdir -p ~/q2-hdstats-tutorial/data/gut-to-soil
cd ~/q2-hdstats-tutorial/data/gut-to-soil

BASE=https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil
curl -L -O "${BASE}/asv-table-ms2.qza"
curl -L -O "${BASE}/taxonomy.qza"
curl -L -O "${BASE}/sample-metadata.tsv"
```

| File | QIIME 2 type | What it is |
|---|---|---|
| `asv-table-ms2.qza` | `FeatureTable[Frequency]` | the 335 × 99 count table used on these pages |
| `taxonomy.qza` | `FeatureData[Taxonomy]` | taxonomic assignments for the ASVs |
| `sample-metadata.tsv` | metadata | full-study sample metadata (1,660 rows); only 99 rows match `asv-table-ms2.qza` |

Two further files are available from the same location and are **not** used here:
`asv-table.qza`, a different tutorial table of 1069 features × 104 samples, and
`asv-seqs.qza`, the representative sequences. The two tables are not
interchangeable; if a command on the next pages behaves strangely, check first
that you passed `asv-table-ms2.qza`.

```{note}
The tier 3 metadata file is called `sample-metadata.tsv`, and so is the tier 2
one. Keep them in separate directories, as above. Every command in this tier
writes and reads under `data/gut-to-soil/`.
```

## The shape of the problem

Two numbers follow from 335 features and 99 samples, and both are arithmetic
rather than results:

- The precision matrix has $335 \times 334 / 2 = 55{,}945$ off-diagonal entries.
  That is the number of edges the graphical lasso is deciding about, and it is
  the reason $\lambda_1$ has to be chosen with a criterion instead of by eye.
- $p / n \approx 3.4$. The empirical covariance is singular — its rank is at most
  $n - 1 = 98$ — so *some* regularization is not optional, it is what makes the
  problem solvable at all.

That ratio is milder than the Tier 2 Atacama table's $\approx 5.6$, which might
suggest this is the easier problem. It is not, and the reason is depth.

## Why the covariance will be noisy

Per-sample sequencing depth in this table ranges from **3 to 1218 reads, with a
median of 261** — the shallowest samples carry almost no information at all.
Counts that small
have three consequences that all land on the covariance matrix:

**Most entries are zero.** In a table this shallow, the majority of feature-sample
cells are structural or sampling zeros, and there is no way to tell which is
which. Every log-ratio transform has to do *something* with them, and the
something you choose changes the geometry more than any penalty parameter will
later.

**The pseudo-count is not a technicality.** `clr` replaces every zero with
`--p-pseudo-count` and then rescales each sample back to its original total. When
a typical non-zero count is in the single or low double digits, putting every zero
cell on a floor of 1 is a large perturbation, and it lands on the cells you know
least about. The alternative, `mclr`
{cite}`yoon2019microbial`, adds nothing and log-transforms only the positive
entries. For a table like this that is usually the more defensible choice, and it
is what [Network Inference](02_network.md) uses throughout, though the plugin
default is `clr`.

**Multinomial noise looks like correlation.** With low depth, the observed
proportions of two features fluctuate together simply because they compete for
the same finite number of reads. A sparse estimator suppresses weak edges, but it
cannot distinguish a weak real edge from a weak sampling artifact — it only
decides how many of each to keep.

```{important}
Read the networks and regressions in this tier as **a demonstration of the
method**, not as findings about composting. If you want findings, start from the
full study (below) at full depth, and expect to spend the effort on filtering and
on stability checks rather than on the solver.
```

## Look at the table before you model it

Do this even though — especially though — the next pages hand you a working
command sequence. Almost every unpleasant surprise in this workflow is visible in
these two visualizations first.

```bash
cd ~/q2-hdstats-tutorial

qiime feature-table summarize \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --m-metadata-file data/gut-to-soil/sample-metadata.tsv \
    --o-feature-frequencies data/gut-to-soil/gts-feature-frequencies.qza \
    --o-sample-frequencies data/gut-to-soil/gts-sample-frequencies.qza \
    --o-summary data/gut-to-soil/gts-table-summary.qzv

qiime metadata tabulate \
    --m-input-file data/gut-to-soil/sample-metadata.tsv \
    --o-visualization data/gut-to-soil/gts-metadata-summary.qzv
```

Open both at [QIIME 2 View](https://view.qiime2.org/). From the table summary,
note the per-sample frequency distribution — particularly the minimum, since a
sample with a handful of reads contributes a row of near-noise to the
covariance — and the feature detail table, which tells you how many features are
observed in only one or two samples.

From the metadata summary, note **which columns exist, which are numeric and
which are categorical, and how many levels the categorical ones have**. You need
this for three separate decisions later: the grouping variable for the
group-graphical-lasso split in [Network Inference](02_network.md), the outcome
column in [Log-Contrast Regression](03_regression.md), and any covariates you add
to the design.

```{note}
The columns present in this metadata file, their levels and their group sizes are
**not restated in this book**, with one exception: the three columns used by
[Log-Contrast Regression](03_regression.md) were read off the
`sample-metadata.tsv` served by the tutorial site and are named there. Everything
else belongs to the upstream study and has not been checked, and printing an
unverified column name would be worse than making you look. Read them off
`gts-metadata-summary.qzv`.

Read them with the row count in mind: `sample-metadata.tsv` is the **full-study**
metadata (1,660 rows), so `gts-metadata-summary.qzv` tabulates every study
sample, not the 99 in `asv-table-ms2.qza`. The level counts and group sizes shown
there are study-wide and are **not** the group sizes you will get when you split
this table — for that, read the sample list off `gts-table-summary.qzv` instead,
which is generated from the table itself.
```

```{note}
The per-sample depths quoted above are verified against `asv-table-ms2.qza`. The
remaining contents of `gts-table-summary.qzv` — feature prevalences and sparsity —
are pending verification against QIIME 2 2026.7 and are not reproduced here.
```

## Filtering, and what it costs

Rare features are where shallow tables do the most damage: a feature seen in
three samples contributes almost no information about a conditional dependence
but still gets a row and a column in a 335 × 335 matrix. Raising the prevalence
floor is the standard remedy.

```bash
qiime feature-table filter-features \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --p-min-samples 10 \
    --o-filtered-table data/gut-to-soil/gts-table-prev10.qza
```

`--p-min-samples 10` keeps features observed in at least ten samples, roughly a
tenth of the design. Ten is a starting point, not a recommendation: re-run
`qiime feature-table summarize` on the filtered table and check how much of the
total count you discarded before accepting it.

```{note}
Filtering features is **not** a neutral preprocessing step for compositional
data. The CLR transform divides each sample by the geometric mean of the features
*present in the table*, so removing features moves the reference and changes every
transformed value, every covariance entry and therefore every edge. A network
built on a filtered table is a statement about that sub-composition and nothing
wider.

The practical rule is: decide the filter once, before any modelling, on grounds
that have nothing to do with the outcome you will later predict; then keep it
fixed. Changing the filter after seeing a network is how you talk yourself into
one.
```

The pages that follow use the **unfiltered** `asv-table-ms2.qza`, so that the
commands match the file you downloaded and any difference you see is not a
filtering difference. If you filter, substitute `gts-table-prev10.qza`
consistently in every subsequent command, including the per-group tables.

## Taxonomy

`taxonomy.qza` is used in two genuinely different ways in this tier, and it is
worth keeping them apart.

In the network chapter it is **interpretation only** — it puts a name on a node.
It is also passed to `qiime gglasso transform-features`, but only because that
action declares `--i-taxonomy` as a required input and then never reads it; the
transformed table is identical whatever taxonomy you pass. See
[Troubleshooting](../../90_reference/04_troubleshooting.md).

In the regression chapter it is **part of the model**. `qiime classo add-taxa`
turns the taxonomic hierarchy into an aggregation matrix and derives penalty
weights from it, so the taxonomy actively determines which clades the model can
select {cite}`bien2021tree`. A wrong or mismatched taxonomy changes the answer
there, whereas in the network chapter it changes only the labels.

## The full study, if you want to scale out

The complete dataset is Zenodo record **15390940**
{cite}`caporaso2025guttosoil` — **1,660 samples**. That is the right starting
point for an actual analysis, and it is a different computational problem. The
scaling below is written conditionally on a feature count, because the feature
count of the full table is not restated here: at $p = 30{,}000$ features there
are about 450 million candidate edges, and a dense
$30{,}000 \times 30{,}000$ float64 covariance matrix is around 7.2 GB before the
solver allocates anything of its own.

In practice you would not hand that matrix to `solve-problem` directly. You would
agglomerate to a taxonomic rank or filter hard on prevalence first, which brings
$p$ down to something a dense solver can hold, and you would do it for
statistical reasons as well as memory ones — at 1,660 samples and $p = 30{,}000$
the ratio is still about 18 to 1.

```{note}
The number of features in the full-study table is **pending verification** against
the deposited artifacts, so no ASV count is quoted here; $p = 30{,}000$ above is a
round number used to show how the arithmetic scales. The sample count, 1,660, is
from the record. Runtimes and memory figures for either table on QIIME 2 2026.7
are likewise pending verification and are not quoted anywhere in this tier.
```

With the files in place and the metadata read, continue to
[Network Inference](02_network.md).
