# The 300-ASV Dataset

This page is bookkeeping: where the Tier 2 table comes from, how the taxonomy
and the transformed artifacts were made, and which file is which. No model is
fitted here — that starts in [Selecting lambda](02_model_selection.md).

It is worth reading rather than skipping, because two of the choices below
(which ASVs were kept, and which log-ratio transform was applied) change the
network more than any penalty parameter does.

## Provenance

The starting point is the same Atacama soil study used throughout this tutorial
{cite}`neilson2017significant`, processed exactly as in the
[QIIME 2 Atacama soils tutorial](https://amplicon-docs.qiime2.org/en/latest/tutorials/atacama-soils.html):
demultiplexing with `q2-demux`, denoising with DADA2 {cite}`callahan2016dada2`
via `q2-dada2`. That is the published, reproducible pipeline; nothing about it
is specific to this tutorial, and the ASV table it produces is the one we
subset.

From that table, two reductions were applied:

1. **Sample matching.** Only samples with a complete row in the sample metadata
   were kept, leaving **$n = 54$ samples**.
2. **Abundance filtering.** The **300 ASVs with the greatest total abundance
   across those samples** were retained, giving **$p = 300$ features**.

```{note}
The exact list of the 54 sample IDs and the 300 feature IDs is fixed by the
artifact bundle rather than restated here, so that a reader can check membership
against the file they actually downloaded instead of against prose. The
per-sample and per-feature summary statistics of the subset are **pending
verification against QIIME 2 2026.7**.
```

### Why 300, and why by total abundance

Two things are being traded off. Keeping more ASVs makes the network more
complete but pushes $p/n$ further into a regime where the empirical covariance is
badly conditioned and the graphical lasso is doing more regularizing than
estimating. Keeping fewer makes the estimate more stable but starts deleting the
organisms you are trying to say something about. At $p = 300$, $n = 54$ the ratio
is already about 5.6 features per sample — high, but within the range where a
sparse estimator is the accepted tool rather than a fig leaf.

Filtering by **total abundance** rather than by prevalence, variance, or
association with an outcome is deliberate. It is a marginal, outcome-blind
criterion: it cannot leak information about any of the environmental variables
that are later used as regression targets in
[Log-Contrast Models at Scale](05_classo_cv.md). A filter chosen with the
outcome in view would invalidate every cross-validated $R^2$ on that page.

```{warning}
Abundance filtering is not a neutral pre-processing step for compositional data.
The CLR transform divides each sample by the geometric mean **of the features
that are present in the table**, so removing features changes the reference and
therefore changes every transformed value — and with it the covariance and the
edges {cite}`aitchison1982statistical,gloor2017microbiome`. A network estimated
on the top 300 ASVs is a statement about that sub-composition, not about the full
community. This is the single most important caveat to carry into
[Interpretation](06_interpretation.md).
```

## Taxonomy

Taxonomic labels were assigned with `q2-feature-classifier`
{cite}`bokulich2018q2` using the naive Bayes classifier trained on
SILVA 138 {cite}`quast2012silva`, and ship as
`atacama-taxonomy-silva138.qza` (`FeatureData[Taxonomy]`).

The taxonomy is used for **interpretation only** — putting a genus name on a
node or on a selected coefficient. It plays no part in estimating the network.

```{note}
`qiime gglasso transform-features` requires `--i-taxonomy`, but the function body
never reads it. This is a registration wart, not a real dependency: the artifact
must be a valid `FeatureData[Taxonomy]`, and its contents do not affect the
transformed table. See
[Troubleshooting](../90_reference/04_troubleshooting.md).

Because the input is required anyway, pass the *correct* taxonomy rather than a
convenient one — the same artifact is what you will need later to map feature IDs
back to names.
```

## The transformed table

`transform-features` moves the counts out of the simplex before any covariance
is computed. Two transforms are available, and they differ in how they handle
zeros:

- **`clr`** replaces every zero with `--p-pseudo-count`, renormalizes each sample
  to the simplex, then takes the centred log-ratio. Only zeros are imputed, but
  in a table where ~94% of entries are zero the pseudo-count still sets the floor
  of the whole log-ratio geometry.
- **`mclr`** {cite}`yoon2019microbial` skips zero imputation entirely — no
  pseudo-count is added, and `--p-pseudo-count` is ignored — and log-transforms
  only the positive entries. For sparse tables this avoids inventing counts that
  were never observed.

```bash
qiime gglasso transform-features \
    --i-table atacama-top-300-table.qza \
    --i-taxonomy atacama-taxonomy-silva138.qza \
    --p-transformation clr \
    --p-pseudo-count 1 \
    --p-no-add-metadata \
    --p-keep-original-id \
    --o-transformed-table atacama-top-300-clr.qza
```

Swap `--p-transformation clr` for `--p-transformation mclr` to produce the
modified variant; drop `--p-pseudo-count`, which that branch does not use.

**Explanation:**

- `--p-transformation`: `clr` or `mclr`. Neither plugin declares `Choices()`, so a
  typo here is a runtime `ValueError`, not a command-line error.
- `--p-pseudo-count`: substituted for zeros before normalization, on the `clr`
  branch only.
- `--p-no-add-metadata`: do **not** append the numeric metadata columns to the
  table. Appending them is the adaptive-model workflow of
  [Adaptive Graphical Lasso](../02_lowdim_gglasso/04_adaptive_glasso.md); for the
  Tier 2 network the covariates are deliberately left *out* of the table, so that
  the latent block has something to find.
- `--m-sample-metadata-file` and `--p-scale-metadata` are omitted deliberately.
  Both matter only when metadata is appended, so with `--p-no-add-metadata` they
  are inert; leaving them out keeps the command honest about what actually
  affects the output.
- `--p-keep-original-id`: keep the real feature IDs instead of relabelling to
  sequential `ASV-n` names. At $p = 300$ this matters — without it you cannot
  trace a node in the network or a coefficient in a log-contrast model back to
  its taxonomy.

```{note}
The output is registered as `FeatureTable[Frequency]` even though its values are
log-ratios and can be negative. That is a typing convenience, not a claim about
the contents — but it does mean the QIIME 2 type system will happily let you feed
a transformed table to an action expecting counts, or transform an
already-transformed table a second time. Keep the naming unambiguous.
```

The transformed table is not a dead end after the covariance is computed:
`qiime gglasso pca` takes it as `--i-table` and multiplies it against the
eigenvectors of the low-rank block, so the feature set and its order must still
match the covariance matrix. That is worked through in
[Latent Components & Covariates](04_latent_pca.md).

## The covariance matrix

```bash
qiime gglasso calculate-covariance \
    --i-table atacama-top-300-clr.qza \
    --p-method scaled \
    --p-bias \
    --o-covariance-matrix atacama-top-300-correlation.qza
```

`--p-method scaled` divides out the diagonal, so the result is a **correlation**
matrix. That is what the graphical lasso in this tutorial is applied to, and it
is why a single $\lambda_1$ is meaningful across all 300 features: on the
unscaled covariance, features with larger variance would effectively be penalized
less than features with smaller variance, and a scalar penalty would mean
different things in different rows.

`--p-bias` (the default) uses the $1/N$ normalization that falls out of the
Gaussian log-likelihood, which is the quantity the solver's objective is written
against.

```{important}
**The shipped correlation matrix is `clr`-backed.** This was an open question in
earlier drafts, because some of the reference text describes an mclr-transformed
table. The provenance of the shipped `atacama-top-300-correlation.qza` settles
it: it records `calculate_covariance` with `method: scaled`, `bias: true`, run on
`atacama-top-300-clr.qza`, whose own provenance records `transformation: clr`.
The `calculate-covariance` command above is therefore exactly the one that
produced it; the `transform-features` command above is not (see the next box).

That does not make the question uninteresting: `clr` and `mclr` do **not** give
the same network. If you regenerate the correlation matrix yourself, regenerate
it from the transformed table you intend to report, and say which one it was.
```

```{warning}
**The shipped `atacama-top-300-clr.qza` predates `--p-keep-original-id`.** Its
provenance records no `keep_original_id` parameter at all, and its features are
sequential `ASV-1` … `ASV-300` names rather than the original feature IDs. Every
downstream artifact in the bundle inherits those names.

Those names cannot be mapped back to feature IDs after the fact. `ASV-k` is
assigned by **position** in the abundance ranking, and 209 of these 300 features
share a total-abundance value with another feature, so the ranking does not
determine which organism got which number — see the warning in
[Interpretation](06_interpretation.md). The only way to get a working taxonomy
join is to rebuild the transformed table with `--p-keep-original-id`, as the
command above shows.
```

## Metadata

`sample-metadata.tsv` is the standard Atacama sample metadata. It carries the
continuous environmental variables used as regression outcomes in
[Log-Contrast Models at Scale](05_classo_cv.md) — elevation, pH, soil
temperature and humidity summaries, and several assay measurements — plus two
categorical variables that are the natural grouping factors for this study:

| Column | Levels | Rows in the metadata file |
|---|---|---|
| `transect-name` | Baquedano / Yungay | 32 / 43 |
| `vegetation` | no / yes | 40 / 35 |

```{warning}
Those counts describe **the metadata file**, which covers more samples than the
modelled table: 75 rows against the $n = 54$ samples that survived matching and
entered the network. The per-group counts *within the 54 modelled samples* are a
subset of these and are **pending confirmation** — do not quote 32/43 or 40/35 as
group sizes for the Tier 2 analysis.
```

Two practical consequences of the metadata layout:

- `qiime gglasso pca` filters the metadata to **numeric** columns before it does
  anything. `transect-name` and `vegetation` are dropped, and naming either in
  `--p-color-by` fails. Recode them numerically first if you want to see them.
  This is worked through in
  [Latent Components & Covariates](04_latent_pca.md).
- The same numeric filter governs which columns `transform-features
  --p-add-metadata` would append, and therefore what `summarize --p-n-cov` must
  be set to. See [Summarizing a Solution](../02_lowdim_gglasso/08_summarize.md).

## Artifact inventory

| File | QIIME 2 type | What it is |
|---|---|---|
| `atacama-top-300-table.qza` | `FeatureTable[Frequency]` | the 300 × 54 count table |
| `atacama-top-300-clr.qza` | `FeatureTable[Frequency]` | log-ratio transformed counts, input to `calculate-covariance` and to `pca` |
| `atacama-top-300-correlation.qza` | `PairwiseFeatureData` | 300 × 300 correlation matrix, input to `solve-problem` |
| `atacama-taxonomy-silva138.qza` | `FeatureData[Taxonomy]` | SILVA 138 assignments |
| `sample-metadata.tsv` | metadata | sample metadata, including the regression outcomes |

The three derived artifacts are shipped so that the modelling pages can be run
without re-deriving them. The commands above describe the **intended recompute**,
not byte-for-byte how the current bundle was built: the shipped
`atacama-top-300-clr.qza` was made before `--p-keep-original-id` existed, as the
warning further up records. Check each artifact's own provenance
(`qiime tools peek` / the `provenance/` directory of the unzipped `.qza`) rather
than assuming the prose and the bundle agree.

```{note}
File sizes, checksums and the exact contents of each artifact are **pending
verification against QIIME 2 2026.7**; they will come from the pinned bundle
manifest rather than being transcribed here.
```

## Getting the files

All Tier 2 artifacts are part of the tutorial data bundle described in
[Download the Tutorial Data](../00_getting_started/03_download_data.md).

```{warning}
The bundle will be published as the Zenodo record *q2-hdstats-tutorial-data*
v1.0.0. **The DOI has not been minted yet**, so the download URL on that page is
currently the placeholder `ZENODO_DOI_PENDING`.
```

With the files in place, continue to [Selecting lambda](02_model_selection.md).
