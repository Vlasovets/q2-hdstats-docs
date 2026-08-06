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
Both commands above are therefore exactly the ones that produced the shipped
artifacts. That was not true of releases before the tier-2 regeneration, whose
`clr` table was built without `--p-keep-original-id` — see the next box.

That does not make the question uninteresting: `clr` and `mclr` do **not** give
the same network. If you regenerate the correlation matrix yourself, regenerate
it from the transformed table you intend to report, and say which one it was.
```

```{note}
**The bundle ships real feature IDs, and this was not always true.**
`atacama-top-300-clr.qza` is built with `--p-keep-original-id`, so its features
carry the original 32-character hexadecimal IDs;
`atacama-top-300-correlation.qza` inherits them from it (`calculate-covariance`
has no such parameter — it takes whatever labels its input carries). Both join
directly against `atacama-taxonomy-silva138.qza`: all 300 features resolve, with
no mapping step.

Earlier releases of this bundle were produced *without* that flag, so their
features were sequential `ASV-1` … `ASV-300` names. If you have an older copy,
regenerate it rather than trying to map the names back. `ASV-k` is assigned by
**position** in the abundance ranking, and 209 of these 300 features share a
total-abundance value with another feature, so the ranking does not determine
which organism got which number — see the warning in
[Interpretation](06_interpretation.md).

Relabelling changes nothing numerically. The graphical-lasso objective is
invariant under simultaneous row/column permutation, so λ = 0.8, 216 edges and
eBIC 16130.0988 are identical either way.

It does change one visible thing: a 32-character ID is unusable as a network
node label or a heatmap tick. Figures should render a short display name —
the deepest informative rank plus a slice of the ID, e.g.
`Rubrobacter (a7b877)` — while keeping the full ID as the key. A 5-character
prefix is already unique across these 300 features. `scripts/export_network.py`
emits exactly this as a `display` column and refuses to run if two display names
would collide.
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
without re-deriving them. The commands above are the ones that produced them, so
running the chain yourself reproduces the bundle. Still, check each artifact's
own provenance (`qiime tools peek` / the `provenance/` directory of the unzipped
`.qza`) rather than assuming the prose and the bundle agree — that assumption is
what let an earlier release ship with `ASV-k` labels while the prose described a
chain that would have produced feature IDs.

```{note}
File sizes and checksums are not transcribed here — they live in the bundle
manifest, which is generated from the published files rather than by hand and is
authoritative if the two ever disagree. See
[Download the tutorial data](../00_getting_started/03_download_data.md).
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
