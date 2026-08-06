# Multiple Graphical Lasso (GGL / FGL)

Everything so far has estimated **one** network from **one** covariance matrix.
The Multiple Graphical Lasso (MGL) estimates $K$ networks at once from $K$
covariance matrices and couples them, so that an edge supported in several
groups is easier to keep than an edge supported in only one. You use it when
your samples fall into groups that you expect to share most of their structure —
two sampling transects, treated and untreated, vegetated and bare soil — and you
want the differences between the groups rather than $K$ independently noisy
networks.

Two couplings are available through `--p-reg`:

- **`GGL`** — Group Graphical Lasso {cite}`danaher2014joint`. The $\lambda_2$
  penalty is a group-lasso penalty on each edge *across* groups. It pushes edges
  to be present in all groups or absent in all groups, but lets the surviving
  edge weights differ freely. This is the default.
- **`FGL`** — Fused Graphical Lasso. The $\lambda_2$ penalty is on the
  *difference* between corresponding edge weights across groups. It pushes edge
  weights to be numerically equal, which is a stronger and more specific
  assumption.

In both cases $\lambda_1$ still controls sparsity within each group, exactly as
in [Single Graphical Lasso](02_sgl.md), and $\lambda_2$ controls how strongly the
groups are tied together. Setting $\lambda_2 = 0$ recovers $K$ independent
single graphical lassos.

```{warning}
**Read the "Known gaps" section in Step 3 before you invest time here.** The
flags documented on this page all exist and all validate, but in the current
release the MGL path does not close end-to-end through the QIIME 2 CLI. Nothing
on this page has been run against QIIME 2 2026.7, and the gaps below were read
out of the plugin source rather than reproduced.
```

## Step 1: Split the toy table into two instances

The Atacama study samples two transects. We use `transect-name` as the grouping
variable and build one feature table per transect with the standard
`q2-feature-table` filter.

```bash
qiime feature-table filter-samples \
    --i-table data/atacama-counts.qza \
    --m-metadata-file data/sample-metadata.tsv \
    --p-where "[transect-name]='Baquedano'" \
    --o-filtered-table data/atacama-counts-baquedano.qza

qiime feature-table filter-samples \
    --i-table data/atacama-counts.qza \
    --m-metadata-file data/sample-metadata.tsv \
    --p-where "[transect-name]='Yungay'" \
    --o-filtered-table data/atacama-counts-yungay.qza
```

```{note}
`--p-where` is a SQLite `WHERE` clause over the metadata, so a column name
containing a hyphen must be bracketed.

Note **which** metadata file the filters read. The tier 1
`selected-atacama-sample-metadata.tsv` carries only `ph`,
`average-soil-relative-humidity`, `elevation` and `average-soil-temperature` —
there is no `transect-name` column in it, and passing it here fails with
`ValueError: Selection of IDs failed with query: ...`. `transect-name` lives in
the tier 2 `sample-metadata.tsv` (see
[Download the Tutorial Data](../00_getting_started/03_download_data.md)), which
covers all 50 samples of the 13-ASV table. On that table the split is even:
**25 samples in Baquedano and 25 in Yungay**.
```

```bash
# list the metadata columns available for grouping
qiime metadata tabulate \
    --m-input-file data/selected-atacama-sample-metadata.tsv \
    --o-visualization data/metadata-summary.qzv
```

Each instance then gets its own transformation and its own covariance matrix.
Transform and covariance are per-instance operations: the whole point of MGL is
that the $K$ covariance matrices are estimated separately and only the
*precision* matrices are coupled.

```bash
qiime gglasso transform-features \
    --i-table data/atacama-counts-baquedano.qza \
    --i-taxonomy data/classification.qza \
    --m-sample-metadata-file data/selected-atacama-sample-metadata.tsv \
    --p-transformation mclr \
    --p-add-metadata False \
    --p-scale-metadata False \
    --o-transformed-table data/atacama-mclr-baquedano.qza

qiime gglasso transform-features \
    --i-table data/atacama-counts-yungay.qza \
    --i-taxonomy data/classification.qza \
    --m-sample-metadata-file data/selected-atacama-sample-metadata.tsv \
    --p-transformation mclr \
    --p-add-metadata False \
    --p-scale-metadata False \
    --o-transformed-table data/atacama-mclr-yungay.qza

qiime gglasso calculate-covariance \
    --i-table data/atacama-mclr-baquedano.qza \
    --p-method scaled \
    --o-covariance-matrix data/atacama-corr-baquedano.qza

qiime gglasso calculate-covariance \
    --i-table data/atacama-mclr-yungay.qza \
    --p-method scaled \
    --o-covariance-matrix data/atacama-corr-yungay.qza
```

```{tip}
`--i-taxonomy` is a required input that `transform-features` never reads. Pass
any valid `FeatureData[Taxonomy]` artifact; its contents do not affect the
result. See [Troubleshooting](../90_reference/04_troubleshooting.md).
```

## Step 2: Build the bookkeeping array

When the $K$ instances do not contain exactly the same features, the solver
needs a map telling it where feature pair $(i, j)$ lives in each instance. That
map is the *bookkeeping array* $G$, a `(2, L, K)` integer array with one slice
per group of overlapping features. `build-groups` constructs it from the list of
tables:

```bash
qiime gglasso build-groups \
    --i-tables data/atacama-counts-baquedano.qza data/atacama-counts-yungay.qza \
    --p-check-groups True \
    --o-group-array data/atacama-groups-transect.qza \
    --verbose
```

`--p-check-groups True` (the default) runs GGLasso's `check_G` validator on the
result and prints the per-instance dimensions $p_k$, the per-instance sample
sizes $N_k$, and the number of groups found. Set `--p-check-groups False` only if
you have already validated the array once and want the noise gone — the check is
cheap and catches malformed input.

```{note}
The printed dimensions, sample sizes and group count are pending verification
against QIIME 2 2026.7 and are not reproduced here.
```

```{warning}
`build-groups` only returns an array when it detects that the instances differ.
If every table carries the same labels it prints *"All datasets have exactly the
same number of features."* and returns nothing, which QIIME 2 cannot turn into a
`TensorData` artifact — the action fails instead of producing an empty result.

There is also a discrepancy worth knowing about: the difference check is applied
to the **column** labels of each `biom.Table.to_dataframe()`, which are sample
identifiers, while the reported $p_k$ counts rows, which are features. Two
tables produced by splitting one table on a metadata column always have disjoint
sample sets, so the check will report them as differing even when their feature
sets are identical. Inspect the exported array (next section) before relying on
it. Both behaviours are read from the plugin source and are pending confirmation
against a live run.
```

## Step 3: The chaining gap, and the export workaround

`build-groups` emits `group_array` as a **`TensorData` artifact**. `solve-problem`
accepts `group_array` as a **`List[Int]` parameter**. These are different kinds
of thing in the QIIME 2 type system — one is `--o-`/`--i-`, the other is `--p-` —
so the two actions **do not chain**. There is no pipeline, no transformer and no
`--i-group-array` input that connects them.

The workaround is to export the artifact and pass the values on the command
line. `TensorData` is a single-file directory format holding a zarr `ZipStore`
called `tensor.zip`, with the array stored under the key `tensor`:

```bash
qiime tools export \
    --input-path data/atacama-groups-transect.qza \
    --output-path data/exported-groups-transect
```

```python
import numpy as np
import zarr

store = zarr.ZipStore("data/exported-groups-transect/tensor.zip", mode="r")
G = np.array(zarr.open(store=store)["tensor"])
store.close()

print(G.shape)                                  # (2, L, K)
print(" ".join(str(int(v)) for v in G.ravel())) # paste into --p-group-array
```

The printed integers go straight into the solver call. Note that only one
covariance artifact appears here, for the reason spelled out in gap 3 below:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-corr-baquedano.qza \
    --p-n-samples 25 25 \
    --p-non-conforming True \
    --p-group-array 0 1 2 0 1 2 \
    --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 10 \
    --p-lambda2-min 0.001 --p-lambda2-max 0.1 --p-n-lambda2 5 \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-nonconforming.qza \
    --verbose
```

```{warning}
`--p-group-array 0 1 2 0 1 2` is a syntactically valid six-integer list, not the
array your data produces. Substitute the integers your own export prints. The
values shown here exist only to make the command shape unambiguous; they are
**not** a result.

`--p-n-samples 25 25` are the Baquedano and Yungay sample counts of the 13-ASV
tier 1 table (50 samples, split evenly); the full 75-sample Atacama metadata
splits 32/43, which is a different table. Replace them with the counts your own
filtered tables report if you are working from a different subset.
```

### Known gaps

Three separate things are broken or missing on this path. All three were read
out of the plugin source; none has been reproduced against QIIME 2 2026.7.

1. **The artifact does not chain to the parameter.** As above. Documented in
   [Troubleshooting](../90_reference/04_troubleshooting.md); the clean fix is
   either a `--i-group-array` input on `solve-problem` or a QIIME 2 pipeline that
   does the export internally.
2. **The flattened list is not reshaped.** `group_array` arrives as a flat
   `List[Int]` and is handed to the GGLasso solver unchanged — nothing in
   q2-gglasso restores the `(2, L, K)` shape that `build-groups` produced and
   that the solver expects. Exporting the array therefore recovers the *values*
   but not the *structure*.
3. **There is no way to supply $K$ covariance matrices.** `solve-problem` takes a
   single `--i-covariance-matrix`, and `PairwiseFeatureData` stores exactly one
   $p \times p$ table. The multi-instance branches of the solver are only entered
   when the covariance array is three-dimensional, i.e. a stack of $K$ matrices.
   With one two-dimensional artifact the single-graphical-lasso branch is taken
   regardless of `--p-reg`, `--p-lambda2-*`, `--p-non-conforming` or
   `--p-group-array`.

Gap 3 is the blocking one: until a semantic type exists that can carry a stack of
covariance matrices, MGL is reachable from the GGLasso Python API but not from
the QIIME 2 interface. The parameters remain documented here because they are
registered, they will be the interface once the input type lands, and their
meaning does not change.

## Step 4: GGL versus FGL, and the $\lambda_2$ grid

`--p-reg` chooses the coupling and `--p-lambda2-min` / `--p-lambda2-max` /
`--p-n-lambda2` build the grid over its strength. The $\lambda_2$ grid obeys the
same rules as the $\lambda_1$ grid — spacing comes from `--p-path-scale`, an
unset pair of bounds falls back to `np.logspace(-1, -4, 5)`, and there is no
`lambda2_path`. See
[Regularization Paths & Model Selection](05_lambda_paths.md) for the full rule.

Both commands below name only `atacama-corr-baquedano.qza`, because
`--i-covariance-matrix` accepts exactly one artifact — the Yungay matrix has
nowhere to go. That is gap 3 in plain sight. Treat the two commands as the
documented flag shape for MGL rather than as a runnable two-group fit.

```bash
# Group Graphical Lasso: shared support, free edge weights
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-corr-baquedano.qza \
    --p-n-samples 25 25 \
    --p-reg GGL \
    --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 10 \
    --p-lambda2-min 0.001 --p-lambda2-max 0.1 --p-n-lambda2 5 \
    --p-path-scale log \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-ggl.qza \
    --verbose

# Fused Graphical Lasso: shared support AND similar edge weights
qiime gglasso solve-problem \
    --i-covariance-matrix data/atacama-corr-baquedano.qza \
    --p-n-samples 25 25 \
    --p-reg FGL \
    --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 10 \
    --p-lambda2-min 0.001 --p-lambda2-max 0.1 --p-n-lambda2 5 \
    --p-path-scale log \
    --p-gamma 0.01 \
    --o-solution data/atacama-solution-fgl.qza \
    --verbose
```

Because both $\lambda_1$ and $\lambda_2$ hold ten and five values respectively,
these runs perform model selection over a 50-point grid and eBIC picks the pair.
Two consequences follow from the rule in the previous chapter. First, the cost is
the *product* of the grids, so a 10 × 5 MGL search is fifty solves of a problem
that is itself $K$ times larger than the SGL equivalent. Second, if you want a
single MGL fit you must pin both grids — `--p-lambda2-min 0.01 --p-lambda2-max
0.01 --p-n-lambda2 1` alongside the equivalent for $\lambda_1$.

The boundary check described in the previous chapter extends to $\lambda_2$ on
multiple-instance problems: a selection at either end of the $\lambda_2$ grid
warns `lambda is on the edge of the interval, try SMALLER lambda2` (or
`try BIGGER lambda2`) — the GGL/FGL branch's general follow-up is shorter than
the one quoted in the previous chapter, just
`The solution might have not reached global minimum!`, without the
`lambda is on the edge of the interval,` prefix used by the single-instance and
non-conforming branches. A $\lambda_2$ selected at the bottom of the range is the
most informative version of that warning — it is the estimator telling you it
would rather not couple the groups at all, which is an argument for fitting them
independently.

```{warning}
`--p-reg` accepts any string — it is registered as a bare `Str` with no
`Choices()`, so `ggl` or `GGl` is not rejected by the command line. Today an
invalid value is silently ignored, because gap 3 means the MGL branch is never
entered and `reg` is never passed to a solver. Once a multi-instance covariance
input exists, only `GGL` and `FGL` will be valid and they are case-sensitive.
```

```{note}
Selected $(\lambda_1, \lambda_2)$ pairs, per-group edge counts and the
GGL-versus-FGL difference for this dataset are pending verification against
QIIME 2 2026.7.
```

## Step 5: Non-conforming groups

The runs above assume both instances measure the same features. That assumption
breaks as soon as you filter each group independently — a taxon observed in
vegetated soil may be absent from bare soil entirely. `--p-non-conforming True`
switches the solver to the variant that handles unequal feature sets, applying
the group penalty only to pairs of variables that actually exist in more than one
instance. It always uses the `GGL` coupling; `--p-reg` is ignored in this mode,
because a fused penalty on an edge that exists in only one group is undefined.

Split on `vegetation` and then drop, per group, the features that group never
observes. `vegetation` is carried by the tier 1
`atacama-selected-covariates-veg.tsv`, so unlike the transect split this one
needs no tier 2 download (the tier 2 `sample-metadata.tsv` works here too):

```bash
qiime feature-table filter-samples \
    --i-table data/atacama-counts.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --p-where "[vegetation]='yes'" \
    --o-filtered-table data/atacama-counts-veg-yes.qza

qiime feature-table filter-samples \
    --i-table data/atacama-counts.qza \
    --m-metadata-file data/atacama-selected-covariates-veg.tsv \
    --p-where "[vegetation]='no'" \
    --o-filtered-table data/atacama-counts-veg-no.qza

qiime feature-table filter-features \
    --i-table data/atacama-counts-veg-yes.qza \
    --p-min-samples 1 \
    --o-filtered-table data/atacama-counts-veg-yes-observed.qza

qiime feature-table filter-features \
    --i-table data/atacama-counts-veg-no.qza \
    --p-min-samples 1 \
    --o-filtered-table data/atacama-counts-veg-no-observed.qza
```

The two tables now genuinely have different feature sets, which is the case
`build-groups` was written for:

```bash
qiime gglasso build-groups \
    --i-tables data/atacama-counts-veg-yes-observed.qza data/atacama-counts-veg-no-observed.qza \
    --p-check-groups True \
    --o-group-array data/atacama-groups-vegetation.qza \
    --verbose
```

Export it exactly as in Step 3, then solve with `--p-non-conforming True` and the
resulting integers in `--p-group-array`. Every caveat listed under "Known gaps"
above applies unchanged. The transform and covariance steps are the same two
commands as in Step 1, run on `atacama-counts-veg-yes-observed.qza` and
`atacama-counts-veg-no-observed.qza`; only the filenames change. Because of
gap 3 only one of the two resulting matrices can be passed to
`--i-covariance-matrix`.

```{important}
On the tier 1 table the vegetation split really does produce different feature
sets: after `--p-min-samples 1` the vegetated subtable keeps 13 features (33
samples) and the bare subtable keeps 9 (17 samples). And in any case
`build-groups` bases its match / no-match decision on the **sample**
identifiers rather than the feature identifiers, so a table split on a metadata
column always yields an array — the "datasets match" path is not reachable this
way. Inspect the exported array (Step 3) before relying on it.

A larger feature space remains the realistic setting for non-conforming MGL: the
300-ASV table in [Tier 2](../04_highdim_atacama/00_index.md) is where the mode
earns its keep.
```

## When GGL beats FGL

The choice is about what you believe is shared between the groups, not about
which fits better.

**Use GGL when you believe the groups share a wiring diagram but not its
strengths.** Two transects of the same desert plausibly host the same
interactions — the same taxa competing for the same nutrients — at different
intensities, because moisture, pH and depth differ. GGL encodes exactly that: an
edge is switched on or off jointly, and once on, each group estimates its own
weight. This is also the safer default when the groups have very different sample
sizes, because a small group borrows *support* from the large one without being
forced toward its coefficient values.

**Use FGL when the groups are ordered or nearly identical replicates.** FGL's
penalty on pairwise *differences* is the right prior for a time course, a dose
series, or technical replicates, where you expect adjacent conditions to have
almost the same numbers and you want the estimator to shrink small differences
away. Applied to two ecologically distinct habitats, that same prior actively
suppresses the between-group differences you were trying to find.

**Use $K$ separate single graphical lassos when the groups may not share
structure at all.** MGL with an aggressive $\lambda_2$ will manufacture agreement
between groups that have none. If you cannot articulate why the groups should
share edges, fit them independently and compare the results — that comparison is
an honest answer, whereas a jointly-estimated pair of near-identical networks is
partly an artifact of the penalty.

A useful diagnostic, once the CLI path is complete, is to fit at
$\lambda_2 \approx 0$ and at the eBIC-selected $\lambda_2$ and compare. If the
networks barely move, the coupling is not doing anything and the simpler
independent fits are preferable. If they move a great deal, check that the shared
edges are ones you can defend on biological grounds rather than ones the penalty
invented.

---

Next: [Latent-Component PCA](07_pca.md) uses the low-rank part of a latent
solution to place samples in a covariate-aware space, and
[Interpretation](09_interpretation.md) compares all the Tier 1 models
side by side.
