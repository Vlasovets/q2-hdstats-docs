# Network Inference

This page estimates one network over all 99 samples, and then sets up the
comparison that this dataset is uniquely suited to: **gut-derived versus
soil-derived communities, estimated jointly**. The first half is the familiar
`transform-features` → `calculate-covariance` → `solve-problem` → `summarize`
chain at a new size. The second half is where the multiple graphical lasso
becomes a scientific question rather than a syntax demonstration — and also where
the current QIIME 2 interface runs out of road, which the page says plainly.

Everything assumes the layout from [The Dataset](01_data.md), with files under
`data/gut-to-soil/`.

## Step 1: Transform

```bash
qiime gglasso transform-features \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --i-taxonomy data/gut-to-soil/taxonomy.qza \
    --p-transformation mclr \
    --p-no-add-metadata \
    --p-no-scale-metadata \
    --p-keep-original-id \
    --o-transformed-table data/gut-to-soil/gts-mclr.qza
```

Three of those choices matter here more than they did on the Atacama tables.

**`--p-transformation mclr`.** The `clr` branch does not add a pseudo-count to
every cell; it *replaces zeros* with `--p-pseudo-count` and then rescales each
sample back to its original total. On a table where most cells are zero, that
puts a large fraction of the features on exactly the same floor value within a
sample, and a block of tied values moves together with that sample's geometric
mean. The resulting apparent co-variation among rare features is an artifact of
the imputation, not of the community. `mclr` {cite}`yoon2019microbial` avoids it
by log-transforming only the positive entries and leaving the zeros as a single
shifted baseline. Use `--p-transformation clr --p-pseudo-count 1` if you want the
comparison, but expect the two networks to differ, and report which one you used.

**`--p-keep-original-id`.** Without it, features are relabelled to sequential
`ASV-n` names ordered by total abundance. At $p = 335$ you cannot recover a node's
identity by eye, and you need the real IDs to join against `taxonomy.qza` later.

**`--p-no-add-metadata`.** The covariates stay out of the table. Appending them is
the adaptive workflow of
[Adaptive Graphical Lasso](../../02_lowdim_gglasso/04_adaptive_glasso.md); here we
want a network of taxa only.

```{tip}
`--i-taxonomy` is required and unread. Passing the correct file costs nothing and
means you have it to hand for interpretation. See
[Troubleshooting](../../90_reference/04_troubleshooting.md).
```

```{important}
`mclr` applies a global shift derived from the **minimum over the whole matrix**
before it fills the zeros back in. That makes the transform depend on the set of
samples in the table: transforming a subset is not the same as transforming
everything and then taking the subset. It matters in Step 6, where the table gets
split into groups, and it is the reason that section transforms each group
separately and says so rather than reusing `gts-mclr.qza`.
```

## Step 2: Covariance

```bash
qiime gglasso calculate-covariance \
    --i-table data/gut-to-soil/gts-mclr.qza \
    --p-method scaled \
    --p-bias \
    --o-covariance-matrix data/gut-to-soil/gts-correlation.qza
```

`--p-method scaled` divides out the diagonal, giving a correlation matrix. That is
what makes a single scalar $\lambda_1$ meaningful across 335 heterogeneous
features: on the unscaled covariance, a high-variance feature would effectively be
penalized less than a low-variance one, and the same $\lambda_1$ would mean
different things in different rows. `--p-bias` (the default) uses the $1/N$
normalization that the solver's Gaussian objective is written against.

## Step 3: Solve, selecting $\lambda_1$ with eBIC

At $p = 335$ there is no eyeballing the penalty. Sweep it and let the extended BIC
{cite}`foygel2010extended` choose:

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/gut-to-soil/gts-correlation.qza \
    --p-n-samples 99 \
    --p-no-latent \
    --p-lambda1-min 0.1 --p-lambda1-max 1.0 --p-n-lambda1 19 \
    --p-path-scale linear \
    --p-gamma 0.3 \
    --o-solution data/gut-to-soil/gts-sgl-path.qza \
    --verbose
```

**Explanation:**

- `--p-n-samples 99`: the number of samples the covariance was computed from.
  eBIC uses it, so it has to be right — if you filtered samples, change it.
- `--p-lambda1-min 0.1 --p-lambda1-max 1.0 --p-n-lambda1 19` with
  `--p-path-scale linear`: the grid `0.10, 0.15, …, 1.00`. A linear path spends
  its points where the edge count actually changes; the default logarithmic
  spacing crowds them at the sparse end.
- `--p-gamma 0.3`: the eBIC's extra penalty on edges. This is the same value
  argued for in [Selecting lambda](../../04_highdim_atacama/02_model_selection.md)
  and it is **not** a constant of nature — $\gamma$ and $\lambda_1$ are chosen
  jointly, so changing $\gamma$ changes which $\lambda_1$ wins. Re-run with
  `--p-gamma 0.5` (the conventional choice) and with the plugin default `0.01`
  and see how far the selection moves. If it moves a lot, say so when you report
  the network.
- `--p-no-latent`: sparse only, for now.

Because `--p-n-lambda1` is greater than one, model selection runs and the
solution carries a `modelselect_stats` group. With a single value it would not,
and `summarize` would silently omit the statistics tab.

```{tip}
`--p-lambda1-path` takes an explicit grid if you would rather name the points
than describe them. There is no `lambda2_path`.
```

```{note}
The selected $\lambda_1$, the eBIC curve and the resulting edge count for this
table are **pending verification against QIIME 2 2026.7**. No numbers are quoted
here, and any you find elsewhere for this dataset have not been produced by the
recompute.
```

## Step 4: Summarize

```bash
qiime gglasso summarize \
    --i-solution data/gut-to-soil/gts-sgl-path.qza \
    --p-width 1800 \
    --p-height 1800 \
    --p-label-size 4pt \
    --o-visualization data/gut-to-soil/gts-sgl-summary.qzv
```

`--p-label-size 4pt` is small on purpose: 335 node labels at the tier 1 setting of
`25pt` overlap into a solid block. `--p-width` and `--p-height` trade file size
for legibility. `--p-n-cov` is left unset because no metadata columns were
appended to the table in Step 1; set it only if you used `--p-add-metadata`.

Open the result at [QIIME 2 View](https://view.qiime2.org/).

## Step 5 (optional): a latent block for the stage gradient

There is a reason to expect this table to have a strong low-rank structure. The
composting process is a single dominant gradient that moves nearly every taxon at
once, and a driver like that shows up in the correlation matrix as dense,
diffuse, non-sparse structure. A sparse-only estimator has to spend edges
representing it. The sparse-plus-low-rank formulation gives it somewhere else to
go: the low-rank block absorbs the shared driver, and the sparse block is left
describing conditional dependence *given* the stage.

```bash
qiime gglasso solve-problem \
    --i-covariance-matrix data/gut-to-soil/gts-correlation.qza \
    --p-n-samples 99 \
    --p-latent \
    --p-lambda1-min 0.4 --p-lambda1-max 0.4 --p-n-lambda1 1 \
    --p-mu1-min 5 --p-mu1-max 20 --p-n-mu1 6 \
    --p-path-scale linear \
    --p-gamma 0.3 \
    --o-solution data/gut-to-soil/gts-slr-mu-scout.qza \
    --verbose
```

That is a **scouting** run: fix $\lambda_1$ at whatever Step 3 selected — `0.4`
above is a stand-in so the command is unambiguous, not a selected value — and
sweep $\mu_1$ to see what ranks are reachable. Larger $\mu_1$ gives a smaller
rank.

```{important}
`--p-rank` is registered and **always raises** — `ValueError` if `--p-latent` is
not set, `NotImplementedError` otherwise. No released
GGLasso can fix the rank of the low-rank block; it exposes only the continuous
$\mu_1$ penalty and reports the achieved rank as an output. Size the block through
`--p-mu1-min` / `--p-mu1-max` and read the rank back out of the solution. The
procedure is worked through in
[Choosing the Latent Rank](../../04_highdim_atacama/03_slr_ranks.md).
```

```{important}
For a latent problem, `lambda1`, `lambda2` **and** `mu1` must all collapse to a
single value before the solver treats the run as a single fit. If you pin
$\lambda_1$ and forget $\mu_1$, the unset $\mu_1$ falls back to a default
ten-point grid and you have silently asked for ten solves; leave $\lambda_1$
unpinned as well and the default fifteen-point $\lambda_1$ grid makes it 150. At
$p = 335$ that is the difference between a coffee and an afternoon.
```

```{note}
Achieved ranks, runtimes and the sparse edge count at each $\mu_1$ are pending
verification against QIIME 2 2026.7.
```

Note that `qiime gglasso pca` needs a **latent** solution — it reads
`solution/lowrank_`, which a sparse-only solution does not have — and it also
needs `--m-sample-metadata-file` despite that being optional in the signature. If
you want the sample projection, this is the solution to give it.

## Step 6: Gut versus soil, jointly

Now the part this dataset is here for.

### Splitting the design

The composting series has a natural two-group structure: samples from the
gut-derived end of the process and samples from the mature compost or soil end.
Read the metadata first and decide which column encodes it and which two levels
you want — [The Dataset](01_data.md) shows the `qiime metadata tabulate` command.

```bash
# Edit these three lines to match your own metadata before running anything below.
GROUP_COLUMN=sample-type
GROUP_A=gut
GROUP_B=soil
```

```{note}
The three values above are placeholders chosen to make the command shape
unambiguous. **They are not a claim about the gut-to-soil metadata.** No column
named `sample-type` was checked for; the one categorical column this book did
read off the served `sample-metadata.tsv` is `SampleType`, and it has **15
levels**, not two (see [Log-Contrast Regression](03_regression.md)), so a
gut-versus-soil split has to be chosen deliberately rather than assumed. Read the
real column and levels off `gts-metadata-summary.qzv` and substitute them. If your
column name contains a space or a hyphen it must stay quoted and bracketed inside
`--p-where`, as written below.
```

```bash
qiime feature-table filter-samples \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --m-metadata-file data/gut-to-soil/sample-metadata.tsv \
    --p-where "[${GROUP_COLUMN}]='${GROUP_A}'" \
    --o-filtered-table data/gut-to-soil/gts-table-a.qza

qiime feature-table filter-samples \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --m-metadata-file data/gut-to-soil/sample-metadata.tsv \
    --p-where "[${GROUP_COLUMN}]='${GROUP_B}'" \
    --o-filtered-table data/gut-to-soil/gts-table-b.qza
```

Each group then drops the features it never observed, and gets its own transform
and its own covariance. The `--p-min-samples 1` filter is what makes the two
instances genuinely *non-conforming* — a taxon absent from the soil-derived end
should be an unmeasured variable there, not a measured zero. Estimating the two
covariances separately is likewise not bureaucracy: the whole premise of the
multiple graphical lasso is that the $K$ covariance matrices are estimated
independently and only the *precision* matrices are coupled.

```bash
for g in a b; do
    qiime feature-table filter-features \
        --i-table "data/gut-to-soil/gts-table-${g}.qza" \
        --p-min-samples 1 \
        --o-filtered-table "data/gut-to-soil/gts-table-${g}-observed.qza"

    qiime gglasso transform-features \
        --i-table "data/gut-to-soil/gts-table-${g}-observed.qza" \
        --i-taxonomy data/gut-to-soil/taxonomy.qza \
        --p-transformation mclr \
        --p-no-add-metadata \
        --p-no-scale-metadata \
        --p-keep-original-id \
        --o-transformed-table "data/gut-to-soil/gts-mclr-${g}.qza"

    qiime gglasso calculate-covariance \
        --i-table "data/gut-to-soil/gts-mclr-${g}.qza" \
        --p-method scaled \
        --p-bias \
        --o-covariance-matrix "data/gut-to-soil/gts-correlation-${g}.qza"
done
```

### The bookkeeping array

If the two groups no longer measure the same features — which they will not, once
you drop features a group never observed — the solver needs a map telling it
where feature pair $(i,j)$ lives in each instance. `build-groups` constructs it:

```bash
qiime gglasso build-groups \
    --i-tables data/gut-to-soil/gts-mclr-a.qza data/gut-to-soil/gts-mclr-b.qza \
    --p-check-groups True \
    --o-group-array data/gut-to-soil/gts-groups.qza \
    --verbose
```

```{note}
The printed per-instance dimensions, sample sizes and group count are pending
verification against QIIME 2 2026.7.

`build-groups` returns an array only when it detects that the instances differ,
and the difference check it applies looks at *sample* labels rather than feature
labels — so two tables produced by splitting one table always register as
differing, whether or not their feature sets match. Inspect the exported array
before relying on it. Both behaviours are described in
[Multiple Graphical Lasso](../../02_lowdim_gglasso/06_multiple_graphical_lasso.md).
```

### Where this stops working

`build-groups` emits `group_array` as a **`TensorData` artifact**;
`solve-problem` accepts `group_array` as a **`List[Int]` parameter**. They do not
chain. The export-and-paste workaround, the fact that the flattened list is never
reshaped back to `(2, L, K)`, and the blocking problem — that `solve-problem`
takes a single `--i-covariance-matrix` and therefore cannot be given $K$
matrices at all — are documented once, in full, in
[Multiple Graphical Lasso](../../02_lowdim_gglasso/06_multiple_graphical_lasso.md).
Read that page before spending time here. It is not repeated because nothing
about it is different at this scale.

```{important}
The consequence for this dataset: **the joint gut-versus-soil fit cannot be run
from the QIIME 2 CLI today.** The commands above are real and produce real
artifacts; the step that would consume them does not exist yet. Two honest ways
forward are to fit the two groups as independent single graphical lassos and
compare them by hand, or to drop to the GGLasso Python API, where the multi-
instance solvers are reachable. The section below is about how to *read* such a
comparison, and it applies either way.
```

## What a shared edge means, and what a differential edge means

This is the part worth thinking about before you run anything, because the
statistics only pay off if the biological question is posed correctly.

The two ends of a composting process are not two samples of one community. A
gut-derived community is anaerobic and host-associated; a mature compost or soil
community is aerobic, has been through a thermophilic phase, and is dominated by
different lineages. Temperature, oxygen, pH, moisture and carbon-to-nitrogen ratio
all change together. Whether that holds in *your* subsample is something to check
against `taxonomy.qza` — but it is the design assumption that makes the joint
estimate interesting.

**A shared edge is one supported in both groups.** Its interpretation is stronger
than it would be within a single group, and the reason is that the confounding is
not shared. Within one group, two taxa can appear conditionally dependent because
both respond to the same unmeasured environmental variable. Across two groups
whose environments differ in almost every respect, that explanation gets harder:
the driver that would have to be doing the work is not the same driver. What is
left is more plausibly an association that travels with the organisms —
cross-feeding, syntrophy, co-aggregation, physical co-transport through the
process — or a pair of taxa that simply survive and disperse together.

Three uninteresting explanations survive that argument and have to be excluded by
hand:

- **The same organism twice.** Two ASVs from different 16S copies of one genome
  are perfectly coupled everywhere. Check whether the taxonomy assigns them to
  the same species.
- **A compositional shadow.** If one taxon dominates a sample, the log-ratio
  values of everything else move against it. CLR-based methods reduce but do not
  eliminate this {cite}`gloor2017microbiome`.
- **A shared technical artifact** — a contaminant pair from the extraction kit,
  or two features that are really one chimera and its parent — which is present
  in both groups precisely because it came from the lab and not from the pile.

**A differential edge is one supported in only one group.** Ecologically this is
the interesting outcome: a coupling that requires a particular regime. An
anaerobic syntrophy that exists in the gut-derived material has no reason to
persist once the pile is aerobic and hot; a thermophile pairing has no reason to
exist before it. That is a hypothesis with a mechanism attached, and it is
testable.

Statistically it is also the **least trustworthy thing the model produces**, for
three compounding reasons, and it is worth being blunt about them:

1. **It is a difference of two noisy estimates.** Splitting 99 samples gives each
   group substantially fewer, and covariance error grows fast as $n$ falls.
2. **The penalty is designed to suppress exactly this.** Both GGL and FGL exist
   to *share* information across groups. An edge that survives as
   group-specific has beaten a prior that was pushing against it — which is
   reassuring — but the same machinery means the size of the difference is
   shrunk, so the estimated contrast is biased toward zero and cannot be read as
   an effect size.
3. **A group-specific edge may be a group-specific *detection*.** If a taxon is
   rare in one group, its edges there are unestimable rather than absent. The
   non-conforming variant (`--p-non-conforming True`) is the principled response:
   it applies the group penalty only to variable pairs that exist in more than
   one instance, instead of pretending an unmeasured pair is a measured zero.

Before believing a differential edge, check that it is stable — that it survives
a change of $\lambda_1$, a change of $\lambda_2$, and resampling of the samples
within each group. A differential edge that appears at one grid point and nowhere
else is a grid point, not a finding.

### GGL or FGL for this design

With two ecologically distinct endpoints, **GGL** is the right coupling. Its
$\lambda_2$ is a group-lasso penalty acting on each edge *across* groups: an edge
tends to be present in all groups or absent in all groups, while the surviving
weights are free to differ. That is exactly the belief "the same organisms
interact at both ends of the process, but not with the same intensity". It is
also the safer choice when the two groups have different sample sizes, because
the smaller group borrows *support* from the larger without being dragged toward
its coefficient values.

**FGL** penalizes the *difference* between corresponding edge weights, which
pushes them to be numerically equal. Applied to gut versus soil, that prior
actively suppresses the between-group differences you were looking for. FGL earns
its place on this dataset in a different configuration: the composting series is a
**time course**, so if you build $K$ instances from ordered timepoints rather than
two endpoints, "adjacent timepoints should have nearly the same edge weights" is
precisely FGL's assumption and precisely what you want. Same data, different
grouping, different penalty.

And if you cannot articulate why the two groups should share edges at all, fit
them separately. A jointly estimated pair of near-identical networks is partly an
artifact of $\lambda_2$, and reporting it as agreement is reporting the penalty.

```{note}
One caveat specific to log-ratio data and joint estimation. `mclr` shifts the
whole matrix by a constant derived from its global minimum before restoring
zeros, so a per-group transform anchors each group differently. Edge $(i,j)$ in
group A and edge $(i,j)$ in group B are therefore not guaranteed to be on exactly
the same scale, and the group penalty assumes they are. With `clr` the
per-sample geometric-mean reference is unaffected by which samples are in the
table, but it *is* affected by which **features** are — so dropping
group-specific features, as non-conforming MGL requires, re-anchors each group
anyway. There is no way to avoid this entirely; the practical response is to note
it as a limitation rather than to pretend the two precision matrices are
measured in identical units.
```

---

Next: [Log-Contrast Regression](03_regression.md) predicts a sample-level outcome
from the same table, which is the other way of asking what the community is
tracking.
