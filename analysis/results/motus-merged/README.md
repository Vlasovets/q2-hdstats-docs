# Merged mOTUs profiles

Everything here is generated. Nothing is edited by hand.

    python analysis/scripts/aggregate_motus_profiles.py

The script merges the 36 per-sample artifacts in `analysis/results/motus-profiles/`,
which are the unmodified output of `qiime motus profile` (stage
`analysis/slurm/35_motus_profile_all.sh`). It reads the BIOM payload out of each `.qza`
with h5py rather than through `qiime2.Artifact.load`, because building a QIIME 2
PluginManager pulls in `q2-composition` and then rpy2 and R, which fails on this
cluster's login node. A `.qza` is a zip with a BIOM inside; the framework is not needed
to read one.

## Files

| file | shape | what it is |
|---|---|---|
| `table-raw-with-unassigned.tsv` | 36 x 379 | the profiler's output as merged, before any decision is taken about it |
| `table.tsv` | 36 x 378 | the same table with the `unassigned` row dropped |
| `taxonomy-378.tsv` | 378 rows | lineage for every merged feature |
| `genus-prev10.tsv` | | genus-level counts at a prevalence-10 floor, kept from an earlier comparison |

Rows are sequencing runs (`run_accession`), columns are features. Counts are marker-gene
inserts, not reads.

## Why there are two tables

`table.tsv` has already had one editorial decision applied to it, and a reader given only
that table cannot see what was removed. mOTUs emits an explicit `unassigned` row for reads
it cannot place. That row is the complement of the profiled fraction, so including it in a
log-ratio would define the composition over a mixture of taxa and not-taxa. It accounts
for **1.8%** of all counts here.

`table-raw-with-unassigned.tsv` exists so the filtering chain can be shown from its true
start rather than from the first step that had already been taken.

## The filtering chain

| step | result |
|---|---|
| merged profile | 36 x **379** |
| drop `unassigned` (1.8% of counts) | 378 features, **93.6%** zeros |
| samples with at least 100 assigned counts | 36 -> **34** samples |
| features present in at least 2 remaining samples | 378 -> **156** features |
| the 100 most abundant of those | **100** features, **85.4%** zeros, median depth **814** |

The last three steps are applied in `analysis/slurm/44_motus_tutorial_run.sh`, which
writes the analysed artifacts to `docs/_data/motus/`. Features are ranked by abundance
rather than prevalence on purpose: on a table this sparse, two features absent from the
same samples correlate strongly through shared zeros alone, and under prevalence
filtering the extended BIC selects the empty graph at every threshold tried.

## On the taxonomy

mOTUs emits lineage **strings**, not a phylogeny, and there is no tree file anywhere in
this analysis. `taxonomy-378.tsv` is `Feature ID` and `Taxon`, in the
`FeatureData[Taxonomy]` format QIIME 2 expects.

Each lineage terminates in `m__<mOTU_id>` rather than a species name. That substitution
is what guarantees unique leaves when two mOTUs share a species name, which
`qiime classo add-taxa` requires: it builds the aggregation matrix from these strings at
run time, and duplicate leaves would silently merge two features into one node.

The published `docs/_data/motus/motus-top100-taxonomy.qza` covers only the 100 analysed
features. Use `taxonomy-378.tsv` if you want to re-filter the table differently.
