# Appendix: Shotgun Metagenomics (MOSHPIT cocoa)

This is a compatibility note, not a tutorial. It records what a shotgun
metagenomics dataset looks like when it reaches q2-gglasso or q2-classo, using the
MOSHPIT cocoa fermentation example as the concrete case — and it states, up front,
that **no feature table ships with this book for it and none of the commands below
were run**.

## The dataset

Cocoa bean fermentation, sampled as **14 shotgun metagenomes**: two seed types
(*Forasteiro* and *Hybrid*) at seven timepoints (0, 24, 48, 72, 96, 120 and
144 hours). Sequencing is Illumina HiSeq X Ten whole-genome shotgun, deposited
under BioProject **PRJNA552479**.

Two features of that design are worth naming immediately, because they determine
everything below: $2 \times 7 = 14$ means **exactly one sample per condition**,
with no within-cell replication, and 14 samples is a small number in an absolute
sense regardless of how the conditions are arranged.

## Where a feature table comes from

MOSHPIT is the QIIME 2 distribution for shotgun metagenomics, and there are three
distinct routes through it that end in a `FeatureTable[Frequency]`:

- **Read-based taxonomic profiling with Kraken2, refined by Bracken.** Kraken2
  assigns reads to taxa against a k-mer database; Bracken redistributes them to a
  chosen rank and produces per-sample abundance estimates. The result is a table
  of taxa × samples.
- **Read-based profiling with Kaiju**, which classifies against a protein
  database instead. Same shape of output, different assignment behaviour.
- **Assembly and binning.** Assemble the reads, bin the contigs into MAGs,
  dereplicate the bins across samples, then estimate per-sample abundance of each
  MAG — typically as TPM. The features are genome bins rather than named taxa.

The action names for those steps belong to the MOSHPIT documentation and are
deliberately not reproduced here; they change between releases, and a wrong
command in this book would be worse than no command.

## What happens when such a table arrives

Nothing special, which is the entire point of this appendix. Both plugins consume
a `FeatureTable` and neither inspects its provenance:

- `qiime gglasso transform-features` accepts
  `FeatureTable[Composition | Frequency | Design]`. A Bracken, Kaiju or MAG table
  registered as `FeatureTable[Frequency]` goes straight in, and the rest of the
  chain — `calculate-covariance`, `solve-problem`, `summarize` — is unchanged from
  [Data Preparation](../02_lowdim_gglasso/01_data_preparation.md).
- `qiime classo transform-features` likewise takes the table as `--i-features` and
  the workflow proceeds as in
  [Data Preparation](../03_lowdim_classo/02_data_preparation.md).

Three practical caveats attach to shotgun-derived tables specifically:

1. **TPM is not a count.** MAG-abundance tables are already normalized for genome
   length and depth. The compositional argument for a log-ratio transform still
   applies — TPM values are relative by construction — but `--p-pseudo-count` and
   `--p-coef`, which are written as count-like floors, no longer have an obvious
   scale. The two are not equally adjustable. `--p-coef` (q2-classo) is registered
   as a `Float`, so it can be set below the smallest non-zero TPM value in your
   table. `--p-pseudo-count` (q2-gglasso) is registered as an `Int`, so it has a
   hard floor of 1 and no fractional value is expressible — which is exactly the
   value a TPM table usually needs. For such tables prefer
   `qiime gglasso transform-features --p-transformation mclr`, which
   log-transforms only the positive entries and needs no floor at all. Note that
   `mclr` is a q2-gglasso option; `qiime classo transform-features` implements
   `clr` only, so on the classo side the answer is to tune `--p-coef`.
2. **`--i-taxonomy` is still required by `qiime gglasso transform-features` and
   still unread.** For a MAG table there may be no natural taxonomy to pass; any
   valid `FeatureData[Taxonomy]` will satisfy the input. See
   [Troubleshooting](../90_reference/04_troubleshooting.md).
3. **`qiime classo add-taxa` does need a real taxonomy**, because trac builds its
   aggregation matrix from the hierarchy. Bracken output at a fixed rank is flat
   and gives the tree nothing to aggregate over; a MAG table needs taxonomic
   assignments attached first.

## Why nothing ships, and why this is not an analysis

**No feature table is published for this study.** Producing one means running
assembly and binning, or holding a local Kraken2 or Kaiju database — large
reference resources that are impractical to bundle with a tutorial and slow to
build from scratch. That cost is the reason this is an appendix rather than a
tier: see [Download the Tutorial Data](../00_getting_started/03_download_data.md),
which lists no files for it.

**Fourteen samples cannot support covariance estimation.** The empirical
covariance matrix of $n$ samples has rank at most $n - 1$, so with $n = 14$ every
estimate lives in at most 13 dimensions while a Bracken table at species level has
features in the thousands. A graphical lasso would still return a network, because
that is what a penalized estimator does; the network would be almost entirely the
penalty. The same applies to the log-contrast models — with one sample per
condition, cross-validation folds contain one or two samples and the resulting
error estimate is not usable.

So: **this appendix demonstrates that shotgun tables are compatible. It does not
demonstrate that they can be analyzed at this sample size.** If you have a shotgun
dataset with tens to hundreds of samples, the workflows in
[Tier 2](../04_highdim_atacama/00_index.md) transfer directly — 16S rather than
shotgun, but identical in every command.

```{note}
No MOSHPIT step, no table, and no model has been run for this appendix against
QIIME 2 2026.7 or any other release. There are no outputs, runtimes or database
sizes quoted anywhere on this page.
```
