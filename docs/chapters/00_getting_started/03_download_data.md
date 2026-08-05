# Download the Tutorial Data

This is the canonical place for download URLs. Every later chapter assumes the
files are already on disk under a `data/` directory and refers to them by path
alone. If a command in a later chapter cannot find its input, the fix is here,
not there.

What each dataset contains and why it is in the book is described alongside its
download block below, tier by tier; the Tier 1 study itself, and the
`ASV-1` … `ASV-13` labelling, are covered in
[Atacama Soil Microbiome](02_datasets.md). Otherwise this page is purely
mechanical: fetch, verify, lay out.

## Where the files go

Pick one project directory and stay in it. Every command in the tutorial is
written relative to that directory, with a `data/` prefix:

```bash
mkdir -p ~/q2-hdstats-tutorial/data
cd ~/q2-hdstats-tutorial
```

```{tip}
A few of the high-dimensional chapters show bare filenames
(`atacama-top-300-correlation.qza` rather than `data/atacama-top-300-correlation.qza`).
Either `cd data` before running those, or add the prefix. Nothing else changes.
```

## The Zenodo record

The Atacama artifacts (tiers 1 and 2) are published as a single Zenodo record,
**`q2-hdstats-tutorial-data` v1.0.0**. Publishing them rather than pointing at
the plugin repositories means the exact bytes used to build this book stay
citable and stay put.

```{warning}
**The DOI has not been minted yet.** Wherever a base URL is needed below, the
literal placeholder `ZENODO_DOI_PENDING` stands in for it. The tier 1 and tier 2
commands on this page therefore cannot be run as written yet — they are the
shape of the final commands, with one substitution outstanding. Search the
sources for `ZENODO_DOI_PENDING` to find every place that has to change when the
record goes live. Tier 3 (below) is downloadable **today**, from its upstream
tutorial site.
```

Set the base URL once and reuse it. Every tier 1 and tier 2 block on this page
assumes this variable is exported in your shell:

```bash
# Substitute the real record URL here once the DOI is minted.
export ZENODO_BASE=ZENODO_DOI_PENDING
```

## The manifest

`docs/_data/manifest.tsv` in this repository is the machine-checkable index of
everything the tutorial downloads. It has five tab-separated columns:

| Column | Meaning |
|---|---|
| `filename` | Basename as it lands on disk |
| `tier` | `1`, `2` or `3` — which part of the book needs it |
| `bytes` | Expected size |
| `sha256` | Expected checksum |
| `url` | Where it came from (Zenodo for tiers 1–2, the gut-to-soil tutorial site for tier 3) |

Prose can drift out of date; the manifest cannot, because it is what the
verification steps below read. Treat it as authoritative when the two disagree,
and open an issue.

Define this helper once — every tier below uses it to turn the manifest into a
checklist that `sha256sum -c` understands. It needs a checkout of this
repository; point `MANIFEST` at wherever yours is:

```bash
export MANIFEST=~/q2-hdstats-docs/docs/_data/manifest.tsv
tier_checklist () {
    awk -F'\t' -v tier="$1" 'NR > 1 && $2 == tier { print $4 "  " $1 }' "$MANIFEST"
}
```

```{note}
**Every `bytes` and `sha256` value in the manifest is real** — all 16 rows — and
each was verified against the file it describes. `tier_checklist N | sha256sum -c`
reports `OK` for every tier, so the verification steps below work today.

The one thing still outstanding is the **`url` column for tiers 1 and 2**, which
reads `ZENODO_DOI_PENDING` until the record is minted. So you can verify files you
already have, but the `curl` commands for those two tiers cannot fetch them yet.
Tier 3 downloads from its upstream tutorial site and works end to end.
```

## Tier 1 — Atacama, 13 ASVs

The reference tier. Every action in both plugins gets its one canonical
demonstration here, on a table small enough to print in full.

```bash
cd ~/q2-hdstats-tutorial/data
curl -L -O "${ZENODO_BASE}/atacama-counts.qza"
curl -L -O "${ZENODO_BASE}/classification.qza"
curl -L -O "${ZENODO_BASE}/selected-atacama-sample-metadata.tsv"
curl -L -O "${ZENODO_BASE}/atacama-selected-covariates-veg.tsv"
```

| File | What it is |
|---|---|
| `atacama-counts.qza` | `FeatureTable[Frequency]`, 13 ASVs × 50 samples. Feature IDs are MD5 hashes; the `ASV-1` … `ASV-13` labels used throughout the book are *produced* by `transform-features --p-keep-original-id False`, and the mapping is in [Atacama Soil Microbiome](02_datasets.md) |
| `classification.qza` | `FeatureData[Taxonomy]` for those 13 ASVs |
| `selected-atacama-sample-metadata.tsv` | Sample metadata passed to `qiime gglasso transform-features` |
| `atacama-selected-covariates-veg.tsv` | The five covariates used by q2-classo: `ph`, `elevation`, `average-soil-relative-humidity`, `average-soil-temperature`, `vegetation` |

Verify:

```bash
tier_checklist 1 > ~/q2-hdstats-tutorial/data/SHA256SUMS.tier1
cd ~/q2-hdstats-tutorial/data && sha256sum --ignore-missing -c SHA256SUMS.tier1
```

`--ignore-missing` matters: it lets you verify a partial download against a
checklist that also names files you deliberately skipped, instead of failing on
their absence.

You need `classification.qza` even for the pure network chapters:
`qiime gglasso transform-features` declares `--i-taxonomy` as a required input
although its function body never reads it. See
[Troubleshooting](../90_reference/04_troubleshooting.md) for that and the other
registration warts.

```{note}
**On the sample count.** The tier 1 table is **13 features × 50 samples**, read
directly from the artifact. Earlier drafts of this book said 49 in prose while
every command passed `--p-n-samples 50`; the commands were right and the prose
was wrong. The metadata TSV has 75 rows because it covers the full Atacama
tutorial, of which 50 samples appear in this table.
```

## Tier 2 — Atacama, 300 ASVs

The same study, scaled up: the 300 most abundant ASVs, used to show model
selection and latent-rank choice in a regime where you cannot eyeball the matrix.

```bash
cd ~/q2-hdstats-tutorial/data
curl -L -O "${ZENODO_BASE}/atacama-top-300-table.qza"
curl -L -O "${ZENODO_BASE}/atacama-top-300-clr.qza"
curl -L -O "${ZENODO_BASE}/atacama-top-300-correlation.qza"
curl -L -O "${ZENODO_BASE}/atacama-taxonomy-silva138.qza"
curl -L -O "${ZENODO_BASE}/sample-metadata.tsv"
```

| File | What it is |
|---|---|
| `atacama-top-300-table.qza` | `FeatureTable[Frequency]`, the top-300-ASV counts |
| `atacama-top-300-clr.qza` | The transformed table, shipped so you can skip the transform step |
| `atacama-top-300-correlation.qza` | `PairwiseFeatureData` — the direct input to `qiime gglasso solve-problem` |
| `atacama-taxonomy-silva138.qza` | `FeatureData[Taxonomy]` for the 300 ASVs |
| `sample-metadata.tsv` | Full sample metadata, including `transect-name` (Baquedano, Yungay) and `vegetation` (yes/no) — the two natural grouping variables for the multiple-graphical-lasso and PCA chapters |

Verify:

```bash
tier_checklist 2 > ~/q2-hdstats-tutorial/data/SHA256SUMS.tier2
cd ~/q2-hdstats-tutorial/data && sha256sum --ignore-missing -c SHA256SUMS.tier2
```

The transformed table and the correlation matrix are both included on purpose.
Recomputing them from `atacama-top-300-table.qza` is the first exercise of the
tier 2 chapters; having the reference versions on disk means a mismatch is
visible immediately rather than propagating into the network.

```{warning}
`sample-metadata.tsv` is a generic name and the tier 3 dataset ships a file with
exactly the same name. Keep tier 3 in its own subdirectory, as below, or one will
overwrite the other.
```

## Tier 3 — Gut-to-soil (16S)

Tier 3 is not on the Zenodo record. Its files belong to the gut-to-soil
composting tutorial and are served from that tutorial's own site, so we link them
at the source rather than mirroring them. **These URLs work today.**

```bash
mkdir -p ~/q2-hdstats-tutorial/data/gut-to-soil
cd ~/q2-hdstats-tutorial/data/gut-to-soil

BASE=https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil
curl -L -O "${BASE}/asv-table-ms2.qza"
curl -L -O "${BASE}/taxonomy.qza"
curl -L -O "${BASE}/sample-metadata.tsv"
```

Those three are what the tier 3 chapters use. Two more are available from the
same location if you want to work outside the book:

```bash
cd ~/q2-hdstats-tutorial/data/gut-to-soil
BASE=https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil
curl -L -O "${BASE}/asv-table.qza"
curl -L -O "${BASE}/asv-seqs.qza"
```

| File | What it is |
|---|---|
| `asv-table-ms2.qza` | The tutorial subsample used here: 335 ASVs × 99 samples |
| `taxonomy.qza` | `FeatureData[Taxonomy]` for the ASVs |
| `sample-metadata.tsv` | Full-study sample metadata (1,660 rows); only 99 rows match `asv-table-ms2.qza` |
| `asv-table.qza` | The other tutorial table, 1069 features × 104 samples (not used here) |
| `asv-seqs.qza` | Representative sequences (not used here) |

Verify:

```bash
tier_checklist 3 > ~/q2-hdstats-tutorial/data/gut-to-soil/SHA256SUMS.tier3
cd ~/q2-hdstats-tutorial/data/gut-to-soil && sha256sum --ignore-missing -c SHA256SUMS.tier3
```

Because these are upstream artifacts we do not control, the manifest records
their `sha256` **as observed when the book was built**. A mismatch does not
necessarily mean a corrupt download — it may mean the upstream tutorial was
regenerated. Check the upstream site before assuming the worst.

```{important}
`asv-table-ms2.qza` is a **10% subsample** of a much larger study and it is
shallow: per-sample depth ranges from **3 to 1218 reads, with a median of 261**. A CLR-transformed covariance
estimated from counts that small is noisy, and the resulting network should be
read as a demonstration of the workflow, not as a biological result. The tier 3
chapters say the same thing where it matters.
```

The full study — 1,660 samples, from the gut-to-soil composting work of Meilander
and colleagues {cite}`meilander2025upcycling` — is Zenodo record **15390940**
{cite}`caporaso2025guttosoil`. It is far too large for a tutorial run on a laptop,
but it is the right starting point if you want to take the workflow somewhere
real. Its feature count is pending verification and is not quoted here; see
[Gut-to-Soil: The Dataset](../05_metagenomics/01_gut_to_soil/01_data.md) for how
the scaling arithmetic is written conditionally on it.

## The cocoa appendix has no download

The MOSHPIT cocoa fermentation example (14 shotgun metagenomes, BioProject
PRJNA552479) is an **appendix**, not a tier, precisely because there is nothing
to download: no feature table is published for it. Reproducing one means running
assembly and binning and holding local Kraken2/Kaiju databases. See
[Appendix: Shotgun Metagenomics (MOSHPIT cocoa)](../99_appendix/02_moshpit_cocoa_note.md)
for what that involves before you commit compute to it.

## Why the checksums are not ceremony

A truncated `.qza` is still a valid ZIP prefix. QIIME 2 will often load it
without complaint and the failure will surface much later, inside a solver, as
something that looks like a numerical problem rather than a broken file. Ten
seconds of `sha256sum` saves an afternoon of debugging the wrong layer.

Once the checksums pass, confirm that QIIME 2 agrees the artifacts are what the
tutorial thinks they are:

```bash
cd ~/q2-hdstats-tutorial
qiime tools peek data/atacama-counts.qza
qiime tools peek data/classification.qza
```

```{note}
The expected `qiime tools peek` output — UUID, semantic type and format — is
pending verification against QIIME 2 2026.7. Most of this book has not yet been
re-run against that release, and the UUID is per-download in any case, so no
captured output is shown here.
```

## The resulting tree

After all three tiers:

```text
~/q2-hdstats-tutorial/
└── data/
    ├── atacama-counts.qza
    ├── classification.qza
    ├── selected-atacama-sample-metadata.tsv
    ├── atacama-selected-covariates-veg.tsv
    ├── atacama-top-300-table.qza
    ├── atacama-top-300-clr.qza
    ├── atacama-top-300-correlation.qza
    ├── atacama-taxonomy-silva138.qza
    ├── sample-metadata.tsv
    └── gut-to-soil/
        ├── asv-table-ms2.qza
        ├── taxonomy.qza
        └── sample-metadata.tsv
```

Everything else the tutorial uses is **derived**: each chapter writes its
transformed tables, covariance matrices, solutions and visualizations back into
`data/`, so the directory grows as you work through the book. Only the files
above have to be fetched.

```{note}
The tier 2 chapters also build a design-matrix table and a mean-imputed outcomes
TSV for the q2-classo section. Whether those ship on the Zenodo record or are
produced by the chapter commands is being settled as the record is assembled; the
manifest is the place that will say.
```

## Next

With the data in place, continue to
[Prerequisites & Installation](../01_installation/01_prerequisites.md), then
[Verifying Your Installation](../01_installation/04_verify.md) — the last step of
which reads one of the artifacts you just downloaded. If you want to see which
chapter demonstrates which command before you start, the
[Command Coverage Matrix](../90_reference/01_command_coverage.md) is the map.
