# Reproducibility and Source Versions

The supplement-derived pages retain the audited calculations. No models were
refitted during the website conversion. Figure PDFs are copied without changing
their contents; PNGs are rendered from those PDFs for browser display.

## Audited Environment

| Component | Source |
| --- | --- |
| QIIME 2 | 2026.7.0 |
| c-lasso | `d22da8d5ea54eb499de2fce82fa6c0e1f0c2dcc0` (GitHub, PR60) |
| q2-classo | `ae574af47b7d79761c01bb8ffcea161716387d11` |
| q2-gglasso | `3f586aa1e074727f497e59481ec1565228042cd3` |
| GGLasso | 0.3.1 |

The c-lasso PyPI release is not a substitute for this audited GitHub revision.
The numerical audit was an established-environment replay, not a clean-machine
or cross-platform installation test.

## Reading the Commands

The tutorial and supplement show non-default parameter settings only.
For example, the basic regression call uses the default R3 formulation;
the R1 comparison explicitly switches concomitant estimation off.
The scientific choices are explained beside the commands rather than
repeating default flags. Version pinning is essential: a different
plugin release may have different defaults.

The [registered action defaults](/downloads/qiime-interface-defaults.json)
record the exact interface used for this edition. The command-harmonization
report retains the previous explicit commands and verifies that omitting
default-valued options leaves the resolved parameter settings unchanged.

## Source Correspondence

The repository retains the authoritative supplement in
`analysis/supplement-v2/source/`, its conversion manifest in
`analysis/supplement-v2/conversion-manifest.json`, and the converter in
`analysis/scripts/migrate_supplement_book.py`.

The initial conversion uses Pandoc and Poppler. Once a page has been adapted for
step-by-step use, keep its supplement-derived baseline alongside the edited
version. `editorial-manifest.json` records both hashes and any explicit command
filename corrections. The checks require the source analysis commands, figures,
citations, and labels to remain available. Added preparation, prediction, and
visualization commands are recorded in the reviewed tutorial revision.

To check the current tutorial, run from the repository root:

```bash
python analysis/scripts/check_supplement_migration.py .
jupyter-book build docs --warningiserror
```

The checks retain the 30 source command blocks in the conversion baseline, all 26
original figure PDFs, and all 10 tables. Each figure has a PDF download
and a provenance JSON. The manifest maps supplement sections to web pages.
Historical tutorial files remain in Git and are excluded from the new navigation;
their URLs redirect to the corresponding new subject areas.

Do not run the initial converter over the edited book. To update it from a new
supplement, convert into a clean staging copy, compare with the retained
baselines, then review the tutorial adaptations.

## Optional Python Figures

For the 15-ASV regression and classification examples, download the plotting ZIP
linked in [Work with 15 ASVs](tutorial/s2-01.md). Extract it in your working
directory and follow the collapsed **Reproduce the figures with Python** steps.
The script reads the model and prediction artifacts created by the QIIME commands;
it does not refit any model. It writes result tables, matched-style figures, and
new provenance records in a separate output directory.

The adapter reuses the established plotting functions. The original scripts are
also retained in the reproduction archive. Portable adapters for the covariate,
data heatmap, and network figures are not yet included in that download.

## Input Distribution

The 15-ASV walkthrough now supplies a prepared-input ZIP with hashes. This
revision is not ready for public release until all other required tutorial inputs
have stable download locations and the workflows have been checked with portable
paths. The supplied Bracken artifacts also require confirmation of permission for
public distribution. Local paths in historical provenance identify the audited
sources; they are not public download URLs.

The statistical-background command examples are schematic, as in the supplement.
Case-study commands require the matching prepared inputs and helper scripts;
do not replace them with historical 13-ASV, top-100 mOTU, or old solver outputs.
