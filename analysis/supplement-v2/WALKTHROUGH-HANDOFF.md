# Atacama Top-15 Online Walkthrough

The website now leads readers through prepared inputs, a fixed 39/14 sample
split, transformations, R3 elevation prediction, vegetation classification,
and three graphical-model fits. The full-data networks use all 53 samples.

## What Was Validated

- All 24 QIIME commands extracted from the page ran sequentially in a clean
  working directory in the established environment. This includes the optional
  R1 comparison used only in the model-selection figure.
- Six input/transform comparisons and thirteen model/prediction comparisons
  match the audited QIIME artifacts: 19 passing numerical checks.
- Both portable plotting modes ran on the newly generated artifacts.
- The four PNG figures and rasterized model-selection PDF reproduce the
  existing tutorial images pixel-for-pixel in the tested environment.
- The original drawing/extraction function bodies are preserved. Only path
  setup, extraction configuration, and a Figure S3 input-directory argument were
  adapted; original scripts remain in the reproduction archive.
- The strict documentation build and all 851 local links/images pass.

These are established-environment replay checks, not a clean installation or
cross-platform test. Rendering can vary between fonts/platforms.

## Downloads and Provenance

`docs/downloads/atacama-top15-inputs.zip` contains 11 prepared input files plus
their manifest. It contains no fitted models. The plotting ZIP contains the
artifact reader, existing plotters, helper files, instructions, and a hashed
source map. CI checks the ZIP bytes against the bundled portable source copies.

`editorial-manifest.json` records the adapted pages and their conversion
baselines. It explicitly permits the mCLR output filename correction, matching
the later network commands. It does not permit changes to numerical analysis
parameters. The original source supplement and all 26 source figure PDFs are
unchanged. The ASV table now includes all 15 original feature IDs; the unsupported
raisebox macro had dropped them during initial conversion. The incorrect
geographic Figure S1 sentence is removed only from the web version.

## Remaining Work

- Package and test a complete student-facing QIIME installation recipe.
- Supply and adapt remaining high-dimensional and metagenomics workflows.
- Add portable covariate/data-heatmap and network plotting adapters.
- Complete browser-based visual review. The available browser tool blocked
  file-URL access; no browser-rendered website QA is claimed.
- Resolve public Bracken input distribution before a production release.

Development-branch CI provides a `tutorial-preview` HTML artifact. Production
deployment remains gated on a push to `main`.
