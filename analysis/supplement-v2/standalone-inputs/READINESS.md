# Standalone tutorial: verified existing assets

The Supplement's scientific workflow already exists. This inventory and input
packaging do not refit models, change filtering, alter splits, or replace figures.

## Verified in this stage

- Existing top-15 bundle: retained, previously replayed sequentially through QIIME.
- Top-300 Atacama bundle: all 11 external files named in the page, original source
  count/metadata tables, split, normalization/ranking records, label map, and two
  portable preparation helpers. The helpers differ only in path assignments;
  the remaining AST nodes match their original generators.
- mOTU bundle: all six distinct external filenames required across regression
  and networks, plus the complete profiled table and three upstream source files.
  The same taxonomy is deliberately provided under both webpage aliases.
- ZIP member hashes verified against source copies; no fitted models bundled.
- Portable Atacama preparation rerun in a temporary clean directory: ranking,
  training-support whitelist, covariate scaling, and normalization summary match
  audited data to numerical tolerance 1e-12. Constraint has shape 1x285; all
  prediction weights have shape 285. Full network table remains 300x53.
- mOTU inputs agree exactly with the full 60x34 table, use a 26/8 split, have
  training support for all 60 features, and share no infant between train/test.

See input-inventory.json, bundle-manifest.json, and input-validation.json.

## Concrete remaining boundaries

- The existing scientific environment and source pins are recorded. Its full
  replay was tested previously; installation on a new machine was not. The
  repository's historical Linux environment lock references an older q2-classo
  build and must not be advertised as the audited PR60 installation recipe.
- The prepared mOTU aliases existed only in the prose before this packaging.
  They are byte-identical copies of the existing final source artifacts.
- Case-study snippets are not yet complete plotting runners: the high-dimensional
  joint regression still needs held-out covariate construction/prediction steps,
  TRAC needs its held-out aggregation/prediction steps displayed, and several
  figures need portable entry points around existing archived generators.
  Original scripts and source artifacts exist; no new statistical method is needed.
- Cocoa source artifacts are available privately. They are intentionally excluded
  from public bundles pending confirmation of redistribution permission.

These are interface, packaging, and installation-validation tasks, not missing
scientific results. Do not replace the existing audit with a fresh model design.
