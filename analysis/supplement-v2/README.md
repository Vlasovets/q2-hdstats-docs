# Supplement-to-Tutorial V2

The retained LaTeX supplement is the scientific source of truth. `source/`
contains that snapshot; `conversion-manifest.json` identifies generated pages,
exact command blocks and source figure hashes. Branding hashes are recorded in
`branding-manifest.json`.

`reproduction/` preserves the original plotting and workflow scripts without
silently rewriting historical paths. `script-path-map.json` maps each original
file to its bundled copy and hash. These snapshots include historical diagnostic
scripts: inclusion here is preservation, not endorsement as a final workflow.
The figure registry identifies the active generator for each included figure.

Before public release, package final workflow inputs and adapt runner paths in
a separately reviewed change. The Bracken data are not bundled or uploaded by
this migration. Do not treat older Zenodo placeholder URLs as usable downloads.

Validate the edited tutorial from the repository root with:

```bash
python analysis/scripts/check_supplement_migration.py .
jupyter-book build docs --warningiserror
```

Building requires the repository's Jupyter Book 1.x requirements. Statistical
model fitting is not part of these documentation commands. The original
converter requires Pandoc and Poppler; run it only in a clean staging copy when
updating from a new supplement. It refuses to overwrite reviewed tutorial edits.

`baselines/` preserves the initial converted Atacama pages. The
`editorial-manifest.json` records their reviewed online adaptations, including
the mCLR filename correction. Source analysis commands remain checked against
the baselines. Prepared top-15 inputs and portable prediction-plot adapters are
available through `docs/downloads/`; other case-study downloads are still pending.

## Source issue awaiting editorial review

The source Atacama introduction refers in prose to a geographic Figure S1,
while the active Figure S1 is the covariate overview. The reviewed online
adaptation removes that incorrect sentence; the original supplement and baseline
remain unchanged. The web adaptation also restores the 15 feature IDs lost by
the initial conversion's unsupported `raisebox` macro.
