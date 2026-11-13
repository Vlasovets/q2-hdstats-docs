# Optional Atacama Figures

Run in the audited QIIME 2 environment (NumPy, pandas, Matplotlib, zarr 2).
Do not install a different statistical solver to generate these plots.

After running the corresponding QIIME commands in the tutorial:

```bash
python plotting/plot_atacama.py regression --input-dir . --output-dir figures/elevation
python plotting/plot_atacama.py classification --input-dir . --output-dir figures/vegetation
```

Regression requires the optional R1 comparison artifact for its CV panel, plus
the R3 CV/stability models and both training/held-out prediction artifacts.
The script lists missing files before plotting. It copies the required artifacts
to the separate output directory, extracts their stored results, and reuses the
established Atacama plotting functions. It does not fit models or choose tuning
parameters. Each run writes provenance with source-artifact and script hashes.

The archived original plotters remain unchanged in the reproduction bundle.
The copies here have portable path setup; extraction globals are configured by
the wrapper and temporary Zarr archives are cleaned up after reading.
