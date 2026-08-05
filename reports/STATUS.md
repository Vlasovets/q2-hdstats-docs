# Recompute status — 2026-08-05T10:36:59+02:00

Written by `slurm/19_final_verify.sh` (job 39111921) on `cpusrv49.scidom.de`.
Everything below was re-run just now; nothing is carried over from an earlier report.

## Environment

```
qiime2     2026.7.0
numpy      2.4.2
pandas     2.3.3
scipy      1.17.1
numba      0.66.0
gglasso    0.3.0
bokeh      3.9.2
zarr       2.18.7
c-lasso    1.0.11
```

## Test suites

- q2-gglasso: `44 passed, 53 warnings in 168.27s (0:02:48)`
- q2-classo:  `22 passed in 4.48s`

## Documentation

- build with `--warningiserror`: **FAIL** — see below
```
EOFError: EOF when reading a line
```

## Generated tables (prose cannot drift from these)

- `atacama-lambda-path.tsv` — 15 rows
- `atacama-mu-rank-map.tsv` — 3 rows
- `atacama-classo-cv.tsv` — 15 rows
- `manifest.tsv` — 14 rows

## Recompute artifacts

| stage | output | state |
|---|---|---|
| 1 lambda path | `atacama-top-300-sgl-linear-path.qza` | present |
| 3 mu=15 rank2 | `atacama-top-300-slr-lambda0.8-mu15.qza` | present |
| 5 comparison | `atacama-top-300-rank0-rank2-comparison-summary.tsv` | present |
| 6 classo CV | regress fits | 15/15 |

## Gates

- **Gate C1 PASS** — the CLI reproduces the reference selection (lambda 0.8, 216 edges).
- Orientation probe: see `reports/ORIENTATION_FINDING.md`.

## Still needs you

- `reports/DECISIONS_NEEDED.md` — open decisions, evidence attached.
- Nothing is committed. `git status` in each of the three repos shows the change set.

## Job log index

```
    39075333_9       q2-classo-cv     FAILED   00:05:47 
   39075333_10       q2-classo-cv     FAILED   00:05:47 
   39075333_11       q2-classo-cv     FAILED   00:05:47 
   39075333_12       q2-classo-cv     FAILED   00:05:47 
   39075333_13       q2-classo-cv     FAILED   00:05:47 
   39075333_14       q2-classo-cv     FAILED   00:05:47 
   39075333_15       q2-classo-cv     FAILED   00:05:47 
    39111491_1        q2-cv-smoke  COMPLETED   00:07:55 
      39111920      q2-classo-sum  COMPLETED   00:00:19 
      39111921           q2-final    RUNNING   00:03:53 
    39111919_1       q2-classo-cv  COMPLETED   00:00:11 
    39111919_2       q2-classo-cv  COMPLETED   00:02:00 
    39111919_3       q2-classo-cv  COMPLETED   00:02:00 
    39111919_4       q2-classo-cv  COMPLETED   00:02:00 
    39111919_5       q2-classo-cv  COMPLETED   00:02:00 
    39111919_6       q2-classo-cv  COMPLETED   00:02:00 
    39111919_7       q2-classo-cv  COMPLETED   00:02:00 
    39111919_8       q2-classo-cv  COMPLETED   00:02:00 
    39111919_9       q2-classo-cv  COMPLETED   00:02:00 
   39111919_10       q2-classo-cv  COMPLETED   00:02:00 
   39111919_11       q2-classo-cv  COMPLETED   00:02:00 
   39111919_12       q2-classo-cv  COMPLETED   00:02:00 
   39111919_13       q2-classo-cv  COMPLETED   00:02:00 
   39111919_14       q2-classo-cv  COMPLETED   00:02:00 
   39111919_15       q2-classo-cv  COMPLETED   00:02:00 
```
