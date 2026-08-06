# Recompute status — 2026-08-06T14:38:17+02:00

Written by `slurm/19_final_verify.sh` (job 39164400) on `cpusrv46.scidom.de`.
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

- q2-gglasso: `45 passed, 53 warnings in 504.82s (0:08:24)`
- q2-classo:  `22 passed in 123.51s (0:02:03)`

## Documentation

- build with `--warningiserror`: **PASS** (warning-clean)

## Generated tables (prose cannot drift from these)

- `atacama-lambda-path.tsv` — 15 rows
- `atacama-mu-rank-map.tsv` — 3 rows
- `atacama-classo-cv.tsv` — 15 rows
- `manifest.tsv` — 16 rows

### Manifest agrees with the published files


**11 verified, 0 mismatched, 0 missing.**

### Release tarballs are not stale

- `q2-hdstats-tutorial-data-tier1.tar.gz` — current
- `q2-hdstats-tutorial-data-tier2.tar.gz` — current

## Recompute artifacts

| stage | output | state |
|---|---|---|
| 1 lambda path | `atacama-top-300-sgl-linear-path.qza` | present |
| 3 mu=15 rank2 | `atacama-top-300-slr-lambda0.8-mu15.qza` | present |
| 5 comparison | `atacama-top-300-rank0-rank2-comparison-summary.tsv` | present |
| 6 classo CV | regress fits | 15/15 |

## Gates

- Gate C1 — no PASS line found in `slurm/logs/lambda_path_*.out`.
- Orientation probe: see `reports/ORIENTATION_FINDING.md`.

## Still needs you

- `reports/DECISIONS_NEEDED.md` — open decisions, evidence attached.
- Nothing is committed. `git status` in each of the three repos shows the change set.

## Job log index

```
   39111919_12       q2-classo-cv  COMPLETED   00:02:00 
   39111919_13       q2-classo-cv  COMPLETED   00:02:00 
   39111919_14       q2-classo-cv  COMPLETED   00:02:00 
   39111919_15       q2-classo-cv  COMPLETED   00:02:00 
      39112941         q2-docs-ci  COMPLETED   00:00:50 
      39113762      q2-orient-fix     FAILED   00:07:47 
      39114457      q2-orient-fix     FAILED   00:03:13 
      39115012        q2-oldvsnew  COMPLETED   00:05:01 
      39115317           q2-tests  COMPLETED   00:11:43 
      39116187        q2-regen-t2     FAILED   00:05:48 
      39116619        q2-regen-t2  COMPLETED   00:06:58 
      39117225        q2-corrdiff  COMPLETED   00:06:16 
      39123093      q2-sidebyside  COMPLETED   00:11:49 
      39123763            q2-perm  COMPLETED   00:08:27 
      39124213         q2-recover  COMPLETED   00:05:20 
      39124999           q2-tests  COMPLETED   00:16:06 
      39125194         q2-docs-ci  COMPLETED   00:01:06 
      39125825          q2-t2-ids  COMPLETED   00:11:16 
      39126519          q2-t2-ids  COMPLETED   00:10:25 
      39133812         q2-docs-ci  COMPLETED   00:01:00 
      39135665         q2-docs-ci  COMPLETED   00:01:00 
      39135666           q2-final  COMPLETED   00:07:30 
      39136342         q2-docs-ci  COMPLETED   00:00:52 
      39164399         q2-docs-ci  COMPLETED   00:01:08 
      39164400           q2-final    RUNNING   00:14:25 
```
