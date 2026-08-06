# q2-hdstats recompute

Regenerates the tier-2 numbers and figures in
[q2-hdstats-docs](https://github.com/Vlasovets/q2-hdstats-docs) from the pinned
Atacama artifacts, on QIIME 2 2026.7.

It exists because the published chapters had drifted from their artifacts — in
one case printing an eBIC table underneath a command that was a single fit and
could not have produced it. Everything the book now states as a number is
written here, by a script, into `docs/_data/*.tsv`, and rendered with
`{csv-table}`. Prose cannot disagree with it.

## Layout

```
config/     params.yaml — the canonical parameterisation, with the reasoning
data/       input artifacts (gitignored; unpacked from the Atacama zips)
envs/       conda environment + the generated lockfile
scripts/    reusable steps (env yml generation, network export, packaging)
slurm/      one script per stage, plus submit_all.sh
reports/    findings and decisions that need a human
results/    solver output (gitignored — regenerates)
publish/    assembled Zenodo release (gitignored — regenerates)
```

## Running it

```bash
bash slurm/submit_all.sh     # then log out
cat reports/STATUS.md        # when you come back
```

`submit_all.sh` preflights the inputs, submits a dependency chain, and ends in a
job that writes `reports/STATUS.md` — environment versions, both test suites, the
docs build, every generated table, which artifacts exist, gate results, and what
still needs a decision. That job is chained with `afterany` on purpose, so it
still reports when an upstream stage fails instead of vanishing with it.

Every stage is idempotent: rerunning skips work that already has output.

## Stages

| stage | script | what it settles |
|---|---|---|
| 1 | `01_lambda_path.sh` | **Gate C1** — does the CLI reproduce the reference λ path? |
| 3 | `03_mu_rank_map.sh` | μ₁ → achieved rank, at the selected λ |
| 5 | `compare_rank0_rank5.py` | rank-0 vs rank-2 network comparison |
| 6 | `11_classo_cv.sh` | 15 cross-validated log-contrast fits (job array) |
| 7 | `12_classo_summary.sh` | aggregates 6 into `docs/_data/atacama-classo-cv.tsv` |
| 8 | `package_release.py` | assembles `publish/`, fills the manifest checksums |
| — | `10_orientation_probe.sh` | settles the `transform_features` axis question |
| — | `19_final_verify.sh` | rebuilds everything, writes `STATUS.md` |
| — | `20_docs_ci_parity.sh` | builds the book from a *throwaway* venv, as CI does |

## Gates

Both passed. **Gate A1**: numba 0.66 compiles GGLasso's JIT kernels under
numpy 2.4.2 with the distribution pins intact. **Gate C1**: `solve-problem`
reproduces the reference selection — λ = 0.8, **216 edges**, min eBIC 16130.0988
against a reference 16130.099454. The full path matches at 11 of 14 grid points;
the three dense-end points differing by one or two edges are recorded in the
chapter rather than smoothed over.

## Read next

- `reports/DECISIONS_NEEDED.md` — open decisions, evidence attached
- `reports/ORIENTATION_FINDING.md` — why `transform_features` **used to** emit
  swapped axes. Fixed; kept for the record, with the proof of neutrality in
  `reports/ORIENTATION_FIX_VERIFICATION.md`
- `reports/PERMUTATION_RECOVERED.md` — why the shipped correlation matrix looked
  unreproducible and was not. Read this before trusting any comparison of two
  artifacts made by label
- `reports/TIER2_KEEP_IDS.md` — the tier-2 rebuild with real feature IDs
- `reports/STATUS.md` — written by the last job of the most recent run
