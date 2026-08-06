> **SUPERSEDED 2026-08-05.** Its "red flag" conclusion — that `solve-problem` may not consume its covariance argument — is wrong; the solver is fine.
>
> The shipped and regenerated correlation matrices are the **same matrix,
> reordered**: spectrum identical to 2.1e-14, and the recovered permutation
> reproduces it exactly (residual 0.000e+00). The apparent discrepancy came
> from aligning the matrices **by label** when `ASV-k` is assigned **by
> position**. See `PERMUTATION_RECOVERED.md` and `ASV_MAPPING_VALIDATION.md`.
>
> Kept for the record. Do not act on the verdict below.

---

# Side by side: same solver, two correlation matrices

## Inputs

- off-diagonal pairs compared : 44,850
- median |difference| : 2.638e-02
- pairs differing by >= 0.1 : 8,054 (18.0%)

## Outputs

| | shipped matrix | regenerated matrix |
|---|---|---|
| selected lambda | 0.8 | 0.8 |
| edges | 216 | 216 |
| min eBIC (gamma=0.3) | 16130.0988 | 16130.0988 |
| eBIC at lambda=0.8 | 16130.0988 | 16130.0988 |

## The whole eBIC curve, not just its minimum

A single matching minimum can be a coincidence. Fifteen matching grid points
cannot be. If every point agrees the two solves saw the same matrix.

| lambda | eBIC shipped | eBIC regenerated | delta |
|---|---|---|---|
| 0.300 | 16165.1289 | 16165.1289 | 1.819e-12 |
| 0.350 | 16935.4070 | 16935.4070 | 7.276e-12 |
| 0.400 | 17183.2847 | 17183.2847 | 7.276e-12 |
| 0.450 | 17285.2228 | 17285.2228 | 1.091e-11 |
| 0.500 | 17160.6331 | 17160.6331 | 3.638e-12 |
| 0.550 | 17120.0045 | 17120.0045 | 3.638e-12 |
| 0.600 | 17029.6856 | 17029.6856 | 3.638e-12 |
| 0.650 | 16897.2111 | 16897.2111 | 3.638e-12 |
| 0.700 | 16719.5750 | 16719.5750 | 3.638e-12 |
| 0.750 | 16594.2176 | 16594.2176 | 3.638e-12 |
| 0.800 | 16130.0988 | 16130.0988 | 3.638e-12 |
| 0.850 | 16313.5846 | 16313.5846 | 3.638e-12 |
| 0.900 | 16619.9001 | 16619.9001 | 0.000e+00 |
| 0.950 | 17117.1005 | 17117.1005 | 0.000e+00 |
| 1.000 | 16200.0000 | 16200.0000 | 0.000e+00 |

Grid points agreeing to 1e-6: **15 / 15**

**Red flag.** The inputs differ across 18% of pairs, yet every single
eBIC grid point agrees to 1e-6. A penalised likelihood cannot be that
insensitive to its own covariance argument. The likely explanation is
that `solve-problem` is not actually consuming the matrix it is handed
-- e.g. it recomputes or re-reads a covariance internally, making
`--i-covariance-matrix` partly decorative. That would be a plugin bug
with consequences well beyond this chapter, and it must be run down
before any number here is published.

Next probe: solve on a deliberately corrupted covariance (add 0.5 to a
block of off-diagonals). If the eBIC path does not move, the input is
being ignored and the bug is confirmed.

## Verdict

**Regeneration preserves the published result**
