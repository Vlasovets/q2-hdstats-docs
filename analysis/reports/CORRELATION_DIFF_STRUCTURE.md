> **SUPERSEDED 2026-08-05.** Its verdict — "differences are widespread, the two matrices are genuinely different objects, INVESTIGATE BEFORE REPUBLISHING" — is wrong.
>
> The shipped and regenerated correlation matrices are the **same matrix,
> reordered**: spectrum identical to 2.1e-14, and the recovered permutation
> reproduces it exactly (residual 0.000e+00). The apparent discrepancy came
> from aligning the matrices **by label** when `ASV-k` is assigned **by
> position**. See `PERMUTATION_RECOVERED.md` and `ASV_MAPPING_VALIDATION.md`.
>
> Kept for the record. Do not act on the verdict below.

---

# Where the shipped and regenerated correlation matrices differ

Comparing the 300 x 300 matrices over their 44,850 off-diagonal pairs.

## Distribution of |difference|

| bucket | pairs | share |
|---|---|---|
| < 1e-9 (identical) | 10,220 | 22.7871% |
| < 1e-6 | 2 | 0.0045% |
| < 1e-3 | 596 | 1.3289% |
| < 0.1 | 25,978 | 57.9220% |
| >= 0.1 | 8,054 | 17.9576% |

- median |difference| : **2.638e-02**
- mean   |difference| : 9.281e-02
- max    |difference| : 1.147e+00

## Are the large differences on the degenerate pairs?

- pairs with |difference| >= 0.1 : **8054**
- of those, pairs where the shipped value is >= 0.99 : **108** (1%)

**Differences are widespread**, so the two matrices are genuinely different objects. Do not republish without deciding which is canonical.

## Note for the tier-2 chapter

Independently of the regeneration question, this table contains near-duplicate features: 2 pairs at exactly r = 1.0 and 130 at r >= 0.99. Perfectly correlated CLR features are proportional in the raw counts, which for rare ASVs usually means an identical presence/absence pattern. They make the empirical correlation matrix singular, and any edge the network draws between such a pair is an artefact of that degeneracy rather than a finding. Worth stating in the chapter.

## Verdict

**INVESTIGATE BEFORE REPUBLISHING**
