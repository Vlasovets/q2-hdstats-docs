> **RESOLVED 2026-08-05.** The orientation bug described here is FIXED
> (q2-gglasso merge `1c295b7`: `transform_features` no longer transposes back,
> and `calculate_covariance` takes `rowvar=False` — the two must change
> together, since either alone is wrong). Verified exactly neutral: old code
> vs new code on identical inputs, max abs difference **0.000e+00**. See
> `ORIENTATION_FIX_VERIFICATION.md`.
>
> This document records the original diagnosis and is kept for the record.
> Read it in the past tense.

---

# transform_features orientation — empirical finding

## 1. The transformer's convention

Imported a 3x2 DataFrame indexed `SAMPLE_a..c` with columns `FEATURE_x/y`.

- stored biom **sample** ids: `[np.str_('SAMPLE_a'), np.str_('SAMPLE_b'), np.str_('SAMPLE_c')]`
- stored biom **observation** ids: `[np.str_('FEATURE_x'), np.str_('FEATURE_y')]`

**The DataFrame index is treated as SAMPLES** — the documented convention holds.

## 2. Round-trip through the real action

- `atacama-counts.qza` -> `to_dataframe()` = **(13, 50)**, index[0] = `409faa5f5353e543bf6d99125c7c0e83`
- `transform-features` output -> `to_dataframe()` = **(50, 13)**, index[0] = `BAQ2420.1.1`

The output's biom **observations are sample IDs** — the artifact's axes are swapped relative to the input.

## 3. Conclusion

The transformer follows the samples-x-features convention, and `transform_features` still emits a swapped artifact. So the swap is introduced by the function itself: it reads `(p, N)`, transposes to `(N, p)` at `_func.py:88`, and something after that — most likely the metadata/relabelling path — leaves the frame indexed by features again, or the transpose is applied to a frame that was already `(N, p)`.

**Actionable:** the bug is inside `transform_features`, not in QIIME 2. Fixing it is safe for consumers that are orientation-agnostic (`pca` already is) but WILL change every stored clr/mclr artifact, so it must be paired with regenerating the published bundle.
