# Orientation fix — controlled verification

Old code vs new code, **same inputs, same parameters, run back to back**.
The earlier comparison against the shipped matrix had an uncontrolled
variable: the parameters that produced it are not recorded anywhere.

| comparison | matched labels | max abs difference |
|---|---|---|
| old code vs new code | 300 | **0.000e+00** |
| shipped vs new code | 300 | 1.147e+00 |
| shipped vs old code | 300 | 1.147e+00 |

**The fix does not change the covariance.** Old and new agree to 0.0e+00, which is exact. That is the expected result: `np.cov(A, rowvar=True)` for A of shape (p, N) and `np.cov(A.T, rowvar=False)` are the same computation, so removing the transpose and flipping `rowvar` cancel by construction.

**The shipped matrix differs from BOTH** by the same margin, so it was not produced by this chain at these parameters — a difference in the transform (clr vs mclr, pseudo-count, metadata columns), not something the fix introduced. Gate C1 consumes the shipped matrix directly and is therefore unaffected either way; but the tier-2 chapter should stop implying that `transform-features -> calculate-covariance` regenerates it until the original parameters are known.

## Verdict

**PASS — safe to merge**
