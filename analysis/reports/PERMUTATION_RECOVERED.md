# The permutation, recovered and verified

Job 26 showed the shipped and regenerated correlation matrices are
spectrally identical (max |eigenvalue difference| 2.1e-14) yet disagree
entrywise by 1.147 -- a reordering, not a numerical defect. This pins the
reordering down exactly.

## Matching

- features                       : 300
- colour classes of size > 1     : 2 (sizes [2, 2])
- classes that failed to match   : 0

Classes of size > 1 are genuinely duplicate rows -- feature pairs at
exactly r = 1.0. Any bijection inside such a class reproduces the matrix,
which the residual below confirms.

## Verification

- `A[perm][:,perm] == B` residual : **0.000e+00**
- features whose position changes : **158 / 300**

**Exact to machine precision.** The shipped matrix and the
regenerated matrix contain the same numbers. There is no numerical
discrepancy between them -- only a relabelling.

## Is the reordering confined to abundance ties?

`rename_index_with_sum` ranks features by ascending total abundance.
Commit `9a4c08b` (2026-07-19) replaced `df.sort_index()` -- quicksort,
**unstable** -- with `sort_values(kind="stable")`. The shipped matrix
was written 2026-06-27, before that change. So the reordering should
move features only WITHIN groups of equal total abundance.

- distinct total-abundance values : 152 / 300
- features sharing a value with another : **209**
- moved features that stay at the same abundance level : **158 / 158**
- features moving ACROSS abundance levels : **0**

**Confirmed.** Every displaced feature keeps its total abundance,
so the two orderings differ only in how ties are broken. That is
exactly the predicted effect of the stable-sort change, and it
fully accounts for the apparent 1.147 discrepancy.

## What this settles

- The documented chain DOES regenerate the shipped correlation matrix.
  The "unknown provenance" finding is retracted.
- Gate C1, lambda = 0.8, 216 edges and eBIC 16130.0988 stand unchanged;
  the graphical-lasso objective is permutation-invariant.
- `ASV-k` is **not** a stable identifier. It is an abundance rank assigned
  by position, so it names a different organism depending on which plugin
  version built the artifact.

## Verdict

**Same matrix, reordered by the tie-break change — no numerical defect**
