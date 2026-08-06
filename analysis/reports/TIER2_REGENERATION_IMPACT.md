> **SUPERSEDED 2026-08-05.** Its closing claim that the shipped matrix's provenance is unknown is wrong; the documented chain does reproduce it.
>
> The shipped and regenerated correlation matrices are the **same matrix,
> reordered**: spectrum identical to 2.1e-14, and the recovered permutation
> reproduces it exactly (residual 0.000e+00). The apparent discrepancy came
> from aligning the matrices **by label** when `ASV-k` is assigned **by
> position**. See `PERMUTATION_RECOVERED.md` and `ASV_MAPPING_VALIDATION.md`.
>
> Kept for the record. Do not act on the verdict below.

---

# Tier-2 regeneration — impact on the published numbers

Regenerating the correlation matrix from the raw counts with the DOCUMENTED
commands, then re-running the lambda path on it.

## Does the headline result survive?

| quantity | published | regenerated | |
|---|---|---|---|
| selected lambda | 0.8 | 0.8 | same |
| edges | 216 | 216 | same |
| min eBIC | 16130.0988 | 16130.0988 | ~same |

**The published result survives regeneration.** Selected lambda and edge
count are unchanged, so switching to a fully reproducible chain costs
nothing scientifically. Regenerate, republish the artifacts, and the
chapter's commands then genuinely reproduce its table.

Closest transform setting to the shipped matrix: `corr-clr-pc2.qza`.
None of clr at pseudo-count 1 or 2, nor mclr, reproduces it exactly unless the
difference above is ~0, so the shipped matrix's provenance remains unknown.

## Verdict

**SAFE TO REGENERATE**
