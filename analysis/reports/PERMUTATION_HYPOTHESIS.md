> **SUPERSEDED 2026-08-05.** Its "could not recover a clean permutation" line was a tie-breaking bug in the recovery, not evidence against the hypothesis; the permutation was later recovered exactly.
>
> The shipped and regenerated correlation matrices are the **same matrix,
> reordered**: spectrum identical to 2.1e-14, and the recovered permutation
> reproduces it exactly (residual 0.000e+00). The apparent discrepancy came
> from aligning the matrices **by label** when `ASV-k` is assigned **by
> position**. See `PERMUTATION_RECOVERED.md` and `ASV_MAPPING_VALIDATION.md`.
>
> Kept for the record. Do not act on the verdict below.

---

# Are the two correlation matrices the same matrix, reordered?

- shipped     : (300, 300), labels ['ASV-1', 'ASV-2', 'ASV-3'] ... ['ASV-299', 'ASV-300']
- regenerated : (300, 300), labels ['ASV-1', 'ASV-2', 'ASV-3'] ... ['ASV-299', 'ASV-300']
- label sets equal : True
- label order equal : True

## Entrywise

- POSITIONAL max|diff| (labels ignored) : **1.146726e+00**
- LABEL-ALIGNED max|diff|               : **1.146726e+00**

## Spectrum (invariant under simultaneous row/column permutation)

- max |eigenvalue difference| : **2.131628e-14**
- trace: 300.000000 vs 300.000000
- Frobenius norm: 51.059104 vs 51.059104

## Permutation recovery

- rows of the regenerated matrix matched into the shipped one : **300 / 300**
- could not recover a clean permutation, so the two are not a pure reordering of each other.

## Feature order in the inputs

- raw table observations[:3] : ['409faa5f5353e543bf6d99125c7c0e83', '1237d5925a7176fced9dda961a86c684', 'a7b877ae6d2f079a15b6b192a4425620']
- clr table observations[:3] : ['ASV-1', 'ASV-2', 'ASV-3']
- clr order == regenerated correlation order : True

## Verdict

**Same spectrum, but no clean permutation recovered.** The matrices are spectrally identical yet not a simple reordering. Investigate before republishing.

