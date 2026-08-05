# Tier 2 regenerated with real feature IDs

Rebuilt with `--p-keep-original-id` so features carry their 32-character
IDs instead of positional `ASV-k` labels (DECISIONS_NEEDED item 7, option
(a)).

## Gate C1 — do the published numbers survive the relabelling?

They must: the graphical-lasso objective is invariant under simultaneous
row/column permutation, and changing `keep_original_id` changes only the
labels and the tie order.

| quantity | reference | regenerated | |
|---|---|---|---|
| selected lambda | 0.8 | 0.8 | same |
| edges | 216 | 216 | same |
| min eBIC (gamma=0.3) | 16130.0995 | 16130.0988 | ~same |

## Does the taxonomy join work?

A failed join in QIIME 2 returns all-`NaN` rather than raising, so this
is asserted rather than assumed.

- features in the transformed table : 300
- that look like real feature IDs (32 hex chars) : **300**
- unresolved against the taxonomy : **0**

**The join resolves every feature.** This is what option (a) was
for: `tax.loc[feature, "Taxon"]` now works directly, with no
mapping step and no reliance on abundance rank.

Example rows:

| feature | taxon |
|---|---|
| `1630a900d1b0…` | d__Bacteria; p__Chloroflexi; c__KD4-96; o__KD4-96; f__KD4-96; g__KD4-96; s__un |
| `3ca3a86dad4a…` | d__Bacteria; p__Gemmatimonadota; c__AKAU4049; o__AKAU4049; f__AKAU4049; g__AKA |
| `21ed96b7c0f4…` | d__Bacteria; p__Gemmatimonadota; c__Longimicrobia; o__Longimicrobiales; f__Lon |
| `c86139f5a29d…` | d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f |

## Artifacts written

To `results/tier2-regen/` only — nothing under `publish/` or `data/` was
touched. Promotion is a separate step.

- `clr.qza` (236 KB)
- `correlation.qza` (698 KB)
- `sgl-lambda08.qza` (1513 KB)
- `sgl-linear-path.qzv` (2384 KB)
- `sgl-path.qza` (1514 KB)
- `slr-mu10p0-pca.qzv` (1436 KB)
- `slr-mu10p0.qza` (2163 KB)
- `slr-mu15p0-pca.qzv` (1433 KB)
- `slr-mu15p0.qza` (2163 KB)
- `slr-mu7p5-pca.qzv` (1439 KB)
- `slr-mu7p5.qza` (2162 KB)

## Verdict

**READY TO PROMOTE — Gate C1 holds and the taxonomy join resolves**
