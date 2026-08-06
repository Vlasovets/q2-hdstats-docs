# Does the documented `ASV-k` recovery actually work?

`06_interpretation.md` route 2 maps `ASV-n` to `abundance-rank = 301 - n`
using `top-300-asvs.tsv`. Testing that against the artifacts.

## Result

- TSV total-abundance matches the table : True
- permuting the regenerated matrix by route 2 reproduces the shipped matrix to : **1.137e+00**
- positions where route 2 differs from the current plugin order : 146 / 300
- features with a UNIQUE total abundance (route 2 is safe for these) : **91 / 300**
- features sharing an abundance value (route 2 is a guess for these) : **209 / 300**

**Route 2 does NOT recover the shipped ordering.** Readers following the documented procedure silently attach the wrong taxon to tied features — up to 209 of 300. The join raises nothing, so the error is invisible.

## Recommendation

Route 1 (`--p-keep-original-id`, now the default) is the only reliable path and should be presented as such rather than as one of two options. Route 2 should carry an explicit warning that it is exact only for the 91 features whose total abundance is unique.

The durable fix is in the plugin: break ties on the feature ID so that `ASV-k` becomes a pure function of the data rather than of the input row order — `df.sum(axis=1).sort_values(kind='stable')` replaced by a sort on `(row_sum, feature_id)`.

## Verdict

**Route 2 is broken — silently mislabels tied features**
