# MGL chapter: source-read claims, reproduced against QIIME 2 2026.7

**Run:** `analysis/slurm/31_mgl_verify.sh`, SLURM job 39369241, 2026-08-09.
**Env:** `q2-2026.7-slr` (QIIME 2 2026.7, gglasso 0.3.0, numpy 2.4.2, zarr 2.18.7).
**Data:** `q2-gglasso/data/atacama-counts.qza` — measured at **50 samples × 13 features**.

`docs/chapters/02_lowdim_gglasso/06_multiple_graphical_lasso.md` said of itself:

> Nothing on this page has been run against QIIME 2 2026.7, and the gaps below
> were read out of the plugin source rather than reproduced.

Source-reading is how the 1.147 discrepancy survived as long as it did. Every
claim on that page is now executed. **Six of seven confirmed; one found to be a
defect in the chapter itself.**

## Results

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| C1 | The chapter's Step 1 filter runs as written | **REFUTED — chapter defect** | shipped tier-1 metadata has no `transect-name`; command exits 1 |
| C2 | The tier-1 table splits 25/25 across transects | **CONFIRMED** | Baquedano 25, Yungay 25 |
| C3 | `build-groups` prints `p_k`, `N_k`, group count | **CONFIRMED** | `p_k [13, 13]`, `N_k [25, 25]`, 78 groups |
| C4 | Identical tables → returns nothing, action fails | **CONFIRMED** | `Expected output view type 'ndarray', received 'NoneType'`, exit 1 |
| C5 | The difference check compares samples, not features | **CONFIRMED** | both instances have `p_k = 13`, yet the array was still produced |
| G1 | `TensorData` does not chain to the `List[Int]` parameter | **CONFIRMED** | 0 occurrences of `--i-group-array`; only `--p-group-array INTEGERS...` |
| G3 | One 2-D covariance ⇒ SGL branch regardless of flags | **CONFIRMED** | `precision_` shape `(13, 13)`, i.e. 2-D |

## C1 — the chapter's own command does not run

Step 1 filters on `[transect-name]`, but the metadata shipped with the plugin
(`q2-gglasso/data/selected-atacama-sample-metadata.tsv`) has exactly five
columns:

    sample-id, ph, average-soil-relative-humidity, elevation, average-soil-temperature

There is no `transect-name`, and the command fails:

    Plugin error from feature-table:
      Selection of IDs failed with query:

The column does exist in the **75-sample** Atacama metadata used by tier 2
(`analysis/data/sample-metadata.tsv`: Baquedano 32, Yungay 43) — the chapter
took the grouping variable from one dataset and the table from another.

The transect is still recoverable for tier 1, because it is encoded in the
sample identifier prefix (`BAQ…` / `YUN…`). That is what this stage does, and
what the chapter should show.

## C2 — `--p-n-samples 25 25` is correct

Splitting the 50-sample table on the derived transect gives **25 and 25**, so
the chapter's sample counts were right. Worth stating plainly: this was one of
the numbers flagged as unverified, and it survived.

## C3 / C5 — `build-groups` works, and the check really is on samples

    Dimensions p_k:  [13, 13]
    Sample sizes N_k:  [25, 25]
    Number of groups found:  78
    Saved TensorData to: groups.qza

`p_k` is `[13, 13]` — the two instances have **identical feature sets** — and
yet `build-groups` still treated them as differing and produced an array. That
is the predicted consequence of comparing `biom.Table.to_dataframe()` **column**
labels, which are sample IDs, and it is now observed rather than inferred.

78 groups is every one of the `13 × 12 / 2 = 78` feature pairs, which is what an
identical-feature pair of instances should yield.

## The export workaround produces the documented shape

    exported array shape: (2, 78, 2)      # (2, L, K) with L = 78, K = 2
    flattened length: 312

So the chapter's `(2, L, K)` description is right, and gap 2 follows
arithmetically: 312 integers arrive as a flat `List[Int]` and nothing restores
the `(2, 78, 2)` structure.

## G1 / G3 — both gaps are real

`qiime gglasso solve-problem --help` offers only `--p-group-array INTEGERS...`;
there is no `--i-group-array` input, so the artifact cannot chain.

Passing `--p-non-conforming True` together with `--p-group-array` and a single
2-D covariance **succeeds** — and silently solves the wrong problem:

    solution/precision_ shape: (13, 13)

Two dimensions, so the single-graphical-lasso branch was taken. The
multi-instance branches need a 3-D covariance stack, which
`PairwiseFeatureData` cannot carry. The command does not error; it returns an
SGL solution while the flags suggest MGL. That is the most dangerous of the
three gaps and it is confirmed.

## What should change in the chapter

1. **Fix Step 1.** Either derive the transect from the sample-ID prefix, or say
   explicitly that the grouping column comes from the tier-2 metadata and the
   tier-1 example cannot use it as written.
2. **Replace the three "pending verification" notes** with the measured values
   above (`p_k [13, 13]`, `N_k [25, 25]`, 78 groups, shape `(2, 78, 2)`).
3. **Strengthen the gap-3 warning.** The chapter says the SGL branch is taken;
   it should say the command *exits 0* while doing so, because a silent wrong
   answer needs a louder caveat than a failure would.
4. **Keep `--p-n-samples 25 25`** — verified correct.

## Reproducing

    sbatch analysis/slurm/31_mgl_verify.sh

The stage is deliberately non-fatal on the claims themselves: several assert
that something is broken, so a non-zero exit is the expected result and is
recorded rather than treated as a stage error.
