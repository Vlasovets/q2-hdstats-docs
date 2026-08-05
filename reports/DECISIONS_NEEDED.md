# Decisions needed

Seven items surfaced during the 2026.7 migration that change **behaviour, artifact
semantics, or a release commitment** rather than prose. Evidence is recorded so
each decision can be made without re-deriving it. Items 1 and 7 have been acted
on; the rest are deliberately left unfixed.

Status of everything else: Gate A1 and Gate C1 both PASS, the plugin tests pass,
the book builds with zero warnings. See `project_q2_hdstats_migration` in the
Claude memory index for the full picture.

**Read item 7 first if you are about to mint the Zenodo DOI** — it is the only
one that becomes irreversible at that point.

---

## 1. `transform_features` swapped axes — **DECIDED AND FIXED** (2026-08-05)

**Fixed**, merged to `christian-review-fixes` as `1c295b7`. Two coupled edits:
`transform_features` no longer transposes back to `(p, N)` before returning, and
`calculate_covariance` takes `rowvar=False` (it had been relying on the swap —
`np.cov` defaults to rows-as-variables, so the two bugs cancelled).

**Why now rather than later:** the Zenodo DOI is not yet minted, so no published,
citable artifact depended on the old orientation. That stops being true the
moment the record goes live, which made this the cheapest possible moment.

**Verified neutral, exactly.** Old code vs new code on identical inputs and
parameters, run back to back: **max abs difference 0.000e+00**. This is the
expected result — `np.cov(A, rowvar=True)` on `(p, N)` and
`np.cov(A.T, rowvar=False)` are the same computation. See
`ORIENTATION_FIX_VERIFICATION.md`.

**~~New finding, still open.~~ RETRACTED 2026-08-05 — the discrepancy was an
artefact of my own comparison.** For two days the shipped
`atacama-top-300-correlation.qza` appeared not to be reproducible by
`transform-features -> calculate-covariance`: max|diff| 1.147, later refined to
"18% of pairs differ by >= 0.1, median 2.6e-02". All of those numbers came from a
comparison that aligned the two matrices **by label** (`df.loc[common, common]`)
when the labels are assigned **by position**. It was comparing organism *i*'s row
against organism *j*'s.

The two matrices are the **same matrix, reordered**. Verified three ways:
spectrum identical to 2.1e-14 with identical trace and Frobenius norm; the
recovered permutation reproduces the matrix exactly (`A[perm][:,perm] == B`,
residual **0.000e+00**); and the lambda path, eBIC at all 15 grid points and the
216 edges are bit-identical, because the graphical-lasso objective is invariant
under simultaneous row/column permutation.

Cause: `rename_index_with_sum` ranks by ascending total abundance and
`--p-no-keep-original-id` assigns `ASV-k` by position. Commit `9a4c08b`
(2026-07-19) replaced `df.sort_index()` — quicksort, **unstable** — with a stable
sort. The shipped matrix was written **2026-06-27**. 209 of 300 features share a
total-abundance value, and exactly **158 features moved, all of them within their
own abundance level** (0 crossed levels).

**The documented chain does reproduce the shipped matrix.** Gate C1 stands. No
question for C. Müller here; only the λ=0.8 item below remains.

**What this did surface, and it is worse:** `ASV-k` is not a stable identifier,
and the recovery procedure the tutorial documented is broken. See item 7.

See `PERMUTATION_RECOVERED.md` and `ASV_MAPPING_VALIDATION.md`.

<details><summary>Original write-up, kept for the record</summary>

### `transform_features` stores its output with swapped axes  — HIGHEST IMPACT

**What happens.** `transform_features` returns a `(p, N)` DataFrame — features as
rows, per the `# p, N` comment at `q2_gglasso/_func.py:68`. QIIME 2's
DataFrame→FeatureTable transformer interprets any DataFrame as *samples × features*.
So the stored artifact has the two axes reversed.

**Verified directly:**

| artifact | `.view(biom.Table).to_dataframe()` | index holds |
|---|---|---|
| `atacama-counts.qza` (raw) | (13, 50) | MD5 feature IDs |
| `mclr.qza` (transform output) | (50, 13) | **sample IDs** |

**Why nobody noticed.** `calculate-covariance` consumes the transformed table in a
way that compensates, so the SGL/SLR pipeline produces correct 13×13 results. The
swap only becomes visible when something else reads the artifact.

**What it currently breaks.**
- `qiime gglasso pca` — worked around: `pca` now matches the table against the
  low-rank matrix's shape instead of assuming an orientation, so both work.
- `qiime gglasso build-groups` on mclr tables — **still broken**. Because the
  column labels are feature IDs rather than sample IDs, both transects appear to
  have identical features, `non_conforming_problem` stays False, the function
  returns `None`, and the CLI fails with
  `TypeError: Expected output view type 'ndarray', received 'NoneType'`.
  This makes Step 2 of `02_lowdim_gglasso/06_multiple_graphical_lasso.md`
  non-runnable as written.

**Options.**
- **(a) Fix it** — transpose before returning. Correct, but every existing
  `*-mclr.qza` / `*-clr.qza` on disk becomes wrong-way-round relative to new ones,
  including artifacts in the shared tutorial bundle. Needs a version bump and a
  migration note.
- **(b) Leave it, document it** — keep compensating at each consumer, as `pca` now
  does. Cheap, but every future reader of these artifacts hits the same trap.
- **(c) Fix it and re-generate** the published bundle as part of the recompute,
  so old and new artifacts are never mixed.

Recommendation: **(c)** if the recompute is happening anyway, since it removes the
mixed-state risk that makes (a) unattractive.

</details>

---

## 2. `--p-stabsel-true-lam` / `--p-lamfixed-true-lam` — behaviour change already applied

**Applied, flagging because it changes results.** Both parameters were **inert**:
q2-classo assigned `param.true_lam`, but c-lasso reads `rescaled_lam`. Assigning an
undefined attribute silently created a new one that nothing consulted.

The correct mapping is an **inversion**. c-lasso, `solver.py:643`:
`rescaled_lam (bool) : (only used if method = 'lam') False if lam = lambda, [True]
if lam = lambda/lambdamax`. q2's `true_lam` means "the lambda given is the real
lambda". So `rescaled_lam = not true_lam`, now applied at all four call sites.

**Blast radius:** only consulted when `stabsel_method == 'lam'`; the default is
`'first'`. Two `.qzv` templates that rendered `{{ ....true_lam }}` were updated to
`rescaled_lam` — otherwise that row would have rendered *blank* rather than
erroring, because q2templates uses jinja2's default `Undefined`.

**Confirm:** is `rescaled_lam = not true_lam` the intended semantics, or was
`true_lam` meant to be passed through unchanged under a different name?

---

## 3. `lambda*_min/max`, `mu1_min/max` registered as `List[Float]`

`_dict.py:11-16` registers all six as `List[Float]`; `_func.py:485-493` annotates
them `float`. `get_range` then builds an `(n, 1)` grid rather than `(n,)`:
`get_range([0.001],[1.0],15)` → shape (15, 1); `get_range(0.001,1.0,15)` → (15,).

Singleton pins collapse via `.item()`, so only multi-point paths are affected.
**Numerically benign for the verified case** — Gate C1 reproduced the reference
exactly through this path. But it is a CLI-contract mismatch: the CLI advertises a
list where the implementation wants a scalar.

Fixing changes the registered signature. Low urgency, non-zero risk.

---

## 4. Is `selected-atacama-sample-metadata.tsv` going on the Zenodo record? — **RESOLVED** (2026-08-05)

**Yes, and it already is.** `package_release.py` sources it from
`q2-gglasso/data/` (the tier-1 `GG` entry) and it is present in
`publish/tier1/`. Its manifest row carries real values, verified against the
shipped copy rather than transcribed:

| | |
|---|---|
| bytes | 2599 (matches) |
| sha256 | `d0e1421c51097405…` (matches) |
| shape | 76 rows × 5 columns — `sample-id`, `ph`, `average-soil-relative-humidity`, `elevation`, `average-soil-temperature` |
| url | `ZENODO_DOI_PENDING` |

Six chapters depend on it: `02_lowdim_gglasso/{01_data_preparation, 04_adaptive_glasso,
06_multiple_graphical_lasso, 07_pca, 08_summarize}` and
`00_getting_started/03_download_data`.

The only thing outstanding is the URL, and that is true of **all 11** tier-1 and
tier-2 rows — it is filled when the DOI is minted, not a decision. Nothing here
blocks the release.

---

## 5. λ = 0.8 vs λ = 0.95 — confirm with C. Müller

Now stronger than when the plan was written. Gate C1 reproduced the **whole** λ
path, and both values are points on it:

| λ | edges | source that uses it |
|---|---|---|
| 0.95 | 145 | `atacama-q2-gglasso-share.zip` |
| 0.80 | 216 | the tutorial chapters, `tutorial-notes.md`, `run_q2_gglasso_lambda08_models.sh` |

So the share bundle is not wrong, just an **earlier point on the same path**.
Recommendation stands: λ=0.8 canonical, regenerate the bundle. One email to
confirm before the recompute overwrites it.

---

## 6. MOSHPIT cocoa — scope confirmed, restating the consequence

Decision taken: demote to an appendix note, no data shipped. Consequence worth
re-stating: tier 3 therefore rests entirely on gut-to-soil, which is **16S
amplicon, not shotgun metagenomics**. If the F1000 reviewer specifically wants
shotgun data, that is not yet answered. The cocoa dataset (n=14) is too small for
covariance estimation regardless, so the real alternative is asking Evan for a
larger mOTU table.

---

## 7. `ASV-k` is a position, not an identifier — **DECIDED AND FIXED** (2026-08-05)

Surfaced by running down the "1.147 discrepancy" in item 1. `ASV-k` is assigned
by **position** in an abundance ranking, so it names a different organism
depending on which build produced the artifact. Two consequences, one already
handled and one not.

**Handled — the plugin.** `rename_index_with_sum` now breaks ties on the feature
ID, so the ordering is a pure function of the table's contents rather than of the
row order it happened to arrive in. Regression test added
(`test_ties_break_on_feature_id_not_input_order`); confirmed to fail against the
previous logic. Ties are not an edge case here — **209 of the 300** Atacama
features share a total-abundance value (61 groups, largest 13).

**Not handled — the published bundle, and this is the user-facing part.**
`06_interpretation.md` told readers to recover feature identity from
`top-300-asvs.tsv` via `ASV-n <-> abundance-rank = 301 - n`. **That procedure is
broken.** Tested directly: permuting the correlation matrix by that mapping fails
to reproduce the shipped matrix, off by **1.137**. Only the **91** features with a
unique total abundance land correctly; the other 209 receive a neighbour's
taxonomy and nothing raises. The chapters have been corrected to remove the
procedure and to state that `--p-keep-original-id` is the only reliable route.

**The decision.** Every tier-2 artifact in the bundle carries `ASV-k` labels from
the pre-fix ordering. Options:

- **(a) Regenerate the whole tier-2 chain with `--p-keep-original-id`.** Features
  carry their real IDs, the taxonomy join works directly, and the labels stop
  being version-dependent. Costs: node labels in every figure become 32-character
  hashes unless shortened for display, and every tier-2 artifact is rebuilt.
- **(b) Regenerate with `--p-no-keep-original-id` under the fixed tie-break.**
  Keeps the readable `ASV-k` labels and they are now deterministic, but they still
  cannot be joined to taxonomy without shipping the mapping alongside, and they
  differ from the currently published ones.
- **(c) Ship as is, document the hazard.** Cheapest, but leaves a published bundle
  whose feature identities cannot be recovered by any documented means.

Recommendation: **(a)**, with a short display-name helper for the figures. It is
the only option that makes the taxonomy join work, and the recompute is happening
anyway. Note the numbers do not change under any option — the objective is
permutation-invariant, so λ = 0.8, 216 edges and eBIC 16130.0988 hold throughout.
Only feature identity is at stake.

**~~Must be decided before the Zenodo DOI is minted~~ — DECIDED AND EXECUTED
2026-08-05: option (a).**

The tier-2 chain was rebuilt with `--p-keep-original-id` (`slurm/28_regenerate_
tier2_keep_ids.sh`) and promoted to `publish/tier2/`.

Verified, not assumed:

| check | result |
|---|---|
| Gate C1 | λ = 0.8, **216 edges**, eBIC 16130.0988 vs reference 16130.0995 |
| feature IDs | **300/300** are real 32-hex IDs |
| taxonomy join | **0 unresolved** against `atacama-taxonomy-silva138.qza` |
| PCA outputs | 3/3, after deriving `n_components` from each achieved rank |
| manifest vs files | 11 rows verified, 0 mismatched, 0 missing |

Two things the pre-commit review caught that would have undone it:

- `package_release.py` sourced the two derived tier-2 artifacts from `data/`,
  which still holds the June `ASV-k` copies. As stage 8 of the documented
  pipeline it would have reverted the bundle, rewritten the manifest to the old
  checksums, and exited 0. Now sourced from `results/tier2-regen/`; `data/` is
  deliberately left as the historical record.
- Nothing verified the manifest against the files it describes. `19_final_verify`
  now recomputes every tier-1/tier-2 size and sha256 and flags stale tarballs.

Figure readability, the known cost of (a), is handled by a `display` column in
`export_network.py` — deepest informative rank plus a short ID slice, e.g.
`Rubrobacter (a7b877)` — with the full ID kept as the key and a hard failure if
two display names would collide.

Residual follow-up, not blocking: no command in the book demonstrates
`--p-keep-original-id False`, so that half of the parameter is now undocumented
by example. Logged as coverage debt in `90_reference/01_command_coverage.md`.
