Verified both files against the plugin source, the registered parameter dicts, the q2cli parser, the actual Atacama metadata, and the installed `q2-2026.7-slr` environment. **Flags are clean** (see note at end); the problems are factual/consistency ones.

## Findings

### 1. `06_multiple_graphical_lasso.md`, lines 172–182 — the export workaround does not run against the documented environment
The Python snippet uses `zarr.ZipStore(...)`, which is the **zarr v2** API. The environment this tutorial targets (`/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr`, QIIME 2 2026.7, q2-gglasso installed editable from `/home/itg/oleg.vlasovets/slr_example/q2-gglasso`) has **zarr 3.1.5**, where `zarr.ZipStore` does not exist (`hasattr(zarr,'ZipStore') == False`; it moved to `zarr.storage.ZipStore`). Verified by running it: `import q2_gglasso` itself dies with `AttributeError: module 'zarr' has no attribute 'hierarchy'` at `q2_gglasso/_pca/_visualizer.py:330`.
**Fix:** either write `from zarr.storage import ZipStore` (v3) or state explicitly that the snippet — like the plugin's own `_transformer.py` — assumes `zarr<3`, and add that q2-gglasso currently fails to import under the zarr 3 shipped in QIIME 2 2026.7. This is arguably a fourth "known gap" and it blocks Step 3 entirely.

### 2. `05_lambda_paths.md`, lines 174–179 — false justification in the "pending verification" note
> "No environment with q2-gglasso installed exists yet, so nothing in this chapter has been re-run"

This is wrong. `q2_gglasso-0.0+259.gea3e63d.dirty` and `q2_classo-0.0.1+148.g5923fb3.dirty` are both installed (editable, `direct_url.json` points at the local checkouts) alongside `q2cli 2026.7.0` and `gglasso 0.3.0` in `q2-2026.7-slr`. The environment exists; it is *broken* (finding 1), which is a different claim.
**Fix:** replace with the accurate reason, e.g. "these values have not been recomputed; the `q2-2026.7-slr` environment currently cannot import q2-gglasso (zarr 3 incompatibility)." Withholding the numbers is right; the stated reason is not.

### 3. `06_multiple_graphical_lasso.md`, lines 191, 206–208 (repeated at 257, 269) — `--p-n-samples 32 43` are the wrong table's counts
I read `share/sample-metadata.tsv` from `atacama-q2-gglasso-share.zip`: 75 rows, `transect-name` = `{Yungay: 43, Baquedano: 32}`. So 32/43 are the counts of the **full** Atacama metadata. But every command on this page filters `data/atacama-counts.qza` using `data/selected-atacama-sample-metadata.tsv`, and Tier 1 is documented as **N = 50** (`01_data_preparation.md:4`; `--p-n-samples 50` in `02_sgl.md`, `03_slr.md`, `04_adaptive_glasso.md`, and in `05_lambda_paths.md` itself). 32 + 43 = 75 ≠ 50, so these cannot be the group sizes of the table being filtered. The chapter also contradicts itself: line 60 says "the per-transect sample counts of the 13-ASV subset are being confirmed and are not quoted here", then line 206 quotes them.
**Fix:** drop the numbers, use a placeholder (`--p-n-samples <N_baquedano> <N_yungay>`), and delete the sentence attributing them to "the Atacama sample metadata" — or state explicitly that they come from the 75-sample full metadata and are *not* the Tier 1 split.

### 4. `06_multiple_graphical_lasso.md`, lines 359–366 vs 138–152 — the two admonitions contradict each other
The Step 5 `{important}` says the two vegetation subtables "may well end up with identical feature sets after filtering, in which case `build-groups` reports that the datasets match and produces no array." But the Step 2 `{warning}` (correctly) explains that `build_groups` compares `df.columns` of `biom.Table.to_dataframe()`, i.e. **sample IDs**, not features (confirmed at `q2_gglasso/_func.py:191`, `columns_dict[i] = df.columns.values.tolist()`). Two tables produced by splitting on a metadata column always have disjoint sample sets, so `non_conforming_problem` is always `True` and an array is *always* returned — the "datasets match" outcome the `{important}` predicts is unreachable for this workflow.
**Fix:** rewrite the `{important}` to say the opposite: because the check is on sample IDs, `build-groups` will return an array even when the feature sets are identical, so the array must be inspected before use (which is what the Step 2 warning already tells the reader).

### 5. `06_multiple_graphical_lasso.md`, lines 295–299 — the `--p-reg` warning contradicts gap 3
> "`--p-reg ggl` (lower case) ... is not rejected by the command line; the failure surfaces inside the solver."

With the only reachable input (a 2-D `PairwiseFeatureData`), `solve_problem` takes the `S.ndim == 2` branch and calls `solve_SGL`, which never passes `reg` to `glasso_problem` (`_func.py:272, 282`). So an invalid `--p-reg` is **silently ignored**, not a failure — exactly what gap 3 (lines 226–231) says. It is `reg` has no `Choices()` that is true (verified: `Choices` is imported but never used in either plugin's `_dict.py`/`plugin_setup.py`).
**Fix:** "`--p-reg` accepts any string — there is no `Choices()`. Today an invalid value is silently ignored, because gap 3 means the MGL branch is never entered. Once the multi-instance input type lands, an invalid value will fail inside the solver."

### 6. `06_multiple_graphical_lasso.md`, lines 287–292 — "the boundary check ... extends" overstates the shared behaviour
The directional λ2 texts are exact (`utils.py:453, 456`). But `05_lambda_paths.md:212–231` says solve-problem emits the directional warning *followed by* `lambda is on the edge of the interval, the solution might have not reached global minimum!`. For MGL that second warning is different: `solve_MGL` emits only `"The solution might have not reached global minimum!"` (`_func.py:355`), with no `lambda is on the edge of the interval,` prefix — unlike `solve_SGL` (`:278`) and `solve_non_conforming` (`:435`).
**Fix:** note the MGL variant's shorter second warning, or drop the "extends ... as described in the previous chapter" framing.

### 7. `06_multiple_graphical_lasso.md`, lines 320–357 — Step 5 has no covariance matrices
Step 5 produces `atacama-counts-veg-{yes,no}-observed.qza` and a group array, then says "then solve with `--p-non-conforming True` and the resulting integers in `--p-group-array`." No `transform-features` / `calculate-covariance` steps are given for the vegetation split, so a reader following literally has no artifact for the required `--i-covariance-matrix`.
**Fix:** either add the two transform/covariance calls (as in Step 1) or say explicitly that the solve step reuses the Step 3 command shape and there is nothing new to run.

### 8. `06_multiple_graphical_lasso.md`, line 349 vs line 121 — inconsistent inputs to `build-groups`
Step 2 passes mclr-transformed tables; Step 5 passes raw count tables. Harmless (only labels matter) but it reads as an inconsistency and invites a reader to wonder which is required.
**Fix:** one sentence saying `build-groups` only reads labels, so either table works — or use the transformed tables in both places.

### 9. `05_lambda_paths.md`, lines 315–317 — `--p-rank` failure mode is only half stated
Correct that it always raises, and `_gglasso_supports_rank()` returns `False` against the installed gglasso 0.3.0 (verified: `fix_latent_rank` not in `ADMM_SGL` signature). But the raise is `NotImplementedError` only when `latent` is truthy; with `latent` unset/False, `solve_problem` raises `ValueError("The 'rank' parameter is only meaningful for the sparse + low-rank model; set latent=True.")` first (`_func.py:572–585`).
**Fix:** "raises `ValueError` without `--p-latent`, and `NotImplementedError` with it — there is no configuration in which it works."

### 10. Low / cross-cutting — `lambda*_min`/`max` are registered as `List[Float]`
`_dict.py` declares `lambda1_min`, `lambda1_max`, `lambda2_min`, `lambda2_max`, `mu1_min`, `mu1_max` as `List[Float]`, while `solve_problem` annotates them `float`. `get_range` then calls `np.logspace(np.log10([x]), np.log10([y]), n)`, which yields an `(n, 1)` grid rather than `(n,)` for every multi-point path in these chapters. Singleton pins are unaffected (`.item()` collapses them). Not introduced by these chapters (`02_sgl.md` has it too) and unverified end-to-end because the plugin will not import, but it means the "three ways" commands may not behave as narrated.
**Fix:** out of scope for these files; worth an issue against `_dict.py`.

## Checked and clean
- **Flags.** Every `--p-`/`--i-`/`--o-`/`--m-` in both files resolves to a registered name: all of `05`'s `solve-problem` params are in `glasso_parameters`; `summarize` `--p-label-size`/`--i-solution`/`--o-visualization`, `transform-features`, `calculate-covariance`, `build-groups` (`--i-tables`, `--p-check-groups`, `--o-group-array`), and the `q2-feature-table` / `metadata tabulate` / `tools export` calls all match `plugin_setup.py`. `--i-group-array` appears only in the prose asserting it does *not* exist — correct.
- **No `qiime classo` command appears in either file**, so `--p-concomitant` / `--p-do-yshift` on `classify` cannot occur (grep for `classo|concomitant|do-yshift` returns nothing).
- **CLI syntax.** `--p-latent True` and `--p-n-samples 32 43` / `--p-lambda1-path 1.0 0.5 …` are valid: q2cli 2026.7 uses `store_maybe` for `Bool` and `append_greedy` for `List` (`q2cli/click/{option,parser}.py`). Note `chapters/90_reference/02_gglasso_parameters.md:111` says `--i-tables` must be repeated per instance — both forms work, so that reference row is merely imprecise, not a defect in these two files.
- **Source-derived claims all verified exactly:** default paths `logspace(0,-4,15)` / `logspace(-1,-4,5)` / `logspace(2,-1,10)` and the `"Default values for lambda1 have been used."` text (`utils.py:310–321`); missing-bound fallbacks 1e-3 / 1 (`utils.py:245–248`); `ValueError: Unknown scale …, use 'log' or 'linear'.` (`utils.py:253`); no `lambda2_path`; the model-selection singleton rule incl. the latent three-way case (`utils.py:318–336`); directional boundary warnings and the absence of any μ1 check (`utils.py:438–457`); `summarize` guarding on `modelselect_stats` (`_summarize/_visualizer.py:234`); `check_G` printing p_k / N_k / group count; `build_groups` returning `None` on match; the flat-`List[Int]` non-reshape; single `--i-covariance-matrix` and the `S.ndim == 2/3` dispatch; `non_conforming` hardcoding `reg="GGL"`; `taxonomy` required-but-unread.
- **Arithmetic.** 1.64× ratio (10^(3/14)=1.638), 0.071 linear step (0.9990/14=0.0714), "nine of fifteen below ≈0.072" (9 exactly), "two candidates above 0.6" (0.6105, 1.0), 15×10=150, 10×5=50 — all correct, all input-derived.
- **Links.** All eleven relative links resolve: `../04_highdim_atacama/{00_index,02_model_selection,03_slr_ranks}.md`, `../90_reference/{02_gglasso_parameters,04_troubleshooting}.md`, `{02_sgl,03_slr,05_lambda_paths,07_pca,08_summarize,09_interpretation}.md`.
- **Citations.** `foygel2010extended` and `danaher2014joint` both present in `docs/references.bib`.
- **Metadata columns.** `transect-name` and `vegetation` (yes/no) both exist, so the `--p-where` clauses are valid.
- **MyST.** Fences balanced in both files (12 pairs in `05`, 19 in `06`); exactly one real H1 each — the extra `^# ` hits are shell comments inside bash blocks.
- **Traps.** `build-groups` TensorData/`List[Int]` chaining gap: present and correct (Step 3 + gap 1). `--p-rank` always-raises: present, see finding 9 for the incomplete half. Neither chapter documents `pca`, so the `--p-latent` / `--m-sample-metadata-file` requirement does not apply here.