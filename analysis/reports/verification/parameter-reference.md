## Verdict: neither file is clean — 12 findings (2 high, 4 medium-high, 6 medium/low)

Neither file contains a single bash block, so check #1 reduces to the tables. **Good news first, since it is load-bearing:** every `--p-` / `--i-` / `--o-` / `--m-` name in both files matches the registration (`glasso_parameters` = 20 params, all present; `regress_parameters` = 36, `classify_parameters` = 34, counts in the prose are right); **no `classify` row uses `--p-concomitant` or `--p-do-yshift`** — both are correctly documented as nonexistent, with a correct explanation (`_func.py:426` forces `problem.formulation.concomitant = False`). **No fabricated results:** no edge counts, eBIC values, R², runtimes, or selected lambdas anywhere. **MyST is valid:** 12 and 14 fence lines, all paired; exactly one H1 each; all 27 relative links resolve (`01_command_coverage.md` now exists). All defaults match the function signatures.

The problems are elsewhere.

### HIGH

**1. `chapters/90_reference/03_classo_parameters.md:213` — factually wrong, and contradicts a sibling chapter.**
Row claims of `cv__nlam`: "Passing both with different values is an error." `_resolve_cv_nlam` (`q2_classo/_func.py:30-45`) only raises when `cv_nlam != 100`:
```python
if cv_nlam != 100 and cv_nlam != cv__nlam:
    raise ValueError(...)
```
So `--p-cv-nlam 100 --p-cv--nlam 50` silently uses 50. `chapters/03_lowdim_classo/05_advanced/02_model_selection.md:165-167` already documents this correctly ("does **not** raise — it silently uses 50"). **Fix:** "Passing both raises only if `cv_nlam` was changed from its default of 100; `--p-cv-nlam 100 --p-cv--nlam 50` silently uses 50."

**2. `chapters/90_reference/02_gglasso_parameters.md:162` and `:209-219` — the model-selection rule is unreachable as written for SGL.**
`n_lambda1` row says leaving it at 1 "is the usual reason model selection silently does not run", and the section says a single fit is what you get when every grid collapses. But in `utils.py:302-336`, `lambda2 = get_range(None, None, n)` returns `np.array([None])` whenever both lambda2 bounds are unset, which is then **unconditionally** replaced by `np.logspace(-1, -4, 5)` (size 5) — outside the `if latent:` branch. The non-latent test is `if lambda1.size == 1 and lambda2.size == 1`, so any SGL run that does not explicitly pin `--p-lambda2-min`/`--p-lambda2-max` gets `model_selection=True` and a `modelselect_stats` group, regardless of `n_lambda1`. (Confirmed by source read; I could not execute — no numpy in this interpreter.) **Fix:** state that a non-latent single fit additionally requires pinning lambda2 to one value, or verify against a real run before publishing. This also weakens the "look for `modelselect_stats`" diagnostic at L218.

### MEDIUM-HIGH

**3. Both files — 10 "Demonstrated in" links point at chapters that never mention the flag.** Verified by grepping each link target:

| File:line | Flag | Claimed chapter | Actually appears in |
|---|---|---|---|
| 02:82 | `--p-pseudo-count` | `02_lowdim_gglasso/01_data_preparation.md` | `04_highdim_atacama/01_data.md`, `05_metagenomics/01_gut_to_soil/01_data.md`, `.../02_network.md` |
| 02:85 | `--p-keep-original-id` | `02_lowdim_gglasso/01_data_preparation.md` | `02_lowdim_gglasso/09_interpretation.md`, `04_highdim_atacama/01_data.md` |
| 02:170 | `--p-mu1-path` | `04_highdim_atacama/03_slr_ranks.md` | `02_lowdim_gglasso/05_lambda_paths.md` |
| 02:255 | `--p-n-cov` | `02_lowdim_gglasso/09_interpretation.md` | `02_lowdim_gglasso/08_summarize.md`, `04_adaptive_glasso.md` |
| 03:135 | `--p-rescale` | `03_lowdim_classo/02_data_preparation.md` | `04_highdim_atacama/05_classo_cv.md`, `05_metagenomics/.../03_regression.md` |
| 03:182, 298 | `--p-intercept` | both `01_logcontrast.md` files | `04_highdim_atacama/05_classo_cv.md`, `06_predict_and_summarize.md`, `05_advanced/02_model_selection.md` |
| 03:212 | `--p-cv-nlam` | `04_highdim_atacama/05_classo_cv.md` | `05_advanced/02_model_selection.md`, `05_metagenomics/.../03_regression.md` |
| 03:297 | `--p-rho` | `04_classification/01_logcontrast.md` | `05_advanced/01_concomitant_formulation.md`, `02_model_selection.md` |
| 03:339 | `--p-maxplot` | `03_lowdim_classo/07_interpretation.md` | `06_predict_and_summarize.md`, `01_generate_data.md` |

**Fix:** repoint each to a chapter listed in the right-hand column. Note 03:212 is especially misleading: `05_classo_cv.md` uses the **deprecated** `--p-cv--nlam` (lines 23, 85), not the current spelling the row advertises.

**4. Both files — 9 parameters labelled `not demonstrated` are used in real bash commands elsewhere.** The page defines the label as "registered but never used in this tutorial", so these are false:
`--p-bias` (02:132 → `04_highdim_atacama/01_data.md`), `--p-width` / `--p-height` (02:252-253 → `02_lowdim_gglasso/08_summarize.md`, `04_highdim_atacama/03_slr_ranks.md`), `--p-do-yshift` (03:178 → `04_highdim_atacama/05_classo_cv.md:21`, `05_metagenomics/.../03_regression.md`), `--p-path-n-active` (03:197), `--p-cv--nlam` (03:213 → `05_classo_cv.md:23`), `--p-cv-logscale` (03:215), `--p-stabsel-percent-ns` (03:240), `--p-lamfixed-true-lam` (03:259) → all `05_advanced/02_model_selection.md` or `05_classo_cv.md`.
A further 7 (`--p-path-numerical-method` 03:196, `--p-cv-numerical-method` 03:208, `--p-stabsel-numerical-method` 03:233, `--p-stabsel-lam` 03:236, `--p-stabsel-true-lam` 03:237, `--p-stabsel-lamin` 03:241, `--p-lamfixed-numerical-method` 03:257) are discussed in tables/prose in `05_advanced/02_model_selection.md` but never run — defensible only if "demonstrated" strictly means "invoked". **Fix:** repoint the first 9; decide and state the convention for the other 7.

**5. `03_classo_parameters.md:268-274` — "`--help` misreports three of these defaults" undercounts, and omits the one that matters most.** The three named (`path`, `cv`, `lamfixed`) are right, but the registered descriptions also misreport: `cv_seed` ("Default value : None", signature `= 1`), all four `*_numerical_method` ("Default value : 'choose'", actual `"not specified"`), `path_n_active` ("Dafault value : False", actual `0`), and — critically — **`classify`'s `rho` description says "Default value = 1.345" (`_dict.py:286`) while the signature is `0.0`**. That last one directly undercuts your own advice at 03:311-312 that the user must set `rho` explicitly. **Fix:** drop the count ("misreports several"), and add the `classify`/`rho` case to the warning at 03:302-314.

**6. `03_classo_parameters.md:130-131` — `--i-features` and `--i-c` are documented as `required` under a **Default** column, but both default to `None`** (`_func.py:152-153`), so QIIME 2 registers them as optional inputs. Omitting `--i-features` does not produce a clean CLI error: `d = len(features.columns)` (`_func.py:158`) raises `AttributeError: 'NoneType' object has no attribute 'columns'`. **Fix:** mark both `None` / optional and add a note that omitting `--i-features` crashes with an `AttributeError` rather than a usage message.

### MEDIUM / LOW

**7. `02_gglasso_parameters.md:85` — `keep_original_id` row omits two things.** `rename_index_with_sum` (`utils.py:548-556`) runs **unconditionally**, so `transform-features` always reorders rows by total abundance even at the default `True`; and the sort is **ascending**, so `ASV-1` is the *least* abundant feature, not the most. **Fix:** "Rows are always reordered by ascending total abundance; with `False` the index is additionally replaced by `ASV-1…ASV-p`, where `ASV-1` is the least abundant."

**8. `02_gglasso_parameters.md:84` — `add_metadata` row omits that only numeric columns are used.** `_func.py:91-92` calls `sample_metadata.filter_columns(column_type="numeric")`, so categorical metadata is silently dropped, and missing values are filled with `0` behind a `warnings.warn("Missing values are imputed with 0!")`. Both are traps for anyone expecting a categorical covariate to become a node. **Fix:** add both facts to the row or the surrounding paragraph.

**9. `02_gglasso_parameters.md:42-44` — internal contradiction on `n_samples`.** The intro says the list-typed params including `n_samples` are "consumed as scalars, so in practice you pass them exactly once", but L159 correctly says "one value per instance", and `list_to_array` (`utils.py:45-50`) only unwraps to a scalar when `len == 1`. **Fix:** remove `n_samples` from that sentence.

**10. `02_gglasso_parameters.md:269` — "`--p-path-scale` → `ValueError` inside grid construction" is conditionally false.** `get_range` (`utils.py:243-244`) returns `np.array([None])` *before* the scale check when both bounds are unset, so a misspelled `--p-path-scale` is silently ignored on any grid you did not bound. **Fix:** "raises `ValueError` — but only when at least one bound is set; otherwise silently ignored."

**11. `03_classo_parameters.md:97` — "returns a sample-major table" overstates.** `transform_features` (`_func.py:112-122`) never transposes; it centres along rows and preserves the caller's index/columns. It *assumes* samples-in-rows input. **Fix:** "preserves the input orientation and centres along rows, so it assumes samples in rows (which is what `regress` expects)."

**12. `02_gglasso_parameters.md:172` — the gglasso page has no equivalent `--help` caveat.** `latent`'s registered description says "The default is False" (`plugin_setup.py:197-200`) while the signature is `None`, mirroring the classo prose/signature drift you do warn about. **Fix:** either add a short warning or accept the asymmetry deliberately.