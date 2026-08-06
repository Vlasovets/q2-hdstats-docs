# Adversarial verification — group "classo-tier1"

## Checks that came back CLEAN

- **Nonexistent flags (item 1): clean.** Every `--p-`, `--i-`, `--o-`, `--m-` flag in every bash block of all three files resolves to a registered name. `qiime classo classify` (`01_generate_data.md:230-242`) uses only `path/cv/cv_seed/stabsel/stabsel_seed/lamfixed` — no `--p-concomitant`, no `--p-do-yshift`. Both are correctly described in prose as *absent* from `classify`, which matches `classify_parameters` in `/home/itg/oleg.vlasovets/slr_example/q2-classo/q2_classo/_dict.py:234-277`. Inputs/outputs (`--i-taxa/--i-features/--i-c/--i-problem/--i-predictions`, `--o-x/--o-c/--o-result/--o-predictions/--o-visualization`) all match `plugin_setup.py`.
- **Broken links (item 3): clean.** All 28 relative links in the three files resolve to existing files (verified by path resolution, not assumption). The `#numerical-methods` anchor matches `## Numerical methods`.
- **Fabricated results (item 2): clean.** No concrete observed number (edge count, eBIC, R², runtime, selected λ) appears in any of the three files. Every number is an input or a registered default, and all defaults I spot-checked match the function signatures in `_func.py` (`n=100, d=80, d_nonzero=5`; `cv_seed=1`; `stabsel_seed=None`; `rho=1.345` regress / `0.0` classify; `path_nlam_log=40`, `path_lamin_log=1e-2`, `cv_nlam=100`, `cv_lamin=1e-3`, `stabsel_b=50/q=10/percent_ns=0.5/threshold=0.7/threshold_label=0.4`, `maxplot=200`). But see finding 4 about the disclaimers.
- **Traps omitted (item 4): N/A for `pca`, `build-groups`, `--p-rank`** — none of the three files mentions them (all three are q2-gglasso). Other trap findings below.
- **MyST validity (item 5): clean.** Exactly one H1 per file, all code/admonition fences balanced (30, 20, 22 fence lines), all directive names valid.

---

## FINDINGS

### 1. `01_generate_data.md:262-271` — the `--i-taxa` example crashes, and violates the file's own warning
The command uses `--p-d 20` with `../data/classification.qza`, which is the **13-ASV** Atacama taxonomy (`chapters/00_getting_started/03_download_data.md:106-107`). That is exactly the regime the warning 10 lines below (`:279-286`) forbids. I reproduced it against the real code:

```
label dtype: <U2 ;  n tips: 13
label after: ['AS','AS',...(x13),'B3','B4','B5','B6','B7','B8','B9']
RAISED: DuplicateNodeError  Duplicate tip name 'AS' found.
```

Two defects: (a) the worked command is broken; (b) the warning claims *"There is no error and no warning"* — false; `tree_to_matrix` → `skbio.TreeNode.find` raises `DuplicateNodeError` as soon as the truncated names collide.

**Fix:** change `--p-d 20` to `--p-d 13` (at `d <= len(tips)` `_func.py:71-72` takes the `label = label2[:d]` branch, which replaces the array wholesale and is safe) and drop `--p-d-nonzero` to ≤ 13; rewrite the warning to say the run *raises* `DuplicateNodeError` rather than silently truncating.

### 2. `05_advanced/02_model_selection.md:35, 100, 112` — `--p-cv-subsets` is a silent no-op; documented as functional
`q2_classo/_func.py:313` does `param.Nsubsets = cv_subsets`, but c-lasso's `CVparameters` attribute is **`Nsubset`** (singular; `classo/solver.py:599`) and `solution_CV` reads `param.Nsubset` (`classo/solver.py:986`). Verified at runtime in the 2026.7 env:
```
after q2 assignment -> Nsubset = 5 | Nsubsets = 10
```
CV always uses 5 folds. The chapter presents it as "number of CV folds" (`:112`), puts it in the worked command (`:100`), and bases the cost column on it (`:35`, "`cv_subsets` path solves"). `06_predict_and_summarize.md:222` is also affected: the CV tab renders `dicocv.Nsubset`, i.e. always 5.

**Fix:** add a gotcha stating `--p-cv-subsets` is currently inert (upstream attribute-name mismatch, `Nsubsets` vs `Nsubset`) and CV is always 5-fold; remove it from the worked command or annotate it. `90_reference/03_classo_parameters.md:211` carries the same wrong claim.

### 3. `05_advanced/02_model_selection.md:217, 243-244, 254, 376, 389` — `--p-stabsel-true-lam` / `--p-lamfixed-true-lam` are inert AND the documented semantics are inverted
`_func.py:323` and `:340` set `param.true_lam`. c-lasso never reads that attribute; it reads **`rescaled_lam`** (`classo/solver.py:691, 745`) and passes `true_lam = not param.rescaled_lam` to the solver (`solver.py:857, 1183, 1343, 1539`). Verified at runtime:
```
after q2 assignment -> rescaled_lam = True | true_lam = False
```
So the flag does nothing, **and** the effective default (`rescaled_lam=True`) means the λ you supply is interpreted as **λ/λ_max** — the opposite of the chapter's "on = the number given is λ" (`:254`, `:217`). This also invalidates the gotcha at `:389` ("the same `--p-lamfixed-lam 0.1` with `--p-no-lamfixed-true-lam` is a different penalty") and makes `--p-lamfixed-true-lam` in the worked command (`:244`) misleading. `06_predict_and_summarize.md:261-263` compounds it: the LAM-fixed tab prints `dicolam.true_lam`, i.e. the value the solver ignored.

**Fix:** state that both `*-true-lam` flags currently have no effect, that λ is always interpreted as λ/λ_max, and that the `.qzv` "Real given lambda" row reports the requested-but-unused value. Same correction needed in `90_reference/03_classo_parameters.md:237, 259`.

### 4. All three files — the "pending verification" disclaimer rests on a false premise
`01_generate_data.md:194`, `05_advanced/02_model_selection.md:345`, `06_predict_and_summarize.md:311` all assert *"No QIIME 2 2026.7 environment exists yet."* It does: `/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr` is `qiime2-2026.7.0` with `q2_classo` and `q2_gglasso` editable-installed pointing at `/ictstr01/home/itg/oleg.vlasovets/slr_example/q2-classo/q2_classo` — the very tree under review. The real blocker is different: `q2_classo` **fails to import** there —
```
_func.py:522:  problem: zarr.hierarchy.Group
AttributeError: module 'zarr' has no attribute 'hierarchy'   (zarr 3.1.5)
```
**Fix:** replace the sentence in all three notes with the accurate reason (q2-classo does not import under the env's zarr 3.1.5, which removed `zarr.hierarchy`, so no command has been executed). Leaving a false justification in place makes the disclaimer un-auditable.

### 5. `05_advanced/02_model_selection.md:215, 377-379` — `--p-stabsel-threshold-label` has *no* effect on the `.qzv`, not merely a cosmetic one
The chapter says it "controls when a bar gets a text label in the profile plot". The QIIME 2 visualizer's `plot_stability` (`q2_classo/_summarize/_visualizer.py:615-661`) is called without `threshold_label` (`:313-322`) and draws no text labels at all — `express.bar(..., hover_data=[...])` only. The labelling behaviour exists solely in c-lasso's own matplotlib `__repr__` (`classo/solver.py:1389`), which QIIME 2 never invokes.

**Fix:** say the parameter is stored in the artifact and displayed nowhere; it changes nothing in the visualization. (`90_reference/03_classo_parameters.md:243` repeats the wrong advice "set it below `stabsel_threshold` or selected features go unlabelled".)

### 6. `06_predict_and_summarize.md:22-27` — command re-creates an artifact the workflow already produced
The `predict` block is byte-identical (same three paths) to `03_regression/01_logcontrast.md:86-91`, which the chapter itself names as the source of its inputs (`:13-16`). A reader following the book in order hits q2cli's "Output path already exists" error. Same hazard, weaker, for `--i-predictions data/regress-predictions_trac.qza` at `:90` (created in `03_regression/02_trac.md:107-108`).

**Fix:** either state "you already built this in Log-Contrast Regression Step 5 — skip to Summarizing if it exists", or write to a distinct name (e.g. `data/regress-predictions_lc_ch06.qza`).

### 7. `06_predict_and_summarize.md:22-92` (minor) — working directory never stated
All paths are `data/...`, implying the tutorial root, whereas `01_generate_data.md:51-60` puts the reader in `smoke-test/` and uses `../data/`. Chapter 06 never says where to stand.
**Fix:** one line stating the commands run from the tutorial root (the directory containing `data/`).

---

**Files inspected (authoritative):** `/home/itg/oleg.vlasovets/slr_example/q2-classo/q2_classo/_dict.py`, `.../plugin_setup.py`, `.../_func.py`, `.../_tree.py`, `.../_summarize/_visualizer.py`, `.../_summarize/assets/*.html`, and the installed `classo` package at `/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr/lib/python3.12/site-packages/classo/`. No files were edited.