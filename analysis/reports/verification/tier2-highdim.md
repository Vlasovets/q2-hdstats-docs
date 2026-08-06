## Findings — Tier 2 (`chapters/04_highdim_atacama/`)

Verified against source at `/home/itg/oleg.vlasovets/slr_example/q2-gglasso/` and against the shipped artifacts in `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/atacama-q2-gglasso-share.zip`.

---

### 1. `03_slr_ranks.md` — lines 131–151 and 321–330: the "single fit" commands are not single fits (HIGH)

Both the Procedure B loop and the final canonical-fit command omit the `lambda2` grid. In `/home/itg/oleg.vlasovets/slr_example/q2-gglasso/q2_gglasso/utils.py`:

```python
lambda2 = get_range(lambda2_min, lambda2_max, n_lambda2, scale=path_scale)  # -> np.array([None])
if lambda2.size == 1 and lambda2[0] is None:
    lambda2 = np.logspace(-1, -4, 5)          # size 5
...
if lambda1.size == 1 and lambda2.size == 1 and mu1.size == 1:   # False
    model_selection = False
else:
    model_selection = True
```

So every command in the `for MU in …` loop (lines 136–147) and the "What to do next" fit (lines 321–330) runs `P.model_selection()`, not `P.solve()`. Three statements become false: the heading "**one single fit per $\mu_1$**", "every grid collapsed to a single value" (line 132), and "A single fit writes **no** `modelselect_stats` group, so `summarize` will show a reduced Statistics tab with no `rank` column" (lines 149–151) — the artifacts these commands produce *will* have `modelselect_stats` and *will* have a `rank` column, which also removes the stated motivation for the zarr-export detour at lines 153–176. The page contradicts its own trap at lines 206–214, which states this rule correctly.

**Fix:** add `--p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1` to both commands. The value must be **positive** — `get_range` applies `np.log10` under the default `path_scale: log`, so `0` is unusable. (Same defect exists in the out-of-scope `02_model_selection.md` SGL command: the non-latent branch also tests `lambda2.size == 1`.)

### 2. `06_interpretation.md` — lines 116–118: "the solution artifact stores matrices, not labels" is factually wrong (HIGH)

`solve_problem` in `q2_gglasso/_func.py` ends with `P.__dict__["labels"] = dict(zip(labels_range, labels))` (taken from `covariance_matrix.columns`), and `_transformer._2` serialises the whole `__dict__`. Confirmed empirically: `problem.zip` inside the shipped `atacama-top-300-slr-lambda0.95-mu10.5.qza` contains a `labels/` group with `labels/0 … labels/299`.

**Fix:** correct the note, and replace the fragile ID recovery (lines 93–97) with a read of the stored labels, e.g. `ids = [root[f"labels/{i}"][()] for i in range(Theta.shape[0])]`, keeping the table-order check only as a cross-check.

### 3. `06_interpretation.md` — lines 93–105: the taxonomy join silently fails with the shipped artifacts (HIGH)

The shipped `atacama-top-300-clr.qza` carries feature IDs `ASV-1 … ASV-300` (its biom contains 300 `ASV-n` strings and **zero** 32-hex IDs), while `atacama-taxonomy-silva138.qza` is keyed on 32-char hex IDs (`taxonomy.tsv` header `Feature ID / Taxon / Consensus`, first key `0006edaf32a13056b89b015df9cb42dd`). So `tax.loc[hub, "Taxon"]` raises `KeyError` and `.join(tax["Taxon"])` returns all-NaN.

**Fix:** state that the snippet requires a clr table built with `--p-keep-original-id`, or map `ASV-n` → feature-id first (the bundle's `top-300-asvs.tsv` carries `feature-id / total-abundance / abundance-rank`, which is that mapping). Note also the third taxonomy column is `Consensus`, not `Confidence`.

### 4. `01_data.md` — lines 92–95 and 118: the CLR pseudo-count description is wrong (HIGH)

The page says clr "**adds** `--p-pseudo-count` to every entry … it is applied to zeros and non-zeros alike". `utils.zero_imputation` does `X = X.replace(0, pseudo_count)` and then rescales each sample back to its original column sum. The pseudo-count **replaces zeros only**; non-zero entries change only through the rescale.

**Fix:** "replaces every zero with `--p-pseudo-count`, then rescales each sample back to its original total, so only zeros are imputed." Same correction for the bullet at line 118.

### 5. `01_data.md` — lines 101–111 vs. 220–222: documented command does not reproduce the shipped artifact (MEDIUM)

Line 220–222 claims "the commands above are given so that you can check that they were derived the way this page claims", and the command passes `--p-keep-original-id`. The shipped `atacama-top-300-clr.qza` was made **without** it (ASV-n IDs, see #3). The `{important}` block at 168–176 hedges only the clr-vs-mclr question, not the ID relabelling — and line 128–131 stresses that without original IDs "you cannot trace a node … back to its taxonomy".

**Fix:** extend the pending-confirmation admonition to cover the feature-ID labelling, or state that the recompute regenerates the artifact with `--p-keep-original-id`.

### 6. `00_index.md` — lines 37, 42–45, 49–52, 56–57: result claims not covered by the page's own warning (MEDIUM)

The `{warning}` at lines 99–110 scopes itself to "the **following** pages". The index itself then asserts un-recomputed outputs as fact:
- line 37: `--p-mu1-min/--p-mu1-max = 15`, "**giving rank 2**" — an achieved rank is an output;
- lines 42–43: "$\lambda_1 = 0.8$ … is the minimiser of the extended BIC along a linear path" — a selected λ;
- lines 44–45: "can multiply the edge set several-fold" — an edge-count claim;
- lines 49–52: "$\gamma = 0.5$ … and the tier-1 default of `0.01` both select a different network on this table" — a model-selection result (it comes from the pending table in `02_model_selection.md`);
- lines 56–57: "The value 15 was found by scouting, at $\lambda_1 = 0.8$".

**Fix:** reword the warning to "everything in this tier, **including this page**", or tag each of these as expected/pending.

### 7. `03_slr_ranks.md` — line 137 vs. lines 221–229: scout grid and loop disagree (LOW)

Procedure A sweeps `5, 7.5, 10, 12.5, 15, 17.5, 20` and the "Building the map" table has all seven rows, but the Procedure B loop is `for MU in 5 7.5 10 12.5 15 20` — **17.5 is missing**. **Fix:** add `17.5` to the loop, or drop the row.

### 8. `04_latent_pca.md` — lines 162–164: panel count contradicts the table (LOW)

"a lot of grids, at **four** panels each for the rank-5 fit" conflicts with the table at line 53 ("6 panels") and with `make_plots`, which builds `itertools.combinations(range(4), 2)` = 6 plots in a 4-wide grid. **Fix:** "at six panels each (laid out on a 4×4 grid)".

### 9. `04_latent_pca.md` — lines 250, 273–275: `r` may exceed `proj`'s column count (LOW, robustness)

`r = np.linalg.matrix_rank(L)` but `utils.PCA` builds `loadings` from `np.argwhere(sig > 1e-9)`. `matrix_rank`'s tolerance (`max(shape) * eps * σ_max`) is far tighter than `1e-9`, so `r` can exceed `proj.shape[1]` and `proj[keep, j]` then raises `IndexError`. **Fix:** `for j in range(proj.shape[1])`. (The page's separate claim that `matrix_rank` is what the visualizer uses is accurate — `pair_plot` line 298 — so leave that alone.)

---

## Checks that came back clean

- **Flags.** Every `--i-`, `--o-`, `--p-`, `--m-` flag in all five files is registered for the action it is used with — `transform-features` (`table`, `taxonomy`, `sample_metadata`, `transformation`, `pseudo_count`, `add_metadata`, `scale_metadata`, `keep_original_id`, `transformed_table`), `calculate-covariance` (`table`, `method`, `bias`, `covariance_matrix`), `solve-problem` (`covariance_matrix`, `n_samples`, `latent`, `lambda1_min/max`, `n_lambda1`, `mu1_min/max`, `n_mu1`, `path_scale`, `gamma`, `solution`), `summarize` (`solution`, `width`, `height`, `label_size`, `visualization`), `pca` (`table`, `solution`, `sample_metadata`, `n_components`, `color_by`, `visualization`). No fabricated parameters.
- **`qiime classo`.** No classo command appears in any of the five files, so the `classify` + `--p-concomitant` / `--p-do-yshift` trap cannot occur here. (Aside, out of scope: `--p-cv--nlam` in `05_classo_cv.md` is **legitimate** — `cv__nlam` is a registered deprecated spelling in `q2_classo/_dict.py`.)
- **Links.** All relative markdown links in all five files resolve to real files under `docs/`, including `../../references.md`, `../90_reference/0{2,3,4}_*.md`, `../99_appendix/01_math.md`, `../05_metagenomics/00_index.md` and all `../02_lowdim_gglasso/*` targets. `atacama-top-300-sgl-lambda0.8.qza`, `…-mu10.qza`, `…-mu7.5.qza`, `…-mu15.qza`, `…-rank2.qza` are all produced by commands in this tier.
- **Citations.** All `{cite}` keys resolve against `docs/references.bib`.
- **MyST.** Exactly one H1 per file (apparent extras at `04_latent_pca.md:258` and `06_interpretation.md:87,92` are Python comments inside fences); all code and admonition fences balanced, including the correct ```` outer fence at `03_slr_ranks.md:24–43`.
- **Traps.** `pca`'s two requirements (`--p-latent` solution, `--m-sample-metadata-file`) are both stated at `04_latent_pca.md:29–37`; "`--p-rank` always raises" is stated at `03_slr_ranks.md:24–43` and `00_index.md:55`; `build-groups` is not documented in any of these five files, so the TensorData/`List[Int]` chaining gap is N/A. (Judgement call: `01_data.md:141–145, 199–205` discusses `pca` behaviour without restating the two requirements, but shows no `pca` command and cross-references `04_latent_pca.md` — acceptable.)
- **Verified as accurate, not fabricated:** metadata table at `01_data.md:186–197` (`sample-metadata.tsv` has exactly 75 data rows; `transect-name` Baquedano 32 / Yungay 43; `vegetation` no 40 / yes 35); `n = 54`, `p = 300`, 44,850 possible edges, ≈5.6 features/sample; the linear μ1 grid `5 … 20` at `n_mu1 = 7`; pair-plot counts 6 and 36 for `n_components` 4 and 9; the quoted `NotImplementedError` text and the `ValueError`-first ordering for `--p-rank` + `--p-no-latent`; default grids `np.logspace(0,-4,15)` and `np.logspace(2,-1,10)`; the strict `n_components < rank(L)` assertion and its message; the unconditional `sample_metadata` dereference; the `seq-depth` fallback and its shape-based axis inference; `filter_columns(column_type="numeric")`; the `correlated_PC` dict-overwrite caveat; and all four artifact semantic types in the `01_data.md` inventory (checked against each `.qza`'s `metadata.yaml`).