## Adversarial review — tier3-metagenomics

All six files reviewed against the live plugin registrations **and** against a real QIIME 2 2026.7 environment (`/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr`, q2cli 2026.7.0 with q2-gglasso + q2-classo installed editable). Every flag below was checked with `--help`, not inferred.

### Blocking — commands that fail as written

**1. `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/docs/chapters/05_metagenomics/01_gut_to_soil/01_data.md`, lines 119–122 — two nonexistent options on `qiime feature-table summarize`.**
In 2026.7 this action is a pipeline. Verified signature: parameter `--m-metadata-file`; outputs `--o-feature-frequencies` (ImmutableMetadata, required), `--o-sample-frequencies` (required), `--o-summary` (visualization). The page passes `--m-sample-metadata-file` and `--o-visualization`; both are rejected (`Error: no such option: --m-sample-metadata-file`), and the two frequency outputs are missing. Fix:
```bash
qiime feature-table summarize \
    --i-table data/gut-to-soil/asv-table-ms2.qza \
    --m-metadata-file data/gut-to-soil/sample-metadata.tsv \
    --o-feature-frequencies data/gut-to-soil/gts-feature-frequencies.qza \
    --o-sample-frequencies data/gut-to-soil/gts-sample-frequencies.qza \
    --o-summary data/gut-to-soil/gts-table-summary.qzv
```
(Downstream references to `gts-table-summary.qzv` at 01_data.md:151 and 03_regression.md:327 survive this fix.) The same wrong form appears at 01_data.md:172 in prose ("re-run `qiime feature-table summarize`") — harmless, but the command it points at must be the corrected one.

**2. `…/05_metagenomics/01_gut_to_soil/02_network.md`, line 287 — `--p-check-groups True` is not valid syntax.**
`check_groups` is plain `Bool`; q2cli renders it as a value-less flag pair. Verified: `qiime gglasso build-groups --help` → `--p-check-groups / --p-no-check-groups [default: True]`. Passing `True` makes click raise `Got unexpected extra argument (True)`. Fix: bare `--p-check-groups` (or omit — it is the default). Note the identical error is pre-existing at `chapters/02_lowdim_gglasso/06_multiple_graphical_lasso.md:122` and `:350`, which this page cites as the canonical treatment.

**3. `…/05_metagenomics/01_gut_to_soil/03_regression.md`, line 132 — `--p-rescale True` is not valid syntax.**
`rescale` is `List[Bool]`; q2cli still renders it as a flag: verified `qiime classo add-covariates --help` → `--p-rescale / --p-no-rescale ...`. Fix: bare `--p-rescale` (one occurrence per entry of `--p-to-add`). Tier 2 already does this correctly (`04_highdim_atacama/05_classo_cv.md:72`); tier 1 has the mirror-image bug (`--p-stratify False` at `03_lowdim_classo/03_regression/02_trac.md:64`), so do not copy tier 1 here.

### Substantive — wrong statements about behaviour

**4. `…/01_gut_to_soil/01_data.md`, lines 89–92 — the CLR pseudo-count description is wrong and contradicts 02_network.md.**
Text: "`clr` adds `--p-pseudo-count` to every entry before normalizing … it is applied to observed and unobserved cells alike". The implementation (`q2_gglasso/utils.py:496` `zero_imputation`) does `X = X.replace(0, pseudo_count)` and then rescales each sample back to its original total — zeros only, non-zeros untouched. 02_network.md:28–31 states this correctly, so the two pages disagree. Fix: replace lines 89–92 with the 02_network.md wording ("replaces zeros with `--p-pseudo-count`, then rescales the sample to its original total"); the "large perturbation" argument still holds, but for the zero cells only.

**5. `…/01_gut_to_soil/02_network.md`, line 59 — wrong step cross-reference.**
The `mclr` warning says the split into groups happens "in Step 5". Step 5 (line 149) is the latent block; the group split is Step 6 (line 205). Fix: "It matters in Step 6".

**6. `…/05_metagenomics/00_index.md`, lines 95–101 — the stated reason for the missing recompute is false.**
"Nothing in this tier has been run against QIIME 2 2026.7, **because no such environment exists yet**." The environment does exist (`~/.conda/envs/q2-2026.7-slr`, q2cli 2026.7.0, both plugins installed) — I used it to verify every signature in this review. Fix: keep the "no captured output" claim, change the justification to "because the recompute has not been run yet". Same wording risk at `01_gut_to_soil/00_index.md:63–67`, which only says "has not been run" — that one is fine.

**7. `…/01_gut_to_soil/03_regression.md`, lines 139–150 — the justification for putting the outcome in the design is fabricated, and the worked example leaks the target.**
The admonition claims the outcome must be added "because `add-covariates` is also how the outcome gets carried into the design table that `split-table` will divide". That is not true: `qiime sample-classifier split-table` takes the target through `--m-metadata-file` / `--m-metadata-column` (verified), and `regress` takes it through `--m-y-file` / `--m-y-column`. Nothing requires the outcome to be a design column. As written, Step 3 inserts `${OUTCOME_COLUMN}` into X with penalty weight 0.1 and Step 5 then regresses y on y, so every downstream number would be meaningless. Fix: delete Step 3 from the main path (it is optional — `qiime classo regress` defaults `c` to the zero-sum constraint when `--i-c` is omitted, and `--i-weights` can take `gts-w-taxa.qza` straight from `add-taxa`), and keep `add-covariates` only as a sidebar with a genuine covariate named in `--p-to-add`.

**8. `…/01_gut_to_soil/01_data.md`, lines 79–80 vs 150–154 — an unverified observed number, contradicted three paragraphs later.**
Line 79–80 quotes "a sampling depth of **260**, and per-sample frequencies sit in the low hundreds"; line 150–154 then says "sample depths, feature prevalences, sparsity — are pending verification … and are not reproduced here". Under the tier's own rule the 260 is a defect: it is a concrete observation about the artifact, not arithmetic, and it is not marked pending. Fix: attribute it explicitly ("the upstream tutorial's own rarefaction figure, not re-measured for this book") or drop it and let the pending note stand.

**9. `…/99_appendix/02_moshpit_cocoa_note.md`, lines 56–61 — advice that cannot be followed.**
Caveat 1 tells the reader to choose `--p-pseudo-count` and `--p-coef` "relative to the smallest non-zero value in your own table". `pseudo_count` is registered as `Int` (`q2_gglasso/plugin_setup.py:70`; CLI shows `--p-pseudo-count INTEGER`), so on a TPM/MAG table it cannot be set below 1 — the exact case the caveat is about. Only `--p-coef` (Float) is tunable. Fix: state the Int restriction and recommend `mclr` for TPM tables, which needs no floor at all.

### Nits

- **`00_index.md:77`** — "No new action. Every command in this tier has already been demonstrated in Tier 1." `qiime feature-table summarize` is *run* only in tier 3; tiers 1–2 only mention it in prose (`02_lowdim_gglasso/05_lambda_paths.md`). Reword to "no new *plugin* action".
- **`02_network.md:103`** — the grid is described as "1.00, 0.95, …, 0.10". `get_range` returns `np.linspace(min, max, n)`, i.e. ascending 0.10 … 1.00. Same 19 values, reversed presentation.
- **`03_regression.md:298`** — `--p-maxplot 200` is the registered default (`[default: 200]`), so the flag is a no-op while the prose implies a deliberate cap for this dataset.
- **`01_data.md:10`** — "(Meilander et al., 2024)" is a plain-text citation with no entry in `docs/references.bib`. Every `{cite}` key used across the six files does resolve (checked all ten).
- **`01_data.md:95`** — "it is what Network Inference uses **by default**" reads as if `mclr` were the plugin default; the plugin default is `clr`.

### Clean

- **Flags/inputs/outputs**: every `--p-`, `--i-`, `--o-`, `--m-` in all bash blocks was matched against the registered signature. Apart from items 1–3, all are valid. Specifically confirmed: no `qiime classo classify` command uses `--p-concomitant` or `--p-do-yshift` (neither is in `classify_parameters`, and `qiime classo classify --help` shows only `--p-huber` / `--p-rho`); `--p-cv-nlam`, `--p-cv--nlam`, `--p-do-yshift`, `--p-concomitant`, `--p-no-lamfixed`, `--p-lambda1-path`, `--p-mu1-*`, `--p-path-scale`, `--p-non-conforming`, `--p-n-cov`, `--p-label-size` all exist as claimed; `split-table` accepts and returns `FeatureTable[Design]` (TypeMatch includes `Design`) so Step 4 → Step 5 chains.
- **Traps**: `pca` (02_network.md:200–203) states both the latent-solution requirement (verified — `_pca/_visualizer.py:359` reads `solution/lowrank_`) and the `--m-sample-metadata-file` requirement (verified — line 361 calls `.filter_columns` unguarded). `build-groups` (02_network.md:304–314) states the `TensorData` → `List[Int]` chaining gap plus the single-`--i-covariance-matrix` blocker. `--p-rank` (02_network.md:177–184) states it always raises `NotImplementedError`, matching `_func.py:578`.
- **Verified-correct numeric claims**: 15-point default λ1 grid, 10-point default μ1 grid, 150 combinations (`utils.py:311,320`); gamma default 0.01 (`_func.py:501`); tier-1 `25pt`; tier-2 p/n ≈ 5.6 (300/54); 55,945 / 98 / 3.4 / 450 M / 7.2 GB / 18:1 all check out as arithmetic. No fabricated edge counts, eBIC values, R², selected λ or runtimes anywhere — every such slot carries a pending-verification note.
- **Links**: all 40-odd relative links in the six files resolve to real files under `docs/`; all six pages are registered in `_toc.yml`.
- **MyST**: exactly one H1 per file (the `# …` lines at 02_network.md:217 and 03_regression.md:32 are shell comments inside ```bash fences, not headings); all code and admonition fences balanced; no nested-fence problems.
- **`…/01_gut_to_soil/00_index.md`** is clean.