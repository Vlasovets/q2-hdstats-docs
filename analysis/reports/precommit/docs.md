> **HISTORICAL SNAPSHOT — 2026-08-04 pre-commit review of the docs repo.**
>
> Two of its findings turned on repository structure that has since changed:
> **docs-13** noted the generator scripts sat in a directory that was "not a
> git repository at all, so nothing is committed anywhere". Those scripts were
> subsequently versioned as `q2-hdstats-recompute`, and on 2026-08-05 that repo
> was merged into this one as `analysis/`. They are now committed alongside the
> book, and the chapters point at `analysis/slurm/`.
>
> Absolute paths of the form `.../q2-hdstats-recompute/...` in the text below
> now correspond to `analysis/...` in this repository. Read as a record of what
> was found, not as a description of the current tree.

---

## CONFIRMED (13)

**[docs-01] `docs/_config.yml` + `requirements.txt` vs `.github/workflows/ci.yml` — published book build will fail (HIGH)**
Verified: `git diff docs/_config.yml` adds `sphinxext.rediraffe` under `sphinx: extra_extensions:` plus 19 `rediraffe_redirects:` entries; `requirements.txt` gains `sphinxext-rediraffe`. `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/.github/workflows/ci.yml` is untouched by this change set (`git status --short` lists no `.github` entry) and its only install step is `mamba install -n conda-env jupyter-book`. `sphinxext-rediraffe` is not a jupyter-book dependency (`sphinxcontrib-bibtex`, the pre-existing extension, is — which is why nobody noticed). `.readthedocs.yml` does `python: install: - requirements: requirements.txt`, so RTD builds fine. I also confirmed rediraffe genuinely ran in the local build — `docs/_build/html/chapters/02_graphical_models/02_sgl.html` is a redirect stub pointing at `../02_lowdim_gglasso/02_sgl.html` — so the config is correct and the failure is CI-only: Sphinx will raise `ExtensionError: Could not import extension sphinxext.rediraffe` on push to `main`.
Fix: add `sphinxext-rediraffe` to the `mamba install` line in `.github/workflows/ci.yml` (or replace it with `pip install -r requirements.txt`).

**[docs-02] `docs/chapters/04_highdim_atacama/06_interpretation.md:197` — consumes an artifact this change set deleted the producer of (MEDIUM)**
Verified: `git diff -M docs/chapters/04_highdim_atacama/02_model_selection.md` shows `-    --o-solution atacama-top-300-sgl-lambda0.8.qza` / `+    --o-solution atacama-top-300-sgl-linear-path.qza` (line 36). `grep -rn "atacama-top-300-sgl" docs/chapters` returns only those two hits: produced as `…-linear-path.qza`, consumed at 06_interpretation.md:197 as `…-lambda0.8.qza`. The `qiime tools export` there fails with a missing-file error, so the SGL-vs-SLR edge comparison cannot be run.
Fix: change line 197 to `atacama-top-300-sgl-linear-path.qza`. This is semantically correct — 03_slr_ranks.md itself states a model-selection run returns one solution (the eBIC-best grid point), i.e. λ = 0.8.

**[docs-03] `docs/chapters/01_installation/04_verify.md:28-33` — the "must list" help output does not match reality (MEDIUM)**
Verified by running it in the env (not inferred). Actual `qiime gglasso --help`:
```
  build-groups          build-groups
  calculate-covariance  calculate_covariance
  pca                   Principal component analysis (PCA)
  solve-problem         solve_problem
  summarize             Summary table
  transform-features    transform-features
```
The page claims `Group array` / `Calculate covariance matrix` / `Solve graphical lasso problem` / `Transform features` for four of six. Cause confirmed at `q2cli/commands.py:331` (`short_help=action['name']`) against `q2_gglasso/plugin_setup.py` lines 99/123/152/258 (`name="transform-features"`, `"build-groups"`, `"calculate_covariance"`, `"solve_problem"`). A reader following this chapter concludes their install is broken. Note the adjacent claims *are* right: `qiime classo --help` does list exactly the eight named actions, and `classify` is no longer mislabelled `regress`.
Fix: paste the real output, or drop the description column.

**[docs-04] `02_lowdim_gglasso/07_pca.md` (~158-169) and `04_highdim_atacama/04_latent_pca.md` (~128-133) — document a heuristic the sibling change set removed (MEDIUM)**
Verified: `git diff q2_gglasso/utils.py` in `/home/itg/oleg.vlasovets/slr_example/q2-gglasso` deletes the `p, n = counts.shape; if p >= n: sum(axis=0) else sum(axis=1)` guess and replaces it with an unconditional `depth = counts.sum(axis=1)` (plus a zero-spread guard); `q2_gglasso/_pca/_visualizer.py` now orients the table against `L.shape[0]` (`if df.shape[1] != p: … df = df.T`) *before* the `get_seq_depth(df)` call at line 417. Both doc passages still say the helper "infers orientation from the table's shape" and instruct Tier 1 readers to always pass `--p-color-by` because the fallback allegedly returns an ASV-indexed column. After the plugin change that is false at any p/n ratio.
Fix: rewrite both notes to state that `get_seq_depth` sums over features unconditionally and `pca` orients the table first. (Cross-repo: only a defect once the q2-gglasso change lands, which is the same batch.)

**[docs-05] `docs/chapters/04_highdim_atacama/05_classo_cv.md` — only Tier 2 page not brought into the pending/verified convention, and its inputs don't exist (MEDIUM)**
Verified: `grep -in "pending\|verif"` over the 113-line file returns nothing. It states an ASV-only CV-R² table (lines ~36-52) and a joint/filtered table (~97-113) and names *Pseudarthrobacter* as fact. `04_highdim_atacama/00_index.md:106-115` explicitly lists "the cross-validated $R^2$ values and selected features ([Log-Contrast Models at Scale](05_classo_cv.md))" under **Still pending**, and `06_interpretation.md` marks the same Pseudarthrobacter reading pending. `git diff -M` on this file shows only two link-path rewrites, so the restructure moved it without applying the convention. Separately, `grep -rn` shows `atacama-top-300-clr-design.qza` and `atacama-classo-outcomes-mean-imputed.tsv` appear *only* on this page — no chapter produces them and `docs/_data/manifest.tsv` does not list them, so none of the three commands is runnable as written.
Fix: add the standard "pending verification against QIIME 2 2026.7" admonition, and either add the two artifacts to `manifest.tsv` or document the commands that build them.

**[docs-06] `docs/chapters/04_highdim_atacama/02_model_selection.md:64` — 1403 should be 1405 (MEDIUM)**
Verified against the actual solver output, not just the committed table. `docs/_data/atacama-lambda-path.tsv` last row is `0.3 / 16165.13 / 1405`; the upstream generated file `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/results/tables/lambda-path.tsv` gives sparsity `0.0313266` at λ = 0.3, and `0.0313266 × C(300,2) = 0.0313266 × 44850 = 1405.0`. (Same arithmetic gives 216 at λ = 0.8, matching the prose and Gate C1.) So **1405 is right and the prose is wrong**, and it is wrong on the one page whose note at lines 54-60 asserts "the prose above cannot drift away from the numbers below."
Fix: change 1403 → 1405 at line 64.

**[docs-07] `docs/chapters/04_highdim_atacama/03_slr_ranks.md` + `00_index.md` — stale "pending" on a map that has been recomputed (LOW)**
Verified: line 17-18 "the values quoted at the end are pending recompute" and lines ~114-117 "the achieved ranks are **pending verification against QIIME 2 2026.7**" sit in the same file as line ~226 "These are measured values, read out of the fitted solutions at the Gate-C1-selected $\lambda = 0.8$" plus a `{csv-table}` from `docs/_data/atacama-mu-rank-map.tsv` (15→2/202/162, 10→5/158/124, 7.5→10/110/92) and a note at ~250 claiming it is generated. Those are exactly the verified numbers, and the recompute script's own gate (`03_mu_rank_map.sh`, lines 98-110) checks them against the asserted table. `00_index.md:106-109` still lists the μ₁→rank map under **Still pending**.
Fix: drop the two stale notes in 03_slr_ranks.md and move the μ₁→rank map into the "Recomputed and confirmed" block of 00_index.md.

**[docs-13] `02_model_selection.md:54-60` and `03_slr_ranks.md:249-257` — the named generator scripts are not in any of the three repos, and one does not write the file it is said to write (LOW)**
Partly re-verified against the reviewer's evidence, which was incomplete but the finding stands and is *worse* than stated. `ls slurm` in the docs repo → no such directory, not gitignored, and absent from q2-gglasso and q2-classo. The scripts do exist, but at `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/{01_lambda_path.sh,03_mu_rank_map.sh}` — a directory that is **not a git repository at all** (`git rev-parse` → "not a git repository"), so nothing is committed anywhere. Worse, `01_lambda_path.sh:101,148` writes `$ROOT/results/tables/lambda-path.tsv` with columns `lambda1 / sparsity / ebic`, whereas the committed `docs/_data/atacama-lambda-path.tsv` has `lambda / eBIC (gamma=0.3) / edges` — a different name, path and schema, i.e. a hand step sits between them. That hand step is precisely where docs-06's 1403/1405 drift came from, which disproves the "generated, not transcribed … cannot drift" claim empirically. (`03_mu_rank_map.sh:94` does emit the exact four committed columns.)
Fix: do **not** just "commit the scripts under `slurm/`" as proposed — they belong to the separate recompute tree. Either vendor them with the transformation step included, or reword both notes to name where the generator actually lives and drop the "cannot drift" guarantee.

**[docs-08] Tier 3 — three mutually contradictory statements about whether the metadata was inspected (LOW)**
Verified all five passages. `05_metagenomics/00_index.md:96-102` "**Nothing in this tier has been run against QIIME 2 2026.7** … none shows captured output"; `01_gut_to_soil/00_index.md:35-42` and `01_data.md:148-153` say column names are not restated because unverified. But `03_regression.md:31-54` says "# These are real column names, read off the downloaded sample-metadata.tsv", sets `OUTCOME_COLUMN="Composting Time Point"` / `CLASS_COLUMN="SampleType"` / `COVARIATE_COLUMN="Compost pH"`, and quotes 1267/1209 row coverage and 15 `SampleType` levels with sizes 799/453. `01_data.md` also contradicts itself in-file: one note says nothing is verified, another says "The per-sample depths quoted above are verified against `asv-table-ms2.qza`" (3-1218, median 261). Meanwhile `02_network.md:219-231` still uses `sample-type`/`gut`/`soil` placeholders that cannot exist in a 15-level `SampleType`.
Fix: pick one story — propagate the verified names into `02_network.md`'s `GROUP_COLUMN`/`GROUP_A`/`GROUP_B` and relax the three "not verified" admonitions, or strip the verified block from 03_regression.md.

**[docs-09] `docs/_toc.yml:9-10` vs `docs/chapters/00_getting_started/02_datasets.md` (LOW)**
Verified: ToC entry is `file: chapters/00_getting_started/02_datasets` / `title: The Four Example Datasets`, but the page (43 lines) opens `# Atacama Soil Microbiome` and its only dataset section lists "50 samples … 13 microbial taxa (ASVs)". `03_download_data.md:9` sends readers there for "The datasets themselves — what they contain and why each one is in the book", which the page does not deliver. One correction to the reviewer's evidence: the second inbound link (`03_download_data.md:108`) points there for the `ASV-1…ASV-13` ↔ MD5 mapping, which the page *does* contain — that link is fine.
Fix: either write the four-dataset overview into the page and update its H1, or revert the ToC title and the line-9 link text.

**[docs-10] Installation chapters disagree on the rebrand release (LOW)**
Verified: `01_prerequisites.md:14`, `02_q2_gglasso.md:8` and `03_q2_classo.md:10` all attribute the `amplicon`→`qiime2` and `qiime2`→`rachis` renames to **2026.4**; `04_verify.md:13` says "Since the **2026.1** rebrand the framework package is called `rachis`". Same event, two dates — one is wrong. I could not determine which from this environment (`qiime info` reports only the end state, `rachis release: 2026.7`), and every install URL in the book pins 2026.7 regardless, so nothing is broken today; it is a factual inconsistency the author must settle.
Fix: pick one release number and use it in all four files.

**[docs-11] `docs/chapters/90_reference/04_troubleshooting.md:133` — states a rule the code does not implement (LOW)**
Verified against source. `q2_classo/_func.py:30-38`:
```python
if cv__nlam is None:
    return cv_nlam

if cv_nlam != 100 and cv_nlam != cv__nlam:
    raise ValueError(...)
```
`--p-cv-nlam 100 --p-cv--nlam 50` → `cv_nlam == 100`, so the guard short-circuits, no raise, and the function returns `50`. The page's "Passing both with *different* values is an error." is therefore wrong, and it contradicts `90_reference/03_classo_parameters.md:225` and `03_lowdim_classo/05_advanced/02_model_selection.md:162-168`, which both document the hole correctly ("does **not** raise — it silently uses 50"). Note the same overstatement is in the function's own docstring at `_func.py:27-28` ("Passing both is an error rather than a silent precedence rule").
Fix: add the "100 is treated as unset" caveat to the troubleshooting entry and to the docstring.

**[docs-12] `docs/chapters/90_reference/03_classo_parameters.md:343-344` — advice contradicts every worked command and the same page's own table (LOW)**
Verified: registered default is `rho: float = 0.0` at `q2_classo/_func.py:410`, wired to `problem.formulation.rho_classification` at :438 (c-lasso's own default for that field is `-1.`, per its docstring: "it has to be strictly smaller then 1", so `0.0` is a legal knee). The page says "you must set it explicitly for `huber` to do anything" — yet every worked `classify` command passes exactly the registered default: `03_lowdim_classo/05_advanced/01_concomitant_formulation.md:130` and `:155` (`--p-rho 0.0`) and `05_metagenomics/01_gut_to_soil/03_regression.md:265` (`--p-huber --p-rho 0.0`). By the page's own logic Huber does nothing in all of them. Line 330 of the same file, and `01_concomitant_formulation.md:114`, describe `0.0` as a working default with no such caveat.
Fix: drop the "must set it explicitly" sentence (keeping the true part — that `classify`'s registered default is `0.0`, not the `1.345` its `--help` text claims via `_dict.py`), or switch the examples to a non-default rho.

**[docs-14] `docs/chapters/02_lowdim_gglasso/05_lambda_paths.md:101-105` — leftover uncertainty note contradicting the 49→50 correction (LOW)**
Verified: the note says "the exact sample count of the 13-ASV subset is being confirmed, so check it against your own `qiime feature-table summarize` output before trusting a published number." The same change set settles it: `03_download_data.md:130-136` says "The tier 1 table is **13 features × 50 samples**, read directly from the artifact… the commands were right and the prose was wrong", and `git diff` shows `-49 samples` → `+50 samples` in `02_datasets.md` and `-$N = 49$` → `+$N = 50$` in `02_lowdim_gglasso/01_data_preparation.md`.
Fix: replace the note with the settled statement (50 samples), keeping the generic "check it on your own table" advice if wanted.

## REFUTED (0)

None of the 14 proposals collapsed under verification. The only substantive correction is to **docs-13**, whose evidence was wrong in the reviewer's favour: the two generator scripts *do* exist (in the untracked, non-git `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/`), but `01_lambda_path.sh` writes `results/tables/lambda-path.tsv` with columns `lambda1/sparsity/ebic`, not `docs/_data/atacama-lambda-path.tsv` with `lambda/eBIC (gamma=0.3)/edges` — so the reviewer's proposed fix ("commit the two scripts under `slurm/`") would not make the documented claim true. Minor evidence corrections also apply to **docs-09** (only one of the two inbound links is misleading) and **docs-10** (the wrong release number does not currently affect any distribution URL, all of which pin 2026.7).