# Adversarial verification — group `download-and-coverage`

Both files pass the flag check: **no nonexistent `--p-`/`--i-`/`--o-`/`--m-` flags in either file**, and no `qiime classo classify` usage of `--p-concomitant` or `--p-do-yshift`. Both are MyST-valid. The defects below are factual and structural.

---

## `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/docs/chapters/00_getting_started/03_download_data.md`

**D1 — line 105 — wrong fact about the tier 1 artifact (highest severity).**
The table says `atacama-counts.qza` is "13 ASVs with human-readable labels `ASV-1` … `ASV-13`". I extracted the biom from `/home/itg/oleg.vlasovets/slr_example/q2-classo/data/atacama-counts.qza` (byte-identical to the copy in the tutorial zip and in `q2-gglasso/data/`): shape `(13, 50)`, and `/observation/ids` are 32-char MD5 hashes — `409faa5f5353e543bf6d99125c7c0e83`, `a7b877ae6d2f079a15b6b192a4425620`, … The `ASV-n` labels are *assigned* by `transform-features --p-keep-original-id False` (see `q2_gglasso/plugin_setup.py:90-94`), they are not a property of the download. This also contradicts the sibling chapter this page links to: `chapters/00_getting_started/02_datasets.md:11-23` prints the hash→`ASV-n` mapping table explicitly.
Fix: "13 ASVs; feature IDs are MD5 hashes. The `ASV-1`…`ASV-13` labels used throughout the book are produced by `transform-features --p-keep-original-id False`; the mapping is in [The Four Example Datasets](02_datasets.md)."
Knock-on: this same false premise is the stated justification for a "Justified" row in the coverage matrix (`01_command_coverage.md:354`, "The tier 1 table already carries short `ASV-n` labels").

**D2 — line 207 — transposed dimensions.**
"`asv-table.qza` | The other tutorial table, 104 features × 1069 samples". I downloaded it from the URL this page gives (line 198); the biom `shape` attribute is `1069, 104`, i.e. **1069 features × 104 samples**. Also self-evidently wrong as written: 104 features cannot be the superset of a 335-feature subsample. Same error is propagated to `chapters/05_metagenomics/01_gut_to_soil/01_data.md:51`.
Fix: "1069 features × 104 samples" in both places.

**D3 — lines 55, 70-79, 113-114, 161-162, 213-214 — the manifest does not exist.**
The page calls `docs/_data/manifest.tsv` "the machine-checkable index", "authoritative", and builds `tier_checklist` plus all three verify blocks on it. `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/docs/_data/` **does not exist** (no such directory), and nothing in this workflow is listed as creating it. So `tier_checklist` writes an empty file and every `sha256sum -c` errors — including tier 3, which the page promises is verifiable today. This is not covered by the `ZENODO_DOI_PENDING` warning, which only excuses the URL and the tier-1/2 checksum *values*.
Fix: create `docs/_data/manifest.tsv` with real tier-3 `bytes`/`sha256` rows (tier 3 content-lengths as served today: `asv-table-ms2.qza` 47343, `asv-table.qza` 65124, `taxonomy.qza` 79579, `asv-seqs.qza` 75123, `sample-metadata.tsv` 354759), or add a warning that the manifest is not yet written, in the same style as the DOI warning.

**D4 — line 3 — false exclusivity claim.**
"This is the **only** page in the book that carries download URLs." It is not:
- `chapters/00_getting_started/02_datasets.md:43` — `https://github.com/Vlasovets/q2-gglasso/blob/main/data/atacama-counts.qza` and `https://data.qiime2.org/2026.7/tutorials/atacama-soils/sample_metadata.tsv`
- `chapters/05_metagenomics/01_gut_to_soil/01_data.md:36-41` — repeats the entire tier-3 `BASE=…readthedocs.io/…` curl block.
Fix: strip the URLs from those two pages, or soften to "the canonical place for download URLs".

**D5 — line 84 — understates the failure mode.**
"the tier 1 and 2 verification steps will not do anything useful yet." With `ZENODO_DOI_PENDING` in the `sha256` column, `sha256sum -c` does not no-op — it exits non-zero with `no properly formatted checksum lines found`. A reader following the page will see a hard error on a documented command.
Fix: say the tier-1/2 verify commands will *fail* until the record is minted.

**D6 — lines 231-234 — citation year unsupported, one figure unverified.**
I resolved Zenodo record 15390940: it is real — *"Upcycling Human Excrement: The Gut Microbiome to Soil Microbiome Axis (supporting data)"*, creators Caporaso J G and Meilander J, DOI `10.5281/zenodo.15390940`, **publication date 2025-05-12**. The "1,660 samples" figure is corroborated (the record's `sample-metadata.tsv` has exactly 1,660 sample rows). But "**(2024)**" is not supported by the record, there is no `Meilander` entry in `docs/references.bib` (grep returns nothing), and the "roughly 30,000 ASVs" figure is not corroborated by anything I could check.
Fix: cite `10.5281/zenodo.15390940` (2025), add a bib entry, and either source or drop "30,000 ASVs". (Same numbers repeat at `05_metagenomics/01_gut_to_soil/01_data.md:213-223`.)

**D7 — line 226 — observed number stated loosely.**
"the upstream alpha-rarefaction curve tops out around a depth of 260, and per-sample frequencies sit in the low hundreds." I computed the per-sample totals from `asv-table-ms2.qza`: **min 3, median 261, max 1218**. 260 is the *median*, not where anything "tops out"; the max is nearly 5× that. The qualitative conclusion (shallow, noisy) survives, the number does not.
Fix: "per-sample depth ranges from 3 to 1218, median 261".

**D8 — line 206 — misleading description.**
"`sample-metadata.tsv` | Sample metadata for the composting time series." The file served at that URL is the **full-study** metadata: 1,660 sample rows / 355 KB, not the 99 samples in `asv-table-ms2.qza`. Readers will hit QIIME 2 warnings about metadata rows with no matching sample.
Fix: "Full-study sample metadata (1,660 rows); only 99 of them match `asv-table-ms2.qza`."

**D9 — lines 99, 107 — unverifiable provenance (report, do not necessarily change).**
`selected-atacama-sample-metadata.tsv` is real (`/home/itg/oleg.vlasovets/slr_example/q2-gglasso/data/selected-atacama-sample-metadata.tsv`, 2,599 bytes) and is used by 6 later chapters, but it is **absent from both Atacama bundles at the repo root**, so nothing establishes it will be on the Zenodo record. This is exactly what D3's manifest is supposed to settle.

**Verified clean on this page** (so it is not re-flagged later): the 5 covariate names on line 108 match the header of `atacama-selected-covariates-veg.tsv` exactly; tier-2 `sample-metadata.tsv` does carry `transect-name` ∈ {Baquedano (32), Yungay (43)} and `vegetation` ∈ {no (40), yes (35)}; the `--i-taxonomy`-required-but-unread claim (lines 121-125) is correct (`_func.py:32-41`, `taxonomy` appears only in the signature and docstring); all five tier-3 URLs return HTTP 200 today; the 10% subsample claim is corroborated by provenance (`demux subsample_paired`, `fraction: 0.1`); 335 × 99 is correct; all 6 relative links resolve; fences balanced (20 pairs); one H1 (the `#` at line 49 is inside a bash fence). The hedged note at 127-134 about the tier-1 sample count is appropriate — the answer is **50 samples** (biom shape `13, 50`) against a 75-row metadata TSV.

---

## `/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/docs/chapters/90_reference/01_command_coverage.md`

**C1 — lines 13-16 vs. lines 88-313 — the completeness claim is false for outputs (highest severity on this page).**
"Every one appears below, and so does every parameter group listed in [02_gglasso_parameters.md] and [03_classo_parameters.md]" — and its own CI check #2 (lines 381-383) demands it. But **`--o-group-array` (line 108) is the only output flag anywhere on the page.** Absent entirely, while documented in the two reference pages this sentence points at: `--o-transformed-table`, `--o-covariance-matrix`, `--o-solution`, `--i-covariance-matrix` (gglasso); `--o-x`, `--o-c`, `--o-aweights`, `--o-new-features`, `--o-new-c`, `--o-new-w`, `--o-result`, `--o-predictions` (classo). Note `solve-problem` is the only action whose *input* is also missing, while `pca`, `predict`, `regress` etc. all have explicit input rows — so the omission is inconsistent, not a deliberate policy.
Fix: add an outputs row per action (mirroring the `build-groups` row that already names `--o-group-array`), add `--i-covariance-matrix` to the `solve-problem` table at line 133, or narrow the sentence to "every input and parameter" and correspondingly narrow CI check #2.

**C2 — lines 139, 159-164 — "always raises `NotImplementedError`" is the wrong exception in one branch.**
`q2_gglasso/_func.py:572-585`:
```python
if rank is not None:
    if not latent:
        raise ValueError("The 'rank' parameter is only meaningful for the sparse + low-rank model; set latent=True.")
    if not _gglasso_supports_rank():
        raise NotImplementedError(...)
```
`--p-rank` without `--p-latent` raises `ValueError`, not `NotImplementedError`. A reader who hits the `ValueError` will not match it against the documented trap.
Fix: "always raises — `ValueError` if `--p-latent` is not set, `NotImplementedError` otherwise, on every released GGLasso ≤ 0.3.0."

**C3 — lines 359-363, 396-402 — unverified negative claim about the environment.**
"there is no QIIME 2 2026.7 environment to capture from" is used to excuse the missing `--help` capture. A conda env named `q2-2026.7-slr` exists at `/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr`. Either the claim is stale or that env lacks the plugins — I did not activate it. As written it is an unverified negative load-bearing for the page's central caveat.
Fix: check the env, then either drop the clause or restate as "the `--help` capture step is not wired into the build".

**Verified clean on this page:** action counts (6 gglasso / 8 classo) correct; all 20 `solve-problem` parameters present across exactly 8 rows, matching `glasso_parameters`; the `classify` "does not have `--p-concomitant` / `--p-do-yshift`" warning (lines 286, 289-296) is correct against `classify_parameters`, and `_func.py:425` hard-sets `problem.formulation.concomitant = False`; `rho` default **0.0** on classify vs **1.345** on regress confirmed (`_func.py:400` vs `_func.py:266`); PATH/CV/StabSel/LAMfixed all default `True` on both actions, so "all on by default" is correct; the two cross-cutting claims at lines 148-156 verified against `utils.py:243-253` (silent 1e-3/1 substitution when only one bound is set, no warning) and `utils.py:310-336` (`warnings.warn` on full default, and for `latent` all three of `lambda1`/`lambda2`/`mu1` must be singleton); `--p-n-samples` is genuinely the only `solve-problem` parameter with no default (`_func.py:482-504`); one-hot label format `<name> = <value>` confirmed at `q2_classo/_func.py:197`; `randomy.tsv` CWD side effect confirmed at `q2_classo/_func.py:97-99`; pca traps confirmed (`_pca/_visualizer.py:359` reads `solution/lowrank_`, line 361 dereferences `sample_metadata` unconditionally); the four non-plugin flags at lines 325/327 belong to `sample-classifier split-table` and `feature-table filter-samples` and are labelled as such. All 33 chapter-key links resolve to real files and all appear in `_toc.yml`. One H1, 9 balanced directive fences.