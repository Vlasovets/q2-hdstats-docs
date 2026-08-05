# Decisions needed

Six items surfaced during the 2026.7 migration that change **behaviour, artifact
semantics, or a release commitment** rather than prose. Each is deliberately left
unfixed. Evidence is recorded so the decision can be made without re-deriving it.

Status of everything else: Gate A1 and Gate C1 both PASS, 57 plugin tests pass,
the book builds with zero warnings. See `project_q2_hdstats_migration` in the
Claude memory index for the full picture.

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

**New finding, still open.** In the same run, the shipped
`atacama-top-300-correlation.qza` differs from **both** old and new output by the
same 1.147, so it was not produced by `transform-features -> calculate-covariance`
at the documented parameters (clr, pseudo-count 1, no-keep-original-id). Some
other transform setting — mclr, a different pseudo-count, appended metadata
columns — produced it. Gate C1 consumes that matrix directly and is unaffected,
but **the tier-2 chapter implies the chain regenerates it, and that claim is not
currently true.** Either recover the original parameters from C. Müller or
regenerate the matrix and re-run Gate C1 against the new one.

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

## 4. Is `selected-atacama-sample-metadata.tsv` going on the Zenodo record?

Six chapters depend on it. It exists only at
`q2-gglasso/data/selected-atacama-sample-metadata.tsv` (2599 bytes) — not in
q2-classo's data dir, and not in either Atacama bundle at the docs repo root.
`docs/_data/manifest.tsv` lists it as a tier-1 row with `ZENODO_DOI_PENDING`.

This is a release decision: confirm it is deposited, then fill in bytes + sha256.

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
