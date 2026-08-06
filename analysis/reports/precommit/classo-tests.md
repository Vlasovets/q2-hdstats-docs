Verification complete. Environment used: ruff 0.16.1 + q2lint (master) installed into a throwaway venv at `/localscratch/oleg.vlasovets/q2dbg/venv/lintenv` (left in place for re-checks; node-local scratch only).

# CONFIRMED

### 1. [ci-ruff-gate-red] — HIGH — `/home/itg/oleg.vlasovets/slr_example/q2-classo/.github/workflows/ci.yml:89`
**Defect stands.** The new `lint` job is red on the exact tree that introduces it, and the in-file comment on lines 86–87 ("so the gate is green and meaningful from day one") is factually false.

Verification: installed ruff 0.16.1 and ran the workflow's literal command in the repo root:
`ruff check --select=F .` → `Found 68 errors`, `EXIT=1`.
Reviewer's breakdown reproduced exactly: 57 F401, 5 F405, 5 F841, 1 F403. Confirmed there is no `[tool.ruff]` section in `pyproject.toml` and no `ruff.toml`/`.ruff.toml`, so no config softens the gate. Confirmed none of the 68 are in the three new test files and none are in the gitignored `_version.py` — a clean CI checkout produces the same 68. Anchor violations spot-checked in source: `q2_classo/_func.py:5` is `from classo import *` (F403, unfixable by `--fix`), `_func.py:551` `classification` and `:552` `dico_ms` are dead (F841), `q2_classo/tests/test_name_formulation.py:6` imports `os` which is never used.

Minimal fix: make the gate honest now — either `continue-on-error: true` on the ruff step, or narrow to the rules that actually pass (`--select=F --ignore=F401,F403,F405,F841`), and rewrite the comment to state the real starting point and the tightening path.

Note (corroboration, out of this change set): the sibling `q2-gglasso` workflow's `ruff check .` is *also* red — 286 errors. So this is not a copy/paste slip; the "green from day one" assumption is wrong in both repos.

### 2. [ci-q2lint-red] — MEDIUM — `.github/workflows/ci.yml:92`
**Defect stands.** Ran the real thing: installed `q2lint` from the same URL the workflow uses and executed it in the repo root:

```
Invalid header: q2_classo/_summarize/__init__.py (header missing/too short)
Invalid header: q2_classo/tests/__init__.py (header missing/too short)
EXIT=1
```

Exactly the two files named, no others. Confirmed in source: `q2_classo/_summarize/__init__.py` is a single line (`from ._visualizer import summarize`) and `q2_classo/tests/__init__.py` is empty. This fails the same `lint` job independently of the ruff step, so fixing only defect 1 leaves CI red.

Minimal fix: add the 7-line QIIME 2 header to those two files, or drop the `q2lint` step.

### 3. [predict-path-xgraph-unused] — MEDIUM — `q2_classo/_summarize/_visualizer.py:743, 767`
**Defect stands, and it is stronger than reported: the mislabeling is unconditional, not conditional.**

Verified in source: `plot_predict_path` computes `xGraph` at 740/743, then line 767 plots `graph_objects.Scatter(x=lambdas, y=error, ...)` while line 770 sets the axis title from `textlam`. ruff independently confirms `F841 Local variable 'xGraph' is assigned to but never used --> q2_classo/_summarize/_visualizer.py:743` — the only F841 among the plot helpers. So the `xGrpah` → `xGraph` rename at 743 is inert here.

The reviewer said this bites "for a default c-lasso PATH". It is actually always: `_func.py:303` (`regress`) and `_func.py:452` (`classify`) both hardcode `param.logscale = True` for PATH parameters — there is no `path_logscale` plugin parameter. `dico_path["logscale"]` therefore reads back `True` in every run, so the caller at line 185 always passes `logscale=True`, and the shipped `predict-path.html` always shows raw lambda values under an axis titled `- log10 lambda / lambdamax`.

This is a *new* user-visible defect from this change set: line 775 replaces the commented-out `offline.plot` with a live `fig.write_html`, so the pane goes from blank to rendered-and-mislabeled.

Minimal fix: line 767 → `graph_objects.Scatter(x=xGraph, y=error, name="L2 error over test set")`, matching `plot_path`'s convention.

(One reviewer sub-claim is overstated and should not be repeated: they call the same rename in `plot_path` "load-bearing" because `logscale=False` would previously have raised NameError. Since `param.logscale` is hardcoded `True`, that branch was unreachable from the plugin — the `plot_path` rename is a harmless dead-code cleanup, not a live bug fix.)

### 4. [summarize-assets-name-matching-only] — LOW — `q2_classo/tests/test_summarize_assets.py:1-12`
**Confirmed, but narrowly: it is a docstring accuracy problem, not a broken test.** Both replays reproduce exactly.

Replayed `_written_filenames()`/`_iframe_targets()` against `git show HEAD:q2_classo/_summarize/_visualizer.py` (every `offline.plot` commented out, all 16 panes blank), with both current and HEAD templates: `missing = ['stabsel-tree.html']`, `orphans = ['StabSel-tree.html']` in both cases. So of the 16 blank panes, only the one case mismatch would have been caught.

Mutation check: deleting `plot_beta`'s `fig.write_html(...)` (line 543) from the current source gives `missing = []`, `orphans = []`, `'# offline.plot(' absent`, `'write_html(' present` — all four tests stay green while lam-beta/lam-refit/cv-refit/stabsel-refit go blank.

The docstring's "each iframe the templates reference has something that writes it" therefore overstates what is checked (template ↔ source-literal filename consistency).

Minimal fix: soften that one sentence to describe filename consistency. **Do not** adopt the proposed ast rewrite — the four tests are sound and non-vacuous (they did catch the real `StabSel-tree.html`/`stabsel-tree.html` bug), and adding an ast walker for a low-value guarantee is scope the change set does not need.

### 5. [cv-nlam-sentinel-100] — LOW — `q2_classo/_func.py:18-45`
**Confirmed on the docstring contradiction and the duplicated literal; the "tests are wrong" framing is not.**

Verified in source: the guard is `if cv_nlam != 100 and cv_nlam != cv__nlam:` (`_func.py:33`), while the docstring three lines above says "Passing both is an error rather than a silent precedence rule". Those disagree whenever `cv_nlam == 100`: `_resolve_cv_nlam(100, 120)` returns 120 with only a DeprecationWarning. The docstring never mentions the "value equals the default ⇒ treat as unset" carve-out.

The literal `100` is the de facto sentinel and appears in three unlinked places: the guard (`:33`), `regress` (`:239`), `classify` (`:384`). Change the registered default and every legacy `--p-cv--nlam N` invocation starts raising ValueError. Verified `q2_classo/tests/test_plugin_setup.py:47-58` only asserts both names are *registered* — it never reads defaults — so nothing in the suite would catch that.

Correction to the reviewer's framing: `test_cv_nlam_alias.py:21-27` asserting `_resolve_cv_nlam(100, 120) == 120` is *correct* — that is precisely the legacy-only call path QIIME 2 produces. The test is right; only the docstring and the unlinked literal are the problem.

Minimal fix: introduce `_CV_NLAM_DEFAULT = 100` at module level and use it in all three sites, and amend the docstring to state the actual rule. **Do not** adopt the proposed `cv_nlam: int = None` signature change — with neither spelling given, `param.Nlam` would then be set to `None` (`_func.py:318`, `:467`) and break the ordinary CV path, and `--help` would stop advertising 100.

# REFUTED

None. All five survived; two (4 and 5) survived only in a narrower form than proposed, and their proposed fixes should be replaced with the minimal ones above.