## CONFIRMED

**[classo-01-numpy2-classo-incompatible] — /home/itg/oleg.vlasovets/slr_example/q2-classo/pyproject.toml (+ requirements.txt, environment-files/q2-classo-qiime2-2026.7.yml, .github/workflows/ci.yml) — HIGH. Survived refutation.**

The change set declares `c-lasso==1.0.11` alongside `numpy>=2.0,<3` in all three dependency files, and adds a CI gate that *requires* numpy 2.x — a combination in which `qiime classo regress` cannot run at all.

Verification I did myself (not taken from the reviewer):
- `numpy 2.4.2`, `hasattr(np, 'infty') == False`; `classo 1.0.11` at `.../site-packages/classo/__init__.py`.
- `np.infty` present at `classo/solve_R1.py:212`, `solve_R2.py:239` and `:293`, `solve_R3.py:205`, `solve_R4.py:211`.
- Ran a 40×12 synthetic problem in the env. Results:
  - `PATH+LAMfixed`, concomitant=True → `AttributeError: np.infty was removed in the NumPy 2.0 release`
  - `PATH+LAMfixed`, concomitant=False → same
  - `CV` only → same; `StabSel` only → same
  - classification (`formulation.classification=True`, PATH+LAMfixed) → **OK**. So `classify` genuinely survives; only `regress` is dead.
- No shim exists: no `np.infty` reference anywhere in q2-classo, no `sitecustomize.py`, no patching in `classo/__init__.py`.
- No test can catch it: `q2_classo/tests/` contains only `test_cv_nlam_alias.py`, `test_name_formulation.py`, `test_plugin_setup.py`, `test_summarize_assets.py`; `grep` for `.solve()` / `regress(` / `classo_problem(` in tests returns nothing. `.github/workflows/ci.yml` only runs `qiime classo regress --help` plus `pytest`, and its version-drift step (`{"numpy": "2", ...}`) actively pins the broken configuration.
- Documented commands would error: `qiime classo regress ...` appears as an executable block in the docs book at `docs/chapters/03_lowdim_classo/05_advanced/02_model_selection.md` (lines 50, 94, 173, 237, 318), `docs/chapters/04_highdim_atacama/05_classo_cv.md` (17, 80), `docs/chapters/05_metagenomics/01_gut_to_soil/03_regression.md:210`.

Minimal fix: patch the four `np.infty` → `np.inf` sites in a pinned fork/vendored copy of c-lasso (upstream 1.0.11 is the latest release and has no fix), and add one real solve test (n≈40, d≈12, PATH+LAMfixed only) so CI can detect this class of breakage. Reverting to numpy<2 is not viable against the 2026.7 distribution.

---

**[classo-03-ruff-gate-not-green] — /home/itg/oleg.vlasovets/slr_example/q2-classo/.github/workflows/ci.yml:86-89 — MEDIUM. Survived refutation, verified by running ruff.**

I installed ruff 0.16.1 into scratch (same version the CI would get) and ran exactly the CI command at the repo root:

```
ruff check --select=F --no-cache .   ->  Found 68 errors.   EXIT=1
```

There is no `[tool.ruff]` in `pyproject.toml` and no `ruff.toml`/`.ruff.toml`/`setup.cfg`/`tox.ini`, so nothing relaxes it. Representative hits: `F403` at `q2_classo/_func.py:5` (`from classo import *`) plus 4× `F405`; `F401` at `q2_classo/plugin_setup.py:15` (`Numeric`) and ~45 other unused imports across `__init__.py`, `_dict.py`, `_formats.py`, `_summarize/_visualizer.py`, `tests/test_name_formulation.py:6`; `F841` at `_func.py:551,552,607` and `_visualizer.py:67,743`. The inline comment ("the gate is green and meaningful from day one") is false on the first push.

Minimal fix: either add `[tool.ruff.lint] ignore = ["F403","F405","F401","F841"]` (or the narrower `--select` that actually passes), or set `continue-on-error: true` on the lint step, and correct the comment. 35 of the 68 are `--fix`-able if you prefer to clear the debt instead.

---

**[classo-05-predict-path-axis-still-mislabeled] — /home/itg/oleg.vlasovets/slr_example/q2-classo/q2_classo/_summarize/_visualizer.py:740-772 — LOW/MEDIUM. Survived refutation, and it is worse than the reviewer stated.**

`plot_predict_path` computes `xGraph` at :740 (logscale) / :743 (linear) but builds the trace with `graph_objects.Scatter(x=lambdas, ...)` at :767, while `fig.update_xaxes(title_text=textlam)` at :770 uses the logscale-derived title. Ruff independently confirms `F841 Local variable 'xGraph' is assigned to but never used --> _visualizer.py:743`. The sibling `plot_path` consumes its `xGraph` at :505, so only that rename had an effect.

The logscale branch is the only one reachable: `_func.py:303` sets `param.logscale = True` unconditionally for PATH, and `build_context` passes `logscale=dico_path["logscale"]` (read back from `model_selection/PATHparameters`) at `_visualizer.py:186`.

Escalating factor the reviewer missed: this is not merely cosmetic-and-pre-existing. The same change set *added* `fig.write_html(os.path.join(directory, name), ...)` at :775 — previously that call was commented out, so `predict-path.html` was never written and the pane referenced by `q2_classo/_summarize/assets/path.html:50` was blank. The change set therefore ships a newly-rendered plot whose x-axis title ("- log10 lambda / lambdamax") does not match its x data (raw lambdas), and `test_summarize_assets.py` only checks that *some* writer exists for each iframe, so it passes.

Minimal fix: `_visualizer.py:767` → `graph_objects.Scatter(x=xGraph, y=error, name="L2 error over test set")`.

## REFUTED

**[classo-02-rescaled-lam-flips-default-output] — REFUTED as a defect (the mechanism is accurate; the conclusion is not actionable).**

Every factual step checks out: `_func.py:262`/`:407` `lamfixed_true_lam: bool = True`; `_func.py:349`/`:498` `param.rescaled_lam = not lamfixed_true_lam` → `False`; c-lasso `solver.py:745` and `:691` set `self.rescaled_lam = True` as the library default; `solver.py:178-180` multiplies `theoretical_lam` by `n` when `not rescaled_lam`; `solver.py:1539` `true_lam=not self.rescaled_lam` and `:1551` `self.lamb = self.lambdamax * self.lam` only when rescaled. So yes, the default LAMfixed lambda changes.

But that is the *intent*, and the change does exactly what its comment says. The old `param.true_lam = ...` wrote a dead attribute (no `__slots__`/`__setattr__` guard in `solver.py`, and `true_lam` appears in c-lasso only as the `Classo(..., true_lam=not param.rescaled_lam)` keyword — never as a parameter-object attribute), so `--p-lamfixed-true-lam` and `--p-stabsel-true-lam` were inert. The new behaviour matches the registered documentation verbatim (`_dict.py:225-231`: "If True and lam = -1., then it will takes the value n\*theoretical_lam. Default value : True"). The inversion direction is right (`rescaled_lam` docstring at `solver.py:729`). The reviewer concedes "the fix is right" and proposes only a changelog entry plus a baseline — process advice, not a correctness or regression defect in the code. Also moot in practice until classo-01 is fixed, since no regression path can execute at all.

**[classo-04-cv-nlam-resolver-ambiguous-100] — REFUTED.**

The mechanics are as described (`_func.py:34` guard is `if cv_nlam != 100 and cv_nlam != cv__nlam`, so `_resolve_cv_nlam(100, 120)` returns 120), but this does not meet the reporting bar:

- It is not a regression: `cv_nlam` is brand new in this change set; nothing that worked before changes.
- The test does **not** assert the wrong thing. `test_deprecated_spelling_is_honoured_and_warns` exercises the shim's whole reason for existing — caller leaves `cv_nlam` at its default and passes only the legacy `--p-cv--nlam`. Returning 120 there is the required behaviour, and `test_both_spellings_conflicting_raises` covers the conflict case. The test is neither unfailable nor weakened.
- The only divergence is a docstring sentence over-promising for one input class (`cv_nlam` explicitly set to exactly 100 *and* a different `cv__nlam`) — a user would have to pass both spellings on one command line with one of them equal to the default. There is no wrong answer, no crash, no data loss; the deprecation warning still fires.
- The suggested remedy (`cv_nlam: int = None`) is itself unverified against QIIME 2's signature validation and would be a riskier change than the imprecision it removes.

At most this is a one-line docstring softening; it is not a correctness defect.