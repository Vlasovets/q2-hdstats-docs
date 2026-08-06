## CONFIRMED

**[GG-PKG-01] `/home/itg/oleg.vlasovets/slr_example/q2-gglasso/.github/workflows/ci.yml:69-82` — new CI step fails on the tree that introduces it (high)**

Verified: I ran the exact inline script from the step at the repo root. It prints
`WOULD EXIT NONZERO: stale bokeh 2 / CDN reference in: q2_gglasso/_summarize/assets/index.html, q2_gglasso/_pca/assets/index.html`.
The match is real and is in both templates at line 9, inside the new `{# … #}` Jinja block: `hardcoded cdn.bokeh.org bokeh-2.4.3 <script> tags`. The grep reads raw file text (`tpl.read_text()`), so Jinja comments are not stripped. The step has no `working-directory`, so it runs in `${{ github.workspace }}` and `pathlib.Path("q2_gglasso").rglob(...)` resolves to exactly those two files. It is the last step in the job, so every push/PR to `dev`/`main` goes red *after* the tests pass, which reads as a visualization regression.

Minimal fix: strip `{#…#}` before scanning, or restrict the scan to lines containing `src=`/`href=`. Re-run the script locally until it prints `templates clean`.

---

**[GG-PKG-02] same step, `ci.yml:69-70` comment — claims coverage it does not provide (low; the reviewer's headline claim is partly wrong, see below)**

Verified: the templates contain no literal script/link tags — `q2_gglasso/_pca/assets/index.html:15` and `q2_gglasso/_summarize/assets/index.html:15` are just `{{ bokeh_resources|safe }}`, and the runtime is injected at render time from `q2_gglasso/_pca/_visualizer.py:444` and `q2_gglasso/_summarize/_visualizer.py:533` (`bokeh_resources=INLINE.render()`). The step never renders a visualization or unzips a `.qzv`. I also grepped `q2_gglasso/tests/*.py`: no test asserts anything about the bokeh runtime in a rendered `.qzv` (only `test_heatmap.py:10-11` imports `file_html`/`CDN`). So the failure mode the comment names — "A bokeh version mismatch renders a blank .qzv with no Python traceback, so assert on the generated HTML rather than trusting the visualizer to raise" — is not covered anywhere in CI, and the comment misleads a maintainer into thinking it is.

Correcting the reviewer: the claim that the check "can never fail" is **false**. It fails today (GG-PKG-01), and after a correct GG-PKG-01 fix it still catches a contributor re-adding a hardcoded `cdn.bokeh.org/…/bokeh-2.x.min.js` tag — exactly the regression just removed (see `git show HEAD:q2_gglasso/_pca/assets/index.html:6-8`). It is a narrow but non-vacuous hygiene guard.

Minimal fix: reword the comment to describe what it does (template hygiene), and drop the "generated HTML" claim. Do not build a `.qzv` in CI on the strength of this finding.

---

**[GG-PKG-03] `/home/itg/oleg.vlasovets/slr_example/q2-gglasso/Makefile:16` — `test` lost its install prerequisite (low)**

Verified: `git show HEAD:Makefile` has `test: all` (→ `install` → `pip install -e .`); the working tree has a bare `test:`. I read `rachis/plugin/testing.py` in the env (qiime2 2026.7 aliases `qiime2.*` to `rachis.*`): `TestPluginBase.setUp` iterates `rachis.sdk.PluginManager().plugins` and matches `plugin_.package == package`, calling `self.fail('%s is not a registered Rachis plugin.')` on no match. Discovery is via the `qiime2.plugins` entry point declared in `setup.py`, so with no install there is no match. `environment-files/q2-gglasso-qiime2-2026.7.yml` does **not** install the plugin — its own header says "the plugin is installed separately (`pip install -e .`)". So `conda env create && make test` now fails where it previously self-installed.

Scope is narrower than the reviewer implies: only `q2_gglasso/tests/test_type.py` subclasses `TestPluginBase` (2 tests) — it is the sole `TestPluginBase` user in the suite. Collection itself still works (`q2_gglasso/tests/__init__.py` exists, so pytest prepends the repo root). CI is unaffected (separate `pip install --no-deps -e .` step). The added comment only justifies the explicit test path, not the dropped prerequisite.

Minimal fix: `test: install`. (Note `dev: all` → `dev:` is *not* a defect — `dev` runs `pip install -e ".[dev]"` itself, and the `dev` extra now exists in `setup.py`.)

---

**[GG-PKG-04] `setup.py:16` vs `README.md:39` — half-applied repo URL change (low)**

Verified at file level: `setup.py:16` is now `URL = "https://github.com/Vlasovets/q2-gglasso.git"` (was bio-datascience, per the diff); `git remote -v` origin is Vlasovets; `README.md:5` CI badge is Vlasovets; but `README.md:39` still reads `pip install git+https://github.com/bio-datascience/q2-gglasso.git`. README.md is untouched by this change set (`git status --short` does not list it), so the change set moves setup.py and leaves the documented install command pointing elsewhere.

Caveat on the reviewer's consequence claim: I did **not** verify over the network what code currently sits on bio-datascience. That the fork's pre-change state had `bokeh==2.4.3` / `numpy<=1.27` is verified locally (`git diff requirements.txt`), and the branch history contains merges *from* `bio-datascience/main`, which makes "the README URL yields pre-fix code" plausible but inferred, not proven. The file-level contradiction is proven.

Minimal fix: make `README.md:39` use the same URL as `setup.py:16` (or revert `setup.py` if bio-datascience is meant to stay canonical) — one of the two, consistently.

## REFUTED

None of the four were fully refuted. The only sub-claim refuted is inside GG-PKG-02: "once GG-PKG-01 is fixed the check becomes one that can never fail" is incorrect — the grep remains a working regression guard against reintroduced hardcoded CDN `<script>` tags, and its proposed heavyweight fix (generate/unzip a `.qzv` in CI) is not warranted by the evidence.