# Draft email — C. Müller

Two questions remain that need your answer. A third (where the shipped
correlation matrix came from) is now resolved and is included only as an FYI,
because it changes one of your artifacts' status.

Not sent. Edit freely — the technical claims are all verified and are safe to
keep, but the tone is a first draft.

---

**Subject:** q2-gglasso / q2-classo on QIIME 2 2026.7 — two things to confirm before we mint the DOI

Hi Christian,

Both plugins now install and pass their tests on QIIME 2 2026.7 (they were
pinned to 2025.4 and could not install at all — `numpy<=1.27`, `gglasso<=0.2.1`,
`bokeh==2.4.3` against a distribution shipping numpy 2.4.2 on Python 3.12).
45 q2-gglasso tests and 22 q2-classo tests pass, the CLI round-trip works, and
the book builds warning-clean. I have not pushed anything yet.

Two questions before the Zenodo record goes live, since the DOI freezes both.

**1. Is λ = 0.8 the canonical parameterisation?**

The evidence says yes, and I would like you to confirm before the recompute
overwrites the share bundle. Both values are points on the same path, which the
CLI now reproduces exactly:

| λ | edges | where it is used |
|---|---|---|
| 0.95 | 145 | `atacama-q2-gglasso-share.zip` |
| 0.80 | 216 | the tutorial chapters, `tutorial-notes.md`, `run_q2_gglasso_lambda08_models.sh` |

File timestamps in the tutorial archive put the λ=0.95 fits earlier (06-27
19:28) than the λ=0.8 models (06-28 09:24), so the share bundle looks like the
earlier exploratory run rather than a disagreement. If that is right, λ=0.8 is
canonical and the share bundle is superseded. Say if you intended otherwise.

**2. Did you mean `stabsel_true_lam` to invert?**

q2-classo assigned `param.true_lam`, but c-lasso reads `rescaled_lam` — assigning
an undefined attribute silently created one that nothing consulted, so both
`--p-stabsel-true-lam` and `--p-lamfixed-true-lam` were inert. c-lasso's
`solver.py:643` documents `rescaled_lam` as "False if lam = lambda, True if
lam = lambda/lambdamax", which is the inverse of q2's `true_lam`, so I applied
`rescaled_lam = not true_lam` at all four call sites. Only consulted when
`stabsel_method == 'lam'`, and the default is `'first'`, so the blast radius is
small — but it does change results for anyone who set it. Confirm that is the
semantics you wanted.

**FYI — the shipped correlation matrix is fine, and one of your artifacts is not
what it appears to be.**

`atacama-top-300-correlation.qza` looked unreproducible: regenerating it from
the documented chain gave a matrix differing by up to 1.147. That turned out to
be an artefact of comparing the two **by label**. They are the same matrix with
the features in a different order — identical spectrum to 2.1e-14, and the
recovered permutation reproduces it exactly (residual 0.000e+00).

The cause is worth knowing because it affects how anyone should read these
artifacts. `--p-no-keep-original-id` assigns `ASV-k` by **position** in an
ascending total-abundance ranking, and 209 of those 300 features share a
total-abundance value with another feature. The relabelling helper used
`df.sort_index()`, whose default quicksort is not stable; switching to a stable
sort moved 158 of the 300 features to different `ASV-k` labels. Every published
number was unaffected — the graphical-lasso objective is invariant under
simultaneous row/column permutation, so λ, the eBIC at every grid point and the
216 edges are bit-identical — which is exactly why nobody noticed.

Two consequences:

- **`ASV-k` is not an identifier**, and the recovery procedure the tutorial
  documented (`ASV-n ↔ abundance-rank = 301 − n`) does not work. Measured: it
  places 146 of 300 features wrongly and raises nothing. I have removed it and
  rebuilt the tier-2 bundle with `--p-keep-original-id`, so features now carry
  real IDs and join to taxonomy directly (300/300 resolve). Your review comment
  asking to keep the original ASV identifier was right, and for a stronger
  reason than either of us had at the time.
- Ties now break on the feature ID, so `ASV-k` is a pure function of the data
  rather than of the row order the table happened to arrive in.

Nothing here changes λ = 0.8, 216 edges or eBIC 16130.0988.

Best,
Oleg

---

## Notes for me, not for the email

- Do not send until the branches are pushed, so the links resolve.
- Item 4 in `DECISIONS_NEEDED.md` (is `selected-atacama-sample-metadata.tsv`
  going on the Zenodo record?) is a release decision, not one for Christian.
- The MOSHPIT / shotgun question is for Evan and the F1000 reviewer, not here.
