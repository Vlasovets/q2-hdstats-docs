# Sidenotes to Section S1

Material cut from Section S1 of the supplement during the editorial pass on
`supplement/s1-prose-pass`. Nothing here is wrong; it is secondary to what S1 is
for, which is stating each estimator and the call that runs it. Kept because
several items are worth having in a response letter, a revision, or the
documentation site.

## Related network estimators

These were briefly at risk during the pass and are worth flagging: the single
graphical lasso paragraph cites gCoda (Fang et al., 2017), SparCC (Friedman &
Alm, 2012) and NetCoMi (Peschel et al., 2021) in one clause. If S1 is trimmed
further, these are the first candidates to cut — but dropping them leaves those
three `supplementary.bib` entries uncited, so cut the entries too or accept the
orphans.

## Why the concomitant coefficients are not hyperparameters

S1 states that the `n/2` in (S9) and the `n` in (S10) are fixed constants arising
from the perspective function of the respective loss (Combettes & Müller, 2018,
2020). The fuller version, cut for length: the perspective construction is what
makes the joint objective in (β, σ) convex, which is why estimating σ jointly
costs nothing in tuning burden. Restore if a reviewer questions the claim.

## Symbol choices made to avoid collisions

The notation pass renamed four symbols. S1 no longer explains why, since the
table now fixes each meaning, but the reasoning matters if anyone edits these
sections later:

| Symbol | Was also used for | Renamed to |
|---|---|---|
| `𝕊` | simplex *and* SPD cone | simplex is `Δ` |
| `Λ` | penalty matrix (S16) *and* rPCA eigenvalues | eigenvalues are `D` |
| `γ` | trac node coefficients *and* the eBIC parameter in S2 | coefficients are `α` |
| `r` | residual *and* `rank(L̂)` | rank is `r_L` |

Also: the count matrix is `N`, not `C`, because `C` is the constraint matrix
throughout S1.2. And `r` (residual) and `r_L` (rank) never appear in the same
subsection, so the reminder that they are unrelated was cut from S1.4.

## Interface observations cut as secondary

- **Zero-sum by default.** No option switches the constraint on. Stated in S1.2.1
  as a property of the model; the fact that there is no flag is an interface
  detail.
- **`--p-rank` is registered but guarded.** In the parameter table
  (Table S2) with a pointer to S1.3.4. The longer explanation — that the rank is
  an output of `μ₁` so a requested rank cannot be honoured directly — is in
  S1.3.4 itself.
- **trac design labels.** After `add-taxa` the design columns are tree nodes, so
  coefficient labels are clade names and a planned feature-level join downstream
  no longer applies. Retained in S1.2.4 in shortened form; the fuller warning
  belongs in the tutorial rather than here.

## Open question for the maintainers

The `concomitant` default was read from the registered signature in
`q2_classo/_func.py`:

```python
def regress(..., concomitant: bool = True, huber: bool = False, ...)
```

with `concomitant` passed straight to `problem.formulation.concomitant`. That
makes R3 the no-option default, which is what Table S3 and S1.2.2 now state. This
has not been confirmed by an actual CLI invocation. If the plugin is reached
through a wrapper that overrides the default, the effective behaviour could differ
and both the table footnote and the R1–R4 pseudocode would need revising.

---

# Sidenotes to Sections S2–S4

Second pass, same branch. Material cut from the two data sections and the summary
tables.

## Figure-reading guidance that repeated the caption

Several paragraphs following a figure restated what its caption already said, then
drew a conclusion the figure supports on its own. Cut in each case; the caption
carries the reading. The pattern was: caption describes panels, following
paragraph re-describes the same panels and adds "which is the point of", "is the
argument of Section S1 in one image", "makes a point that no single outcome
shows".

## The stability-selection threshold example

The regression caption (Figure S6) previously named two taxa on either side of the
0.7 threshold — ASV₂ at 0.68 and ASV₁₁ at 0.75 — to show that no theory places the
line at 0.7 rather than 0.65. The general statement survives in the caption; the
worked example is here because it is a good illustration if a reviewer asks how
arbitrary the threshold is.

## The SLR nesting caveat, in full

S2.2.4 states that the latent block absorbs structure rather than unmasking it on
this table, and that nesting should be checked before an SLR result is
interpreted. The mechanism was cut: removing a dense confounding direction can
reveal an association that was previously cancelled by it, so an SLR network need
not be a subset of the corresponding SGL network. Nesting held here at all three
ranks, but it is not a property of the model.

## `toc` under nested selection

S2.2.5 reports that only `toc` moves the other way under nested cross-validation,
from −0.026 to +0.043, and calls it noise at 43 non-missing samples. The reason it
moves at all was cut: `toc` selects no features on the full data, so under nesting
some folds select a non-zero model and the near-null prediction improves
marginally.

## Duplicated dataset description in S3

The metagenomics section had two paragraphs describing the same subset — one under
"How the subset was obtained", one three paragraphs later opening mid-sentence
with `$p = 100$ mOTUs against $n = 34$ faecal samples`. The second was a leftover
from an earlier draft. The two facts it carried and the first paragraph did not —
the singularity of the covariance at p/n = 2.94, and why this study was chosen —
were merged into it; the rest was duplicate.

Also fixed there: a sentence beginning "$\gamma \geq 0.20$ the criterion selects
the empty graph", missing its opening preposition.

## Interface details cut from S3

- **`--p-rank` raises on every released GGLasso version.** S3.1.3 now says the
  rank is an output and must be reached by tuning μ₁. That the registered
  parameter actively raises an exception is in Table S2's notes column.
- **`--p-cv-subsets` and `--p-cv-seed` exist but no grouping option does.** S3.1.5
  now says the plugin offers no grouping option, which is the operative fact; the
  parameter names are in Table S3.
- **`--p-no-cv-one-se`.** S3.1.4 now says which of the two rules is used is
  material to the result, without naming the flag.

## All-caps emphasis in the summary tables

Sixteen words were set in full capitals for emphasis (`NOT`, `OUTPUT`, `FIXED`,
`PER-SAMPLE`, `TRANSFORMED`, and so on). All normalised to lower case. The
distinctions they marked are carried by the surrounding sentence in every case; if
any needs visual emphasis in a future revision, use `\textbf` rather than caps so
it survives copy-paste into plain text.
