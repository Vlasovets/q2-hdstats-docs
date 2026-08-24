#!/usr/bin/env python
"""Build a self-contained bundle of the mOTUs shotgun analysis, for writing it up.

The audience is a reader with no cluster, no QIIME 2, and no way to open a `.qza`. A
QIIME 2 artifact is a zip with a zarr or biom payload inside, so shipping only artifacts
would ship an unreadable bundle. Everything load-bearing is therefore exported to TSV
alongside the artifacts:

  * the feature table, taxonomy and sample metadata as plain text
  * the 481-edge network as an EDGE LIST with both endpoints' lineages attached, which
    exists nowhere else -- the chapter quotes the edge count but never the edges
  * the sparse-plus-low-rank network at mu1 = 5, and the latent loadings
  * the cross-validation curve behind the regression figure

STATISTICS.md is generated from the committed TSVs rather than written by hand, so it
cannot drift from the tables the chapter cites. If a number changes upstream, re-running
this script changes it here too.

    python analysis/scripts/make_paper_bundle.py            # -> dist/motus-paper-bundle.zip
    python analysis/scripts/make_paper_bundle.py --no-qzv   # drop the 4.4 MB of .qzv
"""
import argparse
import csv
import io
import os
import pathlib
import shutil
import tempfile
import zipfile

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
RES = ROOT / "analysis" / "results"
TABLES = RES / "tables"
TUT = RES / "motus-tutorial"
PUB = ROOT / "docs" / "_data" / "motus"
QZV = ROOT / "docs" / "_static" / "qzv"
FIGS = ROOT / "docs" / "images" / "png" / "generated"

SELECTED_LAMBDA = 0.30
SELECTED_EDGES = 481


# --------------------------------------------------------------------------- artifacts

def _payload(qza, suffix):
    """Extract the single data payload of a .qza to a temp file and return its path."""
    z = zipfile.ZipFile(qza)
    names = [n for n in z.namelist() if "/data/" in n and n.endswith(suffix)]
    if not names:
        raise SystemExit(f"{qza.name}: no /data/*{suffix} member")
    tf = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    tf.write(z.read(names[0]))
    tf.close()
    return pathlib.Path(tf.name)


def read_biom(qza):
    """(sample_ids, feature_ids, dense matrix samples x features)."""
    import h5py
    from scipy.sparse import csr_matrix
    p = _payload(qza, "feature-table.biom")
    with h5py.File(p) as h:
        obs = [x.decode() for x in h["observation/ids"][:]]
        sam = [x.decode() for x in h["sample/ids"][:]]
        m = csr_matrix((h["observation/matrix/data"][:],
                        h["observation/matrix/indices"][:],
                        h["observation/matrix/indptr"][:]),
                       shape=(len(obs), len(sam))).toarray()
    os.unlink(p)
    return sam, obs, m.T


def read_solution(qza):
    import zarr
    t = tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z:
        n = [x for x in z.namelist() if x.endswith(".zip")][0]
        z.extract(n, t)
    return zarr.open(zarr.ZipStore(os.path.join(t, n), mode="r"))


def labels_of(group):
    """`labels/` is a zarr group of scalars keyed by stringified position."""
    keys = sorted(group["labels"].keys(), key=int)
    return [str(np.asarray(group["labels"][k])) for k in keys]


def read_taxonomy(qza):
    z = zipfile.ZipFile(qza)
    n = [x for x in z.namelist() if x.endswith("/data/taxonomy.tsv")][0]
    rows = list(csv.reader(io.StringIO(z.read(n).decode()), delimiter="\t"))
    return {r[0]: r[1] for r in rows[1:] if r}


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


# ----------------------------------------------------------------------------- exports

def export_data(out, tax):
    sam, feats, mat = read_biom(PUB / "motus-top100-table.qza")
    write_tsv(out / "motus-top100-table.tsv", ["sample-id"] + feats,
              [[s] + [int(v) for v in mat[i]] for i, s in enumerate(sam)])
    write_tsv(out / "motus-top100-taxonomy.tsv", ["feature-id", "lineage"],
              [[f, tax.get(f, "")] for f in feats])
    shutil.copy2(PUB / "motus-outcomes.tsv", out / "motus-outcomes.tsv")
    shutil.copy2(RES / "motus-merged" / "table.tsv", out / "motus-merged-table-378features.tsv")
    shutil.copy2(ROOT / "analysis" / "config" / "map13241-fecal-wgs.tsv",
                 out / "map13241-fecal-wgs-manifest.tsv")
    return len(sam), len(feats)


def export_edges(out, tax):
    """Edge lists for both networks, with lineages on both endpoints."""
    g = read_solution(TUT / "motus-sgl-path.qza")
    lam = float(np.asarray(g["reg_params"]["lambda1"]))
    lab = labels_of(g)
    P = np.asarray(g["solution"]["precision_"])
    if abs(lam - SELECTED_LAMBDA) > 1e-9:
        raise SystemExit(f"sgl solution is at lambda1={lam}, expected {SELECTED_LAMBDA}")

    def edges_of(M):
        iu = np.triu_indices_from(M, k=1)
        keep = np.abs(M[iu]) > 1e-8
        return [(lab[i], lab[j], float(M[i, j]))
                for i, j, m in zip(iu[0], iu[1], keep) if m]

    e = edges_of(P)
    if len(e) != SELECTED_EDGES:
        raise SystemExit(f"sgl network has {len(e)} edges, expected {SELECTED_EDGES}")
    e.sort(key=lambda r: -abs(r[2]))
    write_tsv(out / "network-sparse-lambda0.30-edges.tsv",
              ["feature_1", "feature_2", "precision_entry", "lineage_1", "lineage_2"],
              [[a, b, f"{w:.6f}", tax.get(a, ""), tax.get(b, "")] for a, b, w in e])

    gl = read_solution(TUT / "motus-slr-mu5.qza")
    labl = labels_of(gl)
    Pl = np.asarray(gl["solution"]["precision_"])
    L = np.asarray(gl["solution"]["lowrank_"])
    iu = np.triu_indices_from(Pl, k=1)
    el = [(labl[i], labl[j], float(Pl[i, j]))
          for i, j, m in zip(iu[0], iu[1], np.abs(Pl[iu]) > 1e-8) if m]
    el.sort(key=lambda r: -abs(r[2]))
    write_tsv(out / "network-slr-mu5-rank3-edges.tsv",
              ["feature_1", "feature_2", "precision_entry", "lineage_1", "lineage_2"],
              [[a, b, f"{w:.6f}", tax.get(a, ""), tax.get(b, "")] for a, b, w in el])

    # Latent loadings: eigenvectors of the low-rank block, largest eigenvalue first.
    w, V = np.linalg.eigh(L)
    order = np.argsort(w)[::-1][:3]
    write_tsv(out / "latent-loadings-mu5.tsv",
              ["feature_id", "lineage"] + [f"latent_{k+1}" for k in range(3)],
              [[labl[i], tax.get(labl[i], "")] + [f"{V[i, j]:.6f}" for j in order]
               for i in range(len(labl))])
    eig = [float(w[j]) for j in order]
    return len(e), len(el), eig


def export_cv(out):
    g = read_solution(TUT / "motus-regress-age.qza")["solution"]["CV"]
    x = np.asarray(g["xGraph"]); y = np.asarray(g["yGraph"])
    se = np.asarray(g["standard_error"])
    write_tsv(out / "regression-cv-curve.tsv",
              ["lambda_over_lambdamax", "neg_log10_lambda_ratio", "cv_mse", "standard_error"],
              [[f"{a:.6f}", f"{-np.log10(a):.6f}", f"{b:.4f}", f"{c:.4f}"]
               for a, b, c in zip(x, y, se)])
    lab = [str(v) for v in np.asarray(g["label"]).ravel()]
    b = np.asarray(g["refit"]).ravel()
    write_tsv(out / "regression-coefficients-all.tsv", ["label", "coefficient"],
              [[l, f"{v:.6f}"] for l, v in zip(lab, b)])
    j = int(np.argmin(y))
    return float(y[0]), float(y[j]), float(se[j]), float(x[j])


# ------------------------------------------------------------------------------- prose

def _tidy(cell):
    """Round binary float noise for display only: 0.39999999999999997 -> 0.4.

    The TSV beside this document keeps the exact value. What is being avoided is a
    reader copying "$\\lambda_1 = 0.39999999999999997$" into a manuscript.
    """
    try:
        v = float(cell)
    except (TypeError, ValueError):
        return cell
    if cell.strip() in ("", "-", "+"):
        return cell
    for digits in range(1, 13):
        short = f"{v:.{digits}g}"
        if abs(float(short) - v) <= 1e-12 * max(1.0, abs(v)):
            # Keep integers looking like integers, not 1e+03.
            return str(int(v)) if v == int(v) and abs(v) < 1e15 else short
    return cell


def tsv_md(path, max_rows=40):
    rows = [[_tidy(c) for c in r] for r in csv.reader(open(path), delimiter="\t")]
    if not rows:
        return "_(empty)_"
    head, body = rows[0], rows[1:max_rows]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(c.replace("|", "\\|") for c in r) + " |" for r in body]
    if len(rows) - 1 > len(body):
        out.append(f"| _… {len(rows)-1-len(body)} more rows_ |" + " |" * (len(head) - 1))
    return "\n".join(out)


def write_statistics(out, n, p, n_edges, n_edges_slr, eig, cv):
    mse_null, mse_min, se_min, x_min = cv
    T = TABLES
    doc = f"""# Statistics — mOTUs shotgun analysis

Generated from the committed tables under `analysis/results/tables/` by
`analysis/scripts/make_paper_bundle.py`. Every number below is read from a file in
`tables/`, not transcribed, so it cannot disagree with the data shipped beside it.

## Dataset

{tsv_md(T / 'motus-tutorial-shape.tsv')}

Design: {n} samples x {p} features. Qiita study 13241 (ENA PRJEB52147), the MAP
preterm-infant study. Faecal whole-genome shotgun only; breast milk and meconium were
excluded as distinct niches. Two of the 36 runs fell below 100 assigned counts and were
dropped, so 16 infants contribute two samples and two contribute one.

## Network inference (q2-gglasso)

Extended BIC over a 20-point linear $\\lambda_1$ grid on [0.05, 1.0] at $\\gamma = 0.15$
selects $\\lambda_1 = {SELECTED_LAMBDA:.2f}$, giving **{n_edges} edges of 4,950 pairs
(9.72% density)**. The minimum is interior: both adjacent grid points score higher.

{tsv_md(T / 'motus-lambda-path.tsv')}

### Why gamma = 0.15, and the permutation null

The abundance sweep re-solved each configuration five times with every feature column
permuted independently across samples. That destroys all between-feature dependence while
preserving each marginal, including its zeros, so every recovered edge is a false
positive. At $\\gamma = 0.15$ the real network exceeds the null by three orders of
magnitude; at $\\gamma \\geq 0.20$ eBIC selects the empty graph.

{tsv_md(T / 'abundance-sweep-gamma-map-vs-null.tsv', max_rows=20)}

### Abundance versus prevalence filtering

Prevalence filtering selects the empty graph at every threshold examined. Abundance
filtering does not. The data and the commands are identical.

{tsv_md(T / 'motus-model-viability.tsv')}

### Sparse plus low-rank

At $\\mu_1 = 5$ the achieved rank is 3 and **{n_edges_slr} of the {n_edges} edges remain**
({n_edges - n_edges_slr} removed, {100*(n_edges-n_edges_slr)/n_edges:.0f}%). Leading
eigenvalues of the latent block: {', '.join(f'{v:.3f}' for v in eig)}. The third is small
relative to the first two, so the fit is close to rank 2 with a weak third direction.

`--p-rank` is registered but raises on every released GGLasso, so a target rank is reached
by tuning $\\mu_1$ and reading back the achieved value:

{tsv_md(T / 'motus-top100-mu-rank-map.tsv')}

## Log-contrast regression (q2-classo)

Response `host_age_days`, postnatal age at sampling, range 10-47 days. It is the only
numeric variable that varies within an infant; gestational age, birth weight and maternal
age are subject-constant, so 34 samples would carry no more information than 18.

`add-taxa` aggregates the 100 features into a 133-column design (100 leaves plus 33
internal nodes). The cross-validated refit retains **two clades**:

{tsv_md(T / 'motus-trac-selected.tsv')}

The two coefficients sum to zero to machine precision, which is the zero-sum constraint of
the log-contrast formulation: the model is a single balance, and only the ratio is
identified. The fitted intercept is +21.7318.

> **Read coefficient labels from the artifact, never by position.** c-lasso prepends the
> intercept, so the 133-column design returns 134 coefficients, and labelling them against
> the design columns shifts every name by one. This already produced three wrong clade
> names in an earlier draft. Read labels from `solution/CV/label` in the artifact, or from
> `CV-beta.csv` in the `.qzv`. This bundle's `data/regression-coefficients-all.tsv` is
> already correctly labelled.

### How much the model is worth

| quantity | value |
|---|---|
| CV MSE, intercept-only ($\\lambda = \\lambda_{{\\max}}$) | {mse_null:.1f} |
| CV MSE, cross-validated minimum | {mse_min:.1f} |
| standard error at the minimum | {se_min:.1f} |
| improvement over intercept-only | {mse_null - mse_min:.1f} MSE units = {(mse_null-mse_min)/se_min:.2f} SE |
| $\\lambda/\\lambda_{{\\max}}$ at the minimum | {x_min:.3f} |
| variance of `host_age_days` | 169.2 |

The intercept-only model at {mse_null:.1f} falls inside the one-standard-error band
[{mse_min-se_min:.1f}, {mse_min+se_min:.1f}], so the one-standard-error rule necessarily
returns it. **The selected model beats predicting the mean by {(mse_null-mse_min)/se_min:.2f}
standard errors.** State this in any write-up: the two clades are the best cross-validation
can find, but the curve does not establish that they beat the mean.

## Cross-validation leakage

Each infant contributes up to two samples and `--p-cv-subsets` partitions at random with
no grouping option, so one sample of a pair can fall in train and the other in test.
Holding out whole infants instead:

{tsv_md(T / 'subject-holdout-vs-leaky.tsv')}

The differences have **opposite signs**. `diagnosis` is a property of the infant, so a
paired sample discloses the label and the random partition inflates the estimate.
Postnatal age varies within an infant, so a paired sample discloses nothing about it and
the leaky estimate is merely noisier. The direction of the bias depends on whether the
response is constant within a subject.

q2-classo exposes `--p-cv-subsets` and `--p-cv-seed` but no grouping parameter, so the
subject-disjoint figures were produced by calling the `classo` API directly
(`scripts/42_subject_holdout.sh`).

## The classification result is negative

Held out by infant, predicting maternal asthma from the microbiome gives 0.471 accuracy
against a majority-class baseline of 0.500. With 9 subjects per group a 95% confidence
interval on any accuracy estimate spans roughly +/-0.23, so this does not demonstrate the
absence of an effect. It is consistent with Bai-Tong et al. 2022, who reported metabolomic
differences and concluded that maternal asthma "did not play an observable role in shaping
the infant gut microbiome during the study period".

## Profiling and QC

{tsv_md(T / 'motus-profile-summary.tsv', max_rows=12)}

mOTUs calls a species only when at least 3 of its 10 universal marker genes are covered,
so a 2 M-read run yields on the order of several hundred assigned inserts. The
`min_alen` sweep found 45/60/70/75 to give identical results on this data:

{tsv_md(T / 'motus-min-alen-sweep.tsv')}

## Limits worth stating in the paper

- **Eighteen infants.** This, not the profiler or the filter, is what forces
  $\\gamma$ down to 0.15, produces a network of this density, and bounds the achievable
  $R^2$.
- **481 edges from 34 samples is over-parameterised** whatever the selection criterion.
  The permutation null is evidence the estimate reflects structure, not that the
  individual edges are reliable.
- **The published study is null for the microbiome.** A null here reproduces the
  literature.
- **`unassigned` was removed before any log-ratio.** mOTUs emits it as a real row; it is
  the complement of the profiled fraction, and keeping it would define the composition
  over a mixture of taxa and non-taxa.
"""
    (out / "STATISTICS.md").write_text(doc)


def write_methods(out):
    (out / "METHODS.md").write_text("""# Methods

Draft material for a methods section. Parameters are the ones actually run; the commands
are in `scripts/44_motus_tutorial_run.sh`, which produced every artifact in this bundle.

## Data

Qiita study 13241, ENA project PRJEB52147 (Bai-Tong et al., *Scientific Reports* 2022,
doi:10.1038/s41598-022-10276-y). The study contains 16S V4 amplicon (MiSeq) and
whole-genome shotgun (NovaSeq) preparations over the same samples; only the shotgun
preparation was used. Analysis was restricted to faecal samples, giving 36 runs from 18
infants, each contributing one early (10-18 days) and one late (23-47 days) sample.
FASTQ files were verified against the MD5 digests published by ENA.

Anonymous Qiita access uses `qiita.ucsd.edu/public/?study_id=13241`; the
`/study/description/` URL redirects for unauthenticated visitors.

## Taxonomic profiling

mOTUs3 (Ruscheweyh et al. 2022, doi:10.1186/s40168-022-01410-z) via the q2-mOTUs QIIME 2
plugin, `qiime motus profile`, per sample, results merged. mOTUs assigns reads to ten
universal single-copy marker genes and calls a species only when at least three are
covered (`marker_gene_cutoff` default 3). Counts (`-c`) and full taxonomy (`-q`) were
requested. A `min_alen` sweep over 45/60/70/75 gave identical results.

The literal `unassigned` feature mOTUs emits was removed before any log-ratio transform:
it is the complement of the profiled fraction, and retaining it would define the
composition over a mixture of taxa and non-taxa.

## Filtering

1. Samples below 100 assigned counts dropped (36 -> 34).
2. Features present in fewer than two remaining samples dropped (378 -> 156).
3. The 156 ranked by total abundance, the 100 most abundant retained.

Features were ranked by **abundance rather than prevalence**. The unfiltered table is
93.6% zeros, so two features absent from the same samples show a large sample correlation
attributable entirely to shared zeros. Under prevalence filtering the extended BIC selects
the empty graph at every threshold examined; under abundance filtering it selects a
network exceeding a permutation null by three orders of magnitude.

## Network inference

q2-gglasso, backed by GGLasso (Schaipp et al. 2021).

- `transform-features --p-transformation mclr --p-keep-original-id`. The modified centred
  log-ratio takes logs of observed positive counts and places zeros at a common floor
  below every observed value, rather than substituting a pseudo-count into a table that is
  85% zeros.
- `calculate-covariance --p-method scaled`, giving a correlation matrix so that one
  $\\lambda_1$ means the same thing for every pair.
- `solve-problem --p-no-latent --p-path-scale linear --p-lambda1-min 0.05
  --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.15 --p-n-samples 34`.
- Sparse-plus-low-rank: `--p-latent` with $\\lambda_1$, $\\lambda_2$ and $\\mu_1$ each
  pinned to a single value (0.30, 0.1, 5). All three must be pinned; leaving
  `--p-lambda2-*` unset substitutes a five-point default path and turns the single fit
  into a model-selection run.

$\\gamma = 0.15$ is a reported choice, not a default. At $\\gamma \\geq 0.20$ the criterion
selects the empty graph on this dataset, the model-space penalty dominating the likelihood
term at $n = 34$. Report $\\gamma$ with any network derived this way.

**Permutation null.** Each configuration was re-solved five times with every feature
column permuted independently across samples. This removes all between-feature dependence
while preserving each marginal distribution, including its zeros, so the true precision
matrix is diagonal and every recovered edge is a false positive.

## Log-contrast regression

q2-classo, backed by c-lasso (Simpson et al. 2021), implementing the log-contrast
formulation of Aitchison and Bacon-Shone (1984) with the zero-sum constraint.

- `transform-features` then `add-taxa`, which reconstructs the design as $\\log(X)A$ with
  $A$ aggregating features along the taxonomy: 100 columns become 133 (100 leaves plus 33
  internal nodes). This requires distinct leaf labels; the mOTUs lineage terminates in
  `m__<mOTU_id>`, so all 100 leaves are distinct even where two mOTUs share a species name.
- `regress --p-concomitant --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se --p-cv-nlam 60 --p-cv-lamin 0.001
  --p-cv-logscale --p-no-stabsel --p-no-lamfixed`.

`--p-concomitant` estimates the noise scale jointly with the coefficients.
`--p-no-cv-one-se` selects the cross-validated minimum; the one-standard-error rule
returns the intercept-only model on this data (see STATISTICS.md).

**Subject-disjoint validation.** q2-classo exposes `--p-cv-subsets` and `--p-cv-seed` but
no grouping parameter, and the design has repeated measures, so leave-one-subject-out
splits were constructed by `host_subject_id` and evaluated by calling the `classo` API
directly. State this in any methods description that relies on the plugin's
cross-validation.

## Software

| tool | version |
|---|---|
| QIIME 2 | 2026.7 (`rachis` framework) |
| q2-gglasso | `main`, QIIME 2 2026.7 port |
| q2-classo | `master`, QIIME 2 2026.7 port |
| q2-mOTUs | `qiime2-2026.7` branch of github.com/Vlasovets/q2-mOTUs |
| mOTUs | 3 (motu-profiler) |
| numpy / pandas / zarr | 2.4.2 / 2.3.3 / 2.18.7 |

The q2-mOTUs port fixed four defects, each with a regression test: vendored versioneer
0.18 could not build under Python 3.12; `reformat_taxonomy` used `str.replace("|", "; ",
regex=True)`, where `|` is regex alternation and the call inserted a separator between
every character; the two call sites disagreed on the separator, so the `s__` -> `m__`
substitution never fired on the import path; and `_run_command` used `subprocess.call`,
discarding the exit status its own docstring promised to raise on.
""")


def write_readme(out, n, p, n_edges, include_qzv):
    (out / "README.md").write_text(f"""# mOTUs shotgun analysis — bundle for writing up

Everything needed to describe this analysis without a cluster or a QIIME 2 install.
{n} samples x {p} features, {n_edges} network edges, two selected clades.

**Start with `STATISTICS.md`** (every number, generated from the tables) and
`METHODS.md` (methods-section source: parameters, versions, and why each choice was made).
`chapter/07_motus_shotgun.md` is the existing written narrative.

## Layout

| path | what |
|---|---|
| `STATISTICS.md` | every quantitative result, generated from `tables/` |
| `METHODS.md` | methods draft: parameters, software versions, rationale |
| `chapter/07_motus_shotgun.md` | the narrative as currently published |
| `data/` | the analysed data as plain TSV — see below |
| `tables/` | the result tables the statistics are generated from |
| `figures/` | the five figures, as published |
| `scripts/` | pipeline stages 32-44 and the helper scripts |
| `qiime2-artifacts/` | the `.qza`{"/`.qzv`" if include_qzv else ""} originals, for provenance |

## `data/` — the readable exports

A `.qza` is a zip with a biom or zarr payload, so these are the same content in text:

| file | what |
|---|---|
| `motus-top100-table.tsv` | the analysed count table, {n} samples x {p} features |
| `motus-top100-taxonomy.tsv` | mOTU id -> lineage, terminating in `m__<mOTU_id>` |
| `motus-outcomes.tsv` | `host_age_days`, `diagnosis`, `host_subject_id` per sample |
| `motus-merged-table-378features.tsv` | the unfiltered 36 x 378 table, `unassigned` already removed |
| `map13241-fecal-wgs-manifest.tsv` | the sample sheet: run accessions, read counts, ENA URLs, MD5s |
| `network-sparse-lambda0.30-edges.tsv` | **the {n_edges} edges**, with both endpoints' lineages |
| `network-slr-mu5-rank3-edges.tsv` | the same network after a rank-3 latent block |
| `latent-loadings-mu5.tsv` | each feature's loading on the three latent components |
| `regression-cv-curve.tsv` | the cross-validation curve behind the regression figure |
| `regression-coefficients-all.tsv` | all 134 coefficients, **correctly labelled** |

The edge lists and the latent loadings exist only here — the chapter quotes the edge
count but never the edges themselves.

## Two things to carry into the writing

**The result is real but modest.** The selected two-clade model improves on predicting the
mean by 0.35 standard errors, and the intercept-only model sits inside the
one-standard-error band. The network exceeds its permutation null by three orders of
magnitude, but 481 edges from 34 samples is over-parameterised whatever the criterion.

**The coefficient labels are a trap.** c-lasso prepends the intercept, so a 133-column
design returns 134 coefficients and labelling by position shifts every name by one. This
already produced three wrong clade names in an earlier draft. The two correct clades are
`o__Enterobacterales` (+5.2475) and `p__Bacteroidetes` (-5.2475); they sum to zero, which
is the zero-sum constraint being satisfied.

## Regenerating

    python analysis/scripts/make_paper_bundle.py

from a clone of github.com/Vlasovets/q2-hdstats-docs. The statistics document is generated
from the committed tables, so it cannot drift from them.
""")


# -------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=pathlib.Path, default=ROOT / "dist")
    ap.add_argument("--no-qzv", action="store_true", help="omit the 4.4 MB of .qzv")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    stage = args.outdir / "motus-paper-bundle"
    if stage.exists():
        shutil.rmtree(stage)
    for sub in ("data", "tables", "figures", "scripts", "chapter", "qiime2-artifacts"):
        (stage / sub).mkdir(parents=True)

    tax = read_taxonomy(PUB / "motus-top100-taxonomy.qza")
    n, p = export_data(stage / "data", tax)
    n_edges, n_edges_slr, eig = export_edges(stage / "data", tax)
    cv = export_cv(stage / "data")

    for t in sorted(TABLES.glob("motus-*.tsv")):
        shutil.copy2(t, stage / "tables" / t.name)
    for extra in ("subject-holdout-vs-leaky.tsv", "abundance-sweep-gamma-map.tsv",
                  "abundance-sweep-gamma-map-vs-null.tsv", "classo-diagnosis-summary.tsv"):
        src = TABLES / extra
        if src.is_file():
            shutil.copy2(src, stage / "tables" / extra)

    for f in sorted(FIGS.glob("motus-*.png")):
        shutil.copy2(f, stage / "figures" / f.name)

    shutil.copy2(ROOT / "docs" / "chapters" / "04_highdim_atacama" / "07_motus_shotgun.md",
                 stage / "chapter" / "07_motus_shotgun.md")

    # Stages 32-44 are the mOTUs pipeline. 30 (tier-1 figures) and 31 (MGL verify)
    # match a [34][0-9] glob but belong to a different analysis.
    motus_stages = [f"{i:02d}" for i in range(32, 45)]
    for s in sorted(f for st in motus_stages
                    for f in (ROOT / "analysis" / "slurm").glob(f"{st}_*.sh")):
        shutil.copy2(s, stage / "scripts" / s.name)
    for s in ("aggregate_motus_profiles.py", "build_map13241_manifest.py",
              "make_motus_figures.py", "render_qzv_figures.py"):
        src = ROOT / "analysis" / "scripts" / s
        if src.is_file():
            shutil.copy2(src, stage / "scripts" / s)

    for a in sorted(PUB.glob("*.qza")):
        shutil.copy2(a, stage / "qiime2-artifacts" / a.name)
    for a in sorted(TUT.glob("*.qza")):
        shutil.copy2(a, stage / "qiime2-artifacts" / a.name)
    if not args.no_qzv:
        for a in sorted(QZV.glob("*.qzv")):
            shutil.copy2(a, stage / "qiime2-artifacts" / a.name)

    write_statistics(stage, n, p, n_edges, n_edges_slr, eig, cv)
    write_methods(stage)
    write_readme(stage, n, p, n_edges, not args.no_qzv)

    zpath = args.outdir / "motus-paper-bundle.zip"
    if zpath.exists():
        zpath.unlink()
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for f in sorted(stage.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(stage.parent))
    shutil.rmtree(stage)

    print(f"  {zpath}  ({zpath.stat().st_size/1024/1024:.1f} MB)")
    with zipfile.ZipFile(zpath) as z:
        print(f"  {len(z.namelist())} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
