#!/bin/bash
#SBATCH --job-name=q2-perm
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/perm_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/perm_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Hypothesis: the shipped and regenerated correlation matrices hold the SAME
# numbers with the features in a DIFFERENT ORDER, and --p-no-keep-original-id
# hands out the labels ASV-1..ASV-300 by position -- so "ASV-17" names a
# different organism in each file.
#
# This reconciles the two results that looked contradictory:
#   - eBIC agrees at all 15 grid points to 1 ULP, and edges = 216 in both.
#     The graphical-lasso objective is INVARIANT under simultaneous row/column
#     permutation, so a reordered matrix gives a bit-identical path.
#   - a LABEL-aligned entrywise comparison differs on 18% of pairs, because it
#     permutes one matrix into the other's label order while the values stay put.
#
# Decisive because eigenvalues are a permutation invariant: if the spectra match
# but the label-aligned entries do not, the matrices are permutation-equivalent
# and the defect is in the LABELS, not the numbers. If the spectra differ too,
# they are genuinely different matrices and the eBIC agreement is the anomaly.
#
# Writes nothing outside reports/.

set -euo pipefail
ROOT=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/PERMUTATION_HYPOTHESIS.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-perm/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-perm/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

qiime gglasso transform-features \
  --i-table "$ROOT/data/atacama-top-300-table.qza" \
  --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
  --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
  --p-transformation clr --p-pseudo-count 1 --p-no-keep-original-id \
  --o-transformed-table clr.qza >/dev/null
qiime gglasso calculate-covariance \
  --i-table clr.qza --p-method scaled --o-covariance-matrix corr-regen.qza >/dev/null
echo "regenerated"

# keep the CLR table too -- if the hypothesis holds, the feature ORDER in the
# transformed table is what determines the ASV-k assignment.
qiime tools export --input-path clr.qza --output-path clr-exp >/dev/null 2>&1 || true

python - "$ROOT/data/atacama-top-300-correlation.qza" corr-regen.qza \
         "$ROOT/data/atacama-top-300-table.qza" clr.qza "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd

ship_q, regen_q, raw_q, clr_q, report = sys.argv[1:6]

def load_tsv(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        n=[x for x in z.namelist() if x.endswith(".tsv") and "/data/" in x][0]
        z.extract(n,t)
        return pd.read_csv(os.path.join(t,n), sep="\t", index_col=0)

A, B = load_tsv(ship_q), load_tsv(regen_q)
a, b = A.to_numpy(float), B.to_numpy(float)

L=["# Are the two correlation matrices the same matrix, reordered?","",
   f"- shipped     : {A.shape}, labels {list(A.index[:3])} ... {list(A.index[-2:])}",
   f"- regenerated : {B.shape}, labels {list(B.index[:3])} ... {list(B.index[-2:])}",
   f"- label sets equal : {set(A.index)==set(B.index)}",
   f"- label order equal : {list(A.index)==list(B.index)}",""]

pos = float(np.nanmax(np.abs(a-b)))
L += ["## Entrywise", "",
      f"- POSITIONAL max|diff| (labels ignored) : **{pos:.6e}**"]
if set(A.index)==set(B.index):
    Bal = B.loc[list(A.index), list(A.index)].to_numpy(float)
    lab = float(np.nanmax(np.abs(a-Bal)))
    L += [f"- LABEL-ALIGNED max|diff|               : **{lab:.6e}**"]
L += [""]

# --- the permutation invariant -------------------------------------------
ea = np.sort(np.linalg.eigvalsh((a+a.T)/2))
eb = np.sort(np.linalg.eigvalsh((b+b.T)/2))
spec = float(np.max(np.abs(ea-eb)))
L += ["## Spectrum (invariant under simultaneous row/column permutation)","",
      f"- max |eigenvalue difference| : **{spec:.6e}**",
      f"- trace: {np.trace(a):.6f} vs {np.trace(b):.6f}",
      f"- Frobenius norm: {np.linalg.norm(a):.6f} vs {np.linalg.norm(b):.6f}",""]

same_spectrum = spec < 1e-8

# --- try to RECOVER the permutation --------------------------------------
# Two matrices related by B = P A P^T have the same multiset of rows (up to the
# same permutation applied within each row). Match rows by their sorted values.
def rowkey(M, i):
    return tuple(np.round(np.sort(M[i, :]), 9))
keyA = {}
for i in range(a.shape[0]):
    keyA.setdefault(rowkey(a, i), []).append(i)
perm, unmatched = [], 0
for j in range(b.shape[0]):
    cands = keyA.get(rowkey(b, j))
    if cands: perm.append(cands[0])
    else: perm.append(None); unmatched += 1

recovered = unmatched == 0 and len(set(x for x in perm if x is not None)) == b.shape[0]
L += ["## Permutation recovery","",
      f"- rows of the regenerated matrix matched into the shipped one : "
      f"**{b.shape[0]-unmatched} / {b.shape[0]}**"]
if recovered:
    p = np.array(perm)
    resid = float(np.nanmax(np.abs(a[np.ix_(p, p)] - b)))
    n_moved = int((p != np.arange(len(p))).sum())
    L += [f"- verified `A[perm][:,perm] == B` to : **{resid:.3e}**",
          f"- features whose position changes : **{n_moved} / {len(p)}**",""]
    ex = [(list(B.index)[j], list(A.index)[p[j]]) for j in range(len(p)) if p[j] != j][:8]
    if ex:
        L += ["| label in regenerated | same numbers live at, in shipped |","|---|---|"]
        L += [f"| `{u}` | `{v}` |" for u, v in ex] + [""]
else:
    L += ["- could not recover a clean permutation, so the two are not a pure "
          "reordering of each other.",""]

# --- what determines the ordering? ---------------------------------------
try:
    import qiime2, biom
    raw = qiime2.Artifact.load(raw_q).view(biom.Table).to_dataframe()
    clr = qiime2.Artifact.load(clr_q).view(biom.Table).to_dataframe()
    L += ["## Feature order in the inputs","",
          f"- raw table observations[:3] : {list(raw.index[:3])}",
          f"- clr table observations[:3] : {list(clr.index[:3])}",
          f"- clr order == regenerated correlation order : "
          f"{list(clr.index)==list(B.index)}",""]
except Exception as e:
    L += [f"(input-order check skipped: {e})",""]

# --- verdict --------------------------------------------------------------
if same_spectrum and recovered:
    L += ["## Verdict","",
          "**The numbers are identical; the LABELS are permuted.**","",
          "The two files hold the same correlation matrix with the features in a "
          "different order. That is why the lambda path, the eBIC at every grid "
          "point and the edge count are bit-identical: the graphical-lasso "
          "objective does not depend on the ordering of the variables. It is also "
          "why a label-aligned entrywise comparison looked catastrophic -- it was "
          "comparing organism i's row against organism j's row.","",
          "Consequences, in order of severity:","",
          "1. **The 1.147 discrepancy was never real.** Every conclusion drawn "
          "   from it -- that the shipped matrix could not be regenerated, that "
          "   its provenance was unknown -- was an artefact of comparing by label. "
          "   The documented chain DOES reproduce the shipped matrix.",
          "2. **`ASV-k` is not a stable identifier.** `--p-no-keep-original-id` "
          "   assigns it by position, so the same name denotes different organisms "
          "   in two artifacts built from differently-ordered tables. Any chapter "
          "   that maps `ASV-k` to a taxon is only valid for one specific artifact.",
          "3. Gate C1 and the published lambda/edge numbers stand unchanged.",""]
elif same_spectrum:
    L += ["## Verdict","",
          "**Same spectrum, but no clean permutation recovered.** The matrices are "
          "spectrally identical yet not a simple reordering. Investigate before "
          "republishing.",""]
else:
    L += ["## Verdict","",
          "**Genuinely different matrices** -- the spectra differ. The identical "
          "eBIC path at all 15 grid points is then unexplained and points at "
          "`solve-problem` not consuming its covariance argument. Run the "
          "corrupted-covariance probe next.",""]

open(report,"w").write("\n".join(L)+"\n"); print("\n".join(L))
PY
echo "report -> $REPORT"
