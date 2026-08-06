#!/bin/bash
#SBATCH --job-name=q2-classo-cv
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classo_cv_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classo_cv_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --array=1-15
#
# Stage 6 -- log-contrast regression with cross-validation, one array task per
# numeric outcome (15 of them) against the tier-2 300-ASV clr table.
#
# This is the first thing in the whole recompute that exercises `qiime classo
# regress` on real data. Until the np.infty shim landed it could not have run at
# all: c-lasso 1.0.11 calls a NumPy-2-removed alias in five solver files.
#
# Idempotent: an existing output is skipped, so a partial array can be resumed
# by resubmitting.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda

# NOTE the input: the RAW count table, CLR-transformed by q2-classo's own action
# rather than the ready-made atacama-top-300-clr.qza from `qiime gglasso
# transform-features`. Two reasons, one current and one historical.
#
# Current, and the reason to keep it: this is exactly what the tutorial chapter
# documents, so the stage runs the commands a reader would run. It also keeps the
# classo chain independent of the gglasso one -- notably, the 2026-08-05 tier-2
# regeneration that switched features to real IDs does not touch these results.
#
# Historical: this used to be forced. transform_features transposed to (N, p) and
# straight back to (p, N) before returning, and QIIME 2 stores a DataFrame's index
# as SAMPLES, so the stored table's "samples" were feature IDs. q2-classo aligned
# features against y on the sample index, found no overlap, and c-lasso died with
#     IndexError: index 0 is out of bounds for axis 0 with size 0
# on an X of shape (0, 54) -- the first run of this array failed 15/15 that way.
# That bug is FIXED (q2-gglasso merge 1c295b7, coupled with rowvar=False in
# calculate_covariance), so the gglasso artifact would work now. Do not "restore"
# it on that basis: the reason above still stands.
RAW="$ROOT/data/atacama-top-300-table.qza"
OUTCOMES="$ROOT/data/atacama-classo-outcomes-mean-imputed.tsv"
OUT="$ROOT/results/classo"
DESIGN="$OUT/atacama-top-300-classo-clr.qza"
mkdir -p "$OUT" "$ROOT/slurm/logs"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

# Outcome for this array index: header, minus the leading sample-id column.
IDX=${SLURM_ARRAY_TASK_ID:-1}
OUTCOME=$(head -1 "$OUTCOMES" | tr '\t' '\n' | tail -n +2 | sed -n "${IDX}p")
[[ -n "$OUTCOME" ]] || { echo "no outcome at index $IDX" >&2; exit 2; }
# printf, not echo: echo appends a newline and `tr -c` turns it into a trailing
# dash, so filenames came out as "...-elevation-.qza".
SAFE=$(printf '%s' "$OUTCOME" | tr -c 'A-Za-z0-9._-' '-')
RESULT="$OUT/atacama-top-300-regress-${SAFE}.qza"

echo "=== [$IDX/15] outcome: $OUTCOME"
if [[ -f "$RESULT" ]]; then echo "[skip] $RESULT already exists"; exit 0; fi

# Build the CLR design once. Array tasks race here, so write to a task-private
# path and move into place atomically; whoever wins, all tasks then read the
# same artifact.
if [[ ! -f "$DESIGN" ]]; then
  TMPD="$SCRATCH/design-${IDX}.qza"
  qiime classo transform-features --i-features "$RAW" --o-x "$TMPD"
  mv -n "$TMPD" "$DESIGN" 2>/dev/null || true
fi
[[ -f "$DESIGN" ]] || { echo "could not build $DESIGN" >&2; exit 3; }

# Fail loudly rather than handing c-lasso an empty matrix: confirm the design and
# the outcomes actually share sample IDs.
python - "$DESIGN" "$OUTCOMES" <<'PY'
import sys, zipfile, tempfile, os, biom, pandas as pd
qza, tsv = sys.argv[1], sys.argv[2]
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(qza) as z:
    i = [n for n in z.namelist() if n.endswith("feature-table.biom")][0]
    z.extract(i, t); tab = biom.load_table(os.path.join(t, i))
ids = set(tab.ids("sample"))
md = pd.read_csv(tsv, sep="\t", index_col=0)
md = md[~md.index.astype(str).str.startswith("#")]
overlap = ids & set(md.index.astype(str))
print(f"  design samples={len(ids)}  outcomes rows={len(md)}  overlap={len(overlap)}")
if not overlap:
    sys.exit("design and outcomes share NO sample IDs -- check table orientation")
PY

# PATH + CV only. StabSel is deliberately off: it refits B=50 subsamples per
# lambda and would dominate runtime without changing the CV numbers the chapter
# reports. --p-no-lamfixed for the same reason.
qiime classo regress \
  --i-features "$DESIGN" \
  --m-y-file "$OUTCOMES" \
  --m-y-column "$OUTCOME" \
  --p-concomitant \
  --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
  --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale \
  --p-no-stabsel \
  --p-no-lamfixed \
  --o-result "$RESULT" \
  --verbose

echo "wrote $RESULT"
