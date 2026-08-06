#!/bin/bash
#SBATCH --job-name=q2-tests
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/tests_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/tests_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Post-migration verification: run both plugins' test suites plus a CLI
# round-trip that exercises the things unit tests cannot -- artifact round-trip,
# visualizer rendering, and the bokeh-3 template fix.
#
# The bokeh check is the important one: a version mismatch produces a blank .qzv
# with NO Python traceback, so it can only be caught by inspecting the generated
# HTML.

set -euo pipefail

PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
GG=/home/itg/oleg.vlasovets/slr_example/q2-gglasso
CL=/home/itg/oleg.vlasovets/slr_example/q2-classo

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break
  fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; export TMP="$TMPDIR"; export TEMP="$TMPDIR"
WORK="$SCRATCH/work"
mkdir -p "$TMPDIR" "$WORK"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"
conda activate "$PREFIX"

rc=0

echo "############ [1/4] q2-gglasso unit tests ############"
( cd "$GG" && pytest q2_gglasso/tests -v --durations=15 ) || rc=1

echo "############ [2/4] q2-classo unit tests ############"
( cd "$CL" && pytest q2_classo/tests -v --durations=15 ) || rc=1

echo "############ [3/4] CLI round-trip (q2-gglasso) ############"
cd "$WORK"
set -x
qiime gglasso transform-features \
  --i-table "$GG/data/atacama-counts.qza" \
  --i-taxonomy "$GG/data/classification.qza" \
  --m-sample-metadata-file "$GG/data/selected-atacama-sample-metadata.tsv" \
  --p-transformation mclr \
  --o-transformed-table mclr.qza

qiime gglasso calculate-covariance \
  --i-table mclr.qza --p-method scaled --o-covariance-matrix S.qza

# model-selection path (multi-value grid)
qiime gglasso solve-problem \
  --i-covariance-matrix S.qza --p-n-samples 50 \
  --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 5 \
  --p-no-latent --o-solution sol-sgl.qza

# Latent fit with a SMALL mu1. A large mu1 shrinks the low-rank block to rank 0,
# which is a legitimate outcome but leaves pca nothing to project onto -- that
# case is asserted separately below.
qiime gglasso solve-problem \
  --i-covariance-matrix S.qza --p-n-samples 50 \
  --p-lambda1-min 0.05 --p-lambda1-max 0.05 --p-n-lambda1 1 \
  --p-latent --p-mu1-min 0.05 --p-mu1-max 0.05 --p-n-mu1 1 \
  --o-solution sol-slr.qza

# Large mu1 -> rank 0, for the negative test.
qiime gglasso solve-problem \
  --i-covariance-matrix S.qza --p-n-samples 50 \
  --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 3 \
  --p-latent --p-mu1-min 1 --p-mu1-max 10 --p-n-mu1 3 \
  --o-solution sol-rank0.qza

qiime gglasso summarize --i-solution sol-sgl.qza --o-visualization sum-sgl.qzv
qiime gglasso summarize --i-solution sol-slr.qza --o-visualization sum-slr.qzv

# pca must accept BOTH table orientations: a raw FeatureTable is
# features x samples, while transform-features output is samples x features.
qiime gglasso pca \
  --i-table "$GG/data/atacama-counts.qza" --i-solution sol-slr.qza \
  --m-sample-metadata-file "$GG/data/selected-atacama-sample-metadata.tsv" \
  --p-n-components 2 --o-visualization pca.qzv

qiime gglasso pca \
  --i-table mclr.qza --i-solution sol-slr.qza \
  --m-sample-metadata-file "$GG/data/selected-atacama-sample-metadata.tsv" \
  --p-n-components 2 --o-visualization pca-mclr.qzv
set +x

echo "--- pca on a rank-0 solution must fail with a clear message, not matmul ---"
if out=$(qiime gglasso pca \
      --i-table "$GG/data/atacama-counts.qza" --i-solution sol-rank0.qza \
      --m-sample-metadata-file "$GG/data/selected-atacama-sample-metadata.tsv" \
      --p-n-components 2 --o-visualization pca-rank0.qzv 2>&1); then
  echo "FAIL: pca succeeded on a rank-0 solution"; rc=1
elif grep -q "rank 0" <<<"$out"; then
  echo "  OK: clear rank-0 error"
else
  echo "FAIL: unclear error on rank-0 solution:"; echo "$out" | tail -3; rc=1
fi

echo "############ [4/4] assert the bokeh-3 fix actually landed ############"
python - sum-sgl.qzv sum-slr.qzv pca.qzv pca-mclr.qzv <<'PY'
import sys, zipfile, re
bad = []
for qzv in sys.argv[1:]:
    with zipfile.ZipFile(qzv) as z:
        pages = [n for n in z.namelist() if n.endswith("/data/index.html")]
        if not pages:
            bad.append(f"{qzv}: no index.html in artifact"); continue
        for p in pages:
            html = z.read(p).decode("utf-8", "replace")
            if "bokeh-2." in html:
                bad.append(f"{qzv}: stale bokeh 2 reference")
            if "cdn.bokeh.org" in html:
                bad.append(f"{qzv}: still loads bokeh from a CDN (not self-contained)")
            # An inlined bokeh 3 bundle must actually be present, not just absent
            # of the old one -- a blank page would otherwise pass.
            if "Bokeh" not in html:
                bad.append(f"{qzv}: no Bokeh runtime embedded at all")
            m = re.search(r"Bokeh\.version\s*=\s*['\"]([0-9.]+)", html)
            if m:
                print(f"  {qzv}: embedded Bokeh {m.group(1)}")
if bad:
    print("\nFAIL:"); [print("  " + b) for b in bad]; sys.exit(1)
print("  all visualizations carry a self-contained bokeh 3 runtime")
PY
[[ $? -eq 0 ]] || rc=1

echo
echo "############ RESULT: $([[ $rc -eq 0 ]] && echo PASS || echo FAIL) ############"
exit $rc
