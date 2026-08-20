#!/bin/bash
#SBATCH --job-name=q2-env-build
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/env_build_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/env_build_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Phase A1 / Gate A1 -- build the QIIME 2 2026.7 environment and immediately test
# whether numba 0.66 can compile GGLasso's JIT kernels under numpy 2.4.
#
# Extras (gglasso, bokeh, zarr, numcodecs, plotly) are appended to a COPY of the
# distribution lockfile so they solve in one transaction. A second `conda install`
# is what lets the solver drift the distribution's pinned numpy/pandas.
#
# c-lasso is genuinely not on conda-forge, so it goes in by pip --no-deps afterwards.

set -euo pipefail

REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
GGLASSO=/home/itg/oleg.vlasovets/slr_example/q2-gglasso
CLASSO=/home/itg/oleg.vlasovets/slr_example/q2-classo

# ---- scratch -----------------------------------------------------------------
# Never let conda or pip spill into the node's shared /tmp. Unpacking ~680
# packages plus pip wheel builds is several GB of churn, and a full /tmp takes
# down everything else running on the node.
#
# Prefer node-local disk (/localscratch is ~330 GB on these nodes and is real
# local storage, not Lustre), then $LOCAL_SCRATCH if the site ever sets it, then
# per-job scratch on Lustre as a last resort.
#
# This matters: unpacking ~680 conda packages is tens of thousands of small file
# creates, which is exactly the access pattern Lustre is worst at.
for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"
    break
  fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
mkdir -p "$SCRATCH"

export TMPDIR="$SCRATCH/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
export PIP_CACHE_DIR="$SCRATCH/pip-cache"
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"

# NB: CONDA_PKGS_DIRS is deliberately NOT redirected to scratch. ~/.condarc
# already points it at persistent storage (/lustre/groups/shared/apps/conda/pkgs,
# then ~/.conda/pkgs); moving it into a scratch tree we delete on exit would
# force a full re-download of ~680 packages on every run.

# Remove the scratch tree on exit unless KEEP_SCRATCH=1 (for debugging).
cleanup() {
  if [[ "${KEEP_SCRATCH:-0}" != "1" ]]; then
    rm -rf "$SCRATCH"
  else
    echo "KEEP_SCRATCH=1 -- leaving $SCRATCH in place"
  fi
}
trap cleanup EXIT

echo "scratch:  $SCRATCH"
echo "TMPDIR:   $TMPDIR"
df -h "$SCRATCH" | tail -1

mkdir -p "$ROOT/slurm/logs"

if [[ -d "$PREFIX" ]]; then
  echo "Environment already exists at $PREFIX -- remove it first to rebuild." >&2
  exit 1
fi

echo "=== [1/5] conda env create ==="
"$CONDA" env create -p "$PREFIX" -f "$ROOT/envs/q2-slr-qiime2-2026.7.yml"

# The conda env is not the env that runs unless we say so.
#
# ~/.local/lib/python3.12/site-packages shadows this env for 19 packages, and three
# of them are pins the migration declared load-bearing:
#     zarr      2.18.7 -> 3.1.5   (zarr.hierarchy.Group STOPS RESOLVING, and both
#                                  plugins annotate their transformers with it)
#     numpy     2.4.2  -> 2.2.6   (Gate A1 was validated on 2.4.2)
#     numcodecs 0.15.1 -> 0.16.5  (the pin is <0.16)
#
# ~/.local is not ours to clean -- other projects on this account depend on it -- so the
# stages opt out of user site-packages instead. Must be exported BEFORE activation so
# the interpreter never builds a user-site path.
export PYTHONNOUSERSITE=1

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"
conda activate "$PREFIX"

# Assert we got the env the lockfile describes.
#
# This check exists because the shadowing went unnoticed: a silent version swap raises
# no error, it just changes the numbers, while the lockfile keeps saying the right
# thing. Runs AFTER activation, so it tests the env's interpreter and not the login
# node's python3.
python - <<'ENVCHECK' || exit 9
import sys
import numpy, zarr
bad = []
if not zarr.__version__.startswith("2.18"):
    bad.append("zarr %s, expected 2.18.x" % zarr.__version__)
if not numpy.__version__.startswith("2.4"):
    bad.append("numpy %s, expected 2.4.x" % numpy.__version__)
try:
    import zarr.hierarchy  # noqa: F401
except ImportError:
    bad.append("zarr.hierarchy missing -- the plugin transformers cannot load")
if bad:
    sys.stderr.write("ENV CHECK FAILED: %s\n" % "; ".join(bad))
    sys.stderr.write("  PYTHONNOUSERSITE=1 should pin this env. If it did not, something\n")
    sys.stderr.write("  in ~/.local or PYTHONPATH still wins. Do not trust these results.\n")
    sys.exit(9)
ENVCHECK


echo "=== [2/5] verify the solver did NOT drift the distribution pins ==="
python - <<'PY'
import sys, importlib.metadata as md
expect = {"numpy": "2.4.2", "pandas": "2.3.3", "scipy": "1.17.1",
          "numba": "0.66.0", "scikit-learn": "1.7.1"}
bad = []
for pkg, want in expect.items():
    got = md.version(pkg)
    print(f"  {pkg}: {got} (expected {want})")
    if got != want:
        bad.append(f"{pkg} {got} != {want}")
if bad:
    sys.exit("PIN DRIFT: " + "; ".join(bad))
print("  pins intact")
PY

echo "=== [3/5] pip extras (--no-deps: pip cannot see conda pins) ==="
pip install --no-deps 'c-lasso==1.0.11'
pip install --no-deps 'versioneer==0.29' 'versioningit==3.3.0'

echo "=== [4/5] GATE A1: does numba 0.66 compile GGLasso's JIT kernels under numpy 2.4? ==="
python - <<'PY'
import sys, time, traceback
import numpy as np
from gglasso.problem import glasso_problem
from gglasso.helper.data_generation import generate_precision_matrix, sample_covariance_matrix

print("  numpy", np.__version__)
import gglasso, numba
print("  gglasso", gglasso.__version__, "| numba", numba.__version__)

p, N = 30, 100
Sigma, Theta = generate_precision_matrix(p=p, M=2, style="erdos", prob=0.1, seed=1234)
S, _ = sample_covariance_matrix(Sigma, N)

try:
    t0 = time.time()
    P = glasso_problem(S, N=N, reg_params={"lambda1": 0.05}, latent=False)
    P.solve()
    print(f"  SGL solve OK in {time.time()-t0:.1f}s "
          f"(first call includes JIT compile)")

    t0 = time.time()
    P = glasso_problem(S, N=N, reg_params={"lambda1": 0.05, "mu1": 5.0}, latent=True)
    P.solve()
    print(f"  latent SGL solve OK in {time.time()-t0:.1f}s")
    print("  GATE A1: PASS")
except Exception:
    traceback.print_exc()
    print("\n  GATE A1: FAIL -- numba/numpy incompatibility.")
    print("  Next step: re-run with NUMBA_DISABLE_JIT=1 to confirm it is numba,")
    print("  not numpy. The kernels are valid pure Python, just slow.")
    sys.exit(2)
PY

echo "=== [5/5] editable plugin installs + freeze ==="
pip install --no-deps --no-build-isolation -e "$GGLASSO"
pip install --no-deps --no-build-isolation -e "$CLASSO"
qiime dev refresh-cache
qiime info

"$CONDA" env export -p "$PREFIX" --no-builds \
  > "$ROOT/envs/q2-slr-qiime2-2026.7-lock.yml"

echo "DONE. Lockfile: $ROOT/envs/q2-slr-qiime2-2026.7-lock.yml"
