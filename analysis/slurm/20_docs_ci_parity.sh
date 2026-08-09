#!/bin/bash
#SBATCH --job-name=q2-docs-ci
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/docs_ci_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/docs_ci_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Reproduce the docs CI job EXACTLY, from a throwaway environment.
#
# Why this exists: three times in this migration a docs/CI change passed locally
# and would have failed in CI, because the local toolchain was provisioned
# differently from the declared one. Most recently a bare `jupyter-book`
# requirement resolved to 2.1.6 -- a different tool that cannot read this book's
# v1 config -- while the local venv had 1.0.4 from an explicit `<2` install.
#
# The ONLY way to test the CI path is to build the toolchain the way CI does:
# a clean venv, `pip install -r requirements.txt`, nothing pre-existing. Reusing
# a warm venv defeats the entire point, so this deliberately starts from empty.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
DOCS="$REPO"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-docs-ci/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-docs-ci/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# Match CI's interpreter, or say plainly that we could not.
#
# This block used a bare `python3 -m venv`. On this cluster that is Python
# 3.9.21, while CI pins 3.11 through actions/setup-python -- so the stage was
# named "CI parity", printed "mirrors CI's actions/setup-python", and was in
# fact testing a different interpreter than CI on every run. A harness that
# reports parity it does not have is worse than no harness: its FAIL sends you
# hunting a phantom, and its PASS is not evidence CI will pass.
#
# The version is read out of the workflow file so the two cannot drift again.
CI_PY=$(grep -A2 'uses: actions/setup-python' "$REPO/.github/workflows/ci.yml" \
        | grep -oP 'python-version:\s*"\K[0-9.]+' | head -1)
[[ -n "$CI_PY" ]] || { echo "cannot read python-version from ci.yml"; exit 1; }

if command -v "python$CI_PY" >/dev/null 2>&1; then
  PYBIN="python$CI_PY"; PARITY="true"
else
  PYBIN="python3"; PARITY="false"
fi
echo "=== [1/3] clean venv ==="
echo "  CI pins python $CI_PY"
echo "  using: $PYBIN ($($PYBIN -V 2>&1))"
if [[ "$PARITY" != "true" ]]; then
  cat <<EOF

  NOT CI PARITY. python$CI_PY is not installed on this node, so this run uses a
  different interpreter than CI. Treat the result as a smoke test: a FAIL here
  may be an interpreter artefact, and a PASS does not prove CI will pass. The
  authoritative check is the build job on GitHub Actions, which runs on every
  push and every PR to main.

EOF
fi
VENV="$SCRATCH/venv"
"$PYBIN" -m venv "$VENV"
"$VENV/bin/pip" install --upgrade pip -q

echo "=== [2/3] pip install -r requirements.txt (exactly what CI runs) ==="
"$VENV/bin/pip" install -r "$DOCS/requirements.txt" 2>&1 | tail -3
echo "--- resolved toolchain ---"
"$VENV/bin/pip" list 2>/dev/null | grep -iE "^(jupyter-book|sphinx|sphinxcontrib-bibtex|sphinxext-rediraffe) " | sed 's/^/  /'

echo "=== [3/3] jupyter-book build docs --warningiserror ==="
cd "$DOCS"
rm -rf docs/_build
if "$VENV/bin/jupyter-book" build docs --warningiserror > "$SCRATCH/build.log" 2>&1; then
  echo "DOCS CI PARITY: PASS"
  rc=0
else
  echo "DOCS CI PARITY: FAIL"
  echo "--- last 25 lines ---"
  sed -e 's/\x1b\[[0-9;]*m//g' "$SCRATCH/build.log" | tail -25
  rc=1
fi

# Record it where the status report can find it.
{
  echo "# Docs CI parity — $(date -Is)"
  echo
  echo "Built from a THROWAWAY venv via \`pip install -r requirements.txt\`, i.e."
  echo "exactly the toolchain CI resolves — not a warm local environment."
  echo
  echo '```'
  "$VENV/bin/pip" list 2>/dev/null | grep -iE "^(jupyter-book|sphinx|sphinxext-rediraffe) "
  echo '```'
  echo
  echo "Result: **$( [[ $rc -eq 0 ]] && echo PASS || echo FAIL )**"
} > "$ROOT/reports/DOCS_CI_PARITY.md"

exit $rc
