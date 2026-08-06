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

VENV="$SCRATCH/venv"
echo "=== [1/3] clean venv (mirrors CI's actions/setup-python) ==="
python3 -m venv "$VENV"
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
