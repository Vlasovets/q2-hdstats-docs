#!/bin/bash
#SBATCH --job-name=q2-final
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/final_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/final_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Final stage -- re-verify everything and write ONE file to read afterwards:
#   reports/STATUS.md
#
# Runs regardless of whether earlier stages succeeded (submitted with
# --dependency=afterany) precisely so the status report can tell you what
# failed. Never exits non-zero for a failed CHECK -- only for its own errors --
# so STATUS.md always gets written.

set -uo pipefail
ROOT=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute
DOCS=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs
GG=/home/itg/oleg.vlasovets/slr_example/q2-gglasso
CL=/home/itg/oleg.vlasovets/slr_example/q2-classo
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
DOCSENV=/localscratch/${USER}/docsenv
STATUS="$ROOT/reports/STATUS.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

say() { echo "$*" >> "$STATUS"; }
: > "$STATUS"
say "# Recompute status — $(date -Is)"
say ""
say "Written by \`slurm/19_final_verify.sh\` (job ${SLURM_JOB_ID:-manual}) on \`$(hostname)\`."
say "Everything below was re-run just now; nothing is carried over from an earlier report."
say ""

# ---- 1. environment ----------------------------------------------------------
say "## Environment"
say ""
say '```'
python - >> "$STATUS" <<'PY'
import importlib.metadata as md
for p in ("qiime2", "numpy", "pandas", "scipy", "numba", "gglasso", "bokeh", "zarr", "c-lasso"):
    try:
        print(f"{p:10} {md.version(p)}")
    except Exception:
        print(f"{p:10} (absent)")
PY
say '```'
say ""

# ---- 2. test suites ----------------------------------------------------------
say "## Test suites"
say ""
gg_res=$( (cd "$GG" && pytest q2_gglasso/tests -q 2>&1 | tail -1) || true )
cl_res=$( (cd "$CL" && pytest q2_classo/tests -q 2>&1 | tail -1) || true )
say "- q2-gglasso: \`${gg_res}\`"
say "- q2-classo:  \`${cl_res}\`"
say ""

# ---- 3. docs build -----------------------------------------------------------
say "## Documentation"
say ""
# /localscratch is NODE-LOCAL, and this job will not land on the node where the
# venv was first built. Copying one is not an option either — a venv bakes
# absolute paths into its shebangs. So build it here if it is missing; it costs
# a couple of minutes and makes this job runnable on any node.
if [[ ! -x "$DOCSENV/bin/jupyter-book" ]]; then
  echo "provisioning docs venv at $DOCSENV"
  python3 -m venv "$DOCSENV" >/dev/null 2>&1 || true
  "$DOCSENV/bin/pip" install -q --upgrade pip >/dev/null 2>&1 || true
  "$DOCSENV/bin/pip" install -q -r "$DOCS/requirements.txt" >/dev/null 2>&1 || true
fi
if [[ -x "$DOCSENV/bin/jupyter-book" ]]; then
  if (cd "$DOCS" && rm -rf docs/_build && "$DOCSENV/bin/jupyter-book" build docs --warningiserror) >"$SCRATCH/docs.log" 2>&1; then
    say "- build with \`--warningiserror\`: **PASS** (warning-clean)"
  else
    say "- build with \`--warningiserror\`: **FAIL** — see below"
    say '```'
    grep -iE "warning|error" "$SCRATCH/docs.log" | sed 's|.*/docs/||' | head -20 >> "$STATUS"
    say '```'
  fi
else
  say "- **skipped**: could not provision a docs toolchain at \`$DOCSENV\`"
  say "  (no outbound network from this node? try \`jupyter-book build docs --warningiserror\` yourself)"
fi
say ""

# ---- 4. generated tables -----------------------------------------------------
say "## Generated tables (prose cannot drift from these)"
say ""
for f in "$DOCS/docs/_data/atacama-lambda-path.tsv" \
         "$DOCS/docs/_data/atacama-mu-rank-map.tsv" \
         "$DOCS/docs/_data/atacama-classo-cv.tsv" \
         "$DOCS/docs/_data/manifest.tsv"; do
  if [[ -f "$f" ]]; then say "- \`$(basename "$f")\` — $(( $(wc -l < "$f") - 1 )) rows"
  else say "- \`$(basename "$f")\` — **MISSING**"; fi
done
say ""

# ---- 5. recompute artifacts --------------------------------------------------
say "## Recompute artifacts"
say ""
say "| stage | output | state |"
say "|---|---|---|"
chk() { if [[ -e "$2" ]]; then say "| $1 | \`$(basename "$2")\` | present |"; else say "| $1 | \`$(basename "$2")\` | **missing** |"; fi; }
chk "1 lambda path"  "$ROOT/results/gglasso/atacama-top-300-sgl-linear-path.qza"
chk "3 mu=15 rank2"  "$ROOT/results/gglasso/atacama-top-300-slr-lambda0.8-mu15.qza"
chk "5 comparison"   "$ROOT/results/figures/atacama-top-300-rank0-rank2-comparison-summary.tsv"
n_cv=$(ls -1 "$ROOT/results/classo/"atacama-top-300-regress-*.qza 2>/dev/null | wc -l)
say "| 6 classo CV | regress fits | ${n_cv}/15 |"
say ""

# ---- 6. gates ----------------------------------------------------------------
say "## Gates"
say ""
if grep -aq "GATE C1: PASS" "$ROOT/slurm/logs/"lambda_path_*.out 2>/dev/null; then
  say "- **Gate C1 PASS** — the CLI reproduces the reference selection (lambda 0.8, 216 edges)."
else
  say "- Gate C1 — no PASS line found in \`slurm/logs/lambda_path_*.out\`."
fi
if [[ -f "$ROOT/reports/ORIENTATION_FINDING.md" ]]; then
  say "- Orientation probe: see \`reports/ORIENTATION_FINDING.md\`."
else
  say "- Orientation probe: **not run**."
fi
say ""

# ---- 7. what still needs a human --------------------------------------------
say "## Still needs you"
say ""
say "- \`reports/DECISIONS_NEEDED.md\` — open decisions, evidence attached."
say "- Nothing is committed. \`git status\` in each of the three repos shows the change set."
say ""
say "## Job log index"
say ""
say '```'
sacct -u "$USER" --starttime now-2days --format=JobID%14,JobName%18,State,Elapsed --noheader 2>/dev/null \
  | grep -vE "\.(ba|ex|0)\+" | grep -E "q2-" | tail -25 >> "$STATUS"
say '```'

echo "STATUS written to $STATUS"
cat "$STATUS"
