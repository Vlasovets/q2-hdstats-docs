#!/bin/bash
#
# Submit the remaining recompute as one dependency-chained pipeline, then log out.
#
#   bash slurm/submit_all.sh
#
# When you come back, read ONE file:
#   reports/STATUS.md
#
# It is written by the last job and reports the environment, both test suites,
# the docs build, every generated table, which artifacts exist, the gate results
# and what still needs a decision. It is written even if earlier stages fail —
# the final job is chained with `afterany` precisely so it can tell you what
# broke rather than vanishing along with it.
#
# Chain:
#   10 orientation probe  ─┐
#   11 classo CV (array)  ─┼─→ 12 summary ─→ 19 final verify + STATUS.md
#
# 10 and 11 are independent and run concurrently. 12 waits for 11 to succeed
# (it has nothing to summarise otherwise). 19 waits for both, whatever happens.

set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
mkdir -p slurm/logs reports results/{gglasso,classo,figures,tables}

echo "submitting from $ROOT"

# --- preflight: fail fast on missing inputs rather than 15 array tasks later ---
missing=0
for f in data/atacama-top-300-clr.qza \
         data/atacama-classo-outcomes-mean-imputed.tsv \
         data/atacama-top-300-correlation.qza; do
  [[ -f "$f" ]] || { echo "  MISSING INPUT: $f" >&2; missing=1; }
done
[[ -d /home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr ]] || {
  echo "  MISSING ENV: q2-2026.7-slr" >&2; missing=1; }
(( missing == 0 )) || { echo "preflight failed; nothing submitted" >&2; exit 1; }
echo "  preflight OK"

J10=$(sbatch --parsable slurm/10_orientation_probe.sh)
echo "  10 orientation probe : $J10"

J11=$(sbatch --parsable slurm/11_classo_cv.sh)
echo "  11 classo CV (1-15)  : $J11"

# afterok: summarising an array that failed would produce a misleading table.
J12=$(sbatch --parsable --dependency=afterok:"$J11" slurm/12_classo_summary.sh)
echo "  12 classo summary    : $J12  (afterok:$J11)"

# afterany: this one must run even when something upstream failed, so that
# STATUS.md exists and says so.
J19=$(sbatch --parsable --dependency=afterany:"$J10":"$J11":"$J12" slurm/19_final_verify.sh)
echo "  19 final verify      : $J19  (afterany:$J10,$J11,$J12)"

cat <<EOF

Submitted. Safe to log out.

  watch:   squeue -u \$USER
  read:    $ROOT/reports/STATUS.md
  logs:    $ROOT/slurm/logs/

Every stage is idempotent — resubmitting skips work that already has output, so
if a task dies you can just run this again.
EOF
