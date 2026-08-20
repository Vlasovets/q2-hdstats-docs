#!/bin/bash
#SBATCH --job-name=q2-oldvsnew
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/oldvsnew_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/oldvsnew_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Isolate the orientation fix from every other variable.
#
# Job 21 compared the FIXED chain against the SHIPPED correlation matrix and got
# max|diff| = 1.147. That does not implicate the fix: the shipped matrix came out
# of a bundle whose exact transform parameters are unknown, so the comparison had
# an uncontrolled variable.
#
# The controlled experiment is old code vs new code on identical inputs, run back
# to back by checking out each revision. The two paths should be bit-identical:
#   np.cov(A, rowvar=True)  for A of shape (p, N)
#   np.cov(A.T, rowvar=False)
# are the same computation by definition. If they differ, my reasoning about the
# cancellation is wrong and the fix must not ship.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
GG=/home/itg/oleg.vlasovets/slr_example/q2-gglasso
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/ORIENTATION_FIX_VERIFICATION.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-oldvsnew/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-oldvsnew/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

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
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

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


START_REF=$(git -C "$GG" rev-parse --abbrev-ref HEAD)
restore() { git -C "$GG" checkout -q "$START_REF" 2>/dev/null || true; }
trap 'restore; [[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

run_chain () {                      # $1 = label
  local tag="$1"
  ( cd "$SCRATCH"
    qiime gglasso transform-features \
      --i-table "$ROOT/data/atacama-top-300-table.qza" \
      --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
      --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
      --p-transformation clr --p-pseudo-count 1 --p-no-keep-original-id \
      --o-transformed-table "clr-$tag.qza" >/dev/null
    qiime gglasso calculate-covariance \
      --i-table "clr-$tag.qza" --p-method scaled \
      --o-covariance-matrix "corr-$tag.qza" >/dev/null )
  echo "  $tag chain done"
}

echo "=== NEW code (working tree, orientation fixed) ==="
run_chain new

echo "=== OLD code (christian-review-fixes, pre-fix) ==="
git -C "$GG" stash -q -u 2>/dev/null || true
git -C "$GG" checkout -q christian-review-fixes
run_chain old
git -C "$GG" checkout -q "$START_REF"
git -C "$GG" stash pop -q 2>/dev/null || true
echo "  restored $START_REF"

echo "=== compare ==="
python - "$SCRATCH/corr-old.qza" "$SCRATCH/corr-new.qza" \
         "$ROOT/data/atacama-top-300-correlation.qza" "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd

old_q, new_q, ship_q, report = sys.argv[1:5]

def load(qza):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(qza) as z:
        n = [x for x in z.namelist()
             if x.endswith((".csv", ".tsv", ".txt")) and "/data/" in x
             and "MANIFEST" not in x][0]
        z.extract(n, t); p = os.path.join(t, n)
        return pd.read_csv(p, sep="\t" if p.endswith((".tsv", ".txt")) else ",",
                           index_col=0)

old, new, ship = load(old_q), load(new_q), load(ship_q)

def diff(a, b):
    common = [i for i in a.index if i in b.index]
    if not common:
        return None, 0
    A = a.loc[common, common].to_numpy(float)
    B = b.loc[common, common].to_numpy(float)
    return float(np.nanmax(np.abs(A - B))), len(common)

d_on, n_on = diff(old, new)
d_sn, n_sn = diff(ship, new)
d_so, n_so = diff(ship, old)

L = ["# Orientation fix — controlled verification", "",
     "Old code vs new code, **same inputs, same parameters, run back to back**.",
     "The earlier comparison against the shipped matrix had an uncontrolled",
     "variable: the parameters that produced it are not recorded anywhere.", "",
     "| comparison | matched labels | max abs difference |",
     "|---|---|---|",
     f"| old code vs new code | {n_on} | **{d_on:.3e}** |",
     f"| shipped vs new code | {n_sn} | {d_sn:.3e} |",
     f"| shipped vs old code | {n_so} | {d_so:.3e} |", ""]

fix_is_neutral = d_on is not None and d_on < 1e-9
ship_differs_from_both = d_sn > 1e-6 and d_so > 1e-6

if fix_is_neutral:
    L += ["**The fix does not change the covariance.** Old and new agree to "
          f"{d_on:.1e}, which is exact. That is the expected result: "
          "`np.cov(A, rowvar=True)` for A of shape (p, N) and "
          "`np.cov(A.T, rowvar=False)` are the same computation, so removing the "
          "transpose and flipping `rowvar` cancel by construction.", ""]
    if ship_differs_from_both:
        L += ["**The shipped matrix differs from BOTH** by the same margin, so it "
              "was not produced by this chain at these parameters — a difference "
              "in the transform (clr vs mclr, pseudo-count, metadata columns), "
              "not something the fix introduced. Gate C1 consumes the shipped "
              "matrix directly and is therefore unaffected either way; but the "
              "tier-2 chapter should stop implying that "
              "`transform-features -> calculate-covariance` regenerates it until "
              "the original parameters are known.", ""]
else:
    L += [f"**The fix CHANGES the covariance** (max abs {d_on:.3e}). The two "
          "edits do not cancel, so the reasoning behind them is wrong. "
          "Do not merge.", ""]

L += ["## Verdict", "", f"**{'PASS — safe to merge' if fix_is_neutral else 'FAIL — do not merge'}**"]
open(report, "w").write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(0 if fix_is_neutral else 1)
PY
