#!/bin/bash
#SBATCH --job-name=q2-mu-rank
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/mu_rank_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/mu_rank_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Stage 3/4 -- the mu1 -> rank map, at the Gate-C1-selected lambda = 0.8.
#
# WHY: two chapters currently disagree. 04_highdim_atacama/02_model_selection.md
# states the map as established fact (rank 2 @ mu1=15 -> 202 edges / 162 nodes;
# rank 5 @ mu1=10 -> 158 / 124; rank 10 @ mu1=7.5 -> 110 / 92), while
# 03_slr_ranks.md marks every one of those cells "pending". One page is wrong and
# no amount of reading source can say which.
#
# NOTE the --p-lambda2-* triple. Without it lambda2 expands to the 5-point
# default np.logspace(-1,-4,5) and the "single fit" silently becomes a
# model-selection run over lambda2 -- which is exactly the defect found in the
# chapter's own commands.

set -euo pipefail

REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
IN_COV="$ROOT/data/atacama-top-300-correlation.qza"
OUT="$ROOT/results/gglasso"
TAB="$ROOT/results/tables"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break
  fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; export TMP="$TMPDIR"; export TEMP="$TMPDIR"
mkdir -p "$TMPDIR" "$OUT" "$TAB"
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


LAMBDA=0.8
for MU in 15 10 7.5; do
  q="$OUT/atacama-top-300-slr-lambda${LAMBDA}-mu${MU}.qza"
  if [[ -f "$q" ]]; then echo "[skip] $q"; continue; fi
  echo "=== solving lambda=$LAMBDA mu1=$MU ==="
  qiime gglasso solve-problem \
    --i-covariance-matrix "$IN_COV" \
    --p-n-samples 54 \
    --p-latent \
    --p-lambda1-min "$LAMBDA" --p-lambda1-max "$LAMBDA" --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min "$MU" --p-mu1-max "$MU" --p-n-mu1 1 \
    --o-solution "$q"
done

python - "$OUT" "$TAB/mu-rank-map.tsv" "$LAMBDA" <<'PY'
import sys, os, glob, zipfile, tempfile, csv
import numpy as np, zarr

out_dir, out_tsv, lam = sys.argv[1], sys.argv[2], sys.argv[3]

def load(qza):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(qza) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")][0]
        z.extract(inner, t)
        r = zarr.open(zarr.ZipStore(os.path.join(t, inner), mode="r"))
        return (np.asarray(r["solution/precision_"]),
                np.asarray(r["solution/lowrank_"]),
                "modelselect_stats" in r)

rows = []
for mu in ("15", "10", "7.5"):
    q = os.path.join(out_dir, f"atacama-top-300-slr-lambda{lam}-mu{mu}.qza")
    P, L, has_ms = load(q)
    off = np.triu(np.abs(P), 1)
    edges = int((off > 1e-8).sum())
    # a node counts as connected if it has at least one surviving edge
    adj = (np.abs(P) > 1e-8) & ~np.eye(P.shape[0], dtype=bool)
    nodes = int(adj.any(axis=1).sum())
    rank = int(np.linalg.matrix_rank(L))
    rows.append((mu, rank, edges, nodes, has_ms))
    print(f"  mu1={mu:>4}  rank={rank:<3} edges={edges:<5} nodes={nodes:<4} "
          f"single_fit={'no (MODEL SELECTION RAN)' if has_ms else 'yes'}")

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["mu1", "achieved rank", "sparse edges", "connected nodes"])
    for mu, rank, edges, nodes, _ in rows:
        w.writerow([mu, rank, edges, nodes])
print(f"\nwrote {out_tsv}")

# The claim currently asserted as fact in 02_model_selection.md.
claimed = {"15": (2, 202, 162), "10": (5, 158, 124), "7.5": (10, 110, 92)}
print("\n--- against the table asserted in 02_model_selection.md ---")
ok = True
for mu, rank, edges, nodes, _ in rows:
    cr, ce, cn = claimed[mu]
    match = (rank, edges, nodes) == (cr, ce, cn)
    ok &= match
    print(f"  mu1={mu:>4}: got rank/edges/nodes = {rank}/{edges}/{nodes}"
          f"   claimed {cr}/{ce}/{cn}   {'MATCH' if match else 'MISMATCH'}")
print("\nVERDICT:", "the asserted table reproduces" if ok else
      "the asserted table does NOT reproduce -- 02_model_selection.md must be corrected "
      "from the generated mu-rank-map.tsv, and 03_slr_ranks.md's 'pending' was right")
PY
