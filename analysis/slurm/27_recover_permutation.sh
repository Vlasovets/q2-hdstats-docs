#!/bin/bash
#SBATCH --job-name=q2-recover
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/recover_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/recover_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Job 26 established that the shipped and regenerated correlation matrices share
# a spectrum to 2.1e-14, an identical Frobenius norm and an identical trace, and
# that all 300 rows match by value-multiset. They are the same matrix, permuted.
#
# Job 26 nevertheless failed to RECOVER the permutation, because it resolved ties
# greedily (`cands[0]`) and rows with identical value-multisets collided -- there
# are 130 pairs at r >= 0.99 in this matrix, so collisions were guaranteed.
#
# This job recovers the permutation properly using colour refinement (1-WL on a
# weighted graph): start each node coloured by its sorted row, then repeatedly
# recolour by the multiset of (edge weight, neighbour colour). Continuous weights
# individualise every node within a couple of rounds. Then VERIFY the recovered
# permutation reproduces the matrix exactly, and identify which ordering
# principle it corresponds to.
#
# Also saves the regenerated matrix to results/ so this stops being re-derived.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/PERMUTATION_RECOVERED.md"
OUTDIR="$ROOT/results/regen"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-recover/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-recover/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports" "$OUTDIR"
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

cd "$SCRATCH"

if [[ -f "$OUTDIR/atacama-top-300-correlation-regen.qza" ]]; then
  echo "reusing cached regeneration"
  cp "$OUTDIR/atacama-top-300-correlation-regen.qza" corr-regen.qza
  cp "$OUTDIR/atacama-top-300-clr-regen.qza" clr.qza
else
  qiime gglasso transform-features \
    --i-table "$ROOT/data/atacama-top-300-table.qza" \
    --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
    --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
    --p-transformation clr --p-pseudo-count 1 --p-no-keep-original-id \
    --o-transformed-table clr.qza >/dev/null
  qiime gglasso calculate-covariance \
    --i-table clr.qza --p-method scaled --o-covariance-matrix corr-regen.qza >/dev/null
  cp corr-regen.qza "$OUTDIR/atacama-top-300-correlation-regen.qza"
  cp clr.qza        "$OUTDIR/atacama-top-300-clr-regen.qza"
  echo "regenerated and cached to $OUTDIR"
fi

python - "$ROOT/data/atacama-top-300-correlation.qza" corr-regen.qza \
         "$ROOT/data/atacama-top-300-table.qza" "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd

ship_q, regen_q, raw_q, report = sys.argv[1:5]

def load_tsv(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        n=[x for x in z.namelist() if x.endswith(".tsv") and "/data/" in x][0]
        z.extract(n,t)
        return pd.read_csv(os.path.join(t,n), sep="\t", index_col=0)

A = load_tsv(ship_q).to_numpy(float)     # shipped
B = load_tsv(regen_q).to_numpy(float)    # regenerated
n = A.shape[0]

# ---- colour refinement (1-WL on a weighted complete graph) ----------------
# Rounding is load-bearing: without it float noise makes every colour unique in
# round 0 and the refinement is meaningless. 9 dp is far below the 1.147 signal
# and far above float64 noise on values of order 1.
R = 9
def refine(M):
    col = [hash(tuple(np.round(np.sort(M[i, :]), R))) for i in range(n)]
    for _ in range(5):
        new = []
        for i in range(n):
            sig = tuple(sorted((round(float(M[i, j]), R), col[j]) for j in range(n)))
            new.append(hash((col[i], sig)))
        if len(set(new)) == len(set(col)):
            col = new; break
        col = new
    return col

cA, cB = refine(A), refine(B)
uniqA, uniqB = len(set(cA)), len(set(cB))

perm, ambiguous = [None]*n, 0
idxA = {}
for i, c in enumerate(cA):
    idxA.setdefault(c, []).append(i)
for j, c in enumerate(cB):
    cands = idxA.get(c, [])
    if len(cands) == 1:
        perm[j] = cands[0]
    else:
        ambiguous += 1

resolved = all(p is not None for p in perm) and len(set(perm)) == n
resid = None
if resolved:
    p = np.array(perm)
    resid = float(np.max(np.abs(A[np.ix_(p, p)] - B)))

L = ["# Recovering the permutation between the two correlation matrices", "",
     "Job 26 showed the two matrices are spectrally identical "
     "(max |eigenvalue difference| 2.1e-14) but disagree entrywise by 1.147. "
     "That is the signature of a reordering. This recovers it.", "",
     "## Colour refinement", "",
     f"- distinct colours, shipped     : {uniqA} / {n}",
     f"- distinct colours, regenerated : {uniqB} / {n}",
     f"- nodes left ambiguous          : {ambiguous}", ""]

if resolved:
    n_moved = int((p != np.arange(n)).sum())
    L += ["## Verification", "",
          f"- `A[perm][:,perm] == B` to : **{resid:.3e}**",
          f"- features whose position changes : **{n_moved} / {n}**", ""]
    if resid < 1e-9:
        L += ["**Exact.** The shipped matrix and the regenerated matrix are the "
              "same matrix with the features in a different order. There is no "
              "numerical discrepancy between them at all.", ""]

    # ---- which ordering principle does the permutation correspond to? -----
    try:
        import qiime2, biom
        raw = qiime2.Artifact.load(raw_q).view(biom.Table).to_dataframe()  # feat x samp
        tot = raw.sum(axis=1)
        prev = (raw > 0).sum(axis=1)
        feats = list(raw.index)
        cand = {
            "raw table order (identity)": list(range(n)),
            "total abundance, descending": list(np.argsort(-tot.to_numpy(), kind="stable")),
            "total abundance, ascending":  list(np.argsort(tot.to_numpy(), kind="stable")),
            "prevalence, descending":      list(np.argsort(-prev.to_numpy(), kind="stable")),
            "feature ID, lexicographic":   list(np.argsort(np.array(feats), kind="stable")),
        }
        L += ["## What ordering does the shipped matrix use?", "",
              "`perm[j] = i` means regenerated row *j* is shipped row *i*. Testing "
              "the permutation against orderings the original author might have "
              "applied to the table before transforming:", "",
              "| candidate ordering | matches |", "|---|---|"]
        for name, order in cand.items():
            L.append(f"| {name} | {'**yes**' if list(p) == list(order) else 'no'} |")
        L += ["", f"- permutation is its own inverse (a pure swap set): "
                  f"{bool(np.array_equal(p[p], np.arange(n)))}",
              f"- first 8 mappings (regen -> ship): "
              f"{[(j, int(p[j])) for j in range(8)]}", ""]
    except Exception as e:
        L += [f"(ordering-principle check skipped: {e})", ""]
else:
    L += ["## Verification", "",
          f"- permutation NOT fully resolved ({ambiguous} ambiguous nodes). "
          "Colour refinement did not individualise every feature, which happens "
          "when genuinely duplicate rows exist. The spectral evidence for "
          "permutation-equivalence from job 26 still stands.", ""]

L += ["## Consequence", "",
      "The graphical-lasso objective is invariant under simultaneous row/column "
      "permutation, which is why the lambda path, the eBIC at all 15 grid points "
      "and the edge count (216) are bit-identical on the two matrices. The "
      "published *numbers* are unaffected.", "",
      "What IS affected is feature identity: `--p-no-keep-original-id` assigns "
      "`ASV-k` by position, so `ASV-k` denotes different organisms in the two "
      "artifacts. Any statement mapping an `ASV-k` to a taxon is valid only for "
      "the artifact it was derived from.", "",
      "## Verdict", "",
      f"**{'Same matrix, reordered — the 1.147 discrepancy is not a numerical defect' if (resolved and resid is not None and resid < 1e-9) else 'Permutation-equivalent by spectrum; exact permutation not pinned'}**"]

open(report, "w").write("\n".join(L) + "\n"); print("\n".join(L))
PY
echo "report -> $REPORT"
