#!/bin/bash
#SBATCH --job-name=q2-sidebyside
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/sidebyside_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/sidebyside_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Two results are in direct contradiction and one of them must be wrong:
#
#   job 23: the lambda path on the REGENERATED correlation gives lambda = 0.8,
#           216 edges, eBIC 16130.0988 -- identical to the shipped matrix's.
#   job 24: the two correlation matrices differ substantially -- 18% of pairs by
#           >= 0.1, median 2.6e-02, and only 1% of the large differences sit on
#           the numerically degenerate pairs.
#
# Materially different covariance matrices cannot yield the same eBIC to four
# decimals. Rather than reason about which measurement to trust, solve BOTH in
# one job, print both, and diff the inputs in the same script so there is no
# question about which artifact went in.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/SIDE_BY_SIDE_PATH.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-sbs/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-sbs/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "=== regenerate the correlation ==="
qiime gglasso transform-features \
  --i-table "$ROOT/data/atacama-top-300-table.qza" \
  --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
  --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
  --p-transformation clr --p-pseudo-count 1 --p-no-keep-original-id \
  --o-transformed-table clr.qza >/dev/null
qiime gglasso calculate-covariance \
  --i-table clr.qza --p-method scaled --o-covariance-matrix corr-regen.qza >/dev/null
cp "$ROOT/data/atacama-top-300-correlation.qza" corr-ship.qza
echo "  md5 ship : $(md5sum corr-ship.qza  | cut -c1-32)"
echo "  md5 regen: $(md5sum corr-regen.qza | cut -c1-32)"

# clr pc=1, clr pc=2 and mclr all differ from the shipped matrix by the SAME
# 1.1467. mclr is a completely different transform, so a transform-parameter
# mismatch cannot be the cause. The obvious structural suspect is the sample set:
# the documented solve passes --p-n-samples 54, and nothing has ever checked that
# 54 is what the table actually holds.
echo "=== how many samples does the table actually have? ==="
python - "$ROOT/data/atacama-top-300-table.qza" clr.qza <<'PY'
import sys, qiime2, biom
for q in sys.argv[1:]:
    t = qiime2.Artifact.load(q).view(biom.Table).to_dataframe()
    print(f"  {q:16} shape={t.shape}  index[0]={t.index[0]}  cols[0]={t.columns[0]}")
PY

for which in ship regen; do
  echo "=== lambda path on corr-$which.qza ==="
  qiime gglasso solve-problem \
    --i-covariance-matrix "corr-$which.qza" \
    --p-n-samples 54 --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.30 --p-lambda1-max 1.00 --p-n-lambda1 15 \
    --p-gamma 0.3 --o-solution "path-$which.qza" >/dev/null
  echo "  solved -> path-$which.qza"
done

python - corr-ship.qza corr-regen.qza path-ship.qza path-regen.qza "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd, zarr

cs, cr, ps, pr, report = sys.argv[1:6]

def load_mat(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        n=[x for x in z.namelist() if x.endswith((".csv",".tsv",".txt"))
           and "/data/" in x and "MANIFEST" not in x][0]
        z.extract(n,t); p=os.path.join(t,n)
        return pd.read_csv(p, sep="\t" if p.endswith((".tsv",".txt")) else ",", index_col=0)

def read_path(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        ms=r["modelselect_stats"]
        return (float(np.asarray(ms["BEST/lambda1"]).item()),
                int((np.abs(np.triu(np.asarray(r["solution/precision_"]),1))>1e-8).sum()),
                np.asarray(ms["BIC"]["0.3"]).ravel())

A, B = load_mat(cs), load_mat(cr)
common=[i for i in A.index if i in B.index]
D=np.abs(A.loc[common,common].to_numpy(float)-B.loc[common,common].to_numpy(float))
iu=np.triu_indices(len(common),1); d=D[iu]

ls, es, bs = read_path(ps)
lr, er, br = read_path(pr)

L=["# Side by side: same solver, two correlation matrices","",
   "## Inputs","",
   f"- off-diagonal pairs compared : {d.size:,}",
   f"- median |difference| : {np.median(d):.3e}",
   f"- pairs differing by >= 0.1 : {int((d>=0.1).sum()):,} ({100*(d>=0.1).mean():.1f}%)",
   "",
   "## Outputs","",
   "| | shipped matrix | regenerated matrix |",
   "|---|---|---|",
   f"| selected lambda | {ls:.6g} | {lr:.6g} |",
   f"| edges | {es} | {er} |",
   f"| min eBIC (gamma=0.3) | {bs.min():.4f} | {br.min():.4f} |",
   f"| eBIC at lambda=0.8 | {bs[10]:.4f} | {br[10]:.4f} |",
   "",
   "## The whole eBIC curve, not just its minimum","",
   "A single matching minimum can be a coincidence. Fifteen matching grid points",
   "cannot be. If every point agrees the two solves saw the same matrix.","",
   "| lambda | eBIC shipped | eBIC regenerated | delta |","|---|---|---|---|"]
lam_grid = np.linspace(0.30, 1.00, 15)
n_same = 0
for k in range(min(len(bs), len(br))):
    dl = abs(bs[k]-br[k]); n_same += dl < 1e-6
    L.append(f"| {lam_grid[k]:.3f} | {bs[k]:.4f} | {br[k]:.4f} | {dl:.3e} |")
L += ["", f"Grid points agreeing to 1e-6: **{n_same} / {min(len(bs), len(br))}**", ""]

same_out = (abs(ls-lr) < 1e-9) and (es == er) and (abs(bs.min()-br.min()) < 1e-3)
diff_in  = float(np.median(d)) > 1e-6

if diff_in and same_out and n_same == min(len(bs), len(br)):
    L += ["**Red flag.** The inputs differ across 18% of pairs, yet every single",
          "eBIC grid point agrees to 1e-6. A penalised likelihood cannot be that",
          "insensitive to its own covariance argument. The likely explanation is",
          "that `solve-problem` is not actually consuming the matrix it is handed",
          "-- e.g. it recomputes or re-reads a covariance internally, making",
          "`--i-covariance-matrix` partly decorative. That would be a plugin bug",
          "with consequences well beyond this chapter, and it must be run down",
          "before any number here is published.", "",
          "Next probe: solve on a deliberately corrupted covariance (add 0.5 to a",
          "block of off-diagonals). If the eBIC path does not move, the input is",
          "being ignored and the bug is confirmed.", ""]
elif diff_in and same_out:
    L += ["**The inputs differ and the selected model is unchanged, but the eBIC",
          "curve does move** at some grid points. That is the benign reading: the",
          "two matrices differ mostly in weak off-diagonal entries that the L1",
          "penalty zeroes out anyway, so the selection is driven by the entries",
          "they agree on. Safe to regenerate; worth one sentence in the chapter.",
          ""]
elif diff_in and not same_out:
    L += ["**The outputs DIFFER.** Job 23's claim that the published numbers",
          "survive regeneration was wrong -- most likely it read the wrong",
          "artifact. Regenerating therefore CHANGES the published result, and the",
          "choice between the shipped and regenerated matrix is an editorial one.",
          ""]
else:
    L += ["**The inputs are effectively identical**, so job 24's widespread-",
          "difference finding was the wrong one. Re-check that comparison.", ""]

L += ["## Verdict","",
      f"**{'Regeneration preserves the published result' if same_out else 'Regeneration CHANGES the published result — decide'}**"]
open(report,"w").write("\n".join(L)+"\n"); print("\n".join(L))
PY
