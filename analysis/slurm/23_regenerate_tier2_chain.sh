#!/bin/bash
#SBATCH --job-name=q2-regen-t2
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/regen_t2_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/regen_t2_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Does regenerating tier 2 from raw counts change the published numbers?
#
# The shipped atacama-top-300-correlation.qza cannot be reproduced by
# transform-features -> calculate-covariance at the documented parameters: it
# differs from that chain's output by 1.147, and by the same margin before and
# after the orientation fix, so the discrepancy predates it.
#
# That leaves a choice with scientific consequences:
#   (a) keep the shipped matrix canonical -- Gate C1 keeps passing, but the
#       chapter's chain demonstrably does not regenerate it;
#   (b) regenerate everything from raw counts -- reproducible end to end, but the
#       correlation changes, so the lambda path may move.
#
# Neither should be chosen blind. This job runs (b) and reports whether the
# headline result survives: lambda = 0.8, 216 edges, min eBIC 16130.0988.
#
# It writes NOTHING into data/, publish/ or docs/. Purely informational.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/TIER2_REGENERATION_IMPACT.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-regen/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-regen/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

# Try the plausible transform settings. The shipped matrix's parameters are not
# recorded anywhere, so rather than guess once, sweep the small space and report
# which (if any) reproduces it -- that answers the provenance question too.
# pseudo_count is registered as Int, so the reachable space is small: clr at a
# couple of integer pseudo-counts, plus mclr (which ignores it entirely).
# add_metadata is ruled out -- it would append covariate columns and give a
# 305x305 matrix, and the shipped one is 300x300.
for TRANSFORM in clr mclr; do
  for PC in 1 2; do
    [[ "$TRANSFORM" == "mclr" && "$PC" != "1" ]] && continue   # mclr ignores pseudo-count
    tag="${TRANSFORM}-pc${PC}"
    echo "=== transform=$TRANSFORM pseudo-count=$PC ==="
    qiime gglasso transform-features \
      --i-table "$ROOT/data/atacama-top-300-table.qza" \
      --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
      --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
      --p-transformation "$TRANSFORM" --p-pseudo-count "$PC" --p-no-keep-original-id \
      --o-transformed-table "clr-$tag.qza" >/dev/null
    qiime gglasso calculate-covariance \
      --i-table "clr-$tag.qza" --p-method scaled \
      --o-covariance-matrix "corr-$tag.qza" >/dev/null
    echo "  built corr-$tag.qza"
  done
done

echo "=== which candidate matches the shipped matrix? ==="
BEST=$(python - "$ROOT/data/atacama-top-300-correlation.qza" corr-*.qza <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd
def load(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        n=[x for x in z.namelist() if x.endswith((".csv",".tsv",".txt"))
           and "/data/" in x and "MANIFEST" not in x][0]
        z.extract(n,t); p=os.path.join(t,n)
        return pd.read_csv(p, sep="\t" if p.endswith((".tsv",".txt")) else ",", index_col=0)
ship = load(sys.argv[1]); best=(None, float("inf"))
for q in sys.argv[2:]:
    b = load(q)
    common=[i for i in ship.index if i in b.index]
    if not common: continue
    d=float(np.nanmax(np.abs(ship.loc[common,common].to_numpy(float)
                             - b.loc[common,common].to_numpy(float))))
    print(f"  {os.path.basename(q):28} max|diff| = {d:.4e}", file=sys.stderr)
    if d < best[1]: best=(q, d)
print(best[0])
PY
)
echo "  closest: $BEST"

echo "=== lambda path on the REGENERATED correlation (clr, pc=1) ==="
qiime gglasso solve-problem \
  --i-covariance-matrix corr-clr-pc1.qza \
  --p-n-samples 54 --p-no-latent --p-path-scale linear \
  --p-lambda1-min 0.30 --p-lambda1-max 1.00 --p-n-lambda1 15 \
  --p-gamma 0.3 --o-solution path-regen.qza --verbose 2>&1 | tail -2

python - path-regen.qza "$ROOT/data/atacama-top-300-correlation.qza" "$BEST" "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, zarr

path_qza, ship_q, best_q, report = sys.argv[1:5]
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(path_qza) as z:
    i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
    r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
    ms=r["modelselect_stats"]
    lam=np.asarray(ms["LAMBDA"]).ravel()
    ebic=np.asarray(ms["BIC"]["0.3"]).ravel()
    best_lam=float(np.asarray(ms["BEST/lambda1"]).item())
    P=np.asarray(r["solution/precision_"])
edges=int((np.abs(np.triu(P,1))>1e-8).sum())
j=int(np.argmin(ebic))

REF_LAM, REF_EDGES, REF_EBIC = 0.8, 216, 16130.0988
same_lam   = abs(best_lam-REF_LAM) < 1e-9
same_edges = edges == REF_EDGES
survives   = same_lam and same_edges

L=["# Tier-2 regeneration — impact on the published numbers","",
   "Regenerating the correlation matrix from the raw counts with the DOCUMENTED",
   "commands, then re-running the lambda path on it.","",
   "## Does the headline result survive?","",
   "| quantity | published | regenerated | |",
   "|---|---|---|---|",
   f"| selected lambda | {REF_LAM} | {best_lam:.6g} | {'same' if same_lam else '**CHANGED**'} |",
   f"| edges | {REF_EDGES} | {edges} | {'same' if same_edges else '**CHANGED**'} |",
   f"| min eBIC | {REF_EBIC} | {ebic[j]:.4f} | {'~same' if abs(ebic[j]-REF_EBIC)/REF_EBIC < 1e-3 else '**CHANGED**'} |",
   ""]

if survives:
    L += ["**The published result survives regeneration.** Selected lambda and edge",
          "count are unchanged, so switching to a fully reproducible chain costs",
          "nothing scientifically. Regenerate, republish the artifacts, and the",
          "chapter's commands then genuinely reproduce its table.",""]
else:
    L += ["**The published numbers CHANGE under regeneration.** This is now an",
          "editorial decision, not a technical one: either keep the shipped matrix",
          "as canonical and state plainly that the documented chain does not",
          "reproduce it, or adopt the regenerated numbers and update the chapter",
          "and the manuscript. Do not switch silently.",""]

L += [f"Closest transform setting to the shipped matrix: `{os.path.basename(best_q) if best_q else 'none'}`.",
      "None of clr at pseudo-count 1 or 2, nor mclr, reproduces it exactly unless the",
      "difference above is ~0, so the shipped matrix's provenance remains unknown.",
      "", "## Verdict", "",
      f"**{'SAFE TO REGENERATE' if survives else 'REGENERATION CHANGES PUBLISHED NUMBERS — decide'}**"]
open(report,"w").write("\n".join(L)+"\n"); print("\n".join(L))
PY
echo "report -> $REPORT"
