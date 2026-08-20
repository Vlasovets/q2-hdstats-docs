#!/bin/bash
#SBATCH --job-name=q2-corrdiff
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/corrdiff_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/corrdiff_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# WHERE does the shipped correlation matrix differ from the regenerated one?
#
# max|diff| = 1.147 looked alarming until the lambda path came back byte-identical
# on both. Those two facts are only consistent if the difference is concentrated
# in a few cells that the graphical lasso does not care about.
#
# The shipped matrix has 2 off-diagonal entries at exactly r = 1.0 and 130 at
# r >= 0.99 -- feature pairs whose CLR values are proportional, i.e. linearly
# dependent. Correlation there is numerically unstable (near-zero variances in the
# denominator), so a tiny difference upstream can swing those cells wildly while
# leaving everything else alone.
#
# This job tests that explanation: it reports the DISTRIBUTION of differences, not
# just the maximum, and checks whether the large ones sit on the degenerate pairs.
# max is a bad summary and nearly led to the wrong conclusion twice.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/CORRELATION_DIFF_STRUCTURE.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-corrdiff/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-corrdiff/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
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

qiime gglasso transform-features \
  --i-table "$ROOT/data/atacama-top-300-table.qza" \
  --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
  --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
  --p-transformation clr --p-pseudo-count 1 --p-no-keep-original-id \
  --o-transformed-table clr.qza >/dev/null
qiime gglasso calculate-covariance \
  --i-table clr.qza --p-method scaled --o-covariance-matrix corr.qza >/dev/null
echo "regenerated corr.qza"

python - "$ROOT/data/atacama-top-300-correlation.qza" corr.qza "$REPORT" <<'PY'
import sys, zipfile, tempfile, os
import numpy as np, pandas as pd

def load(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        n=[x for x in z.namelist() if x.endswith((".csv",".tsv",".txt"))
           and "/data/" in x and "MANIFEST" not in x][0]
        z.extract(n,t); p=os.path.join(t,n)
        return pd.read_csv(p, sep="\t" if p.endswith((".tsv",".txt")) else ",", index_col=0)

ship, new, report = load(sys.argv[1]), load(sys.argv[2]), sys.argv[3]
common=[i for i in ship.index if i in new.index]
A=ship.loc[common,common].to_numpy(float); B=new.loc[common,common].to_numpy(float)
D=np.abs(A-B)
iu=np.triu_indices(len(common),1)
d=D[iu]; a=A[iu]

n_tot=d.size
tiers=[(1e-9,"< 1e-9 (identical)"),(1e-6,"< 1e-6"),(1e-3,"< 1e-3"),
       (1e-1,"< 0.1"),(np.inf,">= 0.1")]
lines=[]; prev=0.0
for thr,label in tiers:
    m=(d>=prev)&(d<thr) if thr!=np.inf else (d>=prev)
    lines.append(f"| {label} | {int(m.sum()):,} | {100*m.mean():.4f}% |")
    prev=thr

big = d >= 0.1
degenerate = a >= 0.99          # shipped value already at/near perfect correlation
overlap = int((big & degenerate).sum())

L=["# Where the shipped and regenerated correlation matrices differ","",
   f"Comparing the {len(common)} x {len(common)} matrices over their "
   f"{n_tot:,} off-diagonal pairs.","",
   "## Distribution of |difference|","",
   "| bucket | pairs | share |","|---|---|---|", *lines, "",
   f"- median |difference| : **{np.median(d):.3e}**",
   f"- mean   |difference| : {d.mean():.3e}",
   f"- max    |difference| : {d.max():.3e}", "",
   "## Are the large differences on the degenerate pairs?","",
   f"- pairs with |difference| >= 0.1 : **{int(big.sum())}**",
   f"- of those, pairs where the shipped value is >= 0.99 : **{overlap}** "
   f"({100*overlap/max(int(big.sum()),1):.0f}%)", ""]

concentrated = big.sum() < 0.001 * n_tot
explained = overlap >= 0.5 * max(int(big.sum()), 1)

if concentrated and explained:
    L += ["**Explanation confirmed.** The differences are confined to a vanishing "
          "fraction of pairs, and those pairs are the ones the shipped matrix "
          "already reports as near-perfectly correlated. Those are linearly "
          "dependent features -- CLR values proportional across samples -- where "
          "the correlation denominator involves near-zero variances and the value "
          "is numerically unstable. Everything that carries signal agrees.", "",
          "This is why the lambda path is unchanged: lambda = 0.8, 216 edges, "
          "eBIC 16130.0988 on either matrix.", ""]
elif concentrated:
    L += ["**Differences are concentrated but NOT on the degenerate pairs.** "
          "Worth understanding before republishing -- the numerical-instability "
          "explanation does not hold.", ""]
else:
    L += ["**Differences are widespread**, so the two matrices are genuinely "
          "different objects. Do not republish without deciding which is "
          "canonical.", ""]

L += ["## Note for the tier-2 chapter", "",
      "Independently of the regeneration question, this table contains "
      "near-duplicate features: 2 pairs at exactly r = 1.0 and 130 at r >= 0.99. "
      "Perfectly correlated CLR features are proportional in the raw counts, "
      "which for rare ASVs usually means an identical presence/absence pattern. "
      "They make the empirical correlation matrix singular, and any edge the "
      "network draws between such a pair is an artefact of that degeneracy rather "
      "than a finding. Worth stating in the chapter.", "",
      "## Verdict", "",
      f"**{'SAFE TO REGENERATE — differences are numerically meaningless' if (concentrated and explained) else 'INVESTIGATE BEFORE REPUBLISHING'}**"]
open(report,"w").write("\n".join(L)+"\n"); print("\n".join(L))
PY
