#!/bin/bash
#SBATCH --job-name=q2-t2-ids
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/t2ids_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/t2ids_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Rebuild the tier-2 chain with --p-keep-original-id (DECISIONS_NEEDED item 7,
# option (a)). Features then carry their real 32-character IDs instead of
# positional ASV-k labels, so the taxonomy join works directly and the labels
# stop depending on which plugin version built the artifact.
#
# Two things this must prove, not assume:
#
#   1. Gate C1 still passes. The graphical-lasso objective is invariant under
#      simultaneous row/column permutation, so changing the feature ORDER must
#      leave lambda = 0.8, 216 edges and eBIC 16130.0994 untouched. If it does
#      not, something other than the ordering changed and the run is void.
#
#   2. The taxonomy join actually works. That is the entire point of option (a),
#      and it is exactly the kind of claim that reports success while doing
#      nothing -- a failed join returns all-NaN rather than raising. The job
#      asserts every feature resolves to a real taxon.
#
# Writes to results/tier2-regen/ ONLY. Nothing under publish/ or data/ is
# touched; promotion is a separate, deliberate step.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
OUT="$ROOT/results/tier2-regen"
REPORT="$ROOT/reports/TIER2_KEEP_IDS.md"

GAMMA=0.3; LMIN=0.30; LMAX=1.00; NLAM=15; LSEL=0.80; NSAMP=54
MUS=(15.0 10.0 7.5)
SKIPPED_PCA=()   # mu tags whose solution had rank 0, so pca legitimately did not run

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-t2ids/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-t2ids/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$ROOT/reports"
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

echo "############ [1/6] transform-features, keeping real feature IDs ############"
qiime gglasso transform-features \
  --i-table "$ROOT/data/atacama-top-300-table.qza" \
  --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
  --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
  --p-transformation clr --p-pseudo-count 1 --p-keep-original-id \
  --o-transformed-table clr.qza
echo "  -> clr.qza"

echo "############ [2/6] calculate-covariance ############"
qiime gglasso calculate-covariance \
  --i-table clr.qza --p-method scaled --o-covariance-matrix correlation.qza
echo "  -> correlation.qza"

echo "############ [3/6] lambda path (Gate C1 re-verification) ############"
qiime gglasso solve-problem \
  --i-covariance-matrix correlation.qza \
  --p-n-samples "$NSAMP" --p-no-latent --p-path-scale linear \
  --p-lambda1-min "$LMIN" --p-lambda1-max "$LMAX" --p-n-lambda1 "$NLAM" \
  --p-gamma "$GAMMA" --o-solution sgl-path.qza
echo "  -> sgl-path.qza"

echo "############ [4/6] single fit at lambda = $LSEL ############"
qiime gglasso solve-problem \
  --i-covariance-matrix correlation.qza \
  --p-n-samples "$NSAMP" --p-no-latent \
  --p-lambda1-min "$LSEL" --p-lambda1-max "$LSEL" --p-n-lambda1 1 \
  --p-gamma "$GAMMA" --o-solution sgl-lambda08.qza
echo "  -> sgl-lambda08.qza"

echo "############ [5/6] SLR fits + pca ############"
for MU in "${MUS[@]}"; do
  tag="${MU/./p}"
  qiime gglasso solve-problem \
    --i-covariance-matrix correlation.qza \
    --p-n-samples "$NSAMP" --p-latent \
    --p-lambda1-min "$LSEL" --p-lambda1-max "$LSEL" --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min "$MU" --p-mu1-max "$MU" --p-n-mu1 1 \
    --o-solution "slr-mu${tag}.qza"
  echo "  -> slr-mu${tag}.qza"

  # pca defaults to 3 components, but the achieved rank falls as mu rises --
  # mu=15 gives rank 2 here -- and the rank guard correctly refuses rather than
  # silently truncating. Derive the component count from the solution instead of
  # hardcoding it, or the headline fit is the one that fails.
  RANK=$(python - "slr-mu${tag}.qza" <<'PYRANK'
import os
import sys
import tempfile
import zipfile

import numpy as np
import zarr

with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(sys.argv[1]) as z:
    inner = [n for n in z.namelist() if n.endswith("problem.zip")][0]
    z.extract(inner, t)
    r = zarr.open(zarr.ZipStore(os.path.join(t, inner), mode="r"))
    print(int(np.linalg.matrix_rank(np.asarray(r["solution/lowrank_"]), tol=1e-8)))
PYRANK
)
  echo "     achieved rank: $RANK"
  if [[ "$RANK" -lt 1 ]]; then
    echo "     rank 0 — no low-rank component, pca not applicable at this mu"
    SKIPPED_PCA+=("$tag")
    continue
  fi
  NCOMP=$(( RANK < 3 ? RANK : 3 ))
  # pca needs the latent solution AND sample metadata; it fails at
  # _pca/_visualizer.py:352 without the latter even though it is optional.
  qiime gglasso pca \
    --i-solution "slr-mu${tag}.qza" --i-table clr.qza \
    --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
    --p-n-components "$NCOMP" \
    --o-visualization "slr-mu${tag}-pca.qzv"
  echo "     pca OK with $NCOMP components"
done

echo "############ verify every expected pca landed ############"
# The first run of this job printed "READY TO PROMOTE" while silently missing
# slr-mu15p0-pca.qzv -- the headline fit. Assert rather than warn. A mu whose
# solution has rank 0 has no low-rank part to decompose and is legitimately
# skipped above, so exempt exactly those and no others.
for MU in "${MUS[@]}"; do
  tag="${MU/./p}"
  skipped=0
  for s in "${SKIPPED_PCA[@]:-}"; do [[ "$s" == "$tag" ]] && skipped=1; done
  if [[ "$skipped" == "1" ]]; then
    echo "  skipped (rank 0): slr-mu${tag}-pca.qzv"; continue
  fi
  if [[ ! -f "slr-mu${tag}-pca.qzv" ]]; then
    echo "MISSING: slr-mu${tag}-pca.qzv (rank was >= 1, so pca should have run)"
    exit 1
  fi
  echo "  ok: slr-mu${tag}-pca.qzv"
done

echo "############ [6/6] summarize ############"
qiime gglasso summarize --i-solution sgl-path.qza \
  --o-visualization sgl-linear-path.qzv
echo "  -> sgl-linear-path.qzv"

cp -f clr.qza correlation.qza sgl-path.qza sgl-lambda08.qza \
      sgl-linear-path.qzv slr-mu*.qza "$OUT"/ 2>/dev/null || true
cp -f slr-mu*-pca.qzv "$OUT"/ 2>/dev/null || true

python - "$ROOT" "$OUT" sgl-path.qza sgl-lambda08.qza clr.qza "$REPORT" <<'PY'
import os
import sys
import tempfile
import zipfile

import numpy as np
import zarr

root, outdir, path_q, fit_q, clr_q, report = sys.argv[1:7]

REF_LAM, REF_EDGES, REF_EBIC = 0.80, 216, 16130.099453949091


def open_problem(q):
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(q) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")][0]
        z.extract(inner, tmp)
    return zarr.open(zarr.ZipStore(os.path.join(tmp, inner), mode="r"))


r = open_problem(path_q)
ms = r["modelselect_stats"]
ebic = np.asarray(ms["BIC"]["0.3"]).ravel()
best_lam = float(np.asarray(ms["BEST/lambda1"]).item())
P = np.asarray(r["solution/precision_"])
edges = int((np.abs(np.triu(P, 1)) > 1e-8).sum())

same_lam = abs(best_lam - REF_LAM) < 1e-9
same_edges = edges == REF_EDGES
# 1e-3 relative would accept a drift of +/- 16 in eBIC, which is larger than the
# gap between adjacent points on the lambda grid -- it would pass a genuinely
# different model. Relabelling is permutation-only, so the observed residual is
# ~7e-4 absolute (~4e-8 relative); 1e-6 relative still leaves three orders of
# margin over that while catching anything real.
same_ebic = abs(ebic.min() - REF_EBIC) / REF_EBIC < 1e-6
gate = same_lam and same_edges and same_ebic

L = ["# Tier 2 regenerated with real feature IDs", "",
     "Rebuilt with `--p-keep-original-id` so features carry their 32-character",
     "IDs instead of positional `ASV-k` labels (DECISIONS_NEEDED item 7, option",
     "(a)).", "",
     "## Gate C1 — do the published numbers survive the relabelling?", "",
     "They must: the graphical-lasso objective is invariant under simultaneous",
     "row/column permutation, and changing `keep_original_id` changes only the",
     "labels and the tie order.", "",
     "| quantity | reference | regenerated | |", "|---|---|---|---|",
     f"| selected lambda | {REF_LAM} | {best_lam:.6g} | "
     f"{'same' if same_lam else '**CHANGED**'} |",
     f"| edges | {REF_EDGES} | {edges} | {'same' if same_edges else '**CHANGED**'} |",
     f"| min eBIC (gamma=0.3) | {REF_EBIC:.4f} | {ebic.min():.4f} | "
     f"{'~same' if same_ebic else '**CHANGED**'} |", ""]

# ---- the point of the exercise: does the taxonomy join work? ----------------
try:
    import biom
    import pandas as pd
    import qiime2

    clr = qiime2.Artifact.load(os.path.join(outdir, "clr.qza")).view(biom.Table)
    feats = list(clr.to_dataframe().index)
    tax = (qiime2.Artifact.load(os.path.join(root, "data",
                                             "atacama-taxonomy-silva138.qza"))
           .view(pd.DataFrame))

    hexish = sum(1 for f in feats if len(f) == 32 and all(c in "0123456789abcdef" for c in f))
    joined = tax.reindex(feats)["Taxon"]
    unresolved = int(joined.isna().sum())

    L += ["## Does the taxonomy join work?", "",
          "A failed join in QIIME 2 returns all-`NaN` rather than raising, so this",
          "is asserted rather than assumed.", "",
          f"- features in the transformed table : {len(feats)}",
          f"- that look like real feature IDs (32 hex chars) : **{hexish}**",
          f"- unresolved against the taxonomy : **{unresolved}**", ""]
    if unresolved == 0 and hexish == len(feats):
        L += ["**The join resolves every feature.** This is what option (a) was",
              "for: `tax.loc[feature, \"Taxon\"]` now works directly, with no",
              "mapping step and no reliance on abundance rank.", "",
              "Example rows:", "", "| feature | taxon |", "|---|---|"]
        for f in feats[:4]:
            L.append(f"| `{f[:12]}…` | {str(joined[f])[:78]} |")
        L.append("")
    else:
        L += ["**The join does NOT fully resolve.** Investigate before promoting.",
              ""]
    join_ok = unresolved == 0 and hexish == len(feats)
except Exception as exc:  # noqa: BLE001
    L += [f"(taxonomy-join check failed to run: {exc})", ""]
    join_ok = False

L += ["## Artifacts written", "",
      "To `results/tier2-regen/` only — nothing under `publish/` or `data/` was",
      "touched. Promotion is a separate step.", ""]
for f in sorted(os.listdir(outdir)):
    L.append(f"- `{f}` ({os.path.getsize(os.path.join(outdir, f)) // 1024} KB)")

L += ["", "## Verdict", "",
      f"**{'READY TO PROMOTE — Gate C1 holds and the taxonomy join resolves' if (gate and join_ok) else 'DO NOT PROMOTE — see above'}**"]

open(report, "w").write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(0 if (gate and join_ok) else 1)
PY
echo "report -> $REPORT"
