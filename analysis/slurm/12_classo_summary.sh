#!/bin/bash
#SBATCH --job-name=q2-classo-sum
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classo_sum_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classo_sum_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Stage 7 -- aggregate the per-outcome CV fits into one generated table.
#
# Written straight into docs/_data/ so the tier-2 classo chapter renders it with
# {csv-table} instead of transcribing numbers, the same contract the lambda path
# and mu->rank map now use.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
DOCS="$REPO"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/results/tables" "$DOCS/docs/_data"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

python - "$ROOT/results/classo" "$ROOT/results/tables/classo-cv-summary.tsv" \
         "$DOCS/docs/_data/atacama-classo-cv.tsv" <<'PY'
import sys, os, glob, csv, zipfile, tempfile
import numpy as np, zarr

src, out_tsv, docs_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

def open_problem(qza):
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")]
        if not inner:
            return None
        z.extract(inner[0], tmp)
    return zarr.open(zarr.ZipStore(os.path.join(tmp, inner[0]), mode="r"))

rows = []
for qza in sorted(glob.glob(os.path.join(src, "atacama-top-300-regress-*.qza"))):
    outcome = os.path.basename(qza)[len("atacama-top-300-regress-"):-len(".qza")]
    root = open_problem(qza)
    if root is None:
        print(f"  {outcome}: no problem.zip, skipped"); continue
    try:
        refit = np.asarray(root["solution/CV/refit"]).ravel()
        selected = np.asarray(root["solution/CV/selected_param"]).ravel().astype(bool)
        xg = np.asarray(root["solution/CV/xGraph"]).ravel()
        yg = np.asarray(root["solution/CV/yGraph"]).ravel()
        labels = np.asarray(root["data/label"]).ravel()
    except KeyError as e:
        print(f"  {outcome}: missing {e}, skipped"); continue

    n_sel = int(selected.sum())
    j = int(np.argmin(yg)) if yg.size else -1
    rows.append({
        "outcome": outcome,
        "selected features": n_sel,
        "lambda (min CV)": f"{xg[j]:.5g}" if j >= 0 else "n/a",
        "CV error (min)": f"{yg[j]:.5g}" if j >= 0 else "n/a",
        "max |beta|": f"{np.abs(refit).max():.4g}",
    })
    top = sorted(zip(labels[: refit.size], refit), key=lambda t: -abs(t[1]))[:3]
    print(f"  {outcome:38} {n_sel:>4} selected   top: "
          + ", ".join(f"{l}({b:+.3g})" for l, b in top))

if not rows:
    raise SystemExit("no CV results found -- did stage 6 run?")

cols = ["outcome", "selected features", "lambda (min CV)", "CV error (min)", "max |beta|"]
for dest in (out_tsv, docs_tsv):
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["outcome"]):
            w.writerow(r)
    print(f"wrote {dest} ({len(rows)} outcomes)")
PY
