#!/bin/bash
#SBATCH --job-name=q2-lambda-path
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/lambda_path_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/lambda_path_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Stage 1 of the tier-2 recompute -- and GATE C1.
#
# The published tier-2 numbers (lambda=0.8, min eBIC 16130.10, 216 edges) were
# NOT produced by `qiime gglasso solve-problem`. They came from an offline script
# calling raw GGLasso with do_scaling=False and gglasso.helper.model_selection.ebic.
# The chapter nevertheless presents a solve-problem invocation as their source.
#
# This stage runs the linear lambda path THROUGH THE CLI and compares. Nothing
# downstream should be written into the tier-2 chapters until it passes, because
# otherwise the docs keep claiming a command produced numbers it did not.
#
# Exit codes:
#   0  gate passed  -- the CLI reproduces the published selection
#   4  gate FAILED  -- CLI and offline script disagree; see the printed diff
#   other: infrastructure failure

set -euo pipefail

REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda

IN_COV="$ROOT/data/atacama-top-300-correlation.qza"
OUT_DIR="$ROOT/results/gglasso"
TAB_DIR="$ROOT/results/tables"

# ---- scratch (node-local disk; $LOCAL_SCRATCH is not set on this cluster) -----
for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break
  fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; export TMP="$TMPDIR"; export TEMP="$TMPDIR"
mkdir -p "$TMPDIR" "$OUT_DIR" "$TAB_DIR" "$ROOT/slurm/logs"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# ---- environment -------------------------------------------------------------
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


[[ -r "$IN_COV" ]] || { echo "missing input: $IN_COV" >&2; exit 2; }

# ---- read the canonical parameters -------------------------------------------
read -r GAMMA LMIN LMAX NLAM SEL EXP_EBIC EXP_EDGES TOL <<<"$(
python - "$ROOT/config/params.yaml" <<'PY'
import sys, re
# Deliberately not importing yaml: it is not guaranteed in the QIIME 2 env and
# this file is a flat, hand-maintained mapping.
txt = open(sys.argv[1]).read()
def g(key):
    m = re.search(rf"^\s*{key}:\s*([0-9.eE+-]+)\s*(?:#.*)?$", txt, re.M)
    if not m:
        sys.exit(f"could not read '{key}' from params.yaml")
    return m.group(1)
print(g("gamma"), g("lambda1_min"), g("lambda1_max"), g("n_lambda1"),
      g("selected_lambda"), g("expect_min_ebic"), g("expect_n_edges"),
      g("ebic_rel_tol"))
PY
)"

echo "gamma=$GAMMA  lambda in [$LMIN, $LMAX] x $NLAM (linear)"
echo "gate: expect lambda=$SEL, eBIC=$EXP_EBIC, edges=$EXP_EDGES (rel tol $TOL)"

PATH_QZA="$OUT_DIR/atacama-top-300-sgl-linear-path.qza"

# ---- run the CLI-native path -------------------------------------------------
if [[ -f "$PATH_QZA" ]]; then
  echo "[skip] $PATH_QZA already exists"
else
  qiime gglasso solve-problem \
    --i-covariance-matrix "$IN_COV" \
    --p-n-samples 54 \
    --p-no-latent \
    --p-path-scale linear \
    --p-lambda1-min "$LMIN" \
    --p-lambda1-max "$LMAX" \
    --p-n-lambda1 "$NLAM" \
    --p-gamma "$GAMMA" \
    --o-solution "$PATH_QZA" \
    --verbose
fi

qiime gglasso summarize \
  --i-solution "$PATH_QZA" \
  --o-visualization "$OUT_DIR/atacama-top-300-sgl-linear-path.qzv"

# ---- extract the path table and evaluate the gate ----------------------------
python - "$PATH_QZA" "$TAB_DIR/lambda-path.tsv" "$SEL" "$EXP_EBIC" "$EXP_EDGES" "$TOL" "$GAMMA" <<'PY'
import sys, zipfile, tempfile, os, csv
import numpy as np
import zarr

qza, out_tsv, sel, exp_ebic, exp_edges, tol, gamma = sys.argv[1:8]
sel, exp_ebic, exp_edges, tol, gamma = (
    float(sel), float(exp_ebic), int(exp_edges), float(tol), float(gamma)
)

# A .qza is a zip; the GGLassoProblem payload is a single problem.zip (zarr).
with tempfile.TemporaryDirectory() as tmp:
    with zipfile.ZipFile(qza) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")]
        if not inner:
            sys.exit("no problem.zip inside the artifact")
        z.extract(inner[0], tmp)
        root = zarr.open(zarr.ZipStore(os.path.join(tmp, inner[0]), mode="r"))

    if "modelselect_stats" not in root:
        sys.exit("solution carries no modelselect_stats -- the grid collapsed to "
                 "a single fit; check n_lambda1")

    ms = root["modelselect_stats"]
    lam = np.asarray(ms["LAMBDA"]).ravel()
    sp = np.asarray(ms["SP"]).ravel()          # sparsity per grid point
    best_lambda = float(np.asarray(ms["BEST/lambda1"]).item())

    # `modelselect_stats/BIC` is NOT an array of eBIC values -- it is the set of
    # gamma keys, serialised as strings ('0.1', '0.3', '0.5', '0.7'), with the
    # per-gamma curves in a subgroup. Reading it as an array yields numpy.str_
    # and any float formatting then raises
    #   ValueError: Unknown format code 'f' for object of type 'numpy.str_'
    # Select the curve for the gamma we actually solved at.
    ebic = None
    if "BIC" in ms:
        gkey = f"{gamma:g}"
        try:
            ebic = np.asarray(ms["BIC"][gkey]).ravel()
        except (KeyError, IndexError, TypeError):
            available = list(ms["BIC"].keys()) if hasattr(ms["BIC"], "keys") else []
            sys.exit(
                f"no eBIC curve for gamma={gkey}; available: {available}"
            )

    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["lambda1", "sparsity"] + (["ebic"] if ebic is not None else []))
        for i, l in enumerate(lam):
            row = [f"{l:.6g}", f"{sp[i]:.6g}"]
            if ebic is not None:
                row.append(f"{ebic.ravel()[i]:.6f}")
            w.writerow(row)
    print(f"wrote {out_tsv} ({len(lam)} grid points)")

    P = np.asarray(root["solution/precision_"])
    n_edges = int((np.abs(np.triu(P, k=1)) > 1e-8).sum())

print(f"\nCLI selected lambda = {best_lambda:.6g}   (expected {sel})")
print(f"CLI edge count      = {n_edges}            (expected {exp_edges})")
if ebic is not None:
    print(f"CLI min eBIC        = {float(ebic.min()):.4f}    (expected {exp_ebic})")

fail = []
if abs(best_lambda - sel) > 1e-9:
    fail.append(f"lambda {best_lambda} != {sel}")
if n_edges != exp_edges:
    fail.append(f"edges {n_edges} != {exp_edges}")
if ebic is not None and abs(float(ebic.min()) - exp_ebic) > tol * abs(exp_ebic):
    fail.append(f"eBIC {float(ebic.min()):.4f} != {exp_ebic}")

if fail:
    print("\nGATE C1: FAIL -- " + "; ".join(fail))
    print("""
Do NOT paper over this. The offline script used do_scaling=False and GGLasso's
own ebic(); `solve-problem` does its own model selection and may scale the input
covariance. Two acceptable resolutions:
  (a) expose/align the scaling and eBIC convention in q2-gglasso, then re-run;
  (b) keep linear_lambda_model_selection.py as a cited auxiliary script and stop
      presenting solve-problem as the source of these numbers in the chapter.
Pick one before writing any tier-2 chapter text.""")
    sys.exit(4)

print("\nGATE C1: PASS -- the documented command genuinely produces the "
      "documented table.")
PY

echo "stage 1 complete"
