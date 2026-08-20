#!/bin/bash
#SBATCH --job-name=q2-orient-fix
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/orient_fix_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/orient_fix_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Does the orientation fix preserve the science?
#
# transform_features no longer transposes back to (p, N), and calculate_covariance
# now takes rowvar=False. Two changes that must cancel exactly: if they do, the
# chain
#     raw counts -> transform-features -> calculate-covariance
# reproduces the SHIPPED atacama-top-300-correlation.qza, which is the input Gate
# C1 already validated. That is a much stronger check than "it runs".
#
# Also asserts the artifact is now correctly oriented and usable by q2-classo --
# the thing that was impossible before.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$ROOT/reports/ORIENTATION_FIX_VERIFICATION.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-orientfix/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-orientfix/${SLURM_JOB_ID:-manual}}"
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

echo "=== [1/4] transform-features on the raw tier-2 table ==="
qiime gglasso transform-features \
  --i-table "$ROOT/data/atacama-top-300-table.qza" \
  --i-taxonomy "$ROOT/data/atacama-taxonomy-silva138.qza" \
  --m-sample-metadata-file "$ROOT/data/sample-metadata.tsv" \
  --p-transformation clr --p-pseudo-count 1 \
  --p-no-keep-original-id \
  --o-transformed-table clr-new.qza

echo "=== [2/4] calculate-covariance ==="
qiime gglasso calculate-covariance \
  --i-table clr-new.qza --p-method scaled \
  --o-covariance-matrix corr-new.qza

echo "=== [3/4] compare against the shipped correlation matrix ==="
python - "$ROOT/data/atacama-top-300-correlation.qza" corr-new.qza clr-new.qza \
         "$ROOT/data/atacama-classo-outcomes-mean-imputed.tsv" "$REPORT" <<'PY'
import sys, zipfile, tempfile, os, glob
import numpy as np, pandas as pd, biom

ship, new, clr, outcomes, report = sys.argv[1:6]
L = ["# Orientation fix — verification", ""]

def load_pairwise(qza):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(qza) as z:
        names = [n for n in z.namelist() if n.endswith((".csv", ".tsv", ".txt"))
                 and "/data/" in n and "MANIFEST" not in n]
        if not names:
            raise SystemExit(f"{qza}: no data table found; members = {z.namelist()[:8]}")
        z.extract(names[0], t)
        p = os.path.join(t, names[0])
        sep = "\t" if p.endswith((".tsv", ".txt")) else ","
        return pd.read_csv(p, sep=sep, index_col=0)

a, b = load_pairwise(ship), load_pairwise(new)
L += [f"- shipped correlation : {a.shape}", f"- recomputed          : {b.shape}", ""]

same_shape = a.shape == b.shape
if same_shape:
    common = [i for i in a.index if i in b.index]
    if not common:
        L += [f"- shipped labels look like `{list(a.index[:2])}`",
              f"- recomputed labels look like `{list(b.index[:2])}`", "",
              "**No shared labels — the two runs used different "
              "`--p-keep-original-id` settings, so this comparison is "
              "meaningless. Not a statement about the fix.**", ""]
        open(report, "w").write("\n".join(L) + "\n"); print("\n".join(L))
        raise SystemExit(2)
    A = a.loc[common, common].to_numpy(float)
    B = b.loc[common, common].to_numpy(float)
    maxabs = float(np.nanmax(np.abs(A - B)))
    L += [f"- matched labels      : {len(common)}/{len(a.index)}",
          f"- max |difference|    : {maxabs:.3e}", ""]
    ok = maxabs < 1e-6
    L += [("**The fixed chain reproduces the shipped correlation matrix.** The two "
           "changes cancel exactly, so every downstream result -- including Gate "
           "C1 -- is unaffected." if ok else
           "**MISMATCH.** The fix changes the covariance, so the science is NOT "
           "preserved. Do not merge."), ""]
else:
    ok = False
    L += ["**Shape mismatch — the fix did not produce a p x p matrix.**", ""]

# --- orientation + classo usability -----------------------------------------
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(clr) as z:
    i = [n for n in z.namelist() if n.endswith("feature-table.biom")][0]
    z.extract(i, t); tab = biom.load_table(os.path.join(t, i))
md = pd.read_csv(outcomes, sep="\t", index_col=0)
md = md[~md.index.astype(str).str.startswith("#")]
overlap = set(tab.ids("sample")) & set(md.index.astype(str))

L += ["## Artifact orientation", "",
      f"- biom samples      : {len(tab.ids('sample'))}  e.g. `{tab.ids('sample')[0]}`",
      f"- biom observations : {len(tab.ids('observation'))}  e.g. `{tab.ids('observation')[0]}`",
      f"- overlap with the outcomes file : **{len(overlap)}**", ""]
usable = len(overlap) > 0
L += [("**q2-classo can now use this artifact.** Before the fix the same join "
       "matched zero samples and c-lasso died with `IndexError: index 0 is out "
       "of bounds for axis 0 with size 0`." if usable else
       "**Still unusable by q2-classo — the join matches nothing.**"), ""]

L += ["## Verdict", "", f"**{'PASS' if (ok and usable) else 'FAIL'}**"]
open(report, "w").write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(0 if (ok and usable) else 1)
PY

echo "=== [4/4] q2-classo end-to-end on the newly oriented artifact ==="
qiime classo regress \
  --i-features clr-new.qza \
  --m-y-file "$ROOT/data/atacama-classo-outcomes-mean-imputed.tsv" \
  --m-y-column ph \
  --p-concomitant --p-path --p-path-nlam-log 20 \
  --p-cv --p-cv-subsets 5 --p-cv-nlam 20 \
  --p-no-stabsel --p-no-lamfixed \
  --o-result regress-check.qza \
  && echo "q2-classo regress on gglasso output: OK (was impossible before)"
