#!/bin/bash
#SBATCH --job-name=q2-mgl-verify
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/mglverify_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/mglverify_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Reproduce -- or refute -- every claim on the MGL chapter.
#
# 02_lowdim_gglasso/06_multiple_graphical_lasso.md is 434 lines and says of
# itself: "Nothing on this page has been run against QIIME 2 2026.7, and the
# gaps below were read out of the plugin source rather than reproduced." It
# then makes six specific, falsifiable claims. Source-reading is how the
# 1.147 discrepancy survived as long as it did, so this stage runs them.
#
# It is deliberately NON-FATAL on the claims themselves: several are assertions
# that something is BROKEN, so a failure to run is the expected result and must
# be recorded rather than treated as a stage error. The stage fails only if it
# cannot execute the experiment at all.
#
# Claims under test:
#   C1  the chapter's Step 1 filter uses [transect-name], but the tier-1
#       metadata shipped in q2-gglasso/data has no such column
#   C2  the tier-1 table splits 25/25 across the two transects
#   C3  build-groups prints per-instance p_k, N_k and a group count
#   C4  build-groups returns nothing for identically-labelled tables, and the
#       action fails rather than emitting an empty artifact
#   C5  the difference check compares COLUMN labels (sample IDs), not features,
#       so a metadata split always reports the instances as differing
#   C6  gap 1: TensorData does not chain to the List[Int] parameter
#       gap 3: with one 2-D covariance the SGL branch is taken regardless of
#              --p-non-conforming / --p-reg / --p-group-array

set -uo pipefail          # NOT -e: individual claims are allowed to fail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
GG="${Q2_HDSTATS_PLUGINS:-/home/itg/oleg.vlasovets/slr_example}/q2-gglasso/data"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
REPORT="$REPO/analysis/reports/mgl-verification.md"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-mglverify/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-mglverify/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$(dirname "$REPORT")"
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

echo "############ [0/6] the table and its transect split ############"
python - "$GG" <<'PY'
import sys, pandas as pd, qiime2
gg = sys.argv[1]
t = qiime2.Artifact.load(f"{gg}/atacama-counts.qza").view(pd.DataFrame)
print(f"  atacama-counts.qza: {t.shape[0]} samples x {t.shape[1]} features")
pref = pd.Series([s[:3] for s in t.index]).value_counts().to_dict()
print(f"  transect split by sample-id prefix: {pref}")
md = pd.read_csv(f"{gg}/selected-atacama-sample-metadata.tsv", sep="\t")
print(f"  shipped tier-1 metadata columns: {list(md.columns)}")
print(f"  has 'transect-name'? {'transect-name' in md.columns}")
# Emit a metadata file that DOES carry the grouping, derived from the ID prefix.
grp = pd.DataFrame({"sample-id": t.index,
                    "transect-name": ["Baquedano" if s.startswith("BAQ") else "Yungay"
                                      for s in t.index]})
grp.to_csv("transect-md.tsv", sep="\t", index=False)
print(f"  wrote transect-md.tsv for {len(grp)} samples")
PY

echo "############ [1/6] C1: does the chapter's own filter command work? ############"
qiime feature-table filter-samples \
  --i-table "$GG/atacama-counts.qza" \
  --m-metadata-file "$GG/selected-atacama-sample-metadata.tsv" \
  --p-where "[transect-name]='Baquedano'" \
  --o-filtered-table c1-should-fail.qza 2>c1.err
echo "  exit=$? (non-zero CONFIRMS C1: the shipped metadata cannot drive this command)"
head -3 c1.err | sed 's/^/    /'

echo "############ [2/6] C2: the actual split, using derived metadata ############"
for g in Baquedano Yungay; do
  qiime feature-table filter-samples \
    --i-table "$GG/atacama-counts.qza" --m-metadata-file transect-md.tsv \
    --p-where "[transect-name]='$g'" --o-filtered-table "counts-${g}.qza" >/dev/null 2>&1
  n=$(python -c "
import pandas as pd, qiime2
print(qiime2.Artifact.load('counts-${g}.qza').view(pd.DataFrame).shape[0])")
  echo "  $g: $n samples"
done

echo "############ [3/6] C3 + C5: build-groups on the two transects ############"
qiime gglasso build-groups \
  --i-tables counts-Baquedano.qza counts-Yungay.qza \
  --p-check-groups True --o-group-array groups.qza --verbose 2>&1 | tail -20 | sed 's/^/    /'
echo "  exit=$?"

echo "############ [4/6] C4: build-groups on two IDENTICAL tables ############"
cp counts-Baquedano.qza dup.qza
qiime gglasso build-groups \
  --i-tables counts-Baquedano.qza dup.qza \
  --p-check-groups True --o-group-array groups-dup.qza 2>&1 | tail -8 | sed 's/^/    /'
echo "  exit=$? (non-zero CONFIRMS C4)"

echo "############ [5/6] the export workaround: shape of the array ############"
if [[ -f groups.qza ]]; then
  qiime tools export --input-path groups.qza --output-path exported >/dev/null 2>&1
  python - <<'PY'
import numpy as np, zarr, pathlib
p = pathlib.Path("exported/tensor.zip")
if not p.is_file():
    print(f"  no tensor.zip; exported/ holds: {[x.name for x in pathlib.Path('exported').iterdir()]}")
else:
    s = zarr.ZipStore(str(p), mode="r")
    G = np.array(zarr.open(store=s)["tensor"]); s.close()
    print(f"  exported array shape: {G.shape}   (chapter says (2, L, K))")
    print(f"  flattened length: {G.size}")
    print(f"  first 24 values: {' '.join(str(int(v)) for v in G.ravel()[:24])}")
PY
else
  echo "  groups.qza was not produced -- nothing to export"
fi

echo "############ [6/6] C6: gap 1 and gap 3 ############"
echo "  -- gap 1: is there an --i-group-array input on solve-problem?"
qiime gglasso solve-problem --help 2>/dev/null | grep -c 'i-group-array' \
  | xargs -I{} echo "     occurrences of --i-group-array: {} (0 CONFIRMS gap 1)"
qiime gglasso solve-problem --help 2>/dev/null | grep -E 'group-array' | sed 's/^/     /'

echo "  -- gap 3: does --p-non-conforming change the branch with one 2-D covariance?"
qiime gglasso transform-features \
  --i-table counts-Baquedano.qza --i-taxonomy "$GG/classification.qza" \
  --m-sample-metadata-file "$GG/selected-atacama-sample-metadata.tsv" \
  --p-transformation mclr --o-transformed-table mclr-baq.qza >/dev/null 2>&1
qiime gglasso calculate-covariance \
  --i-table mclr-baq.qza --p-method scaled --o-covariance-matrix corr-baq.qza >/dev/null 2>&1
qiime gglasso solve-problem \
  --i-covariance-matrix corr-baq.qza --p-n-samples 32 --p-no-latent \
  --p-non-conforming True --p-group-array 0 1 2 0 1 2 \
  --p-lambda1-min 0.01 --p-lambda1-max 1 --p-n-lambda1 5 \
  --p-gamma 0.01 --o-solution nc.qza 2>&1 | tail -6 | sed 's/^/     /'
echo "     exit=$?"
if [[ -f nc.qza ]]; then
  python - <<'PY'
import os, tempfile, zipfile, zarr
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile("nc.qza") as z:
    i = [n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i, t)
    r = zarr.open(zarr.ZipStore(os.path.join(t, i), mode="r"))
    keys = sorted(r.keys())
    print(f"     solution keys: {keys}")
    if "solution" in r:
        sol = sorted(r["solution"].keys())
        print(f"     solution/ : {sol}")
        import numpy as np
        if "precision_" in sol:
            P = np.asarray(r["solution"]["precision_"])
            print(f"     precision_ shape: {P.shape}  "
                  f"({'2-D => SGL branch, CONFIRMS gap 3' if P.ndim == 2 else '3-D => MGL branch, REFUTES gap 3'})")
PY
fi

echo "############ RESULT: experiment completed -- read the claims above ############"
