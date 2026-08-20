#!/bin/bash
#SBATCH --job-name=q2-woltka-cmp
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/woltkacmp_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/woltkacmp_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Same study, same samples, different profiler. Is the DATASET too thin, or is mOTUs?
#
# Stage 36 established that a mOTUs profile of these 36 fecal metagenomes cannot support
# a network: eBIC prefers the empty graph at every prevalence threshold, with two of the
# three selections pinned to the path boundary. That is a fact about the mOTUs table --
# 84% zeros, 705 assigned marker-gene inserts from a 2 M-read run -- and NOT necessarily
# a fact about the study.
#
# mOTUs counts only reads landing on ten universal marker genes. Woltka maps against
# whole reference genomes, so it uses a far larger fraction of the same reads. Qiita
# already publishes a Woltka WoLr2 profile of this exact preparation:
#
#     artifact 212377, none.biom, 5,452 features x 96 samples, genome level
#
# Running the identical chain on it answers the question that decides whether the tier
# survives with a different table or has to move to a different study. It costs one
# 376 MB download and a few minutes of solver time -- far cheaper than guessing.
#
# NOTE ON THE DOWNLOAD. Qiita ignores HTTP Range headers: every probe pulls the whole
# object, and public_artifact_download for 212377 returns 4.3 GB that is not a zip. The
# study-wide BIOM bundle IS a proper zip and contains the member we want, so we fetch
# that once and cache it.
#
# Taxonomy: none.biom carries NO observation metadata, so a placeholder taxonomy is
# synthesised purely to satisfy `gglasso transform-features`, whose --i-taxonomy is
# required and never read. That is fine for a network comparison; it would NOT be fine
# for classo add-taxa (trac), which needs a real hierarchy -- flagged, not papered over.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/woltka-compare"
TABLES="$ROOT/results/tables"
CACHE="${Q2_MAP_CACHE:-/lustre/scratch/users/oleg.vlasovets/map13241}"
ZIP="$CACHE/qiita-13241-biom.zip"
ZIP_URL="https://qiita.ucsd.edu/public_download/?data=biom&study_id=13241"
MEMBER="BIOM/212377/none.biom"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-woltkacmp/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-woltkacmp/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES" "$CACHE"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
python - <<'ENVCHECK' || exit 9
import sys, numpy, zarr
bad=[]
if not zarr.__version__.startswith("2.18"): bad.append("zarr "+zarr.__version__)
if not numpy.__version__.startswith("2.4"): bad.append("numpy "+numpy.__version__)
if bad: sys.stderr.write("ENV CHECK FAILED: %s\n"%"; ".join(bad)); sys.exit(9)
ENVCHECK
cd "$SCRATCH"

echo "############ [1/4] fetch the Qiita BIOM bundle (cached) ############"
if [[ -f "$ZIP" ]] && [[ $(stat -c%s "$ZIP") -gt 300000000 ]]; then
  echo "  [skip] $ZIP present ($(stat -c%s "$ZIP") bytes)"
else
  curl -L --fail --retry 3 -o "$ZIP" "$ZIP_URL" 2>/dev/null \
    || { echo "  FAIL download"; exit 2; }
  echo "  downloaded $(stat -c%s "$ZIP") bytes"
fi

echo "############ [2/4] extract the Woltka table and align it to our 36 samples ############"
python - "$ZIP" "$MEMBER" "$SHEET" "$SCRATCH" <<'PY' || exit 3
import sys, io, csv, zipfile, h5py, numpy as np, pandas as pd, qiime2, biom
zp, member, sheet, out = sys.argv[1:5]
z = zipfile.ZipFile(zp)
names = [n for n in z.namelist() if n.endswith(member)]
if not names: sys.exit(f"  FAIL member {member} not in bundle")
raw = z.read(names[0]); print(f"  {names[0]}: {len(raw):,} bytes")

with h5py.File(io.BytesIO(raw), "r") as f:
    obs = [i.decode() for i in f["observation/ids"][:]]
    smp = [i.decode() for i in f["sample/ids"][:]]
    data = f["sample/matrix/data"][:]; indices = f["sample/matrix/indices"][:]
    indptr = f["sample/matrix/indptr"][:]
M = np.zeros((len(smp), len(obs)))
for j in range(len(smp)):
    M[j, indices[indptr[j]:indptr[j+1]]] = data[indptr[j]:indptr[j+1]]
w = pd.DataFrame(M, index=smp, columns=obs)
print(f"  woltka table: {w.shape[0]} samples x {w.shape[1]} features")

rows = list(csv.DictReader(open(sheet), delimiter="\t"))
want = {r["sample_name"]: r for r in rows}
hit = [s for s in w.index if s in want]
print(f"  of our {len(want)} fecal WGS samples, {len(hit)} are present in the Woltka table")
if len(hit) < 20: sys.exit("  FAIL too few samples matched -- check the id convention")
w = w.loc[hit]
w = w.loc[:, w.sum(axis=0) > 0]
print(f"  restricted: {w.shape[0]} samples x {w.shape[1]} features with any count")
print(f"  depth: min {w.sum(1).min():.0f} median {w.sum(1).median():.0f} max {w.sum(1).max():.0f}")

md = pd.DataFrame([{"sample-id": s,
                    "host_age_days": want[s]["host_age_days"],
                    "diagnosis": want[s]["diagnosis"],
                    "host_subject_id": want[s]["host_subject_id"]} for s in w.index]
                  ).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")

# Placeholder taxonomy: --i-taxonomy is required by transform-features and never read.
tax = pd.DataFrame({"Taxon": ["Unassigned"] * w.shape[1]},
                   index=pd.Index(w.columns, name="Feature ID"))
qiime2.Artifact.import_data("FeatureData[Taxonomy]", tax).save(f"{out}/taxonomy.qza")

for thr in (3, 5, 10):
    prev = (w > 0).sum(axis=0)
    sub = w.loc[:, prev >= thr]
    sub = sub.loc[:, sub.sum(axis=0) > 0]
    tab = biom.Table(sub.T.values, list(sub.columns), list(sub.index))
    qiime2.Artifact.import_data("FeatureTable[Frequency]", tab).save(f"{out}/wtable-prev{thr}.qza")
    print(f"  prev>={thr:2d}: p={sub.shape[1]:5d} n={sub.shape[0]}  "
          f"zeros={100*float((sub.to_numpy()==0).mean()):.1f}%")
PY

echo "############ [3/4] identical network chain ############"
for thr in 3 5 10; do
  echo "  --- prevalence >= $thr ---"
  T="$SCRATCH/wtable-prev${thr}.qza"
  N=$(python -c "
import qiime2, pandas as pd
print(qiime2.Artifact.load('$T').view(pd.DataFrame).shape[0])")
  qiime gglasso transform-features --i-table "$T" --i-taxonomy "$SCRATCH/taxonomy.qza" \
    --p-transformation mclr --p-keep-original-id \
    --o-transformed-table "$SCRATCH/wmclr-${thr}.qza" >"$SCRATCH/wtf-${thr}.log" 2>&1 \
    || { echo "    FAIL transform:"; tail -5 "$SCRATCH/wtf-${thr}.log" | sed 's/^/      /'; continue; }
  qiime gglasso calculate-covariance --i-table "$SCRATCH/wmclr-${thr}.qza" \
    --p-method scaled --o-covariance-matrix "$SCRATCH/wcorr-${thr}.qza" >/dev/null 2>&1 \
    || { echo "    FAIL covariance"; continue; }
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/wcorr-${thr}.qza" \
    --p-n-samples "$N" --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.3 \
    --o-solution "$OUT/wsgl-prev${thr}.qza" 2>&1 | tail -2 | sed 's/^/    /'
done

echo "############ [4/4] compare against the mOTUs result ############"
python - "$OUT" "$TABLES/woltka-vs-motus-viability.tsv" <<'PY'
import os, sys, glob, zipfile, tempfile, numpy as np, pandas as pd, zarr
out, tsv = sys.argv[1], sys.argv[2]
rows = []
for p in sorted(glob.glob(os.path.join(out, "wsgl-prev*.qza"))):
    thr = int(os.path.basename(p).split("prev")[1].split(".")[0])
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i = [n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i, t)
        r = zarr.open(zarr.ZipStore(os.path.join(t, i), mode="r"))
        ms = r["modelselect_stats"]
        L = np.asarray(ms["LAMBDA"]).ravel(); o = np.argsort(L); L = L[o]
        sp = np.asarray(ms["SP"]).ravel()[o]
        gk = sorted(ms["BIC"].keys(), key=float)
        g = "0.3" if "0.3" in gk else gk[len(gk)//2]
        c = np.asarray(ms["BIC"][g]).ravel()[o]
        P = np.asarray(r["solution"]["precision_"]); pp = P.shape[0]
        pairs = pp*(pp-1)//2
        off = ~np.eye(pp, dtype=bool)
        edges = int((np.abs(P[off]) > 1e-8).sum() // 2)
        j = int(np.argmin(c))
        interior = (0 < j < len(L)-1 and c[j-1] > c[j]+1e-9 and c[j+1] > c[j]+1e-9)
        rows.append(dict(profiler="woltka", prevalence=thr, p=pp, n_pairs=pairs,
                         selected_lambda=float(L[j]), edges_at_selection=edges,
                         sparsity_at_selection=float(sp[j]),
                         genuine_interior_min=interior,
                         max_edges_on_path=int(sp.max()*pairs)))
d = pd.DataFrame(rows)
prev = os.path.join(os.path.dirname(tsv), "motus-model-viability.tsv")
if os.path.isfile(prev):
    m = pd.read_csv(prev, sep="\t"); m.insert(0, "profiler", "motus")
    keep = ["profiler","prevalence","p","edges","selected_lambda"]
    print("\n  mOTUs (stage 36):")
    print("   " + m[keep].to_string(index=False).replace("\n", "\n   "))
print("\n  Woltka (this stage):")
print("   " + d.to_string(index=False).replace("\n", "\n   "))
d.to_csv(tsv, sep="\t", index=False)
print(f"\n  -> {tsv}")
alive = d[(d.edges_at_selection > 0) & d.genuine_interior_min]
print("\n  --- VERDICT ---")
if len(alive):
    print(f"    Woltka DOES support a network: {len(alive)} of {len(d)} thresholds give a")
    print("    non-empty graph at a genuine interior eBIC minimum. The study is usable;")
    print("    mOTUs was the limiting factor.")
else:
    print("    Woltka does NOT rescue it either -- eBIC still prefers the empty graph.")
    print("    That points at the STUDY (34 samples, 18 infants), not the profiler.")
PY
echo "############ RESULT: comparison complete ############"
