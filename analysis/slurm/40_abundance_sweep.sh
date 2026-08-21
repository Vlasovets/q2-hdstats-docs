#!/bin/bash
#SBATCH --job-name=q2-abund-sweep
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/abundsweep_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/abundsweep_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Find a SPECIES-LEVEL, p > n configuration that both plugins can illustrate honestly.
#
# The genus-collapsed rescue (stage 39) works but lands at p=12 < n, which is tier-1
# scale and does not exercise the high-dimensional machinery the book exists to teach.
# The target here is the tier-2 regime: p comfortably above n, a readable sparse network,
# and a log-contrast regression that selects a sensible number of features.
#
# WHY ABUNDANCE RATHER THAN PREVALENCE. The mOTUs table is 84% zeros and 182 of its
# cells are exactly 1. Two features that are absent in the same 28 samples have a large
# sample correlation driven entirely by shared zeros, not by biology, and a prevalence
# filter KEEPS such a pair as long as each clears the threshold. An abundance filter
# removes them, because a feature whose whole signal is a handful of 1s carries almost
# no counts. That is the mechanism behind the 12-40% densities seen in stage 39, and it
# is the thing to attack.
#
# THE CHECK THAT DECIDES EVERYTHING: A PERMUTATION NULL.
# For each surviving configuration, every feature column is permuted INDEPENDENTLY
# across samples. That destroys all between-feature dependence while preserving each
# feature's marginal distribution exactly -- same counts, same number of zeros. The true
# precision matrix of the permuted data is diagonal, so every edge recovered from it is
# a false positive. Comparing real against permuted edge counts at the same lambda is
# the difference between "we found 50 edges" and "we found 50 edges where noise alone
# gives 45". Nothing in stages 36-39 tested this, and the shared-zero mechanism above is
# exactly the kind of artefact it catches.
#
# Both sources are swept, because the earlier comparison showed they behave very
# differently: mOTUs species (378 features, median depth 740) and Woltka genome-level
# (4432 features, median depth 2.7 M).

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
CACHE="${Q2_MAP_CACHE:-/lustre/scratch/users/oleg.vlasovets/map13241}"
ZIP="$CACHE/qiita-13241-biom.zip"
OUT="$ROOT/results/abundance-sweep"
TABLES="$ROOT/results/tables"
NPERM="${Q2_NPERM:-5}"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-abundsweep/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-abundsweep/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ [1/4] build abundance-filtered tables ############"
python - "$MERGED" "$ZIP" "$SHEET" "$SCRATCH" "$NPERM" <<'PY' || exit 3
import sys, os, io, csv, zipfile, h5py, numpy as np, pandas as pd, qiime2, biom
merged, zp, sheet, out, nperm = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
rng = np.random.default_rng(1)

meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
by_name = {r["sample_name"]: r for r in meta.values()}

# --- source A: mOTUs species level, indexed by run accession
m = pd.read_csv(merged, sep="\t", index_col=0)
m = m.loc[m.sum(1) >= 100]
print(f"  mOTUs species : {m.shape[0]} samples x {m.shape[1]} features")

# --- source B: Woltka genome level, indexed by sample_name -> remap to run accession
z = zipfile.ZipFile(zp)
name = [n for n in z.namelist() if n.endswith("BIOM/212377/none.biom")][0]
with h5py.File(io.BytesIO(z.read(name)), "r") as f:
    obs=[i.decode() for i in f["observation/ids"][:]]; smp=[i.decode() for i in f["sample/ids"][:]]
    data=f["sample/matrix/data"][:]; idx=f["sample/matrix/indices"][:]; ptr=f["sample/matrix/indptr"][:]
W=np.zeros((len(smp),len(obs)))
for j in range(len(smp)): W[j, idx[ptr[j]:ptr[j+1]]]=data[ptr[j]:ptr[j+1]]
w=pd.DataFrame(W,index=smp,columns=obs)
name2run={r["sample_name"]: r["run_accession"] for r in meta.values()}
w=w.loc[[s for s in w.index if s in name2run]]
w.index=[name2run[s] for s in w.index]
w=w.loc[[s for s in m.index if s in w.index]]      # same 34 samples as mOTUs
w=w.loc[:, w.sum(0)>0]
print(f"  woltka genome : {w.shape[0]} samples x {w.shape[1]} features")

md = pd.DataFrame([{"sample-id": s,
                    "host_age_days": meta[s]["host_age_days"],
                    "diagnosis": meta[s]["diagnosis"],
                    "host_subject_id": meta[s]["host_subject_id"]} for s in m.index]
                  ).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")

def emit(tag, sub, permute=False, seed=0):
    if permute:
        r = np.random.default_rng(seed)
        a = sub.to_numpy().copy()
        for j in range(a.shape[1]):
            r.shuffle(a[:, j])          # each feature independently -> diagonal truth
        sub = pd.DataFrame(a, index=sub.index, columns=sub.columns)
    tab = biom.Table(sub.T.values, list(sub.columns), list(sub.index))
    qiime2.Artifact.import_data("FeatureTable[Frequency]", tab).save(f"{out}/t-{tag}.qza")
    tx = pd.DataFrame({"Taxon": list(sub.columns)},
                      index=pd.Index(sub.columns, name="Feature ID"))
    qiime2.Artifact.import_data("FeatureData[Taxonomy]", tx).save(f"{out}/x-{tag}.qza")
    return sub

combos = []
SOURCES = {"motus": m, "woltka": w}
TOPNS   = [40, 60, 80, 100, 150]
for src, tbl in SOURCES.items():
    n = tbl.shape[0]
    # A floor of 2 samples removes strict singletons before ranking, so "most abundant"
    # cannot be satisfied by one enormous count in one sample.
    base = tbl.loc[:, (tbl > 0).sum(0) >= 2]
    order = base.sum(0).sort_values(ascending=False)
    for N in TOPNS:
        if N > base.shape[1]:
            continue
        sub = base[list(order.head(N).index)]
        sub = sub.loc[:, sub.sum(0) > 0]
        if sub.shape[1] < 10:
            continue
        tag = f"{src}-top{N}"
        s = emit(tag, sub)
        combos.append(tag)
        print(f"    {tag:16s} p={s.shape[1]:4d} n={s.shape[0]}  "
              f"{'p>n' if s.shape[1] > s.shape[0] else 'p<=n'}  "
              f"zeros={100*float((s.to_numpy()==0).mean()):5.1f}%  "
              f"median depth {s.sum(1).median():9.0f}  "
              f"min feature total {s.sum(0).min():.0f}")
        for k in range(nperm):
            emit(f"{tag}-perm{k}", sub, permute=True, seed=100 + k)
            combos.append(f"{tag}-perm{k}")
open(f"{out}/combos.txt", "w").write("\n".join(combos))
print(f"  {len(combos)} tables written ({nperm} permutations each)")
PY

echo "############ [2/4] solve every table ############"
while read -r tag; do
  [[ -z "$tag" ]] && continue
  [[ -f "$OUT/sgl-${tag}.qza" ]] && { echo "  [skip] $tag"; continue; }
  T="$SCRATCH/t-${tag}.qza"
  N=$(python -c "
import qiime2,pandas as pd; print(qiime2.Artifact.load('$T').view(pd.DataFrame).shape[0])")
  qiime gglasso transform-features --i-table "$T" --i-taxonomy "$SCRATCH/x-${tag}.qza" \
    --p-transformation mclr --p-keep-original-id \
    --o-transformed-table "$SCRATCH/m-${tag}.qza" >/dev/null 2>&1 || { echo "  $tag FAIL transform"; continue; }
  qiime gglasso calculate-covariance --i-table "$SCRATCH/m-${tag}.qza" \
    --p-method scaled --o-covariance-matrix "$SCRATCH/c-${tag}.qza" >/dev/null 2>&1 || { echo "  $tag FAIL cov"; continue; }
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/c-${tag}.qza" \
    --p-n-samples "$N" --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.3 \
    --o-solution "$OUT/sgl-${tag}.qza" >/dev/null 2>&1 || { echo "  $tag FAIL solve"; continue; }
  echo "  $tag solved"
done < "$SCRATCH/combos.txt"

echo "############ [3/4] real vs permuted, across gamma ############"
python - "$OUT" "$TABLES/abundance-sweep-gamma-map.tsv" <<'PY'
import os,sys,glob,re,zipfile,tempfile,numpy as np,pandas as pd,zarr
out,tsv=sys.argv[1],sys.argv[2]
GAMMAS=(0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)
def curve(p):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        ms=r["modelselect_stats"]
        L=np.asarray(ms["LAMBDA"]).ravel(); o=np.argsort(L); L=L[o]
        sp=np.asarray(ms["SP"]).ravel()[o]
        c1=np.asarray(ms["BIC"]["0.1"]).ravel()[o]; c3=np.asarray(ms["BIC"]["0.3"]).ravel()[o]
        pp=np.asarray(r["solution"]["precision_"]).shape[0]
        return L,sp,(c3-c1)/0.2,c1-0.1*((c3-c1)/0.2),pp
rows=[]
for p in sorted(glob.glob(os.path.join(out,"sgl-*.qza"))):
    tag=os.path.basename(p)[4:-4]
    m=re.match(r"^(.*?)(?:-perm(\d+))?$", tag)
    base,perm=m.group(1),m.group(2)
    L,sp,B,A,pp=curve(p); pairs=pp*(pp-1)//2
    for g in GAMMAS:
        c=A+g*B; j=int(np.argmin(c))
        interior=bool(0<j<len(L)-1 and c[j-1]>c[j]+1e-9 and c[j+1]>c[j]+1e-9)
        rows.append(dict(config=base,perm=(perm is not None),p=pp,pairs=pairs,gamma=g,
                         selected_lambda=round(float(L[j]),3),
                         edges=int(round(sp[j]*pairs)),
                         density_pct=round(100*sp[j],3),
                         interior=interior))
d=pd.DataFrame(rows); d.to_csv(tsv,sep="\t",index=False)

real=d[~d.perm]; null=d[d.perm]
agg=null.groupby(["config","gamma"]).edges.agg(["mean","max"]).reset_index() \
        .rename(columns={"mean":"null_mean_edges","max":"null_max_edges"})
cmp=real.merge(agg,on=["config","gamma"],how="left")
cmp["excess"]=cmp.edges-cmp.null_mean_edges
print("\n  REAL vs PERMUTED NULL (null = each feature shuffled independently)")
print("  " + cmp[["config","p","gamma","selected_lambda","edges","null_mean_edges",
                  "null_max_edges","excess","density_pct","interior"]]
      .to_string(index=False).replace("\n","\n  "))
cmp.to_csv(tsv.replace(".tsv","-vs-null.tsv"),sep="\t",index=False)

good=cmp[(cmp.edges>0)&cmp.interior&(cmp.p>34)&
         (cmp.edges>2*cmp.null_mean_edges.fillna(0))&(cmp.edges>cmp.null_max_edges.fillna(0))]
print("\n  --- p>n configurations whose edge count BEATS the null ---")
print("  " + (good.to_string(index=False).replace("\n","\n  ") if len(good) else "none"))
print(f"\n  -> {tsv}")
PY

echo "############ [4/4] classo regression for the surviving configurations ############"
python - "$TABLES/abundance-sweep-gamma-map-vs-null.tsv" "$SCRATCH/pick.txt" <<'PY'
import sys, pandas as pd
try:
    d=pd.read_csv(sys.argv[1],sep="\t")
    g=d[(d.edges>0)&d.interior&(d.p>34)&(d.edges>2*d.null_mean_edges.fillna(0))]
    picks=sorted(set(g.config))[:4] if len(g) else sorted(set(d[d.p>34].config))[:4]
except Exception:
    picks=[]
open(sys.argv[2],"w").write("\n".join(picks))
print("  regressing:", picks)
PY
while read -r tag; do
  [[ -z "$tag" ]] && continue
  qiime classo transform-features --i-features "$SCRATCH/t-${tag}.qza" \
    --o-x "$SCRATCH/cx-${tag}.qza" >/dev/null 2>&1 || { echo "  $tag FAIL classo transform"; continue; }
  qiime classo regress --i-features "$SCRATCH/cx-${tag}.qza" \
    --m-y-file "$SCRATCH/outcomes.tsv" --m-y-column host_age_days \
    --p-concomitant --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale --p-no-stabsel --p-no-lamfixed \
    --o-result "$OUT/reg-${tag}.qza" >/dev/null 2>&1 || { echo "  $tag FAIL regress"; continue; }
  python - "$OUT/reg-${tag}.qza" "$tag" <<'PY'
import sys,os,zipfile,tempfile,numpy as np,zarr
p,tag=sys.argv[1],sys.argv[2]
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
    n=[x for x in z.namelist() if x.endswith(".zip")][0]; z.extract(n,t)
    r=zarr.open(zarr.ZipStore(os.path.join(t,n),mode="r"))
    s=r["solution"]
    if "CV" in s and "refit" in s["CV"]:
        b=np.asarray(s["CV"]["refit"]).ravel()
        nz=int((np.abs(b)>1e-10).sum())
        print(f"  {tag}: CV-refit selects {nz} of {len(b)} coefficients")
    else:
        print(f"  {tag}: no CV refit stored")
PY
done < "$SCRATCH/pick.txt"
echo "############ RESULT: complete ############"
