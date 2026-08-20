#!/bin/bash
#SBATCH --job-name=q2-woltka-300
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/woltka300_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/woltka300_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# The Woltka table at a tutorial-sized p, matching how tier 2 was built.
#
# Stage 37 showed the Woltka profile of these same 36 samples is vastly richer than the
# mOTUs one -- median depth 2,735,878 against 740, and 4,432 features with any count.
# But it then ran for over 90 minutes on a single lambda path at p=3,115, because the
# graphical lasso is roughly cubic in p. That runtime is itself a result: a p of three
# thousand is not a tutorial, it is an overnight job, and a reader cannot follow it.
#
# Tier 2 faced exactly this and answered it by taking the top 300 ASVs by abundance.
# Doing the same here keeps the comparison honest -- same filter rule, same lambda grid,
# same gamma -- and makes the numbers directly comparable to the published tier-2
# result of 216 edges at lambda 0.8.
#
# Two prevalence/abundance variants, because "top 300" is a choice and choices get
# examined here rather than assumed:
#     top 300 by total count, among features present in >= 3 samples
#     top 300 by prevalence,  ties broken by total count
#
# WHAT WOULD COUNT AS SUCCESS: a non-empty graph at a GENUINE interior eBIC minimum --
# both neighbours strictly worse, not the first index of a flat plateau of empty
# solutions, which is the trap stage 36 fell into and now guards against.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/woltka-top300"
TABLES="$ROOT/results/tables"
CACHE="${Q2_MAP_CACHE:-/lustre/scratch/users/oleg.vlasovets/map13241}"
ZIP="$CACHE/qiita-13241-biom.zip"
MEMBER="BIOM/212377/none.biom"
TOPN="${Q2_WOLTKA_TOPN:-300}"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-woltka300/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-woltka300/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"
[[ -f "$ZIP" ]] || { echo "missing $ZIP -- run stage 37 first (it caches the bundle)"; exit 2; }

echo "############ [1/3] top-$TOPN Woltka tables ############"
python - "$ZIP" "$MEMBER" "$SHEET" "$SCRATCH" "$TOPN" <<'PY' || exit 3
import sys, io, csv, zipfile, h5py, numpy as np, pandas as pd, qiime2, biom
zp, member, sheet, out, topn = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
z = zipfile.ZipFile(zp)
name = [n for n in z.namelist() if n.endswith(member)][0]
with h5py.File(io.BytesIO(z.read(name)), "r") as f:
    obs=[i.decode() for i in f["observation/ids"][:]]
    smp=[i.decode() for i in f["sample/ids"][:]]
    data=f["sample/matrix/data"][:]; idx=f["sample/matrix/indices"][:]
    ptr=f["sample/matrix/indptr"][:]
M=np.zeros((len(smp),len(obs)))
for j in range(len(smp)): M[j, idx[ptr[j]:ptr[j+1]]]=data[ptr[j]:ptr[j+1]]
w=pd.DataFrame(M,index=smp,columns=obs)

rows=list(csv.DictReader(open(sheet),delimiter="\t"))
want={r["sample_name"]:r for r in rows}
w=w.loc[[s for s in w.index if s in want]]
w=w.loc[:, w.sum(0)>0]
print(f"  fecal subset: {w.shape[0]} samples x {w.shape[1]} features")

md=pd.DataFrame([{"sample-id":s,"host_age_days":want[s]["host_age_days"],
                  "diagnosis":want[s]["diagnosis"],
                  "host_subject_id":want[s]["host_subject_id"]} for s in w.index]
                ).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")

prev=(w>0).sum(0); tot=w.sum(0)
base=w.loc[:, prev>=3]
sel={
  "abund": tot[base.columns].sort_values(ascending=False).head(topn).index,
  "preval": pd.DataFrame({"p":prev[base.columns],"t":tot[base.columns]})
            .sort_values(["p","t"],ascending=False).head(topn).index,
}
for tag, cols in sel.items():
    sub=base[list(cols)]
    sub=sub.loc[:, sub.sum(0)>0]
    tab=biom.Table(sub.T.values, list(sub.columns), list(sub.index))
    qiime2.Artifact.import_data("FeatureTable[Frequency]", tab).save(f"{out}/w300-{tag}.qza")
    print(f"  {tag:7s}: p={sub.shape[1]} n={sub.shape[0]} "
          f"zeros={100*float((sub.to_numpy()==0).mean()):.1f}% "
          f"median depth {sub.sum(1).median():.0f}")
    t=pd.DataFrame({"Taxon":["Unassigned"]*sub.shape[1]},
                   index=pd.Index(sub.columns,name="Feature ID"))
    qiime2.Artifact.import_data("FeatureData[Taxonomy]", t).save(f"{out}/tax-{tag}.qza")
PY

echo "############ [2/3] network chain on each ############"
for tag in abund preval; do
  echo "  --- top-$TOPN by $tag ---"
  T="$SCRATCH/w300-${tag}.qza"
  N=$(python -c "
import qiime2,pandas as pd; print(qiime2.Artifact.load('$T').view(pd.DataFrame).shape[0])")
  qiime gglasso transform-features --i-table "$T" --i-taxonomy "$SCRATCH/tax-${tag}.qza" \
    --p-transformation mclr --p-keep-original-id \
    --o-transformed-table "$SCRATCH/m-${tag}.qza" >"$SCRATCH/tf-${tag}.log" 2>&1 \
    || { echo "    FAIL transform:"; tail -5 "$SCRATCH/tf-${tag}.log" | sed 's/^/      /'; continue; }
  qiime gglasso calculate-covariance --i-table "$SCRATCH/m-${tag}.qza" \
    --p-method scaled --o-covariance-matrix "$SCRATCH/c-${tag}.qza" >/dev/null 2>&1 \
    || { echo "    FAIL covariance"; continue; }
  t0=$SECONDS
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/c-${tag}.qza" \
    --p-n-samples "$N" --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.3 \
    --o-solution "$OUT/w300-sgl-${tag}.qza" 2>&1 | tail -2 | sed 's/^/    /'
  echo "    wall $((SECONDS-t0))s"
done

echo "############ [3/3] verdict ############"
python - "$OUT" "$TABLES/woltka-top300-viability.tsv" <<'PY'
import os,sys,glob,zipfile,tempfile,numpy as np,pandas as pd,zarr
out,tsv=sys.argv[1],sys.argv[2]
rows=[]
for p in sorted(glob.glob(os.path.join(out,"w300-sgl-*.qza"))):
    tag=os.path.basename(p).split("w300-sgl-")[1].split(".")[0]
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        ms=r["modelselect_stats"]
        L=np.asarray(ms["LAMBDA"]).ravel(); o=np.argsort(L); L=L[o]
        sp=np.asarray(ms["SP"]).ravel()[o]
        gk=sorted(ms["BIC"].keys(),key=float); g="0.3" if "0.3" in gk else gk[len(gk)//2]
        c=np.asarray(ms["BIC"][g]).ravel()[o]
        P=np.asarray(r["solution"]["precision_"]); pp=P.shape[0]
        pairs=pp*(pp-1)//2; off=~np.eye(pp,dtype=bool)
        edges=int((np.abs(P[off])>1e-8).sum()//2)
        j=int(np.argmin(c))
        interior=(0<j<len(L)-1 and c[j-1]>c[j]+1e-9 and c[j+1]>c[j]+1e-9)
        print(f"\n  --- top-300 by {tag}: p={pp} ---")
        print("    lambda  edges     eBIC")
        for Lv,Sv,Cv in zip(L,sp,c):
            mark=" <-- min" if abs(Cv-c[j])<1e-9 else ""
            print(f"    {Lv:6.3f} {Sv*pairs:7.0f} {Cv:10.1f}{mark}")
        rows.append(dict(selection=tag,p=pp,selected_lambda=float(L[j]),
                         edges_at_selection=edges,
                         density_pct=round(100*edges/pairs,3),
                         genuine_interior_min=interior,
                         max_edges_on_path=int(sp.max()*pairs)))
d=pd.DataFrame(rows); d.to_csv(tsv,sep="\t",index=False)
print("\n" + d.to_string(index=False))
print(f"\n  -> {tsv}")
ok=d[(d.edges_at_selection>0)&d.genuine_interior_min]
print("\n  --- VERDICT ---")
if len(ok):
    print(f"    Woltka at p=300 SUPPORTS a network: {len(ok)} of {len(d)} selections give a")
    print("    non-empty graph at a genuine interior eBIC minimum.")
else:
    print("    Still empty at a genuine interior minimum. With 36 samples the study,")
    print("    not the profiler, is the binding constraint.")
PY
echo "############ RESULT: complete ############"
