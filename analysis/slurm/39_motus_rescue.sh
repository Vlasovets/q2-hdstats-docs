#!/bin/bash
#SBATCH --job-name=q2-motus-rescue
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusrescue_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusrescue_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Can the mOTUs table be filtered into something usable, rather than abandoned?
#
# WHAT ACTUALLY FAILED, AND WHY FILTERING MIGHT STILL HELP.
#
# Stages 36-38 established that neither mOTUs (p=93, 84% zeros) nor Woltka (p=300,
# 11.5% zeros, median depth 2.7 M) supports a network at n=36: eBIC selects the empty
# graph, and the gamma map is a cliff from ~3,400 edges to zero with no middle.
#
# The obvious reading is "n is too small, filtering changes p not n, so give up". That
# is right as far as it goes, but it skips two things neither previous stage tried:
#
#   1. TAXONOMIC COLLAPSE POOLS COUNTS. Every filter so far DISCARDED features. Summing
#      mOTUs to genus keeps their reads and puts them on fewer, deeper features. A table
#      that is 84% zeros at species level may be far denser at genus level. This is the
#      one operation that improves depth rather than trading it away.
#
#   2. p < n IS A DIFFERENT PROBLEM. Every configuration tested so far had p >> n
#      (93, 300, 3115 against 34-36), so the sample covariance was massively
#      rank-deficient and the estimate lived entirely on the penalty. At p=12, n=34 the
#      sample covariance is FULL RANK and the graphical lasso is doing ordinary work.
#      The single p<n case that was run (mOTUs prevalence>=10) was only ever scored at
#      gamma=0.3, and gamma=0.3 has since been shown to force the empty graph on this
#      sample size regardless.
#
# So: collapse to genus and family, sweep prevalence hard enough to cross p<n, and score
# every configuration at every gamma. eBIC is linear in gamma at fixed edge count, so
# two stored curves give the exact selection at any gamma with no extra solving.
#
# SUCCESS CRITERION, fixed in advance: a GENUINE interior eBIC minimum (both neighbours
# strictly worse -- not the first index of a flat plateau of empty graphs) with an edge
# count that is neither 0 nor near-saturated, at a gamma that is not knife-edge.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
PROF="$ROOT/results/motus-profiles"
OUT="$ROOT/results/motus-rescue"
TABLES="$ROOT/results/tables"
MIN_DEPTH="${Q2_MOTUS_MIN_DEPTH:-100}"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motusrescue/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motusrescue/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ [1/3] collapse and filter ############"
python - "$MERGED" "$SHEET" "$PROF" "$SCRATCH" "$MIN_DEPTH" <<'PY' || exit 3
import sys, os, io, csv, glob, zipfile, numpy as np, pandas as pd, qiime2, biom
merged, sheet, prof, out, min_depth = sys.argv[1:6]
t = pd.read_csv(merged, sep="\t", index_col=0)
t = t.loc[t.sum(1) >= float(min_depth)]
print(f"  {t.shape[0]} samples after the depth filter, {t.shape[1]} species-level mOTUs")

tax = {}
for tp in sorted(glob.glob(os.path.join(prof, "*-taxonomy.qza"))):
    with zipfile.ZipFile(tp) as z:
        n = [x for x in z.namelist() if x.endswith("/data/taxonomy.tsv")][0]
        for r in csv.reader(io.StringIO(z.read(n).decode()), delimiter="\t"):
            if r and r[0] != "Feature ID":
                tax[r[0]] = r[1]

def rank_key(fid, depth):
    """Lineage truncated to `depth` ranks; '' when the lineage is too short."""
    parts = [p.strip() for p in str(tax.get(fid, "")).split(";") if p.strip()]
    if len(parts) < depth:
        return ""
    return "; ".join(parts[:depth])

meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
md = pd.DataFrame([{"sample-id": s,
                    "host_age_days": meta.get(s, {}).get("host_age_days", ""),
                    "diagnosis": meta.get(s, {}).get("diagnosis", ""),
                    "host_subject_id": meta.get(s, {}).get("host_subject_id", "")}
                   for s in t.index]).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")

# rank depth: 6 = genus (k p c o f g), 5 = family, 4 = order. 7 would be the mOTU
# itself, i.e. no collapse at all.
LEVELS = {"species": 7, "genus": 6, "family": 5, "order": 4}
combos = []
for lname, d in LEVELS.items():
    if d == 7:
        coll = t.copy()
    else:
        keys = {f: rank_key(f, d) for f in t.columns}
        usable = [f for f in t.columns if keys[f]]
        dropped = t.shape[1] - len(usable)
        coll = t[usable].T.groupby([keys[f] for f in usable]).sum().T
        print(f"  {lname:8s}: {t.shape[1]} -> {coll.shape[1]} features "
              f"({dropped} dropped for a lineage shorter than {d} ranks)")
    for thr in (3, 5, 10, 15, 20, 25):
        sub = coll.loc[:, (coll > 0).sum(0) >= thr]
        sub = sub.loc[:, sub.sum(0) > 0]
        if sub.shape[1] < 4:
            continue
        tag = f"{lname}-prev{thr}"
        tab = biom.Table(sub.T.values, list(sub.columns), list(sub.index))
        qiime2.Artifact.import_data("FeatureTable[Frequency]", tab).save(f"{out}/t-{tag}.qza")
        tx = pd.DataFrame({"Taxon": list(sub.columns)},
                          index=pd.Index(sub.columns, name="Feature ID"))
        qiime2.Artifact.import_data("FeatureData[Taxonomy]", tx).save(f"{out}/x-{tag}.qza")
        z = 100 * float((sub.to_numpy() == 0).mean())
        combos.append(tag)
        print(f"    {tag:18s} p={sub.shape[1]:4d} n={sub.shape[0]}  zeros={z:5.1f}%  "
              f"median depth {sub.sum(1).median():8.0f}  "
              f"{'p<n' if sub.shape[1] < sub.shape[0] else 'p>=n'}")
open(f"{out}/combos.txt", "w").write("\n".join(combos))
PY

echo "############ [2/3] network chain on every configuration ############"
while read -r tag; do
  [[ -z "$tag" ]] && continue
  T="$SCRATCH/t-${tag}.qza"
  N=$(python -c "
import qiime2,pandas as pd; print(qiime2.Artifact.load('$T').view(pd.DataFrame).shape[0])")
  qiime gglasso transform-features --i-table "$T" --i-taxonomy "$SCRATCH/x-${tag}.qza" \
    --p-transformation mclr --p-keep-original-id \
    --o-transformed-table "$SCRATCH/m-${tag}.qza" >"$SCRATCH/tf.log" 2>&1 || {
      echo "  $tag: FAIL transform"; continue; }
  qiime gglasso calculate-covariance --i-table "$SCRATCH/m-${tag}.qza" \
    --p-method scaled --o-covariance-matrix "$SCRATCH/c-${tag}.qza" >/dev/null 2>&1 || {
      echo "  $tag: FAIL covariance"; continue; }
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/c-${tag}.qza" \
    --p-n-samples "$N" --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.3 \
    --o-solution "$OUT/sgl-${tag}.qza" >"$SCRATCH/sp.log" 2>&1 || {
      echo "  $tag: FAIL solve"; tail -3 "$SCRATCH/sp.log" | sed 's/^/      /'; continue; }
  echo "  $tag: solved"
done < "$SCRATCH/combos.txt"

echo "############ [3/3] gamma map for every configuration ############"
python - "$OUT" "$TABLES/motus-rescue-gamma-map.tsv" <<'PY'
import os,sys,glob,zipfile,tempfile,numpy as np,pandas as pd,zarr
out,tsv=sys.argv[1],sys.argv[2]
GAMMAS=(0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.70)
rows=[]
for p in sorted(glob.glob(os.path.join(out,"sgl-*.qza"))):
    tag=os.path.basename(p)[4:-4]
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        ms=r["modelselect_stats"]
        L=np.asarray(ms["LAMBDA"]).ravel(); o=np.argsort(L); L=L[o]
        sp=np.asarray(ms["SP"]).ravel()[o]
        pp=np.asarray(r["solution"]["precision_"]).shape[0]
        n=int(np.asarray(r["solution"]["n_samples"])) if "n_samples" in r["solution"] else 34
        pairs=pp*(pp-1)//2
        c1=np.asarray(ms["BIC"]["0.1"]).ravel()[o]; c3=np.asarray(ms["BIC"]["0.3"]).ravel()[o]
        B=(c3-c1)/0.2; A=c1-0.1*B
        for g in GAMMAS:
            c=A+g*B; j=int(np.argmin(c))
            interior=bool(0<j<len(L)-1 and c[j-1]>c[j]+1e-9 and c[j+1]>c[j]+1e-9)
            e=int(round(sp[j]*pairs))
            rows.append(dict(config=tag,p=pp,n=n,p_lt_n=bool(pp<n),gamma=g,
                             selected_lambda=round(float(L[j]),4),edges=e,
                             density_pct=round(100*e/pairs,3) if pairs else 0.0,
                             genuine_interior_min=interior))
d=pd.DataFrame(rows); d.to_csv(tsv,sep="\t",index=False)

good=d[(d.edges>0)&d.genuine_interior_min&(d.density_pct<25)]
print("\n  configurations with a GENUINE interior minimum and a non-saturated graph:")
if len(good):
    print("   "+good.to_string(index=False).replace("\n","\n   "))
else:
    print("    none")
print("\n  --- best per configuration (widest gamma window with edges) ---")
for cfg, sub in d.groupby("config"):
    ok=sub[(sub.edges>0)&sub.genuine_interior_min]
    p0,n0=sub.p.iloc[0],sub.n.iloc[0]
    if len(ok):
        print(f"    {cfg:18s} p={p0:4d} n={n0}  works at gamma "
              f"{sorted(ok.gamma.tolist())}  edges {sorted(set(ok.edges))}")
    else:
        print(f"    {cfg:18s} p={p0:4d} n={n0}  no gamma gives an interior non-empty graph")
print(f"\n  -> {tsv}")
PY
echo "############ RESULT: complete ############"
