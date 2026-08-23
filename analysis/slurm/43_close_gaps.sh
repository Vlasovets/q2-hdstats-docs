#!/bin/bash
#SBATCH --job-name=q2-close-gaps
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/closegaps_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/closegaps_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# The two things the recommended configuration has not been tested on.
#
# The recommendation is mOTUs species level, top-100 by abundance: p=100, n=34,
# gamma=0.15 selecting lambda=0.30 with 481 edges against a permutation null of 0.2.
# Two parts of the book's workflow have never been run against it, and a tier built on
# an untested half is how the previous tier ended up with placeholder numbers.
#
# GAP 1 -- SPARSE + LOW-RANK. Tier 2 features `--p-latent` prominently: the network is
# reported at rank 0 and again with a latent block absorbing the dominant measured
# gradients. Nothing here has used it. The mu1 -> achieved-rank map has to be measured,
# not guessed, because `--p-rank` is registered but always raises on released GGLasso;
# the rank is an OUTPUT reached by tuning mu1.
#
# Two traps, both documented in the book and both easy to trip:
#   * a single fit requires lambda1 AND lambda2 AND mu1 to each hold exactly one value.
#     Leaving lambda2 unset substitutes a five-point default path and silently turns
#     this into a model-selection run that can move away from the lambda1 we fixed.
#   * `qiime gglasso pca` needs --p-n-components <= the ACHIEVED rank, and it fails
#     without --m-sample-metadata-file even though the signature marks it optional.
#     n_components is therefore derived from the fit, never hardcoded.
#
# GAP 2 -- TRAC. `classo add-taxa` rebuilds the design as log(X)A, where A aggregates
# features up a taxonomy, so coefficients attach to clades rather than single mOTUs.
# It is registered as taking FeatureData[Taxonomy] and converts internally to a tree.
# The mOTUs lineage is a 7-rank string whose leaf is m__<mOTU_id>, unique per feature,
# which is exactly the shape trac needs -- but "should work" is not evidence.
#
# Both gaps are cheap at p=100. The point is to find out now rather than after chapters
# quote numbers from them.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
PROF="$ROOT/results/motus-profiles"
OUT="$ROOT/results/close-gaps"
TABLES="$ROOT/results/tables"
LAMBDA1=0.30            # selected at gamma=0.15 by the abundance sweep

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-closegaps/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-closegaps/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ [1/5] rebuild the recommended table ############"
python - "$MERGED" "$SHEET" "$PROF" "$SCRATCH" <<'PY' || exit 3
import sys, os, io, csv, glob, zipfile, pandas as pd, qiime2, biom
merged, sheet, prof, out = sys.argv[1:5]
t = pd.read_csv(merged, sep="\t", index_col=0)
t = t.loc[t.sum(1) >= 100]
base = t.loc[:, (t > 0).sum(0) >= 2]
sub = base[list(base.sum(0).sort_values(ascending=False).head(100).index)]
sub = sub.loc[:, sub.sum(0) > 0]
print(f"  table: p={sub.shape[1]} n={sub.shape[0]}")
qiime2.Artifact.import_data("FeatureTable[Frequency]",
    biom.Table(sub.T.values, list(sub.columns), list(sub.index))).save(f"{out}/top100.qza")

tax = {}
for tp in sorted(glob.glob(os.path.join(prof, "*-taxonomy.qza"))):
    with zipfile.ZipFile(tp) as z:
        n = [x for x in z.namelist() if x.endswith("/data/taxonomy.tsv")][0]
        for r in csv.reader(io.StringIO(z.read(n).decode()), delimiter="\t"):
            if r and r[0] != "Feature ID": tax[r[0]] = r[1]
td = pd.DataFrame({"Taxon": [tax.get(f, "Unassigned") for f in sub.columns]},
                  index=pd.Index(sub.columns, name="Feature ID"))
qiime2.Artifact.import_data("FeatureData[Taxonomy]", td).save(f"{out}/taxonomy.qza")
leaves = td.Taxon.str.split("; ").str[-1]
print(f"  taxonomy: {len(td)} features, {leaves.nunique()} unique leaves, "
      f"{int((td.Taxon=='Unassigned').sum())} unassigned")

meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
md = pd.DataFrame([{"sample-id": s, "host_age_days": meta[s]["host_age_days"],
                    "diagnosis": meta[s]["diagnosis"],
                    "host_subject_id": meta[s]["host_subject_id"]} for s in sub.index]
                  ).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")
PY

qiime gglasso transform-features --i-table "$SCRATCH/top100.qza" \
  --i-taxonomy "$SCRATCH/taxonomy.qza" --p-transformation mclr --p-keep-original-id \
  --o-transformed-table "$SCRATCH/mclr.qza" >/dev/null 2>&1 || { echo "FAIL transform"; exit 4; }
qiime gglasso calculate-covariance --i-table "$SCRATCH/mclr.qza" --p-method scaled \
  --o-covariance-matrix "$SCRATCH/corr.qza" >/dev/null 2>&1 || { echo "FAIL cov"; exit 4; }

echo "############ [2/5] GAP 1: mu1 -> achieved rank ############"
# lambda2 is pinned explicitly. Leaving it unset substitutes a 5-point default path,
# which turns each of these single fits into a model-selection run.
for MU in 30 20 15 10 7.5 5 3 2 1; do
  TAG=$(echo "$MU" | tr -d '.')
  [[ -f "$OUT/slr-mu${TAG}.qza" ]] && { echo "  [skip] mu=$MU"; continue; }
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/corr.qza" \
    --p-n-samples 34 --p-latent \
    --p-lambda1-min $LAMBDA1 --p-lambda1-max $LAMBDA1 --p-n-lambda1 1 \
    --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
    --p-mu1-min $MU --p-mu1-max $MU --p-n-mu1 1 \
    --o-solution "$OUT/slr-mu${TAG}.qza" >"$SCRATCH/slr-$TAG.log" 2>&1 \
    && echo "  mu=$MU solved" \
    || { echo "  mu=$MU FAILED:"; tail -3 "$SCRATCH/slr-$TAG.log" | sed 's/^/      /'; }
done

echo "############ [3/5] read the rank / edge trade-off ############"
python - "$OUT" "$TABLES/motus-top100-mu-rank-map.tsv" <<'PY'
import os,sys,glob,zipfile,tempfile,numpy as np,pandas as pd,zarr
out,tsv=sys.argv[1],sys.argv[2]
rows=[]
for p in sorted(glob.glob(os.path.join(out,"slr-mu*.qza"))):
    tag=os.path.basename(p)[6:-4]
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        s=r["solution"]
        P=np.asarray(s["precision_"]); pp=P.shape[0]
        L=np.asarray(s["lowrank_"]) if "lowrank_" in s else None
        rank=int(np.linalg.matrix_rank(L, tol=1e-8)) if L is not None else None
        ev=np.sort(np.linalg.eigvalsh(L))[::-1][:5] if L is not None else []
        # P[off] with a boolean mask returns a FLAT array, so .sum(1) raises
        # "axis 1 is out of bounds for array of dimension 1". Mask a copy instead and
        # keep it 2-D. This crashed the first run after all nine fits had succeeded.
        Poff = P.copy(); np.fill_diagonal(Poff, 0.0)
        edges=int((np.abs(Poff)>1e-8).sum()//2)
        nodes=int(((np.abs(Poff)>1e-8).sum(axis=1)>0).sum())
        rows.append(dict(mu1=float(tag.replace("75","7.5") if tag=="75" else tag),
                         p=pp, rank=rank, edges=edges, connected_nodes=nodes,
                         top_eigenvalues=", ".join(f"{v:.3f}" for v in ev)))
d=pd.DataFrame(rows).sort_values("mu1",ascending=False)
d.to_csv(tsv,sep="\t",index=False)
print(d.to_string(index=False))
print(f"\n  -> {tsv}")
print("\n  for reference, the rank-0 fit at the same lambda1 had 481 edges")
PY

echo "############ [4/5] pca on the headline latent fit ############"
python - "$OUT" > "$SCRATCH/best.txt" <<'PY'
import os,sys,glob,zipfile,tempfile,numpy as np,zarr
out=sys.argv[1]; best=None
for p in sorted(glob.glob(os.path.join(out,"slr-mu*.qza"))):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))
        L=r["solution"].get("lowrank_")
        if L is None: continue
        rk=int(np.linalg.matrix_rank(np.asarray(L),tol=1e-8))
        # smallest rank >= 2 -- parsimonious but with something for pca to plot
        if rk>=2 and (best is None or rk<best[1]): best=(p,rk)
print(f"{best[0]} {best[1]}" if best else "NONE 0")
PY
read -r BESTQZA BESTRANK < "$SCRATCH/best.txt"
if [[ "$BESTQZA" != "NONE" ]]; then
  NC=$(( BESTRANK < 3 ? BESTRANK : 3 ))
  echo "  using $(basename "$BESTQZA") (rank $BESTRANK) with --p-n-components $NC"
  qiime gglasso pca --i-solution "$BESTQZA" --i-table "$SCRATCH/mclr.qza" \
    --m-sample-metadata-file "$SCRATCH/outcomes.tsv" --p-n-components "$NC" \
    --o-visualization "$OUT/pca.qzv" >"$SCRATCH/pca.log" 2>&1 \
    && echo "  pca.qzv written" \
    || { echo "  pca FAILED:"; tail -5 "$SCRATCH/pca.log" | sed 's/^/      /'; }
else
  echo "  no fit reached rank >= 2 -- pca skipped (it needs a latent block to plot)"
fi

echo "############ [5/5] GAP 2: trac via add-taxa ############"
qiime classo transform-features --i-features "$SCRATCH/top100.qza" \
  --o-x "$SCRATCH/x.qza" >/dev/null 2>&1 || { echo "  FAIL classo transform"; exit 5; }
if qiime classo add-taxa --i-features "$SCRATCH/x.qza" --i-taxa "$SCRATCH/taxonomy.qza" \
     --o-x "$OUT/x-trac.qza" --o-aweights "$OUT/w-trac.qza" >"$SCRATCH/trac.log" 2>&1; then
  echo "  add-taxa: ok"
  python - "$SCRATCH/x.qza" "$OUT/x-trac.qza" "$OUT/w-trac.qza" <<'PY'
import sys, pandas as pd, qiime2, numpy as np
x0,x1,w = sys.argv[1:4]
a=qiime2.Artifact.load(x0).view(pd.DataFrame); b=qiime2.Artifact.load(x1).view(pd.DataFrame)
print(f"    design before add-taxa: {a.shape[0]} x {a.shape[1]}")
print(f"    design after  add-taxa: {b.shape[0]} x {b.shape[1]}   "
      f"(+{b.shape[1]-a.shape[1]} internal nodes)")
print(f"    example aggregated columns: {list(b.columns)[:3]}")
try:
    ww=qiime2.Artifact.load(w).view(pd.Series)
    print(f"    weights: {len(ww)} values, range {float(ww.min()):.4g}-{float(ww.max()):.4g}")
except Exception as e:
    print(f"    weights: could not view as Series ({type(e).__name__})")
PY
  echo "  --- regress the aggregated design on host_age_days ---"
  # Split into two statements: a heredoc cannot be followed by `||`, which is what made
  # the first version of this script unparseable.
  if qiime classo regress --i-features "$OUT/x-trac.qza" --i-weights "$OUT/w-trac.qza" \
       --m-y-file "$SCRATCH/outcomes.tsv" --m-y-column host_age_days \
       --p-concomitant --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
       --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
       --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale --p-no-stabsel --p-no-lamfixed \
       --o-result "$OUT/reg-trac.qza" >"$SCRATCH/regtrac.log" 2>&1; then
    python - "$OUT/reg-trac.qza" "$OUT/x-trac.qza" <<'PY'
import sys,os,zipfile,tempfile,numpy as np,pandas as pd,zarr,qiime2
p,xq=sys.argv[1:3]
cols=list(qiime2.Artifact.load(xq).view(pd.DataFrame).columns)
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
    n=[x for x in z.namelist() if x.endswith(".zip")][0]; z.extract(n,t)
    r=zarr.open(zarr.ZipStore(os.path.join(t,n),mode="r"))
    s=r["solution"]
    if "CV" in s and "refit" in s["CV"]:
        b=np.asarray(s["CV"]["refit"]).ravel()
        nz=np.argwhere(np.abs(b)>1e-10).ravel()
        print(f"    CV refit selects {len(nz)} of {len(b)} clades")
        for i in nz[:10]:
            name=cols[i] if i < len(cols) else f"<coef {i}>"
            print(f"      {b[i]:+.4f}  {str(name)[:80]}")
PY
  else
    echo "    regress FAILED:"; tail -4 "$SCRATCH/regtrac.log" | sed 's/^/      /'
  fi
else
  echo "  add-taxa FAILED:"; tail -8 "$SCRATCH/trac.log" | sed 's/^/    /'
fi
echo "############ RESULT: complete ############"
