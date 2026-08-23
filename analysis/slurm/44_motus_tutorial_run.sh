#!/bin/bash
#SBATCH --job-name=q2-motus-tut
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motustut_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motustut_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# THE CANONICAL TUTORIAL RUN. This stage IS the chapter.
#
# Every command below appears verbatim in the shotgun section of the high-dimensional
# chapter, and every output the chapter quotes is captured here. That coupling is the
# point: the previous metagenomics tier was deleted because it printed commands nobody
# had run and numbers nobody had produced. A chapter whose commands live in a stage
# cannot drift from its own results without the stage failing.
#
# Starts from the MERGED mOTUs table, which is committed (analysis/results/motus-merged/
# table.tsv, 68 KB). Profiling the 36 metagenomes that produced it needs 13 GiB of reads
# and a 2.9 GiB reference database, so the chapter describes that step and does not ask
# the reader to repeat it -- stages 32 to 35 do, for anyone who wants to.
#
# Publishes into docs/_data/motus/ the small artifacts a reader needs to follow along:
# the filtered table, its taxonomy, and the outcomes file.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
PROF="$ROOT/results/motus-profiles"
OUT="$ROOT/results/motus-tutorial"
PUB="$REPO/docs/_data/motus"
TABLES="$ROOT/results/tables"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motustut/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motustut/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$PUB" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$OUT"

echo "############ [1/6] the reader's starting artifacts ############"
python - "$MERGED" "$SHEET" "$PROF" "$PUB" "$TABLES" <<'PY' || exit 3
import sys, os, io, csv, glob, zipfile, pandas as pd, qiime2, biom
merged, sheet, prof, pub, tables = sys.argv[1:6]
t = pd.read_csv(merged, sep="\t", index_col=0)
n0 = t.shape[0]
t = t.loc[t.sum(1) >= 100]
print(f"  depth filter (>=100 assigned counts): {n0} -> {t.shape[0]} samples")
base = t.loc[:, (t > 0).sum(0) >= 2]
print(f"  prevalence floor (>=2 samples): {t.shape[1]} -> {base.shape[1]} features")
sub = base[list(base.sum(0).sort_values(ascending=False).head(100).index)]
sub = sub.loc[:, sub.sum(0) > 0]
zeros = 100 * float((sub.to_numpy() == 0).mean())
print(f"  top-100 by total abundance: p={sub.shape[1]} n={sub.shape[0]}  "
      f"zeros={zeros:.1f}%  median depth {sub.sum(1).median():.0f}")

qiime2.Artifact.import_data("FeatureTable[Frequency]",
    biom.Table(sub.T.values, list(sub.columns), list(sub.index))).save(f"{pub}/motus-top100-table.qza")

tax = {}
for tp in sorted(glob.glob(os.path.join(prof, "*-taxonomy.qza"))):
    with zipfile.ZipFile(tp) as z:
        nm = [x for x in z.namelist() if x.endswith("/data/taxonomy.tsv")][0]
        for r in csv.reader(io.StringIO(z.read(nm).decode()), delimiter="\t"):
            if r and r[0] != "Feature ID": tax[r[0]] = r[1]
td = pd.DataFrame({"Taxon": [tax[f] for f in sub.columns]},
                  index=pd.Index(sub.columns, name="Feature ID"))
qiime2.Artifact.import_data("FeatureData[Taxonomy]", td).save(f"{pub}/motus-top100-taxonomy.qza")
print(f"  taxonomy: {len(td)} lineages, {td.Taxon.str.split('; ').str[-1].nunique()} unique leaves")

meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
md = pd.DataFrame([{"sample-id": s, "host_age_days": meta[s]["host_age_days"],
                    "diagnosis": meta[s]["diagnosis"],
                    "host_subject_id": meta[s]["host_subject_id"]} for s in sub.index]
                  ).set_index("sample-id")
md.to_csv(f"{pub}/motus-outcomes.tsv", sep="\t")
print(f"  outcomes: {md.shape[0]} rows, {md.host_subject_id.nunique()} infants")

# The table the chapter renders as a csv-table.
pd.DataFrame({
    "property": ["samples (n)", "features (p)", "p/n", "zeros", "median assigned depth",
                 "infants", "samples per infant", "diagnosis balance"],
    "value": [sub.shape[0], sub.shape[1], round(sub.shape[1]/sub.shape[0], 1),
              f"{zeros:.1f}%", f"{sub.sum(1).median():.0f}", md.host_subject_id.nunique(),
              2, "17 asthma / 17 control"],
}).to_csv(f"{tables}/motus-tutorial-shape.tsv", sep="\t", index=False)
PY

T="$PUB/motus-top100-table.qza"; X="$PUB/motus-top100-taxonomy.qza"; Y="$PUB/motus-outcomes.tsv"

echo "############ [2/6] gglasso: transform + covariance ############"
echo "\$ qiime gglasso transform-features --p-transformation mclr --p-keep-original-id"
qiime gglasso transform-features --i-table "$T" --i-taxonomy "$X" \
  --p-transformation mclr --p-keep-original-id --o-transformed-table "$OUT/mclr.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime gglasso calculate-covariance --p-method scaled"
qiime gglasso calculate-covariance --i-table "$OUT/mclr.qza" --p-method scaled \
  --o-covariance-matrix "$OUT/correlation.qza" 2>&1 | sed 's/^/  /'

echo "############ [3/6] gglasso: solve the lambda path (gamma = 0.15) ############"
echo "\$ qiime gglasso solve-problem --p-n-samples 34 --p-no-latent --p-gamma 0.15"
qiime gglasso solve-problem --i-covariance-matrix "$OUT/correlation.qza" \
  --p-n-samples 34 --p-no-latent --p-path-scale linear \
  --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.15 \
  --o-solution "$OUT/motus-sgl-path.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime gglasso summarize"
qiime gglasso summarize --i-solution "$OUT/motus-sgl-path.qza" --p-label-size 6pt \
  --o-visualization "$OUT/motus-sgl-summary.qzv" 2>&1 | sed 's/^/  /'

echo "############ [4/6] gglasso: the latent block ############"
echo "\$ qiime gglasso solve-problem --p-latent --p-mu1 5"
qiime gglasso solve-problem --i-covariance-matrix "$OUT/correlation.qza" \
  --p-n-samples 34 --p-latent \
  --p-lambda1-min 0.30 --p-lambda1-max 0.30 --p-n-lambda1 1 \
  --p-lambda2-min 0.1 --p-lambda2-max 0.1 --p-n-lambda2 1 \
  --p-mu1-min 5 --p-mu1-max 5 --p-n-mu1 1 \
  --o-solution "$OUT/motus-slr-mu5.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime gglasso pca --p-n-components 3"
qiime gglasso pca --i-solution "$OUT/motus-slr-mu5.qza" --i-table "$OUT/mclr.qza" \
  --m-sample-metadata-file "$Y" --p-n-components 3 \
  --o-visualization "$OUT/motus-latent-pca.qzv" 2>&1 | sed 's/^/  /'

echo "############ [5/6] classo: trac regression on host_age_days ############"
echo "\$ qiime classo transform-features"
qiime classo transform-features --i-features "$T" --o-x "$OUT/classo-x.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime classo add-taxa"
qiime classo add-taxa --i-features "$OUT/classo-x.qza" --i-taxa "$X" \
  --o-x "$OUT/classo-x-trac.qza" --o-aweights "$OUT/classo-w-trac.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime classo regress --m-y-column host_age_days --p-no-cv-one-se"
qiime classo regress --i-features "$OUT/classo-x-trac.qza" --i-weights "$OUT/classo-w-trac.qza" \
  --m-y-file "$Y" --m-y-column host_age_days \
  --p-concomitant --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
  --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale --p-no-stabsel --p-no-lamfixed \
  --o-result "$OUT/motus-regress-age.qza" 2>&1 | sed 's/^/  /'
echo "\$ qiime classo summarize"
qiime classo summarize --i-problem "$OUT/motus-regress-age.qza" --p-maxplot 40 \
  --o-visualization "$OUT/motus-regress-summary.qzv" 2>&1 | sed 's/^/  /'

echo "############ [6/6] emit the numbers the chapter quotes ############"
python - "$OUT" "$TABLES" <<'PY'
import os,sys,zipfile,tempfile,numpy as np,pandas as pd,zarr,qiime2
out,tables=sys.argv[1:3]
def prob(q):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(q) as z:
        i=[n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i,t)
        return zarr.open(zarr.ZipStore(os.path.join(t,i),mode="r"))

r=prob(f"{out}/motus-sgl-path.qza"); ms=r["modelselect_stats"]
L=np.asarray(ms["LAMBDA"]).ravel(); o=np.argsort(L); L=L[o]
sp=np.asarray(ms["SP"]).ravel()[o]
c1=np.asarray(ms["BIC"]["0.1"]).ravel()[o]; c3=np.asarray(ms["BIC"]["0.3"]).ravel()[o]
B=(c3-c1)/0.2; A=c1-0.1*B; c=A+0.15*B
pp=np.asarray(r["solution"]["precision_"]).shape[0]; pairs=pp*(pp-1)//2
j=int(np.argmin(c))
pd.DataFrame({"lambda1":L,"edges":(sp*pairs).round().astype(int),
              "sparsity":sp.round(6),"ebic_gamma0.15":c.round(3)}
             ).to_csv(f"{tables}/motus-lambda-path.tsv",sep="\t",index=False)
print(f"  eBIC(gamma=0.15) minimum at lambda1={L[j]:.2f}, "
      f"{sp[j]*pairs:.0f} edges of {pairs}, {100*sp[j]:.2f}% density")

rl=prob(f"{out}/motus-slr-mu5.qza")["solution"]
P=np.asarray(rl["precision_"]); Lo=np.asarray(rl["lowrank_"])
Po=P.copy(); np.fill_diagonal(Po,0.0)
ev=np.sort(np.linalg.eigvalsh(Lo))[::-1][:4]
print(f"  latent fit mu1=5: rank {np.linalg.matrix_rank(Lo,tol=1e-8)}, "
      f"{int((np.abs(Po)>1e-8).sum()//2)} edges, "
      f"eigenvalues {', '.join(f'{v:.3f}' for v in ev)}")

cols=list(qiime2.Artifact.load(f"{out}/classo-x-trac.qza").view(pd.DataFrame).columns)
with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(f"{out}/motus-regress-age.qza") as z:
    n=[x for x in z.namelist() if x.endswith(".zip")][0]; z.extract(n,t)
    s=zarr.open(zarr.ZipStore(os.path.join(t,n),mode="r"))["solution"]
    b=np.asarray(s["CV"]["refit"]).ravel()
nz=np.argwhere(np.abs(b)>1e-10).ravel()
sel=pd.DataFrame({"clade":[cols[i] if i<len(cols) else f"coef{i}" for i in nz],
                  "coefficient":[round(float(b[i]),4) for i in nz]})
sel.to_csv(f"{tables}/motus-trac-selected.tsv",sep="\t",index=False)
print(f"  trac selects {len(nz)} of {len(b)} clades:")
print("   "+sel.to_string(index=False).replace("\n","\n   "))
PY

echo "############ artifact sizes (what would be committed) ############"
ls -la "$PUB" | awk '$1 ~ /^-/ {printf "  %-40s %8.1f KB\n", $9, $5/1024}'
ls -la "$OUT"/*.qzv | awk '{printf "  %-40s %8.1f KB\n", $9, $5/1024}'
echo "############ RESULT: complete ############"
