#!/bin/bash
#SBATCH --job-name=q2-classo-diag
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classodiag_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/classodiag_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Why does q2-classo select one coefficient, and is that a setting or the data?
#
# The abundance sweep found a network that works: mOTUs top-100 by abundance, p=100,
# n=34, a genuine interior eBIC minimum at gamma=0.15 giving 481 edges against a
# permutation null of 0.2. But the log-contrast regression on the SAME table returned
# "1 of 101 coefficients" for host_age_days, and a tier has to illustrate both plugins.
#
# "1 of 101" has two very different explanations and they demand opposite responses:
#
#   A SETTING. --p-cv-one-se applies the one-standard-error rule, which deliberately
#   returns the SPARSEST model within one standard error of the cross-validated minimum.
#   When the CV curve is shallow that rule lands on the intercept-only model by design.
#   Every earlier run used it, and no run without it has ever been made on this data.
#
#   THE DATA. If cross-validated error never improves on the intercept-only model at any
#   lambda, then nothing in these 100 features predicts postnatal age and no setting will
#   conjure it.
#
# These are distinguished by looking at the CV CURVE, not by trying settings until one
# gives a nonzero answer -- that would be selecting a hyperparameter on the outcome we
# want, which is how tutorials end up teaching overfitting.
#
# So this stage reports, for each configuration: the CV error against lambda, whether it
# has an interior minimum, and how the selected support grows along the path. Then it
# runs the variants.
#
# ALSO TESTED: diagnosis (classification). Its expected result is a null -- the source
# publication reports no microbiome effect of maternal asthma, and 9 v 9 infants cannot
# detect a small one. It is included so the comparison is on the record, not because a
# positive result is anticipated.
#
# CV LEAKAGE IS REAL HERE AND IS NOT FIXED BY THIS STAGE. Each infant contributes two
# samples, so a random 5-fold split puts the same infant in train and test. q2-classo
# exposes no grouping parameter. Every number below is therefore optimistic; the stage
# also fits a subject-disjoint split so the two can be compared.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/classo-diag"
TABLES="$ROOT/results/tables"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-classodiag/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-classodiag/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ [1/4] rebuild the top-100 table and the outcomes ############"
python - "$MERGED" "$SHEET" "$SCRATCH" <<'PY' || exit 3
import sys, csv, numpy as np, pandas as pd, qiime2, biom
merged, sheet, out = sys.argv[1:4]
t = pd.read_csv(merged, sep="\t", index_col=0)
t = t.loc[t.sum(1) >= 100]
meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
base = t.loc[:, (t > 0).sum(0) >= 2]
sub = base[list(base.sum(0).sort_values(ascending=False).head(100).index)]
sub = sub.loc[:, sub.sum(0) > 0]
print(f"  top-100 by abundance: p={sub.shape[1]} n={sub.shape[0]}")
biom_t = biom.Table(sub.T.values, list(sub.columns), list(sub.index))
qiime2.Artifact.import_data("FeatureTable[Frequency]", biom_t).save(f"{out}/top100.qza")

rows=[]
for s in sub.index:
    m=meta[s]
    rows.append({"sample-id": s, "host_age_days": m["host_age_days"],
                 "diagnosis": m["diagnosis"], "host_subject_id": m["host_subject_id"]})
md=pd.DataFrame(rows).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")
print(f"  outcomes: {md.shape[0]} samples, "
      f"{md.host_subject_id.nunique()} subjects, "
      f"age {md.host_age_days.astype(float).min():.0f}-{md.host_age_days.astype(float).max():.0f} d")

# A subject-disjoint split, so the leaky estimate can be compared against an honest one.
rng=np.random.default_rng(7)
subs=sorted(md.host_subject_id.unique()); rng.shuffle(subs)
test=set(subs[:max(3,len(subs)//4)])
md["split"]=["test" if s in test else "train" for s in md.host_subject_id]
md.to_csv(f"{out}/outcomes-split.tsv", sep="\t")
print(f"  subject-disjoint split: {(md.split=='train').sum()} train / "
      f"{(md.split=='test').sum()} test samples, "
      f"{len(subs)-len(test)}/{len(test)} subjects")
PY

qiime classo transform-features --i-features "$SCRATCH/top100.qza" \
  --o-x "$SCRATCH/x.qza" >/dev/null 2>&1 || { echo "  FAIL classo transform"; exit 4; }

echo "############ [2/4] regression variants ############"
run_reg () {   # $1 tag, $2 outcome, $3... extra flags
  local tag="$1" col="$2"; shift 2
  [[ -f "$OUT/reg-${tag}.qza" ]] && { echo "  [skip] $tag"; return 0; }
  qiime classo regress --i-features "$SCRATCH/x.qza" \
    --m-y-file "$SCRATCH/outcomes.tsv" --m-y-column "$col" \
    --p-concomitant --p-path --p-path-nlam-log 80 --p-path-lamin-log 0.0005 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 \
    --p-cv-nlam 80 --p-cv-lamin 0.0005 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed "$@" \
    --o-result "$OUT/reg-${tag}.qza" >"$SCRATCH/${tag}.log" 2>&1 \
    && echo "  $tag: ok" \
    || { echo "  $tag: FAIL"; tail -4 "$SCRATCH/${tag}.log" | sed 's/^/      /'; }
}
run_reg age-onese   host_age_days --p-cv-one-se
run_reg age-cvmin   host_age_days --p-no-cv-one-se
run_reg age-stabsel host_age_days --p-no-cv-one-se --p-stabsel --p-stabsel-threshold 0.7

echo "############ [3/4] classification on diagnosis ############"
if [[ ! -f "$OUT/cls-diagnosis.qza" ]]; then
  qiime classo classify --i-features "$SCRATCH/x.qza" \
    --m-y-file "$SCRATCH/outcomes.tsv" --m-y-column diagnosis \
    --p-path --p-path-nlam-log 80 --p-path-lamin-log 0.0005 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
    --p-cv-nlam 80 --p-cv-lamin 0.0005 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed \
    --o-result "$OUT/cls-diagnosis.qza" >"$SCRATCH/cls.log" 2>&1 \
    && echo "  diagnosis: ok" \
    || { echo "  diagnosis: FAIL"; tail -5 "$SCRATCH/cls.log" | sed 's/^/      /'; }
fi

echo "############ [4/4] read the CV curves ############"
python - "$OUT" "$TABLES/classo-diagnosis-summary.tsv" <<'PY'
import os,sys,glob,zipfile,tempfile,numpy as np,pandas as pd,zarr
out,tsv=sys.argv[1],sys.argv[2]
rows=[]
for p in sorted(glob.glob(os.path.join(out,"*.qza"))):
    tag=os.path.basename(p)[:-4]
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(p) as z:
        nm=[x for x in z.namelist() if x.endswith(".zip")]
        if not nm: continue
        z.extract(nm[0],t)
        r=zarr.open(zarr.ZipStore(os.path.join(t,nm[0]),mode="r"))
        s=r.get("solution")
        if s is None: continue
        print(f"\n  === {tag} ===  solution keys: {sorted(s.keys())}")
        if "CV" in s:
            cv=s["CV"]; print(f"      CV keys: {sorted(cv.keys())}")
            err=None
            for k in ("xGraph","yGraph","standard_error","MSE","mse"):
                if k in cv:
                    v=np.asarray(cv[k]).ravel()
                    print(f"      {k}: len {len(v)}  min {v.min():.4g}  max {v.max():.4g}")
                    if k in ("yGraph","MSE","mse"): err=v
            if err is not None and len(err)>2:
                j=int(np.argmin(err))
                interior=0<j<len(err)-1
                # A CV curve with no interior minimum means nothing beat the null model.
                print(f"      CV error min at index {j} of {len(err)-1} "
                      f"({'INTERIOR' if interior else 'BOUNDARY'}); "
                      f"err[0]={err[0]:.4g} min={err[j]:.4g} err[-1]={err[-1]:.4g}")
                print(f"      improvement over the most-penalised end: "
                      f"{100*(err[0]-err[j])/err[0]:.2f}%" if err[0] else "")
                rows.append(dict(config=tag,cv_len=len(err),cv_min_index=j,
                                 cv_interior=interior,err_first=float(err[0]),
                                 err_min=float(err[j]),err_last=float(err[-1])))
            if "refit" in cv:
                b=np.asarray(cv["refit"]).ravel()
                print(f"      CV refit: {int((np.abs(b)>1e-10).sum())} of {len(b)} nonzero")
        if "PATH" in s and "BETAS" in s["PATH"]:
            B=np.asarray(s["PATH"]["BETAS"])
            nz=(np.abs(B)>1e-10).sum(1)
            print(f"      PATH support size along lambda: min {nz.min()} max {nz.max()}"
                  f"  (first 8: {list(nz[:8])})")
if rows:
    pd.DataFrame(rows).to_csv(tsv,sep="\t",index=False); print(f"\n  -> {tsv}")
PY
echo "############ RESULT: complete ############"
