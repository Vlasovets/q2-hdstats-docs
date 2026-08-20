#!/bin/bash
#SBATCH --job-name=q2-motus-models
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusmodels_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusmodels_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Does the mOTUs table actually support the two models the tier would be built on?
#
# Gate M3 cleared its stated thresholds (p=94, n=36 after a prevalence-3 filter) but the
# table is 84.6% zeros, 182 cells are exactly 1, and 21 of the 94 features never exceed
# a count of 2. Whether that supports a covariance estimate is an empirical question,
# and arguing about it is not evidence. So: run the chain.
#
# Three prevalence thresholds, because the filter is the one free choice that most
# changes p, and picking it by eye would be exactly the kind of unexamined decision the
# rest of this pipeline exists to avoid:
#     >= 3 samples   p ~ 94   permissive
#     >= 5 samples   p ~ 39
#     >= 10 samples  p ~ 12   conservative
#
# For each: mclr transform -> scaled covariance -> eBIC lambda path -> edge count.
# mclr rather than clr because the table is sparse and mclr is the transform the book
# already argues for in that regime.
#
# Then a log-contrast regression of host_age_days, the outcome chosen because it is the
# only variable that varies WITHIN an infant (the design is 18 subjects x 2 timepoints,
# 13-29 days apart), so the repeated measures are informative rather than
# pseudo-replication.
#
# WHAT WOULD COUNT AS FAILURE. An empty or near-complete graph at every lambda; an eBIC
# curve with no interior minimum; or a regression that selects nothing at any lambda.
# Any of those means the tier cannot carry a network chapter regardless of what p says.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/motus-models"
TABLES="$ROOT/results/tables"
MIN_DEPTH="${Q2_MOTUS_MIN_DEPTH:-100}"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motusmodels/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motusmodels/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
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
[[ -r "$MERGED" ]] || { echo "missing $MERGED"; exit 2; }

echo "############ [1/4] build filtered artifacts + outcome file ############"
python - "$MERGED" "$SHEET" "$SCRATCH" "$MIN_DEPTH" <<'PY' || exit 3
import sys, csv, pandas as pd, numpy as np, qiime2, biom
merged, sheet, out, min_depth = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
t = pd.read_csv(merged, sep="\t", index_col=0)          # samples x features
meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}

depth = t.sum(axis=1)
keep_s = depth[depth >= min_depth].index
print(f"  dropping {len(t) - len(keep_s)} samples below {min_depth:.0f} assigned counts")
t = t.loc[keep_s]

# Outcomes, indexed by the run accession the table uses.
rows = []
for r in t.index:
    m = meta.get(r, {})
    rows.append({"sample-id": r, "host_age_days": m.get("host_age_days", ""),
                 "diagnosis": m.get("diagnosis", ""),
                 "host_subject_id": m.get("host_subject_id", "")})
md = pd.DataFrame(rows).set_index("sample-id")
md.to_csv(f"{out}/outcomes.tsv", sep="\t")
print(f"  outcomes.tsv: {md.shape[0]} samples")

# gglasso transform-features REQUIRES --i-taxonomy even though its body never reads it
# (documented in the book as a known wart). So one merged taxonomy artifact covering the
# union of features has to exist before the chain can run at all.
import glob, csv as _csv, os
tax = {}
for tp in sorted(glob.glob(os.path.join(os.path.dirname(merged), "..",
                                        "motus-profiles", "*-taxonomy.qza"))):
    import zipfile, io
    with zipfile.ZipFile(tp) as z:
        n = [x for x in z.namelist() if x.endswith("/data/taxonomy.tsv")][0]
        for r in _csv.reader(io.StringIO(z.read(n).decode()), delimiter="\t"):
            if r and r[0] != "Feature ID":
                tax[r[0]] = r[1]
taxdf = pd.DataFrame({"Taxon": [tax.get(f, "Unassigned") for f in t.columns]},
                     index=pd.Index(t.columns, name="Feature ID"))
qiime2.Artifact.import_data("FeatureData[Taxonomy]", taxdf).save(f"{out}/taxonomy.qza")
print(f"  taxonomy.qza: {len(taxdf)} features, "
      f"{sum(1 for v in taxdf.Taxon if v == 'Unassigned')} without a lineage")

for thr in (3, 5, 10):
    prev = (t > 0).sum(axis=0)
    sub = t.loc[:, prev >= thr]
    # A feature can survive the prevalence filter yet be all-zero in the retained
    # samples once shallow ones are dropped; drop those too.
    sub = sub.loc[:, sub.sum(axis=0) > 0]
    tab = biom.Table(sub.T.values, list(sub.columns), list(sub.index))
    qiime2.Artifact.import_data("FeatureTable[Frequency]", tab).save(f"{out}/table-prev{thr}.qza")
    z = 100 * float((sub.to_numpy() == 0).mean())
    print(f"  prev>={thr:2d}: p={sub.shape[1]:3d} n={sub.shape[0]}  zeros={z:.1f}%")
PY

echo "############ [2/4] network chain at each threshold ############"
for thr in 3 5 10; do
  echo "  --- prevalence >= $thr ---"
  T="$SCRATCH/table-prev${thr}.qza"
  N=$(python -c "
import qiime2, pandas as pd
print(qiime2.Artifact.load('$T').view(pd.DataFrame).shape[0])")
  qiime gglasso transform-features --i-table "$T" \
    --i-taxonomy "$SCRATCH/taxonomy.qza" \
    --p-transformation mclr --p-keep-original-id \
    --o-transformed-table "$SCRATCH/mclr-${thr}.qza" >"$SCRATCH/tf-${thr}.log" 2>&1 \
    || { echo "    FAIL transform:"; tail -6 "$SCRATCH/tf-${thr}.log" | sed 's/^/      /'; continue; }
  qiime gglasso calculate-covariance --i-table "$SCRATCH/mclr-${thr}.qza" \
    --p-method scaled --o-covariance-matrix "$SCRATCH/corr-${thr}.qza" >/dev/null 2>&1 \
    || { echo "    FAIL covariance"; continue; }
  qiime gglasso solve-problem --i-covariance-matrix "$SCRATCH/corr-${thr}.qza" \
    --p-n-samples "$N" --p-no-latent --p-path-scale linear \
    --p-lambda1-min 0.05 --p-lambda1-max 1.0 --p-n-lambda1 20 --p-gamma 0.3 \
    --o-solution "$OUT/sgl-prev${thr}.qza" 2>&1 | tail -2 | sed 's/^/    /'
done

echo "############ [3/4] log-contrast regression of host_age_days ############"
qiime classo transform-features --i-features "$SCRATCH/table-prev3.qza" \
  --o-x "$SCRATCH/classo-x.qza" >"$SCRATCH/cx.log" 2>&1 \
  || { echo "  FAIL classo transform:"; tail -6 "$SCRATCH/cx.log" | sed 's/^/    /'; }
if [[ -f "$SCRATCH/classo-x.qza" ]]; then
  qiime classo regress --i-features "$SCRATCH/classo-x.qza" \
    --m-y-file "$SCRATCH/outcomes.tsv" --m-y-column host_age_days \
    --p-concomitant --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.001 \
    --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-cv-one-se \
    --p-cv-nlam 60 --p-cv-lamin 0.001 --p-cv-logscale \
    --p-no-stabsel --p-no-lamfixed \
    --o-result "$OUT/regress-age.qza" 2>&1 | tail -3 | sed 's/^/  /'
fi

echo "############ [4/4] read the models out ############"
python - "$OUT" "$TABLES/motus-model-viability.tsv" <<'PY'
import os, sys, glob, zipfile, tempfile, numpy as np, pandas as pd, zarr
out, tsv = sys.argv[1], sys.argv[2]

def problem(qza):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(qza) as z:
        i = [n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i, t)
        r = zarr.open(zarr.ZipStore(os.path.join(t, i), mode="r"))
        return {k: r[k] for k in r.keys()}, r

rows = []
for p in sorted(glob.glob(os.path.join(out, "sgl-prev*.qza"))):
    thr = int(os.path.basename(p).split("prev")[1].split(".")[0])
    try:
        _, r = problem(p)
        sol = r["solution"]
        P = np.asarray(sol["precision_"])
        off = ~np.eye(P.shape[0], dtype=bool)
        edges = int((np.abs(P[off]) > 1e-8).sum() // 2)
        pairs = P.shape[0] * (P.shape[0] - 1) // 2
        ms = r.get("modelselect_stats")
        lam = sel = ebic_min = None
        if ms is not None and "LAMBDA" in ms:
            L = np.asarray(ms["LAMBDA"]).ravel()
            gk = sorted(ms["BIC"].keys(), key=float) if "BIC" in ms else []
            if gk:
                g = "0.3" if "0.3" in gk else gk[len(gk)//2]
                c = np.asarray(ms["BIC"][g]).ravel()
                o = np.argsort(L); L, c = L[o], c[o]
                j = int(np.argmin(c)); sel = float(L[j]); ebic_min = float(c[j])
                lam = f"{L.min():.3g}-{L.max():.3g}"
                # "argmin is not at an endpoint" is NOT sufficient to call a minimum
                # interior. Every empty graph scores the SAME eBIC, so a run of empty
                # solutions is a flat plateau and argmin returns its first index --
                # which looked like an interior optimum at lambda=0.65 when in fact
                # eBIC was 408.0 for every lambda from 0.65 to 1.0. Require the
                # neighbours to be strictly worse.
                interior = (0 < j < len(L) - 1
                            and c[j - 1] > c[j] + 1e-9 and c[j + 1] > c[j] + 1e-9)
        rows.append(dict(prevalence=thr, p=P.shape[0], edges=edges, pairs=pairs,
                         density=round(100*edges/pairs, 2) if pairs else 0,
                         selected_lambda=sel, min_ebic=ebic_min,
                         lambda_range=lam, interior_min=bool(interior)))
    except Exception as e:
        rows.append(dict(prevalence=thr, p=None, edges=None, pairs=None, density=None,
                         selected_lambda=None, min_ebic=None, lambda_range=None,
                         interior_min=None))
        print(f"  prev{thr}: could not read ({type(e).__name__}: {e})")

if rows:
    d = pd.DataFrame(rows)
    d.to_csv(tsv, sep="\t", index=False)
    print(d.to_string(index=False))
    print(f"\n  -> {tsv}")

reg = os.path.join(out, "regress-age.qza")
if os.path.isfile(reg):
    with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(reg) as z:
        names = [n for n in z.namelist() if n.endswith(".zip")]
        if names:
            z.extract(names[0], t)
            r = zarr.open(zarr.ZipStore(os.path.join(t, names[0]), mode="r"))
            print("\n  regression keys:", sorted(r.keys()))
            if "solution" in r:
                s = r["solution"]
                print("  solution/:", sorted(s.keys()))
                if "CV" in s and "refit" in s["CV"]:
                    b = np.asarray(s["CV"]["refit"]).ravel()
                    print(f"  CV refit: {len(b)} coefficients, "
                          f"{int((np.abs(b) > 1e-10).sum())} nonzero")
else:
    print("\n  no regression artifact produced")
PY
echo "############ RESULT: models run, read the numbers above ############"
