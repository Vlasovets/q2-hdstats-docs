#!/bin/bash
#SBATCH --job-name=q2-motus-smoke
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motussmoke_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motussmoke_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# One sample, end to end, and a sweep of the one parameter that could silently empty
# the table.
#
# WHY THIS EXISTS RATHER THAN PROFILING ALL 36 IMMEDIATELY.
#
# mOTUs calls a species only when at least `marker_gene_cutoff` (default 3) of its ten
# universal marker genes are covered, and an alignment only counts when it is at least
# `min_alen` bases long. q2-mOTUs defaults min_alen to 75, matching mOTUs itself, which
# is a sensible default for the 100-150 bp reads mOTUs was designed around.
#
# THESE READS ARE 73.2 TO 75.5 BASES LONG (base_count / read_count / 2 across the 36
# fecal runs). So the default threshold sits on top of, or above, the entire read. If
# that discards most alignments the profile is empty or near-empty -- and because
# mOTUs exits 0 and writes a well-formed file either way, an empty table looks exactly
# like a real result until someone counts the rows. Running 36 samples against an
# untested threshold risks 36 confidently wrong profiles.
#
# So: sweep min_alen on one real sample and read the curve before committing.
#
# Sample choice. ERR9528684 (13241.J13W2.S) is the smallest run in the set at 5.9 MiB
# and 73,108 reads -- it makes the sweep fast, but it is the shallowest sample and will
# under-report. ERR9528734 (13241.J24W2.S, 153 MiB, 2.0 M reads) is near the median and
# is what the decision should actually be read off. Both are run: the small one proves
# the plumbing, the median one answers the question.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
BWA_ENV="${Q2_BWA_ENV:-/home/itg/oleg.vlasovets/.conda/envs/bwa-0.7.19}"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/motus-smoke"
CACHE="${Q2_MAP_FASTQ:-/lustre/scratch/users/oleg.vlasovets/map13241/fastq}"

SMALL=ERR9528684
MEDIAN=ERR9528734
ALENS=(45 60 70 75)

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motussmoke/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motussmoke/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$CACHE"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
export PATH="$BWA_ENV/bin:$PATH"

python - <<'ENVCHECK' || exit 9
import sys, numpy, zarr
bad = []
if not zarr.__version__.startswith("2.18"): bad.append("zarr " + zarr.__version__)
if not numpy.__version__.startswith("2.4"): bad.append("numpy " + numpy.__version__)
if bad:
    sys.stderr.write("ENV CHECK FAILED: %s\n" % "; ".join(bad)); sys.exit(9)
ENVCHECK
command -v bwa >/dev/null || { echo "bwa not on PATH"; exit 9; }
[[ -r "$SHEET" ]] || { echo "missing sample sheet: $SHEET"; exit 2; }

cd "$SCRATCH"
FAIL=0

fetch_run () {   # $1 = run accession; downloads both mates into $CACHE, verifies md5
  local run="$1"
  python - "$SHEET" "$run" "$CACHE" <<'PY'
import csv, hashlib, os, sys, urllib.request
sheet, run, cache = sys.argv[1:4]
row = next(r for r in csv.DictReader(open(sheet), delimiter="\t")
           if r["run_accession"] == run)
for i in (1, 2):
    url, md5 = row[f"fastq_{i}_url"], row[f"fastq_{i}_md5"]
    dest = os.path.join(cache, os.path.basename(url))
    if os.path.isfile(dest):
        h = hashlib.md5(open(dest, "rb").read()).hexdigest()
        if h == md5:
            print(f"    [skip] {os.path.basename(dest)} present, md5 ok"); continue
        print(f"    md5 mismatch on {dest}, refetching")
    print(f"    fetching {os.path.basename(url)}")
    urllib.request.urlretrieve(url, dest)
    h = hashlib.md5(open(dest, "rb").read()).hexdigest()
    if h != md5:
        sys.exit(f"    FAIL md5 {h} != {md5} for {dest}")
    print(f"    ok md5 {h}")
PY
}

profile_at () {  # $1 = run, $2 = min_alen  -> writes $OUT/<run>-alen<N>-{table,taxonomy}.qza
  local run="$1" alen="$2"
  local tag="${run}-alen${alen}"
  local tbl="$OUT/${tag}-table.qza"
  if [[ -f "$tbl" ]]; then echo "    [skip] $tbl exists"; return 0; fi

  local d="$SCRATCH/$tag"; mkdir -p "$d/reads"
  ln -sf "$CACHE/${run}_1.fastq.gz" "$d/reads/${run}_00_L001_R1_001.fastq.gz"
  ln -sf "$CACHE/${run}_2.fastq.gz" "$d/reads/${run}_00_L001_R2_001.fastq.gz"

  qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --input-path "$d/reads" --output-path "$d/reads.qza" >/dev/null 2>&1 \
    || { echo "    FAIL import for $tag"; return 1; }

  # Timing with SECONDS, not /usr/bin/time -- that binary is not present on these
  # compute nodes, and wrapping the real command in it made every profile call fail.
  local t0=$SECONDS
  qiime motus profile --i-samples "$d/reads.qza" \
    --p-threads "${SLURM_CPUS_PER_TASK:-8}" --p-min-alen "$alen" \
    --o-table "$tbl" --o-taxonomy "$OUT/${tag}-taxonomy.qza" 2>&1 \
    | tail -6 | sed 's/^/      /'
  echo "      wall $((SECONDS - t0))s"
  [[ -f "$tbl" ]] || { echo "    FAIL no table for $tag"; return 1; }
}

echo "############ [1/4] fetch the two smoke samples ############"
for r in "$SMALL" "$MEDIAN"; do echo "  $r"; fetch_run "$r" || FAIL=1; done

echo "############ [2/4] sweep min_alen on the SMALL run ($SMALL) ############"
for a in "${ALENS[@]}"; do echo "  min_alen=$a"; profile_at "$SMALL" "$a" || FAIL=1; done

echo "############ [3/4] sweep min_alen on the MEDIAN run ($MEDIAN) ############"
for a in "${ALENS[@]}"; do echo "  min_alen=$a"; profile_at "$MEDIAN" "$a" || FAIL=1; done

echo "############ [4/4] what did each threshold actually yield? ############"
python - "$OUT" "$ROOT/results/tables/motus-min-alen-sweep.tsv" <<'PY'
import os, sys, glob, pandas as pd, qiime2
out, tsv = sys.argv[1], sys.argv[2]
os.makedirs(os.path.dirname(tsv), exist_ok=True)
rows = []
for p in sorted(glob.glob(os.path.join(out, "*-table.qza"))):
    tag = os.path.basename(p)[:-len("-table.qza")]
    run, alen = tag.rsplit("-alen", 1)
    df = qiime2.Artifact.load(p).view(pd.DataFrame)
    s = df.iloc[0] if len(df) else pd.Series(dtype=float)
    nz = s[s > 0]
    unassigned = float(s.get("unassigned", 0))
    total = float(s.sum())
    assigned = total - unassigned
    rows.append(dict(run=run, min_alen=int(alen), features_total=df.shape[1],
                     features_nonzero=int((nz.index != "unassigned").sum()),
                     assigned_counts=assigned, unassigned_counts=unassigned,
                     pct_assigned=round(100 * assigned / total, 3) if total else 0.0))
d = pd.DataFrame(rows).sort_values(["run", "min_alen"])
d.to_csv(tsv, sep="\t", index=False)
print(d.to_string(index=False).replace("\n", "\n    ").rjust(4))
print(f"\n    -> {tsv}")
if not len(d) or d["features_nonzero"].max() == 0:
    print("    FAIL every threshold produced an empty profile"); sys.exit(1)
PY
[[ $? -eq 0 ]] || FAIL=1

if [[ $FAIL -eq 0 ]]; then echo "############ RESULT: PASS ############"
else echo "############ RESULT: FAIL ############"; fi
exit $FAIL
