#!/bin/bash
#SBATCH --job-name=q2-motus-prof
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusprof_%A_%a.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusprof_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --array=1-36%8
#
# Profile all 36 fecal WGS runs. One array task per run.
#
# THIS STAGE EXISTS TO SETTLE GATE M3, NOT TO PRODUCE A FINISHED ANALYSIS.
#
# The smoke test (stage 34) profiled two samples and found the yield alarmingly low:
#
#     ERR9528684    73 k reads  ->    11 assigned counts,  5 nonzero features
#     ERR9528734   2.0 M reads  ->   705 assigned counts, 12 nonzero features
#
# The taxa are biologically right for a preterm NICU gut -- Staphylococcus, Klebsiella,
# Escherichia, Enterococcus, Corynebacterium -- so nothing is broken. mOTUs simply
# counts only reads landing on ten universal marker genes, and infant stool is
# host-dominated, so a 2 M-read run yields a few hundred marker inserts.
#
# The question Gate M3 asks is whether the UNION across all 36 samples, at a usable
# per-sample depth, can support a covariance estimate. Two samples cannot answer that;
# 36 can, and at ~2 minutes each it is cheap to find out properly rather than
# extrapolate. If the answer is no, these numbers are the evidence for saying so -- and
# they are worth reporting either way, because "what mOTUs yields on shallow
# host-dominated metagenomes" is a real and publishable caveat.
#
# min_alen is left at the default 75. Stage 34 swept 45/60/70/75 and every value gave
# IDENTICAL counts, so the concern that a 75-base threshold would reject 73-75 base
# reads was unfounded -- the threshold applies to the alignment after mapping.
#
# Idempotent: an existing table for a run is skipped.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
BWA_ENV="${Q2_BWA_ENV:-/home/itg/oleg.vlasovets/.conda/envs/bwa-0.7.19}"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/motus-profiles"
CACHE="${Q2_MAP_FASTQ:-/lustre/scratch/users/oleg.vlasovets/map13241/fastq}"
MIN_ALEN="${Q2_MOTUS_MIN_ALEN:-75}"
CUTOFF="${Q2_MOTUS_MARKER_CUTOFF:-3}"

IDX="${SLURM_ARRAY_TASK_ID:-1}"
RUN=$(awk -F'\t' -v n="$IDX" 'NR==n+1 {print $11}' "$SHEET")
[[ -n "$RUN" ]] || { echo "no run at row $IDX of $SHEET"; exit 2; }

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motusprof/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motusprof/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$CACHE"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
export PATH="$BWA_ENV/bin:$PATH"
command -v bwa >/dev/null || { echo "bwa not on PATH"; exit 9; }

TBL="$OUT/${RUN}-table.qza"
if [[ -f "$TBL" ]]; then echo "[skip] $TBL exists"; exit 0; fi
echo "task $IDX -> $RUN"

python - "$SHEET" "$RUN" "$CACHE" <<'PY' || exit 3
import csv, hashlib, os, sys, urllib.request
sheet, run, cache = sys.argv[1:4]
row = next(r for r in csv.DictReader(open(sheet), delimiter="\t")
           if r["run_accession"] == run)
for i in (1, 2):
    url, md5 = row[f"fastq_{i}_url"], row[f"fastq_{i}_md5"]
    dest = os.path.join(cache, os.path.basename(url))
    if os.path.isfile(dest):
        if hashlib.md5(open(dest, "rb").read()).hexdigest() == md5:
            print(f"  [skip] {os.path.basename(dest)} md5 ok"); continue
        print(f"  md5 mismatch, refetching {dest}")
    print(f"  fetching {os.path.basename(url)}")
    urllib.request.urlretrieve(url, dest)
    h = hashlib.md5(open(dest, "rb").read()).hexdigest()
    if h != md5: sys.exit(f"  FAIL md5 {h} != {md5}")
    print("  md5 ok")
PY

D="$SCRATCH/$RUN"; mkdir -p "$D/reads"
ln -sf "$CACHE/${RUN}_1.fastq.gz" "$D/reads/${RUN}_00_L001_R1_001.fastq.gz"
ln -sf "$CACHE/${RUN}_2.fastq.gz" "$D/reads/${RUN}_00_L001_R2_001.fastq.gz"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --input-path "$D/reads" --output-path "$D/reads.qza" >/dev/null 2>&1 \
  || { echo "FAIL import"; exit 4; }

t0=$SECONDS
qiime motus profile --i-samples "$D/reads.qza" \
  --p-threads "${SLURM_CPUS_PER_TASK:-8}" \
  --p-min-alen "$MIN_ALEN" --p-marker-gene-cutoff "$CUTOFF" \
  --o-table "$TBL" --o-taxonomy "$OUT/${RUN}-taxonomy.qza" 2>&1 | tail -4
echo "wall $((SECONDS - t0))s"
[[ -f "$TBL" ]] || { echo "FAIL no table for $RUN"; exit 5; }
echo "done $RUN"
