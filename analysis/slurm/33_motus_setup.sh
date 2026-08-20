#!/bin/bash
#SBATCH --job-name=q2-motus-setup
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motussetup_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motussetup_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Everything `qiime motus profile` needs that is not already in the environment:
# the mOTUs profiler, the bwa aligner, and the 2.9 GiB reference database.
#
# THREE DELIBERATE CHOICES, each avoiding a way this could go wrong.
#
# 1. bwa goes in its OWN conda environment, not the QIIME 2 one.
#    A second `conda install` into a solved distribution env is precisely what lets the
#    solver drift a pinned package -- the risk the migration called out and worked hard
#    to avoid. bwa is a standalone C binary with no Python involvement, so it has no
#    business sharing an environment with numpy 2.4.2. The profiling stage puts it on
#    PATH instead.
#    Also NOT using the upstream Makefile's bwa step: it git-clones and compiles into
#    the current directory, then symlinks $(PWD)/bwa/bwa into the conda prefix, leaving
#    a build tree and a link that breaks the moment anyone tidies up.
#
# 2. motu-profiler is pip-installed with --no-deps into the QIIME 2 env.
#    Checked first: motu-profiler 3.1.0 declares zero Python dependencies, so this
#    cannot perturb the solve. --no-deps anyway, because pip cannot see conda's pins.
#
# 3. The database lives on shared storage and is SYMLINKED into site-packages.
#    motus resolves its database as dirname(motus/motus.py) + "/db_mOTU" unless given a
#    global `-db` argument -- and q2-mOTUs never passes one, the command it builds is
#    fixed. So the path is not negotiable, but what sits at it can be a link. Keeping
#    ~13 GiB of reference data out of the conda environment means the env stays
#    rebuildable from its lockfile, and the database survives an env rebuild.
#
# Idempotent throughout: re-running skips anything already in place.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
BWA_ENV="${Q2_BWA_ENV:-/home/itg/oleg.vlasovets/.conda/envs/bwa-0.7.19}"
DB_ROOT="${Q2_MOTUS_DB_ROOT:-/lustre/groups/itg/shared/oleg.vlasovets/motus-db}"

DB_URL="https://zenodo.org/records/7778108/files/db_mOTU_v3.1.0.tar.gz"
DB_TGZ="$DB_ROOT/db_mOTU_v3.1.0.tar.gz"
DB_MD5="f841c36150025af837f7a9a358c9a3c3"     # hardcoded in motus/downloadDB.py

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motussetup/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motussetup/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$DB_ROOT"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

python - <<'ENVCHECK' || exit 9
import sys, numpy, zarr
bad = []
if not zarr.__version__.startswith("2.18"):
    bad.append("zarr %s, expected 2.18.x" % zarr.__version__)
if not numpy.__version__.startswith("2.4"):
    bad.append("numpy %s, expected 2.4.x" % numpy.__version__)
if bad:
    sys.stderr.write("ENV CHECK FAILED: %s\n" % "; ".join(bad)); sys.exit(9)
ENVCHECK

FAIL=0

echo "############ [1/5] bwa, in its own environment ############"
if [[ -x "$BWA_ENV/bin/bwa" ]]; then
  echo "  [skip] $BWA_ENV/bin/bwa already present"
else
  conda create -y -q -p "$BWA_ENV" -c bioconda -c conda-forge bwa=0.7.19 2>&1 | tail -4 | sed 's/^/    /'
fi
if [[ -x "$BWA_ENV/bin/bwa" ]]; then
  # bwa prints its banner to stderr and exits 1 with no args -- that is success here.
  echo "  bwa: $("$BWA_ENV/bin/bwa" 2>&1 | grep -m1 -i version || echo '(no version line)')"
else
  echo "  FAIL bwa not installed"; FAIL=1
fi

echo "############ [2/5] motu-profiler into the QIIME 2 env ############"
if python -c "import motus" 2>/dev/null; then
  echo "  [skip] motus python package already importable"
else
  pip install --no-deps -q "motu-profiler==3.1.0" 2>&1 | tail -3 | sed 's/^/    /'
fi
MOTUS_PKG=$(python -c "import motus, os; print(os.path.dirname(motus.__file__))" 2>/dev/null)
if [[ -z "$MOTUS_PKG" ]]; then echo "  FAIL motus package not importable"; FAIL=1
else echo "  motus package at: $MOTUS_PKG"; fi

echo "############ [3/5] reference database (2.9 GiB) ############"
if [[ -d "$DB_ROOT/db_mOTU" ]]; then
  echo "  [skip] $DB_ROOT/db_mOTU already extracted"
else
  if [[ -f "$DB_TGZ" ]] && [[ "$(md5sum "$DB_TGZ" | cut -d' ' -f1)" == "$DB_MD5" ]]; then
    echo "  [skip] tarball present and md5 matches"
  else
    echo "  downloading $DB_URL"
    curl -L --fail --retry 5 --retry-delay 10 -C - -o "$DB_TGZ" "$DB_URL" \
      || { echo "  FAIL download"; FAIL=1; }
    got=$(md5sum "$DB_TGZ" | cut -d' ' -f1)
    if [[ "$got" != "$DB_MD5" ]]; then
      echo "  FAIL md5 mismatch: got $got expected $DB_MD5"; FAIL=1
    else
      echo "  md5 ok: $got"
    fi
  fi
  if [[ $FAIL -eq 0 ]]; then
    echo "  extracting into $DB_ROOT"
    tar -xzf "$DB_TGZ" -C "$DB_ROOT" || { echo "  FAIL extract"; FAIL=1; }
  fi
fi
[[ -d "$DB_ROOT/db_mOTU" ]] && echo "  db files: $(ls -1 "$DB_ROOT/db_mOTU" | wc -l)"

echo "############ [4/5] link the database where motus will look for it ############"
# Not negotiable: motus computes dirname(motus.py) + "/db_mOTU" and q2-mOTUs never
# passes -db. A symlink satisfies it without putting the data in the env.
if [[ -n "$MOTUS_PKG" ]]; then
  LINK="$MOTUS_PKG/db_mOTU"
  if [[ -L "$LINK" || -d "$LINK" ]]; then
    echo "  [skip] $LINK exists -> $(readlink -f "$LINK" 2>/dev/null || echo 'real dir')"
  else
    ln -s "$DB_ROOT/db_mOTU" "$LINK" && echo "  linked $LINK -> $DB_ROOT/db_mOTU"
  fi
fi

echo "############ [5/5] does motus actually start and find its database? ############"
export PATH="$BWA_ENV/bin:$PATH"
if motus -h >"$SCRATCH/motus-h.txt" 2>&1; then
  echo "  motus -h: ok"
  grep -m1 -i 'version' "$SCRATCH/motus-h.txt" | sed 's/^/    /'
else
  echo "  motus -h exited non-zero:"; head -12 "$SCRATCH/motus-h.txt" | sed 's/^/    /'
fi
# The real question is whether the profiler resolves its DB, which -h does not test.
if motus profile -h >"$SCRATCH/motus-profile-h.txt" 2>&1; then
  echo "  motus profile -h: ok"
else
  echo "  motus profile -h exited non-zero (may be normal for this tool):"
  head -8 "$SCRATCH/motus-profile-h.txt" | sed 's/^/    /'
fi
python - "$DB_ROOT" <<'PY'
import os, sys, motus
pkg = os.path.dirname(motus.__file__)
db = os.path.join(pkg, "db_mOTU")
print(f"  resolved db path: {db}")
print(f"  exists: {os.path.exists(db)}   is symlink: {os.path.islink(db)}")
if os.path.isdir(db):
    names = sorted(os.listdir(db))
    print(f"  entries: {len(names)}; e.g. {names[:4]}")
    ver = os.path.join(db, "db_mOTU_versions")
    if os.path.isfile(ver):
        print("  db version file, first lines:")
        for line in open(ver).read().splitlines()[:3]:
            print(f"    {line}")
else:
    sys.exit(1)
PY
[[ $? -eq 0 ]] || FAIL=1

if [[ $FAIL -eq 0 ]]; then
  echo "############ RESULT: PASS -- profiler, aligner and database are in place ############"
  echo "  Profiling stages must put bwa on PATH:  export PATH=\"$BWA_ENV/bin:\$PATH\""
else
  echo "############ RESULT: FAIL ############"
fi
exit $FAIL
