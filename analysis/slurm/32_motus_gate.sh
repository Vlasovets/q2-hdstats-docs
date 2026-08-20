#!/bin/bash
#SBATCH --job-name=q2-motus-gate
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusgate_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/motusgate_%j.err
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Gate M1 + M2: can q2-mOTUs be used from this environment at all?
#
# WHY A COMPUTE NODE. The login node cannot build a QIIME 2 PluginManager: q2_composition
# imports rpy2, which raises
#     "9 arguments passed to .Internal(list.files) which requires 8"
# -- an R/rpy2 version mismatch, NOT a missing package. Every action that touches
# Artifact.import_data trips it. Compute nodes are unaffected; stages 30 and 31 both ran
# real qiime commands there. So this gate must be submitted, never run interactively.
#
# WHAT IS ALREADY KNOWN (checked on the login node, which can still do plain imports):
#   - every module q2-mOTUs imports resolves on 2026.7; a `qiime2` shim survives
#     alongside `rachis`
#   - the plugin object BUILDS unmodified, registering `profile` and `import_table`
# So this gate is not asking "does it import" -- it is asking whether the registered
# actions actually run and produce correctly shaped artifacts.
#
# M1  `qiime motus --help` lists both actions after a cache refresh.
# M2  `import-table` on the checked-in fixture yields FeatureTable[Frequency] +
#     FeatureData[Taxonomy], with:
#       - the table oriented samples-as-index (see the orientation note below)
#       - feature IDs matching (ref|meta|ext)_mOTU_v3_\d+
#       - the literal `unassigned` feature present (it must be dropped downstream)
#       - a 7-rank lineage suitable for q2-classo add-taxa
#
# ORIENTATION IS TESTED, NOT ASSUMED. _utils.extract_table_tax ends with `.T`, and a
# trailing transpose is EXACTLY the shape of the bug that sat in
# q2-gglasso.transform_features for months -- invisible because calculate_covariance
# compensated for it. The fixture has a known number of samples, so the test is simply
# whether the DataFrame has that many ROWS.
#
# Uses only the checked-in fixture: no reference database, no reads, no network.

set -uo pipefail          # NOT -e: a failing gate is a recorded result, not a crash
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
MOTUS="${Q2_MOTUS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-mOTUs}"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-motusgate/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-motusgate/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# The conda env is not the env that runs unless we say so. See any other stage for the
# full rationale; ~/.local shadows 19 packages including zarr 2.18.7 -> 3.1.5.
export PYTHONNOUSERSITE=1

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"

python - <<'ENVCHECK' || exit 9
import sys
import numpy, zarr
bad = []
if not zarr.__version__.startswith("2.18"):
    bad.append("zarr %s, expected 2.18.x" % zarr.__version__)
if not numpy.__version__.startswith("2.4"):
    bad.append("numpy %s, expected 2.4.x" % numpy.__version__)
try:
    import zarr.hierarchy  # noqa: F401
except ImportError:
    bad.append("zarr.hierarchy missing -- the plugin transformers cannot load")
if bad:
    sys.stderr.write("ENV CHECK FAILED: %s\n" % "; ".join(bad))
    sys.exit(9)
ENVCHECK

cd "$SCRATCH"
FAIL=0

echo "############ [1/4] install q2-mOTUs (editable, --no-deps) ############"
# --no-deps is mandatory: pip cannot see the conda pins and would happily drop a wheel
# over numpy 2.4.2. motu-profiler is NOT installed here -- import-table only parses a
# TSV, so the reference database and the `motus` binary are not needed for this gate.
pip install -e "$MOTUS" --no-deps -q 2>&1 | tail -3 | sed 's/^/    /'
python -c "import q2_motus, sys; print('    q2_motus at', q2_motus.__file__)" || FAIL=1

echo "############ [2/4] M1: does the CLI see it? ############"
qiime dev refresh-cache 2>&1 | tail -2 | sed 's/^/    /'
if qiime motus --help >motus-help.txt 2>&1; then
  echo "    qiime motus --help: ok"
  grep -E '^\s+(profile|import-table)' motus-help.txt | sed 's/^/      /'
  for a in profile import-table; do
    grep -q -- "$a" motus-help.txt \
      && echo "    ok   action present: $a" \
      || { echo "    FAIL action missing: $a"; FAIL=1; }
  done
else
  echo "    FAIL qiime motus --help did not run:"; sed 's/^/      /' motus-help.txt | head -20; FAIL=1
fi

echo "############ [3/4] M2: import-table on the checked-in fixture ############"
FIXTURE="$MOTUS/q2_motus/tests/data/motus-merged-abundance-1.tsv"
echo "    fixture: $FIXTURE"
python - "$FIXTURE" <<'PY'
import sys, re, pandas as pd, qiime2

fx = sys.argv[1]
raw = open(fx).read().splitlines()
# The mOTUs profile format: comment lines, then a header
#   #mOTU <TAB> consensus_taxonomy <TAB> <sample...>
hdr = [l for l in raw if l.startswith("#mOTU")][0].split("\t")
fixture_samples = hdr[2:]
fixture_rows = [l for l in raw if l and not l.startswith("#")]
print(f"    fixture declares {len(fixture_samples)} samples: {fixture_samples}")
print(f"    fixture data rows (mOTUs): {len(fixture_rows)}")

art = qiime2.Artifact.import_data("MotusMergedAbundanceTable", fx)
from qiime2.plugins import motus
res = motus.actions.import_table(motus_table=art)
tbl, tax = res.table, res.taxonomy
print(f"    table type   : {tbl.type}")
print(f"    taxonomy type: {tax.type}")

df = tbl.view(pd.DataFrame)
tx = tax.view(pd.DataFrame)
checks = []
checks.append(("table is FeatureTable[Frequency]", str(tbl.type) == "FeatureTable[Frequency]", str(tbl.type)))
checks.append(("taxonomy is FeatureData[Taxonomy]", str(tax.type) == "FeatureData[Taxonomy]", str(tax.type)))

# ORIENTATION: QIIME 2 stores a DataFrame index as SAMPLES. The fixture has a known
# sample count, so correct orientation means exactly that many rows.
checks.append((f"rows == {len(fixture_samples)} samples (orientation)",
               df.shape[0] == len(fixture_samples), f"shape={df.shape}"))
checks.append(("more features than samples in the fixture",
               df.shape[1] >= df.shape[0], f"shape={df.shape}"))

pat = re.compile(r"^(ref|meta|ext)_mOTU_v3_\d+$")
feat = list(df.columns)
motu_like = [f for f in feat if pat.match(str(f))]
checks.append(("feature IDs look like mOTU v3 ids", len(motu_like) > 0,
               f"{len(motu_like)}/{len(feat)} match; e.g. {motu_like[:3]}"))
checks.append(("literal 'unassigned' feature is present (must be dropped downstream)",
               "unassigned" in [str(f) for f in feat], str("unassigned" in [str(f) for f in feat])))

vals = df.to_numpy()
integral = bool((vals == vals.astype("int64")).all())
print(f"    counts integer-valued? {integral}   (mode=insert.scaled_counts may be fractional)")

# Taxonomy must be a real hierarchy for q2-classo add-taxa (trac) to aggregate over.
#
# Counting semicolons is NOT sufficient and this check previously false-passed on
# garbage. `.str.replace("|", "; ", regex=True)` treats "|" as regex alternation, which
# matches the empty string everywhere and inserts "; " between every CHARACTER. The
# result is a plausible semicolon-delimited string with a huge "rank count" -- it scored
# "rank depth min=12 max=810" and passed a >=6 check while being complete nonsense.
# So: assert the SHAPE of each rank, not the number of separators.
col = tx.columns[0]
lineages = tx[col].astype(str)
rank_re = re.compile(r"^[kpcofgsm]__")
ranks = lineages.iloc[0].split("; ")
n_ranks = len(ranks)
well_formed = all(rank_re.match(r.strip()) for r in ranks)
shredded = any(len(r.strip()) <= 2 for r in ranks)
print(f"    example lineage: {lineages.iloc[0][:160]}")
print(f"    ranks: {n_ranks}")
checks.append(("lineage has 6-8 ranks", 6 <= n_ranks <= 8, f"{n_ranks}"))
checks.append(("every rank carries a k__/p__/.../m__ prefix", well_formed,
               f"offenders={[r for r in ranks if not rank_re.match(r.strip())][:3]}"))
checks.append(("lineage is not character-shredded", not shredded,
               f"shortest rank={min(len(r.strip()) for r in ranks)}"))
has_m = bool(lineages.str.contains(r"m__").any())
print(f"    lineage uses m__<mOTU_id> in place of s__ ? {has_m}")

bad = 0
for name, ok, got in checks:
    print(f"    {'ok  ' if ok else 'FAIL'} {name:<58} {got}")
    bad += (not ok)
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] || FAIL=1

echo "############ [4/4] verdict ############"
if [[ $FAIL -eq 0 ]]; then
  echo "############ RESULT: PASS -- gates M1 and M2 clear ############"
else
  echo "############ RESULT: FAIL -- see the FAIL lines above ############"
  echo "  If orientation failed, _utils.extract_table_tax's trailing .T is wrong for"
  echo "  QIIME 2 and must be fixed as a COUPLED change (see the q2-gglasso axis-swap"
  echo "  fix: removing the transpose alone produced an N x N covariance)."
fi
exit $FAIL
