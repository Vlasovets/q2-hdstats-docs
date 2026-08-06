#!/bin/bash
#SBATCH --job-name=q2-orient
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/orient_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute/slurm/logs/orient_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Settle the transform_features axis question EMPIRICALLY and write a report.
#
# It is still open because two plausible explanations survive: either QIIME 2's
# DataFrame->FeatureTable transformer does not use the samples-x-features
# convention assumed, or something after transform_features' own `X = X.T`
# re-orients the frame. transform_features DOES transpose to (N, p) at
# _func.py:88, yet mclr.qza still reads back with sample IDs as biom
# observations. One of those two must be wrong.
#
# This job decides it with unambiguous labels and writes
# reports/ORIENTATION_FINDING.md. It changes NO source.

set -euo pipefail
ROOT=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-recompute
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
GG=/home/itg/oleg.vlasovets/slr_example/q2-gglasso

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-hdstats/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-hdstats/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$ROOT/reports"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

python - "$ROOT/reports/ORIENTATION_FINDING.md" "$GG" <<'PY'
import sys, pandas as pd, biom, qiime2

out_md, gg = sys.argv[1], sys.argv[2]
lines = ["# transform_features orientation — empirical finding", ""]

# ---- 1. What does the DataFrame -> FeatureTable transformer actually do? ----
df = pd.DataFrame([[1., 2.], [3., 4.], [5., 6.]],
                  index=["SAMPLE_a", "SAMPLE_b", "SAMPLE_c"],
                  columns=["FEATURE_x", "FEATURE_y"])
art = qiime2.Artifact.import_data("FeatureTable[Frequency]", df)
tab = art.view(biom.Table)
samples, obs = set(tab.ids("sample")), set(tab.ids("observation"))
convention_holds = samples == {"SAMPLE_a", "SAMPLE_b", "SAMPLE_c"}

lines += [
    "## 1. The transformer's convention",
    "",
    "Imported a 3x2 DataFrame indexed `SAMPLE_a..c` with columns `FEATURE_x/y`.",
    "",
    f"- stored biom **sample** ids: `{sorted(tab.ids('sample'))}`",
    f"- stored biom **observation** ids: `{sorted(tab.ids('observation'))}`",
    "",
    ("**The DataFrame index is treated as SAMPLES** — the documented convention "
     "holds." if convention_holds else
     "**The DataFrame index is treated as OBSERVATIONS (features)** — the "
     "convention is the opposite of what the plugin assumes."),
    "",
]

# ---- 2. Round-trip the real function -----------------------------------------
counts = qiime2.Artifact.load(f"{gg}/data/atacama-counts.qza")
taxo = qiime2.Artifact.load(f"{gg}/data/classification.qza")
md = qiime2.Metadata.load(f"{gg}/data/selected-atacama-sample-metadata.tsv")
raw = counts.view(biom.Table).to_dataframe()

from qiime2.plugins import gglasso as gg_actions
res = gg_actions.actions.transform_features(
    table=counts, taxonomy=taxo, sample_metadata=md, transformation="mclr")
outdf = res.transformed_table.view(biom.Table).to_dataframe()

raw_obs_are_features = not str(raw.index[0]).startswith("BAQ") and not str(raw.index[0]).startswith("YUN")
out_obs_are_samples = str(outdf.index[0]).startswith(("BAQ", "YUN"))

lines += [
    "## 2. Round-trip through the real action",
    "",
    f"- `atacama-counts.qza` -> `to_dataframe()` = **{raw.shape}**, index[0] = `{raw.index[0]}`",
    f"- `transform-features` output -> `to_dataframe()` = **{outdf.shape}**, index[0] = `{outdf.index[0]}`",
    "",
    (f"The output's biom **observations are sample IDs** — the artifact's axes are "
     f"swapped relative to the input." if out_obs_are_samples else
     "The output's observations are feature IDs — orientation is preserved."),
    "",
    "## 3. Conclusion",
    "",
]

if convention_holds and out_obs_are_samples:
    lines += [
        "The transformer follows the samples-x-features convention, and "
        "`transform_features` still emits a swapped artifact. So the swap is "
        "introduced by the function itself: it reads `(p, N)`, transposes to "
        "`(N, p)` at `_func.py:88`, and something after that — most likely the "
        "metadata/relabelling path — leaves the frame indexed by features again, "
        "or the transpose is applied to a frame that was already `(N, p)`.",
        "",
        "**Actionable:** the bug is inside `transform_features`, not in QIIME 2. "
        "Fixing it is safe for consumers that are orientation-agnostic (`pca` "
        "already is) but WILL change every stored clr/mclr artifact, so it must "
        "be paired with regenerating the published bundle.",
    ]
elif not convention_holds:
    lines += [
        "The transformer does NOT use the samples-x-features convention, so "
        "`transform_features`' `X = X.T` is what creates the mismatch. Removing "
        "that transpose would store the artifact correctly — but every existing "
        "clr/mclr artifact was written with it, so old and new artifacts would "
        "disagree.",
        "",
        "**Actionable:** do NOT change the producer in isolation. Make consumers "
        "orientation-agnostic first (as `pca` now is), then flip the producer and "
        "regenerate the bundle in the same release.",
    ]
else:
    lines += [
        "Orientation is preserved end to end — the earlier observation of swapped "
        "axes did not reproduce here. Re-check how the mclr artifact under test "
        "was produced before acting on the escalation.",
    ]

open(out_md, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
PY

echo "wrote $ROOT/reports/ORIENTATION_FINDING.md"
