#!/bin/bash
#SBATCH --job-name=q2-docs-art
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/docsart_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/docsart_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Everything the book displays, regenerated from the solution artifacts.
#
# WHY THIS STAGE EXISTS. export_network.py and make_docs_figures.py were both
# run by hand. The network exports consequently went stale without anyone
# noticing: they kept positional ASV-k labels from before the tier-2 rebuild
# while publish/tier2 had moved to real 32-hex feature IDs. Nothing re-ran them
# and nothing checked them, so a node table whose `genus` column silently
# described the wrong organisms sat in results/ for as long as it took someone
# to look. That is the same failure mode as package_release.py sourcing tier 2
# from data/ -- a manual step outside the pipeline.
#
# Anything the chapters render must be produced by a stage, not by a person.
#
# It also ASSERTS the numbers the chapters quote. A re-run that changes the
# model should fail here, loudly, rather than silently produce figures that
# disagree with the prose around them.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda

REGEN="$ROOT/results/tier2-regen"
NET="$ROOT/results/gglasso"
TAX="$ROOT/data/atacama-taxonomy-silva138.qza"

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
mkdir -p "$NET" "$ROOT/slurm/logs"

for f in "$REGEN/sgl-lambda08.qza" "$REGEN/slr-mu15p0.qza" "$REGEN/clr.qza" "$TAX"; do
  [[ -f "$f" ]] || { echo "missing input: $f — run 28_regenerate_tier2_keep_ids.sh first"; exit 1; }
done

echo "############ [1/3] network exports ############"
python "$ROOT/scripts/export_network.py" \
  --solution "$REGEN/sgl-lambda08.qza" --taxonomy "$TAX" --clr "$REGEN/clr.qza" \
  --stem "$NET/atacama-top-300-network-sgl-lambda0.8"
python "$ROOT/scripts/export_network.py" \
  --solution "$REGEN/slr-mu15p0.qza" --taxonomy "$TAX" --clr "$REGEN/clr.qza" \
  --stem "$NET/atacama-top-300-network-slr-lambda0.8-rank2"

echo "############ [2/3] assert the exports match what the chapters say ############"
python - "$NET" <<'PY'
import csv, sys
net = sys.argv[1]
def rows(p): return list(csv.DictReader(open(p), delimiter='\t'))
def keys(p): return {tuple(sorted((r['source'], r['target']))) for r in rows(p)}

S = keys(f"{net}/atacama-top-300-network-sgl-lambda0.8-edges.tsv")
L = keys(f"{net}/atacama-top-300-network-slr-lambda0.8-rank2-edges.tsv")
nodes = rows(f"{net}/atacama-top-300-network-sgl-lambda0.8-nodes.tsv")

# Feature identity: the whole reason tier 2 was rebuilt. A positional ASV-k label
# here means the export was built from a pre-rebuild artifact.
hexish = sum(1 for n in nodes
             if len(n['asv']) == 32 and all(c in '0123456789abcdef' for c in n['asv'].lower()))
unresolved = sum(1 for n in nodes if not n.get('genus') or n['genus'].strip().lower()
                 in ('', 'unassigned', 'nan'))

checks = [
    ("sgl edges == 216",            len(S) == 216,        len(S)),
    ("slr rank-2 edges == 202",     len(L) == 202,        len(L)),
    ("removed by rank 2 == 14",     len(S - L) == 14,     len(S - L)),
    ("added by rank 2 == 0",        len(L - S) == 0,      len(L - S)),
    ("all 300 nodes are hex IDs",   hexish == 300,        hexish),
    ("no unresolved genus",         unresolved == 0,      unresolved),
]
bad = 0
for name, ok, got in checks:
    print(f"  {'ok  ' if ok else 'FAIL'} {name:32} got {got}")
    bad += not ok
if bad:
    print(f"\n  {bad} assertion(s) failed. The chapters quote 216/202/14/0 and real "
          "feature IDs; either the model changed or the export is built from the "
          "wrong artifact. Fix one before regenerating the figures.")
    sys.exit(1)
PY

echo "############ [3/3] figures the chapters render ############"
python "$ROOT/scripts/make_docs_figures.py"

echo "############ RESULT: PASS ############"
