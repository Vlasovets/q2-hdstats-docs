#!/usr/bin/env bash
# Final tutorial workflow: 50k filter, temporal holdout, C2 rho=0, no refit.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
: "${QIIME_ENV:=/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7}"
export PATH="$QIIME_ENV/bin:$PATH"
export R_HOME="$QIIME_ENV/lib/R"
export MPLCONFIGDIR="$ROOT/.mplconfig"
export NUMBA_CACHE_DIR="$ROOT/.numba-cache"
export CONDA_PREFIX="$ROOT/.q2cli-env"
mkdir -p "$CONDA_PREFIX/var/q2cli"
ART="$ROOT/artifacts"

qiime feature-table filter-features \
  --i-table "$ART/bracken-species-raw.qza" \
  --p-min-samples 14 --p-min-frequency 50000 \
  --o-filtered-table "$ART/bracken-species-prev14-total50k.qza"

qiime classo transform-features \
  --i-features "$ART/bracken-species-prev14-total50k.qza" \
  --p-transformation clr --p-coef 0.5 \
  --o-x "$ART/bracken-species-prev14-total50k-clr05.qza"

for split in train test; do
  qiime feature-table filter-samples \
    --i-table "$ART/bracken-species-prev14-total50k.qza" \
    --m-metadata-file "$ROOT/metadata.tsv" \
    --p-where "[analysis_set]='$split'" --p-no-filter-empty-features \
    --o-filtered-table "$ART/bracken-species-prev14-total50k-$split.qza"
  qiime classo transform-features \
    --i-features "$ART/bracken-species-prev14-total50k-$split.qza" \
    --p-transformation clr --p-coef 0.5 \
    --o-x "$ART/bracken-species-prev14-total50k-$split-clr05.qza"
done

# The standard numerical method and default CV seed are used without overrides.
qiime classo classify \
  --i-features "$ART/bracken-species-prev14-total50k-train-clr05.qza" \
  --m-y-file "$ROOT/metadata.tsv" --m-y-column bean_type \
  --p-no-path --p-cv --p-no-stabsel \
  --p-cv-subsets 5 --p-cv-one-se --p-cv-lamin 0.001 \
  --p-huber --p-rho 0 --p-no-refit \
  --o-result "$ART/bracken-bean-type-train-clr05-c2-rho0-cv-no-refit-pr60.qza"

for split in train test; do
  qiime classo predict \
    --i-features "$ART/bracken-species-prev14-total50k-$split-clr05.qza" \
    --i-problem "$ART/bracken-bean-type-train-clr05-c2-rho0-cv-no-refit-pr60.qza" \
    --o-predictions "$ART/bracken-bean-type-train-clr05-c2-rho0-cv-no-refit-pr60-$split-predictions.qza"
done

python "$ROOT/scripts/plot_bracken_final_classification.py"
