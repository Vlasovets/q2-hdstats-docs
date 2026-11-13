#!/usr/bin/env bash
# QIIME 2-native S3 regression stage: complete mOTU table, SLC versus TRAC.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/work/micromamba-root/envs/q2-hdstats-qiime2-2026.7"
export PATH="$ENV/bin:$PATH"
export R_HOME="$ENV/lib/R"
export MPLCONFIGDIR="$ROOT/.mplconfig"
export NUMBA_CACHE_DIR="$ROOT/.numba_cache"

ART="$ROOT/artifacts"
EXP="$ROOT/exports"
mkdir -p "$ART" "$EXP"

# Complete source table -> QIIME 2 inputs. The importer explicitly drops the
# mOTUs 'unassigned' complement before this workflow's QIIME actions.
"$ENV/bin/python" "$ROOT/scripts/prepare_complete_motus_inputs.py" "$ROOT"

# Retain samples with at least 100 assigned marker-gene counts: 36 -> 34.
qiime feature-table filter-samples \
  --i-table "$ART/motus-complete-profiled.qza" \
  --p-min-frequency 100 \
  --o-filtered-table "$ART/motus-depth-ge100.qza"

# The single feature filter: prevalence at least 4 of 34 retained samples (10%).
# It retains all eligible mOTUs; there is deliberately no top-abundance filter.
qiime feature-table filter-features \
  --i-table "$ART/motus-depth-ge100.qza" \
  --p-min-samples 4 \
  --o-filtered-table "$ART/motus-prev10pct.qza"

qiime tools export \
  --input-path "$ART/motus-prev10pct.qza" \
  --output-path "$EXP/motus-prev10pct-table"

# Fixed subject-level holdout. The test samples come from four infants absent
# from training; CV and stability selection below operate on training only.
qiime feature-table filter-samples \
  --i-table "$ART/motus-prev10pct.qza" \
  --m-metadata-file "$EXP/motus-outcomes-all36.tsv" \
  --p-where "[split]='train'" \
  --p-no-filter-empty-features \
  --o-filtered-table "$ART/motus-prev10pct-train.qza"

qiime feature-table filter-samples \
  --i-table "$ART/motus-prev10pct.qza" \
  --m-metadata-file "$EXP/motus-outcomes-all36.tsv" \
  --p-where "[split]='test'" \
  --p-no-filter-empty-features \
  --o-filtered-table "$ART/motus-prev10pct-test.qza"

# Both models use CLR with zeros replaced by 0.5 and the same age endpoint.
qiime classo transform-features \
  --i-features "$ART/motus-prev10pct-train.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "$ART/motus-clr-pseudocount0.5-train.qza"

qiime classo transform-features \
  --i-features "$ART/motus-prev10pct-test.qza" \
  --p-transformation clr \
  --p-coef 0.5 \
  --o-x "$ART/motus-clr-pseudocount0.5-test.qza"

# SLC: ordinary sparse log-contrast regression (zero-sum constraint).
qiime classo regress \
  --i-features "$ART/motus-clr-pseudocount0.5-train.qza" \
  --m-y-file "$EXP/motus-outcomes-all36.tsv" \
  --m-y-column host_age_days \
  --p-concomitant \
  --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.05 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
  --p-cv--nlam 60 --p-cv-lamin 0.05 --p-cv-logscale \
  --p-stabsel --p-stabsel-seed 1 --p-stabsel-b 100 --p-stabsel-q 10 \
  --p-stabsel-threshold 0.4 --p-stabsel-threshold-label 0.2 \
  --p-no-lamfixed \
  --o-result "$ART/motus-age-slc-r3.qza"

qiime classo predict \
  --i-features "$ART/motus-clr-pseudocount0.5-train.qza" \
  --i-problem "$ART/motus-age-slc-r3.qza" \
  --o-predictions "$ART/motus-age-slc-r3-train-predictions.qza"

qiime classo predict \
  --i-features "$ART/motus-clr-pseudocount0.5-test.qza" \
  --i-problem "$ART/motus-age-slc-r3.qza" \
  --o-predictions "$ART/motus-age-slc-r3-test-predictions.qza"

# TRAC: same regression settings, with taxonomy-defined aggregation weights.
qiime classo add-taxa \
  --i-features "$ART/motus-clr-pseudocount0.5-train.qza" \
  --i-taxa "$ART/motus-complete-taxonomy.qza" \
  --o-x "$ART/motus-clr-trac-train.qza" \
  --o-aweights "$ART/motus-trac-weights.qza"

qiime classo add-taxa \
  --i-features "$ART/motus-clr-pseudocount0.5-test.qza" \
  --i-taxa "$ART/motus-complete-taxonomy.qza" \
  --o-x "$ART/motus-clr-trac-test.qza" \
  --o-aweights "$ART/motus-trac-weights-test.qza"

qiime classo regress \
  --i-features "$ART/motus-clr-trac-train.qza" \
  --i-weights "$ART/motus-trac-weights.qza" \
  --m-y-file "$EXP/motus-outcomes-all36.tsv" \
  --m-y-column host_age_days \
  --p-concomitant \
  --p-path --p-path-nlam-log 60 --p-path-lamin-log 0.05 \
  --p-cv --p-cv-subsets 5 --p-cv-seed 1 --p-no-cv-one-se \
  --p-cv--nlam 60 --p-cv-lamin 0.05 --p-cv-logscale \
  --p-stabsel --p-stabsel-seed 1 --p-stabsel-b 100 --p-stabsel-q 10 \
  --p-stabsel-threshold 0.4 --p-stabsel-threshold-label 0.2 \
  --p-no-lamfixed \
  --o-result "$ART/motus-age-trac-r3.qza"

qiime classo predict \
  --i-features "$ART/motus-clr-trac-train.qza" \
  --i-problem "$ART/motus-age-trac-r3.qza" \
  --o-predictions "$ART/motus-age-trac-r3-train-predictions.qza"

qiime classo predict \
  --i-features "$ART/motus-clr-trac-test.qza" \
  --i-problem "$ART/motus-age-trac-r3.qza" \
  --o-predictions "$ART/motus-age-trac-r3-test-predictions.qza"

qiime classo summarize \
  --i-problem "$ART/motus-age-slc-r3.qza" \
  --p-maxplot 40 \
  --o-visualization "$ART/motus-age-slc-r3-summary.qzv"

qiime classo summarize \
  --i-problem "$ART/motus-age-trac-r3.qza" \
  --p-maxplot 40 \
  --o-visualization "$ART/motus-age-trac-r3-summary.qzv"

qiime tools export \
  --input-path "$ART/motus-age-slc-r3-summary.qzv" \
  --output-path "$EXP/motus-age-slc-r3-summary"

qiime tools export \
  --input-path "$ART/motus-age-trac-r3-summary.qzv" \
  --output-path "$EXP/motus-age-trac-r3-summary"

for artifact in \
  "$ART/motus-age-slc-r3-train-predictions.qza" \
  "$ART/motus-age-slc-r3-test-predictions.qza" \
  "$ART/motus-age-trac-r3-train-predictions.qza" \
  "$ART/motus-age-trac-r3-test-predictions.qza"; do
  name="$(basename "$artifact" .qza)"
  qiime tools export --input-path "$artifact" --output-path "$EXP/$name"
done
