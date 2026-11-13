from pathlib import Path
import json
import numpy as np
import pandas as pd
from biom import Table
from biom.util import biom_open

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
TOP15 = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'
FULL = ROOT / 'work/atacama-full'
OUT.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(FULL / 'table.tsv', sep='\t', skiprows=[0], index_col=0)
counts.index.name = 'feature-id'
counts = counts.apply(pd.to_numeric, errors='raise')

metadata = pd.read_csv(FULL / 'sample-metadata.tsv', sep='\t')
types = metadata.iloc[0].copy()
metadata = metadata.iloc[1:].copy()
metadata = metadata.set_index('sample-id')

source_samples = list(counts.columns)
split = pd.read_csv(TOP15 / 'train-test-split-53-samples.tsv', sep='\t').set_index('sample-id')
analysis_samples = list(split.index)
train_samples = list(split.index[split['split'] == 'train'])
test_samples = list(split.index[split['split'] == 'test'])

# Selection is defined on the 53 analysis samples used everywhere downstream.
# This avoids all-zero feature columns after dropping YUN1005.1.1 and keeps the
# high-dimensional q2-classo and q2-gglasso examples on the same nondegenerate
# table. Ties are resolved by total relative abundance, then total count, then
# feature ID for determinism.
source_counts = counts[analysis_samples]
lib_sizes = source_counts.sum(axis=0).replace(0, np.nan)
rel = source_counts.div(lib_sizes, axis=1)
ranking = pd.DataFrame({
    'feature_id': source_counts.index,
    'prevalence_samples_53': (source_counts > 0).sum(axis=1).astype(int),
    'total_relative_abundance_53': rel.sum(axis=1),
    'mean_relative_abundance_53': rel.mean(axis=1),
    'total_count_53': source_counts.sum(axis=1),
})
ranking = ranking.sort_values(
    ['prevalence_samples_53', 'total_relative_abundance_53', 'total_count_53', 'feature_id'],
    ascending=[False, False, False, True],
    kind='mergesort',
).reset_index(drop=True)
ranking.insert(0, 'prevalence_rank', np.arange(1, len(ranking) + 1))
top300 = ranking.head(300).copy()
top300.to_csv(OUT / 'top300-prevalence-relative-abundance-ranking-from-53-samples.tsv', sep='	', index=False)

# Both tutorial dimensions use this single deterministic 53-sample ranking.
top15 = ranking.head(15).copy()
top15.to_csv(TOP15 / 'top15-prevalence-ranking-from-53-samples.tsv', sep='\t', index=False)

subset = counts.loc[top300['feature_id'], analysis_samples]
zero_features = subset.index[subset.sum(axis=1) <= 0].tolist()
if zero_features:
    raise ValueError(f'Top-300 53-sample table contains zero-sum features: {zero_features}')
subset.to_csv(OUT / 'top300-count-table-53-samples.tsv', sep='	')

# Learn prediction-feature eligibility using training counts only. The full
# top-300 table remains unchanged for the all-sample network analyses.
training_support = subset[train_samples].sum(axis=1) > 0
model_subset = subset.loc[training_support]
supported_ids = pd.DataFrame(index=model_subset.index)
supported_ids.index.name = 'feature-id'
supported_ids.to_csv(OUT / 'top300-training-supported-feature-ids.tsv', sep='\t')
excluded = top300.loc[~top300['feature_id'].isin(model_subset.index)].copy()
excluded.to_csv(OUT / 'top300-training-absent-features.tsv', sep='\t', index=False)

# BIOM tables for QIIME import: observations/features x samples.
def write_biom(df, path):
    table = Table(df.values.astype(float), observation_ids=list(df.index), sample_ids=list(df.columns))
    with biom_open(str(path), 'w') as fh:
        table.to_hdf5(fh, generated_by='prepare_top300_inputs.py')

write_biom(subset, OUT / 'top300-count-table-53-samples.biom')
write_biom(model_subset[train_samples], OUT / 'top300-count-table-train-53-samples.biom')
write_biom(model_subset[test_samples], OUT / 'top300-count-table-test-53-samples.biom')

top15_subset = counts.loc[top15['feature_id'], analysis_samples]
if (top15_subset.sum(axis=1) <= 0).any():
    raise ValueError('Top-15 53-sample table contains an all-zero feature')
top15_subset.to_csv(TOP15 / 'top15-count-table-53-samples.tsv', sep='\t')
write_biom(top15_subset, TOP15 / 'top15-count-table-53-samples.biom')

# QIIME-compatible split and outcome metadata.
split_q2 = split.copy()
split_q2.index.name = 'sample-id'
split_q2.loc['#q2:types'] = ['categorical', 'numeric', 'categorical', 'categorical', 'numeric', 'numeric']
split_q2 = split_q2.loc[['#q2:types'] + analysis_samples]
split_q2.index.name = 'sample-id'
split_q2.to_csv(OUT / 'train-test-split-53-q2.tsv', sep='\t')

for name, samples in [('train', train_samples), ('test', test_samples), ('all', analysis_samples)]:
    elev = metadata.loc[samples, ['elevation']].copy()
    elev.index.name = 'sample-id'
    elev = pd.concat([pd.DataFrame({'elevation':['numeric']}, index=['#q2:types']), elev])
    elev.index.name = 'sample-id'
    elev.to_csv(OUT / f'elevation-{name}-53-metadata-q2.tsv', sep='\t')

# Covariate panel: interpretable, not near-duplicate by train-set pairwise |r| <= 0.65.
raw_covariates = pd.DataFrame(index=analysis_samples)
raw = metadata.loc[analysis_samples]
for col in ['depth', 'ph', 'ec', 'toc', 'amplicon-concentration', 'extract-concen']:
    raw_covariates[col] = pd.to_numeric(raw[col], errors='coerce')
raw_covariates['log10-amplicon-concentration'] = np.log10(raw_covariates['amplicon-concentration'])
raw_covariates['log10-extract-concen'] = np.log10(raw_covariates['extract-concen'])
selected_covariates = ['log10-extract-concen', 'ph', 'ec', 'log10-amplicon-concentration', 'toc', 'depth']

train_raw = raw_covariates.loc[train_samples, selected_covariates]
means = train_raw.mean(axis=0, skipna=True)
imputed_train = train_raw.fillna(means)
centered_train = imputed_train - means
covariate_train_l2_norms = np.sqrt((centered_train ** 2).sum(axis=0))
if (covariate_train_l2_norms == 0).any():
    raise ValueError(f'Zero norm for covariates: {list(covariate_train_l2_norms[covariate_train_l2_norms == 0].index)}')

# Match covariate scale to the microbial design rather than forcing unit-norm
# covariates. This reproduces q2-classo's CLR transform exactly: zeros are
# replaced by 0.5, log-counts are computed, and each sample row is centered.
train_counts = model_subset[train_samples].T.astype(float)
train_clr = train_counts.mask(train_counts <= 0.0, 0.5)
train_clr = np.log(train_clr)
train_clr = train_clr.sub(train_clr.mean(axis=1), axis=0)
asv_train_l2_norms = np.sqrt((train_clr ** 2).sum(axis=0))
median_asv_train_l2_norm = float(np.median(asv_train_l2_norms))

normalized = raw_covariates[selected_covariates].fillna(means)
normalized = (normalized - means) / covariate_train_l2_norms * median_asv_train_l2_norm
renamed = {c: f'asvscaled-{c}' for c in selected_covariates}
normalized = normalized.rename(columns=renamed)

cov_meta = normalized.copy()
cov_meta.index.name = 'sample-id'
cov_meta = pd.concat([pd.DataFrame({c:['numeric'] for c in cov_meta.columns}, index=['#q2:types']), cov_meta])
cov_meta.index.name = 'sample-id'
cov_meta.to_csv(OUT / 'top300-elevation-augmented-covariates-normalized-from-train.tsv', sep='\t')

for name, samples in [('train', train_samples), ('test', test_samples)]:
    part = normalized.loc[samples].copy()
    part.index.name = 'sample-id'
    part = pd.concat([pd.DataFrame({c:['numeric'] for c in part.columns}, index=['#q2:types']), part])
    part.index.name = 'sample-id'
    part.to_csv(OUT / f'top300-elevation-augmented-covariates-normalized-from-train-{name}.tsv', sep='\t')

summary = pd.DataFrame({
    'covariate': selected_covariates,
    'q2_column': [renamed[c] for c in selected_covariates],
    'train_mean_before_imputation': [means[c] for c in selected_covariates],
    'train_l2_norm_after_centering_before_asv_scaling': [covariate_train_l2_norms[c] for c in selected_covariates],
    'target_median_asv_train_l2_norm': [median_asv_train_l2_norm] * len(selected_covariates),
    'train_l2_norm_after_asv_scaling': [float(np.sqrt((normalized.loc[train_samples, renamed[c]] ** 2).sum())) for c in selected_covariates],
    'missing_train': [int(train_raw[c].isna().sum()) for c in selected_covariates],
    'missing_test': [int(raw_covariates.loc[test_samples, c].isna().sum()) for c in selected_covariates],
})
summary.to_csv(OUT / 'top300-covariate-normalization-summary.tsv', sep='\t', index=False)

corr_input = raw_covariates.loc[train_samples, selected_covariates].fillna(means)
corr = corr_input.corr()
corr.to_csv(OUT / 'top300-selected-covariate-train-correlation.tsv', sep='\t')

manifest = {
    'feature_filter': 'rank by prevalence on the 53 analysis samples, tie-break by total relative abundance, then total count and feature id',
    'n_source_samples_for_filter': len(analysis_samples),
    'n_analysis_samples': len(analysis_samples),
    'n_train_samples': len(train_samples),
    'n_test_samples': len(test_samples),
    'n_network_asvs': len(subset),
    'n_prediction_asvs': len(model_subset),
    'n_asvs_absent_from_training': int((~training_support).sum()),
    'prediction_feature_filter': 'retain only features with positive total counts in training; apply the same feature IDs to held-out samples before CLR',
    'dropped_sample_before_filtering': 'YUN1005.1.1',
    'selected_covariates': selected_covariates,
    'normalized_covariate_columns': [renamed[c] for c in selected_covariates],
    'covariate_scaling': 'training-set mean imputation and centering, then L2 norm matched to median top-300 CLR ASV column norm from the training design',
    'target_median_asv_train_l2_norm': median_asv_train_l2_norm,
    'covariate_weight': 1.0,
    'covariate_weight_note': 'Covariates were scaled to the median training-set CLR-ASV column norm, so covariates and ASV coordinates use equal penalty weights.',
}
(OUT / 'top300-workflow-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps(manifest, indent=2))
