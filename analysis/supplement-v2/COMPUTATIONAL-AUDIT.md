# Final Computational Audit

Audit date: 2026-11-12. Model fitting was performed through the established QIIME 2 environment and its plugins, not direct c-lasso/GGLasso solver calls. This is an established-environment replay, not a newly provisioned cross-platform installation test.

## Source Pins

- c-lasso: GitHub commit `d22da8d5ea54eb499de2fce82fa6c0e1f0c2dcc0` (PR60); installed Python source matches the retained GitHub source, excluding generated version metadata.
- q2-classo: `ae574af47b7d79761c01bb8ffcea161716387d11`.
- q2-gglasso: `3f586aa1e074727f497e59481ec1565228042cd3`.
- QIIME 2: 2026.7.0; GGLasso: 0.3.1.

See `environment.json` and `environment/` for source hashes, snapshots, observed pip packages, and the macOS conda explicit specification. The c-lasso version string alone does not identify this GitHub build; use its commit and source hashes.

## Resolved Training-Support Issue

Atacama retains 53 samples, with the established 39/14 training/test split. The global prevalence-ranked table contains 300 ASVs and no globally empty feature columns. However, 15 of these ASVs are absent from every training sample. The approved correction removes those 15 predictors before the prediction-model CLR transform and uses the same 285 columns for held-out prediction. It does not rerank or rename ASVs. The complete 300-ASV table remains the input to full-data network analyses.

The canonical `outputs/fresh-atacama-workflow/q2-classo-top300/` folder contains the retained feature whitelist, excluded-feature table, supported 53-sample artifact, corrected preparation script, dynamic constraint/weight script, and QIIME workflow. Covariate scaling uses the median training CLR column norm, 5.594393674358808; the joint zero-sum constraint covers only the 285 ASVs and leaves six covariates unconstrained. All 291 predictors have unit penalty weights.

Independent replay of the corrected R3 fits matches. Both corrected classification fits also match identical-input repeats. Selected predictors and rounded prediction metrics remain unchanged.

## Final Prediction Results

| Example | Selection | Held-out result |
| --- | --- | --- |
| Atacama top15 elevation R3 CV | 2 ASVs | r=0.6715, RMSE=629.75 |
| Atacama top15 elevation R3 stability | 4 ASVs; threshold 1/3, q=4, B=100 | r=0.7083, RMSE=590.78 |
| Atacama 285-supported elevation R3 | 5 ASVs; lambda=0.007626985859 | r=0.6650, RMSE=613.73 |
| Atacama 285-supported elevation R3 plus covariates | 13 ASVs + 2 covariates; lambda=0.21268217473 | r=0.8784, RMSE=421.36 |
| Atacama 285-supported vegetation C1 CV | 17 ASVs | accuracy 10/14 |
| Atacama 285-supported vegetation C1 stability | 9 ASVs | accuracy 11/14 |
| mOTU SLC | 3 terms | r=0.72, RMSE=9.3 |
| mOTU TRAC CV | 2 terms | r=0.47, RMSE=11.7 |
| mOTU TRAC stability | 4 terms | r=0.532, RMSE=11.28 |
| Bracken temporal classification C2, rho=0 | 5 species; lambda=0.327454916288 | accuracy 3/4 |

mOTU uses n=34, p=60 with a subject-disjoint 26/8 regression split. All 60 features have training support. The complete upstream filtering and SLC/TRAC workflows reproduce. Stability selection is strictly `probability > threshold`, matching the backend, including plotting selection colors.

The final Bracken example uses the 50,000-count filtering rule, first five time points per variety for training, and final two time points per variety for testing. It uses the standard C2 solver, rho=0, default CV seed behavior, and no refit. A repeat matches. The 144-hour Forasteiro sample SRR9640350 has score +0.01992792962 and remains misclassified. Do not substitute earlier C1 diagnostics or historical scores.

## Network Replays

All eight final network model artifacts reproduce through q2-gglasso with their provenance-derived grids. Matrix labels, covariance normalization/scaling, and numerical validity were checked.

| Data | Model | gamma | lambda | mu | Rank | Edges |
| --- | --- | --- | --- | --- | --- | --- |
| Atacama top15 mCLR | SGL | 0.01 | 0.3 | - | - | 13 |
| Atacama top15 mixed | SGL | 0.01 | 0.3 | - | - | 100 total; 4 ASV-ASV |
| Atacama top15 mCLR | Sparse + low rank | 0.01 | 0.3 | 1.5 | 1 | 6 sparse |
| Atacama top300 CLR | SGL | 0.25 | 0.65 | - | - | 408 |
| Atacama top300 mixed | SGL | 0.25 | 0.7 | - | - | 1282 total; 146 ASV-ASV |
| Atacama top300 CLR | Sparse + low rank | 0.25 | 0.65 fixed | 15 | 2 | 376 sparse |
| mOTU mCLR | SGL | 0.10 | 0.25 | - | - | 276 |
| mOTU CLR | Sparse + low rank | 0.10 | 0.25 | 5 | 1 | 248 sparse |

No selected searched parameter lies on its grid boundary; the fixed top300 latent-model lambda is not a searched parameter. Atacama robust-PCA plots were checked against the six original QIIME visualizations and exported coordinates. Their axis percentages are intentionally omitted.

## Reproduction and Figure Provenance

- `commands.sh` plus stage-specific command records and `logs/` preserve executed QIIME commands. Exact input/model paths are in the stage check JSON files.
- `training-support/checks.json` and `training-support/contracts/checks.json` are authoritative for corrected high-dimensional prediction; original `design-contracts/checks.json` documents the discovery, not an unresolved final failure.
- Other authoritative stage reports: `motus/checks.json`, `networks/checks.json`, `PCA-exports/checks.json`, `metadata-figures/checks.json`, `training-support/classification-repeatability/checks.json`, and `bracken/repeatability-checks.json`.
- `figure-registry.json` and `.tsv` map all 26 included figures to source PDFs, generator scripts, inputs, hashes, and sidecars. Installed figures receive individual provenance JSONs. Existing plotters and historical outputs are preserved.
- The first two figures were rebuilt from original metadata and QIIME transforms, including raw-scale mean imputation, geographic ordering, and average-linkage correlation clustering.
- Audit driver scripts are in the parent `outputs/` directory, named `audit_*`, `repeat_supported_classification.py`, `finalize_tutorial_figure_registry.py`, and `finish_tutorial_audit.py`. The final Bracken runner is `outputs/bracken-beans/scripts/run_bracken_final_classification_qiime.sh`.

Historical comparison reports retain backend-related differences and superseded zero-support fits. Do not combine them with the corrected release artifacts or delete their provenance.

## Publication Boundary

Before publishing the online tutorial, package the required inputs with portable paths and an installation recipe pinned to these GitHub revisions. The supplied Bracken QZAs currently live in project assets; confirm their permitted public distribution and stable download location. No public-data availability or clean-machine rebuild is claimed by this audit.

## Final Installation Check

All 26 included PDFs match their registered sources, with no missing-provenance flags. The corrected-stage reports listed above contain 101 passing checks and no failures. The edited supplement compiles with zero undefined citations and references. Updated prediction, classification, and bean-result pages were visually inspected. The newly introduced long command overflow was corrected; four pre-existing minor overfull-box locations elsewhere in the appendix remain outside this computational update. The prevalence-ranking prose now explicitly matches the 53-sample implementation and its relative-abundance/count/identifier tie breaks.
