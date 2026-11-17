from __future__ import annotations

from pathlib import Path
import io
import json
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd
import zarr

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
TOP15 = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'
SCRIPT = Path(__file__)
sys.path.insert(0, str(ROOT / 'outputs/fresh-atacama-workflow/tutorial-figure-generation'))
from provenance_utils import write_provenance

META = pd.read_csv(TOP15 / 'atacama-53-sample-metadata.tsv', sep='\t', comment='#').set_index('sample-id')
LABELS = pd.read_csv(OUT / 'top300-asv-label-map.tsv', sep='\t').set_index('feature_id')


def open_problem(path: Path):
    with zipfile.ZipFile(path) as artifact:
        name = next(n for n in artifact.namelist() if n.endswith('/data/problem.zip'))
        data = artifact.read(name)
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(tmp)
    return zarr.open(tmp, mode='r')


def artifact_arrays(path: Path) -> dict[str, np.ndarray]:
    with zipfile.ZipFile(path) as artifact:
        name = next(n for n in artifact.namelist() if n.endswith('/data/problem.zip'))
        data = artifact.read(name)
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(tmp)
    root = zarr.open(tmp, mode='r')
    return {k: np.asarray(root[k]) for k in root.array_keys()}


def display_label(feature: str) -> str:
    cov = {
        'asvscaled-log10-extract-concen': 'log10 extract conc.',
        'asvscaled-ph': 'pH',
        'asvscaled-ec': 'EC',
        'asvscaled-log10-amplicon-concentration': 'log10 amplicon conc.',
        'asvscaled-toc': 'TOC',
        'asvscaled-depth': 'depth',
    }
    if feature in cov:
        return cov[feature]
    if feature == 'intercept':
        return 'intercept'
    return LABELS.loc[feature, 'display_label'] if feature in LABELS.index else feature


def source(feature: str) -> str:
    if feature == 'intercept':
        return 'intercept'
    return 'covariate' if feature.startswith('asvscaled-') else 'asv'


def regression_metrics(y: np.ndarray, yhat: np.ndarray) -> dict[str, float]:
    resid = yhat - y
    r = float(np.corrcoef(y, yhat)[0, 1]) if len(y) > 1 else np.nan
    return {
        'n': int(len(y)),
        'pearson_r': r,
        'r2': float(1.0 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)),
        'rmse': float(np.sqrt(np.mean(resid ** 2))),
        'mae': float(np.mean(np.abs(resid))),
    }


def extract_prediction(path: Path, observed_col: str, split: str) -> pd.DataFrame:
    arr = artifact_arrays(path)
    samples = arr['sample_labels'].astype(str)
    pred_key = next(k for k in arr if k.startswith('Yhat'))
    out = pd.DataFrame({'sample-id': samples, 'set': split, 'prediction': arr[pred_key].astype(float)})
    out[f'observed_{observed_col}'] = META.loc[samples, observed_col].to_numpy()
    return out


def solution_table(problem, key: str, model: str, coef_col: str) -> pd.DataFrame:
    g = problem[f'solution/{key}']
    labels = np.asarray(g['label']).astype(str)
    selected = np.asarray(g['selected_param']).astype(bool)
    refit = np.asarray(g['refit']).astype(float)
    beta = np.asarray(g['beta']).astype(float) if 'beta' in g else refit
    return pd.DataFrame({
        'model': model,
        'feature': labels,
        'source': [source(x) for x in labels],
        'display_label': [display_label(x) for x in labels],
        'selected': selected,
        'beta': beta,
        coef_col: refit,
    })


def extract_cv_curve(problem, key: str = 'CV') -> pd.DataFrame:
    g = problem[f'solution/{key}']
    lam = np.asarray(g['xGraph']).astype(float)
    return pd.DataFrame({
        'lambda': lam,
        'cv_error': np.asarray(g['yGraph']).astype(float),
        'standard_error': np.asarray(g['standard_error']).astype(float),
        'is_lambda_min': np.arange(len(lam)) == int(g.attrs['index_min']),
        'is_lambda_1se': np.arange(len(lam)) == int(g.attrs['index_1SE']),
    })


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    vals = sorted(set(y_true))
    recalls = []
    for val in vals:
        mask = y_true == val
        recalls.append(float(np.mean(y_pred[mask] == val)) if mask.any() else np.nan)
    return float(np.nanmean(recalls))


def class_metrics(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split, d in pred.groupby('set', sort=False):
        rows.append({
            'set': split,
            'n': len(d),
            'accuracy': float(d['correct'].mean()),
            'balanced_accuracy': balanced_accuracy(d['observed_vegetation'].to_numpy(), d['predicted_vegetation'].to_numpy()),
        })
    return pd.DataFrame(rows)


def extract_regression():
    specs = [
        ('ASV-only', 'atacama-top300-elevation-asv-only-r3-cv.qza',
         'atacama-top300-elevation-asv-only-r3-train-predictions.qza',
         'atacama-top300-elevation-asv-only-r3-test-predictions.qza'),
        ('ASV + covariates', 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza',
         'atacama-top300-elevation-covariate-augmented-equalweights-r3-train-predictions.qza',
         'atacama-top300-elevation-covariate-augmented-equalweights-r3-test-predictions.qza'),
    ]
    pred_rows, metric_rows, coef_rows, cv_rows = [], [], [], []
    inputs = []
    for model, problem_fn, train_fn, test_fn in specs:
        problem_path = OUT / problem_fn
        problem = open_problem(problem_path)
        inputs.append(problem_path)
        train = extract_prediction(OUT / train_fn, 'elevation', 'train')
        test = extract_prediction(OUT / test_fn, 'elevation', 'test')
        for d in (train, test):
            d['model'] = model
            d['observed_elevation'] = d['observed_elevation'].astype(float)
            d['predicted_elevation'] = d.pop('prediction')
            d['transect-name'] = META.loc[d['sample-id'], 'transect-name'].to_numpy()
            d['site-name'] = META.loc[d['sample-id'], 'site-name'].to_numpy()
            pred_rows.append(d)
        metrics = {'model': model}
        for split, d in [('train', train), ('test', test)]:
            m = regression_metrics(d['observed_elevation'].to_numpy(float), d['predicted_elevation'].to_numpy(float))
            metrics.update({f'{split}_{k}': v for k, v in m.items()})
        sol = problem['solution/CV']
        coef = solution_table(problem, 'CV', model, 'coefficient')
        metrics.update({
            'lambda_min': float(sol.attrs['lambda_min']),
            'lambda_1se': float(sol.attrs['lambda_1SE']),
            'sigma': float(sol.attrs.get('sigma', np.nan)),
            'selected_total': int(coef.query("feature != 'intercept' and selected").shape[0]),
            'selected_asvs': int(coef.query("source == 'asv' and selected").shape[0]),
            'selected_covariates': int(coef.query("source == 'covariate' and selected").shape[0]),
        })
        metric_rows.append(metrics)
        coef_rows.append(coef.rename(columns={'coefficient': 'coefficient'}))
        cv = extract_cv_curve(problem)
        cv.insert(0, 'model', model)
        cv_rows.append(cv)
        inputs.extend([OUT / train_fn, OUT / test_fn])

    pred = pd.concat(pred_rows, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    coefs = pd.concat(coef_rows, ignore_index=True)
    selected = coefs.query("feature != 'intercept' and selected").copy()
    cv = pd.concat(cv_rows, ignore_index=True)
    outputs = [
        OUT / 'top300-elevation-asv-vs-covariate-augmented-predictions.tsv',
        OUT / 'top300-elevation-asv-vs-covariate-augmented-metrics.tsv',
        OUT / 'top300-elevation-asv-vs-covariate-augmented-coefficients.tsv',
        OUT / 'top300-elevation-asv-vs-covariate-augmented-selected-coefficients.tsv',
        OUT / 'top300-elevation-asv-vs-covariate-augmented-cv-curves.tsv',
    ]
    pred.to_csv(outputs[0], sep='\t', index=False)
    metrics.to_csv(outputs[1], sep='\t', index=False)
    coefs.to_csv(outputs[2], sep='\t', index=False)
    selected.to_csv(outputs[3], sep='\t', index=False)
    cv.to_csv(outputs[4], sep='\t', index=False)
    write_provenance(
        OUT / 'top300-elevation-asv-vs-covariate-augmented-qza-extraction',
        description='TSV extraction from regenerated q2-classo top-300 elevation R3 QIIME artifacts',
        project_root=ROOT,
        script=SCRIPT,
        inputs=inputs + [TOP15 / 'atacama-53-sample-metadata.tsv', OUT / 'top300-asv-label-map.tsv'],
        outputs=outputs,
        parameters={'models': ['ASV-only R3', 'ASV + covariates R3'], 'source': 'QIIME artifacts only'},
    )


def extract_classification():
    specs = [
        ('cv', 'atacama-top300-vegetation-classification-cv.qza',
         'atacama-top300-vegetation-classification-cv-train-predictions.qza',
         'atacama-top300-vegetation-classification-cv-test-predictions.qza',
         'refit_beta', 'CV'),
        ('stabsel', 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04.qza',
         'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-train-predictions.qza',
         'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-test-predictions.qza',
         'stabsel_refit_beta', 'StabSel'),
    ]
    inputs = []
    output_paths = []
    for kind, problem_fn, train_fn, test_fn, coef_col, sol_key in specs:
        problem = open_problem(OUT / problem_fn)
        inputs.append(OUT / problem_fn)
        pred_parts = []
        for split, fn in [('train', train_fn), ('test', test_fn)]:
            d = extract_prediction(OUT / fn, 'vegetation', split)
            d['score'] = d.pop('prediction')
            d['predicted_vegetation'] = np.where(d['score'] >= 0, 'no', 'yes')
            d['correct'] = d['predicted_vegetation'].eq(d['observed_vegetation'])
            d['elevation'] = META.loc[d['sample-id'], 'elevation'].astype(float).to_numpy()
            pred_parts.append(d)
            inputs.append(OUT / fn)
        pred = pd.concat(pred_parts, ignore_index=True)
        metrics = class_metrics(pred)
        coefs = solution_table(problem, sol_key, kind, coef_col)

        if kind == 'cv':
            cv = extract_cv_curve(problem)
            cv_path = OUT / 'atacama-top300-vegetation-classification-cv-curve.tsv'
            cv.to_csv(cv_path, sep='\t', index=False)
            output_paths.append(cv_path)
            prefix = OUT / 'atacama-top300-vegetation-classification-cv'
        else:
            g = problem['solution/StabSel']
            profile = pd.DataFrame({
                'feature': np.asarray(g['label']).astype(str),
                'stability': np.asarray(g['distribution']).astype(float),
            })
            profile['source'] = profile['feature'].map(source)
            profile['display_label'] = profile['feature'].map(display_label)
            profile = profile.query("feature != 'intercept'").sort_values('stability', ascending=False)
            prof_path = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-profile.tsv'
            profile.to_csv(prof_path, sep='\t', index=False)
            output_paths.append(prof_path)
            prefix = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04'

        paths = [
            Path(str(prefix) + '-train-test-predictions.tsv'),
            Path(str(prefix) + '-metrics.tsv'),
            Path(str(prefix) + '-coefficients.tsv'),
        ]
        pred.to_csv(paths[0], sep='\t', index=False)
        metrics.to_csv(paths[1], sep='\t', index=False)
        coefs.to_csv(paths[2], sep='\t', index=False)
        output_paths.extend(paths)

    write_provenance(
        OUT / 'atacama-top300-vegetation-classification-qza-extraction',
        description='TSV extraction from regenerated q2-classo top-300 vegetation QIIME artifacts',
        project_root=ROOT,
        script=SCRIPT,
        inputs=inputs + [TOP15 / 'atacama-53-sample-metadata.tsv', OUT / 'top300-asv-label-map.tsv'],
        outputs=output_paths,
        parameters={'models': ['C1 CV', 'C1 stability selection'], 'positive_score_class': 'no', 'source': 'QIIME artifacts only'},
    )


def main():
    extract_regression()
    extract_classification()
    print(json.dumps({'status': 'ok', 'output_dir': str(OUT)}, indent=2))


if __name__ == '__main__':
    main()
