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
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'
SCRIPT = Path(__file__)
sys.path.insert(0, str(ROOT / 'outputs/fresh-atacama-workflow/tutorial-figure-generation'))
from provenance_utils import write_provenance

META = pd.read_csv(OUT / 'atacama-53-sample-metadata.tsv', sep='\t', comment='#').set_index('sample-id')
LABELS = pd.read_csv(OUT / 'top15-asv-label-map.tsv', sep='\t').set_index('feature_id')


def open_problem(path: Path):
    with zipfile.ZipFile(path) as artifact:
        name = next(n for n in artifact.namelist() if n.endswith('/data/problem.zip'))
        data = artifact.read(name)
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(tmp)
    return zarr.open(tmp, mode='r')


def artifact_arrays(path: Path) -> dict[str, np.ndarray]:
    root = open_problem(path)
    return {k: np.asarray(root[k]) for k in root.array_keys()}


def display_label(feature: str) -> str:
    if feature == 'intercept':
        return 'intercept'
    return LABELS.loc[feature, 'display_label'] if feature in LABELS.index else feature


def extract_prediction(path: Path, observed_col: str, split: str) -> pd.DataFrame:
    arr = artifact_arrays(path)
    pred_key = next(k for k in arr if k.startswith('Yhat'))
    samples = arr['sample_labels'].astype(str)
    out = pd.DataFrame({'sample-id': samples, 'set': split, 'prediction': arr[pred_key].astype(float)})
    out[f'observed_{observed_col}'] = META.loc[samples, observed_col].to_numpy()
    return out


def regression_metrics(y: np.ndarray, yhat: np.ndarray) -> dict[str, float]:
    resid = yhat - y
    return {
        'n': int(len(y)),
        'pearson_r': float(np.corrcoef(y, yhat)[0, 1]) if len(y) > 1 else np.nan,
        'r2': float(1.0 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)),
        'rmse': float(np.sqrt(np.mean(resid ** 2))),
        'mae': float(np.mean(np.abs(resid))),
    }


def solution_table(problem, key: str, model: str, coef_col: str) -> pd.DataFrame:
    g = problem[f'solution/{key}']
    labels = np.asarray(g['label']).astype(str)
    selected = np.asarray(g['selected_param']).astype(bool)
    refit = np.asarray(g['refit']).astype(float)
    beta = np.asarray(g['beta']).astype(float) if 'beta' in g else refit
    return pd.DataFrame({
        'model': model,
        'feature': labels,
        'display_label': [display_label(x) for x in labels],
        'selected': selected,
        'beta': beta,
        coef_col: refit,
    })


def cv_curve(problem, model: str) -> pd.DataFrame:
    g = problem['solution/CV']
    lam = np.asarray(g['xGraph']).astype(float)
    return pd.DataFrame({
        'model': model,
        'lambda': lam,
        'cv_error': np.asarray(g['yGraph']).astype(float),
        'standard_error': np.asarray(g['standard_error']).astype(float),
        'is_lambda_min': np.arange(len(lam)) == int(g.attrs['index_min']),
        'is_lambda_1se': np.arange(len(lam)) == int(g.attrs['index_1SE']),
        'lambda_min': float(g.attrs['lambda_min']),
        'lambda_1se': float(g.attrs['lambda_1SE']),
        'sigma': float(g.attrs.get('sigma', np.nan)),
    })


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    vals = sorted(set(y_true))
    return float(np.mean([np.mean(y_pred[y_true == v] == v) for v in vals]))


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
        ('R1 CV', 'atacama-top15-elevation-r1-cv.qza',
         'atacama-top15-elevation-r1-cv-train-predictions.qza',
         'atacama-top15-elevation-r1-test-predictions.qza'),
        ('R3 CV', 'atacama-top15-elevation-r3-cv.qza',
         'atacama-top15-elevation-r3-cv-train-predictions.qza',
         'atacama-top15-elevation-r3-test-predictions.qza'),
    ]
    pred_rows, metric_rows, coef_rows, cv_rows, inputs = [], [], [], [], []
    for model, problem_fn, train_fn, test_fn in specs:
        problem = open_problem(OUT / problem_fn)
        inputs.append(OUT / problem_fn)
        for split, fn in [('train', train_fn), ('test', test_fn)]:
            d = extract_prediction(OUT / fn, 'elevation', split)
            d['model'] = model
            d['observed_elevation'] = d['observed_elevation'].astype(float)
            d['predicted_elevation'] = d.pop('prediction')
            d['transect-name'] = META.loc[d['sample-id'], 'transect-name'].to_numpy()
            d['site-name'] = META.loc[d['sample-id'], 'site-name'].to_numpy()
            pred_rows.append(d)
            inputs.append(OUT / fn)
        metric = {'model': model}
        for split in ['train', 'test']:
            d = pred_rows[-2 if split == 'train' else -1]
            m = regression_metrics(d['observed_elevation'].to_numpy(float), d['predicted_elevation'].to_numpy(float))
            metric.update({f'{split}_{k}': v for k, v in m.items()})
        sol = problem['solution/CV']
        coef = solution_table(problem, 'CV', model, 'coefficient')
        metric.update({
            'lambda_min': float(sol.attrs['lambda_min']),
            'lambda_1se': float(sol.attrs['lambda_1SE']),
            'sigma': float(sol.attrs.get('sigma', np.nan)),
            'selected_total': int(coef.query("feature != 'intercept' and selected").shape[0]),
        })
        metric_rows.append(metric)
        coef_rows.append(coef)
        cv_rows.append(cv_curve(problem, model))

    stab = open_problem(OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza')
    inputs.append(OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza')
    stab_coef = solution_table(stab, 'StabSel', 'R3 stabsel q4 threshold 0.4', 'stabsel_refit_beta')
    dist = pd.DataFrame({
        'feature': np.asarray(stab['solution/StabSel/label']).astype(str),
        'stability': np.asarray(stab['solution/StabSel/distribution']).astype(float),
    })
    dist['display_label'] = dist['feature'].map(display_label)
    dist = dist.query("feature != 'intercept'").sort_values('stability', ascending=False)
    stab_pred_rows = []
    for split, fn in [('train', 'atacama-top15-elevation-r3-stabsel-threshold04-train-predictions.qza'),
                      ('test', 'atacama-top15-elevation-r3-stabsel-threshold04-test-predictions.qza')]:
        d = extract_prediction(OUT / fn, 'elevation', split)
        d['model'] = 'R3 stabsel q4 threshold 0.4'
        d['observed_elevation'] = d['observed_elevation'].astype(float)
        d['predicted_elevation'] = d.pop('prediction')
        d['transect-name'] = META.loc[d['sample-id'], 'transect-name'].to_numpy()
        d['site-name'] = META.loc[d['sample-id'], 'site-name'].to_numpy()
        stab_pred_rows.append(d)
        inputs.append(OUT / fn)
    stab_metrics = {'model': 'R3 stabsel q4 threshold 0.4'}
    for split, d in zip(['train', 'test'], stab_pred_rows):
        m = regression_metrics(d['observed_elevation'].to_numpy(float), d['predicted_elevation'].to_numpy(float))
        stab_metrics.update({f'{split}_{k}': v for k, v in m.items()})
    stab_metrics['selected_total'] = int(stab_coef.query("feature != 'intercept' and selected").shape[0])

    outputs = [
        OUT / 'atacama-top15-elevation-r1-r3-cv-curves.tsv',
        OUT / 'atacama-top15-elevation-r1-r3-cv-coefficients.tsv',
        OUT / 'atacama-top15-elevation-r1-r3-train-test-predictions.tsv',
        OUT / 'atacama-top15-elevation-r1-r3-train-test-metrics.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-coefficients.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-predictions.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-metrics.tsv',
    ]
    pd.concat(cv_rows, ignore_index=True).to_csv(outputs[0], sep='\t', index=False)
    pd.concat(coef_rows, ignore_index=True).to_csv(outputs[1], sep='\t', index=False)
    pd.concat(pred_rows, ignore_index=True).to_csv(outputs[2], sep='\t', index=False)
    pd.DataFrame(metric_rows).to_csv(outputs[3], sep='\t', index=False)
    stab_coef.to_csv(outputs[4], sep='\t', index=False)
    dist.to_csv(outputs[5], sep='\t', index=False)
    pd.concat(stab_pred_rows, ignore_index=True).to_csv(outputs[6], sep='\t', index=False)
    pd.DataFrame([stab_metrics]).to_csv(outputs[7], sep='\t', index=False)
    write_provenance(
        OUT / 'atacama-top15-elevation-qza-extraction',
        description='TSV extraction from q2-classo top-15 elevation QIIME artifacts',
        project_root=ROOT,
        script=SCRIPT,
        inputs=inputs + [OUT / 'atacama-53-sample-metadata.tsv', OUT / 'top15-asv-label-map.tsv'],
        outputs=outputs,
        parameters={'source': 'QIIME artifacts only'},
    )


def extract_classification():
    problem = open_problem(OUT / 'atacama-top15-vegetation-classification-cv.qza')
    pred_parts = []
    inputs = [OUT / 'atacama-top15-vegetation-classification-cv.qza']
    for split, fn in [('train', 'atacama-top15-vegetation-classification-cv-train-predictions.qza'),
                      ('test', 'atacama-top15-vegetation-classification-cv-test-predictions.qza')]:
        d = extract_prediction(OUT / fn, 'vegetation', split)
        d['score'] = d.pop('prediction')
        d['predicted_vegetation'] = np.where(d['score'] >= 0, 'no', 'yes')
        d['correct'] = d['predicted_vegetation'].eq(d['observed_vegetation'])
        d['elevation'] = META.loc[d['sample-id'], 'elevation'].astype(float).to_numpy()
        pred_parts.append(d)
        inputs.append(OUT / fn)
    pred = pd.concat(pred_parts, ignore_index=True)
    metrics = class_metrics(pred)
    coef = solution_table(problem, 'CV', 'C1 CV', 'refit_beta')
    cv = cv_curve(problem, 'C1 CV')
    outputs = [
        OUT / 'atacama-top15-vegetation-classification-cv-curve.tsv',
        OUT / 'atacama-top15-vegetation-classification-cv-coefficients.tsv',
        OUT / 'atacama-top15-vegetation-classification-cv-train-test-predictions.tsv',
        OUT / 'atacama-top15-vegetation-classification-cv-metrics.tsv',
    ]
    cv.to_csv(outputs[0], sep='\t', index=False)
    coef.to_csv(outputs[1], sep='\t', index=False)
    pred.to_csv(outputs[2], sep='\t', index=False)
    metrics.to_csv(outputs[3], sep='\t', index=False)
    write_provenance(
        OUT / 'atacama-top15-vegetation-classification-qza-extraction',
        description='TSV extraction from q2-classo top-15 vegetation QIIME artifacts',
        project_root=ROOT,
        script=SCRIPT,
        inputs=inputs + [OUT / 'atacama-53-sample-metadata.tsv', OUT / 'top15-asv-label-map.tsv'],
        outputs=outputs,
        parameters={'positive_score_class': 'no', 'source': 'QIIME artifacts only'},
    )


def main():
    extract_regression()
    extract_classification()
    print(json.dumps({'status': 'ok', 'output_dir': str(OUT)}, indent=2))


if __name__ == '__main__':
    main()
