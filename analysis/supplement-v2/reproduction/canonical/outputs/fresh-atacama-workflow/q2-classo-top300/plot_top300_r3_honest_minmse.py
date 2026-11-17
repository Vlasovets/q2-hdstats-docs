from pathlib import Path
import io
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import zarr
from scipy.stats import pearsonr
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
TOP15 = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'

metadata = pd.read_csv(TOP15 / 'atacama-53-sample-metadata.tsv', sep='\t', comment='#').set_index('sample-id')
label_map = pd.read_csv(OUT / 'top300-asv-label-map.tsv', sep='\t').set_index('feature_id')
one_se_metrics = pd.read_csv(OUT / 'top300-elevation-r3-asv-vs-covariate-augmented-metrics.tsv', sep='\t').set_index('model')
min_metrics = pd.read_csv(OUT / '_diagnostic-top300-elevation-covariate-augmented-r3-cv-min-metrics.tsv', sep='\t').set_index('split')
min_coefs = pd.read_csv(OUT / '_diagnostic-top300-elevation-covariate-augmented-r3-cv-min-selected-coefficients.tsv', sep='\t')

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7.0,
    'axes.titlesize': 7.8,
    'axes.labelsize': 7.3,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.2,
    'legend.fontsize': 6.4,
    'figure.titlesize': 9.0,
    'axes.linewidth': 0.75,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
blue = '#356F9F'
red = '#C43B4D'
green = '#3E8A5B'
grid = '#E6E6E6'
ref = '#404040'
markers = {'Baquedano': 'o', 'Yungay': 's'}
lims = [850, 4850]


def open_problem(path):
    with zipfile.ZipFile(path) as artifact:
        name = next(n for n in artifact.namelist() if n.endswith('/data/problem.zip'))
        data = artifact.read(name)
    d = tempfile.mkdtemp()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(d)
    return zarr.open(d, mode='r')


def as_text(a):
    x = np.asarray(a)
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in x])


def prediction_table(qza, split, model):
    p = open_problem(qza)
    df = pd.DataFrame({
        'sample-id': as_text(p['sample_labels']),
        'predicted_elevation': np.asarray(p['YhatCV'], dtype=float),
    })
    df['split'] = split
    df['model'] = model
    df['observed_elevation'] = metadata.loc[df['sample-id'], 'elevation'].astype(float).to_numpy()
    df['transect-name'] = metadata.loc[df['sample-id'], 'transect-name'].to_numpy()
    df['site-name'] = metadata.loc[df['sample-id'], 'site-name'].to_numpy()
    return df[['model', 'sample-id', 'split', 'observed_elevation', 'predicted_elevation', 'transect-name', 'site-name']]


def metric(y, yp):
    y = np.asarray(y, dtype=float)
    yp = np.asarray(yp, dtype=float)
    return {
        'n': len(y),
        'pearson_r': float(pearsonr(y, yp).statistic),
        'r2': float(r2_score(y, yp)),
        'rmse': float(np.sqrt(mean_squared_error(y, yp))),
        'mae': float(mean_absolute_error(y, yp)),
    }


min_pred = pd.concat([
    prediction_table(OUT / '_diagnostic-atacama-top300-elevation-covariate-augmented-r3-cv-min-train-predictions.qza', 'train', 'ASV + covariates, min-MSE'),
    prediction_table(OUT / '_diagnostic-atacama-top300-elevation-covariate-augmented-r3-cv-min-test-predictions.qza', 'test', 'ASV + covariates, min-MSE'),
], ignore_index=True)
min_pred.to_csv(OUT / 'top300-elevation-r3-covariate-augmented-minmse-predictions.tsv', sep='\t', index=False)

summary_rows = []
for split, dd in min_pred.groupby('split'):
    row = metric(dd['observed_elevation'], dd['predicted_elevation'])
    row['split'] = split
    summary_rows.append(row)
pd.DataFrame(summary_rows).to_csv(OUT / 'top300-elevation-r3-covariate-augmented-minmse-metrics.tsv', sep='\t', index=False)


def cv_curve_data(path):
    p = open_problem(path)
    cv = p['solution/CV']
    return {
        'x': np.asarray(cv['xGraph'], dtype=float),
        'y': np.asarray(cv['yGraph'], dtype=float),
        'se': np.asarray(cv['standard_error'], dtype=float),
        'lambda_1se': float(cv.attrs['lambda_1SE']),
        'lambda_min': float(cv.attrs['lambda_min']),
    }


def plot_cv_curves():
    fig = plt.figure(figsize=(8.7, 3.15))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.18, top=0.80, wspace=0.30)
    gs = GridSpec(1, 2, figure=fig)
    specs = [
        ('ASV-only', OUT / 'atacama-top300-elevation-asv-only-r3-cv.qza', blue),
        ('ASV + covariates', OUT / 'atacama-top300-elevation-covariate-augmented-r3-cv.qza', red),
    ]
    for i, (model, path, color) in enumerate(specs):
        ax = fig.add_subplot(gs[0, i])
        d = cv_curve_data(path)
        order = np.argsort(d['x'])
        x, y, se = d['x'][order], d['y'][order], d['se'][order]
        ax.plot(x, y, color=color, linewidth=1.05)
        ax.fill_between(x, y - se, y + se, color=color, alpha=0.16, linewidth=0)
        ax.axvline(d['lambda_min'], color=ref, linewidth=0.8, linestyle=':', label='min-MSE')
        ax.axvline(d['lambda_1se'], color=ref, linewidth=0.9, linestyle='--', label='1-SE')
        ax.set_xscale('log')
        ax.invert_xaxis()
        ax.grid(True, color=grid, linewidth=0.55)
        ax.set_title(model, pad=4)
        ax.set_xlabel('lambda / lambda max')
        ax.set_ylabel('CV prediction error')
        if model == 'ASV-only':
            m = one_se_metrics.loc[model]
            txt = f"1-SE selected={int(m['selected_total'])}\ntest R2={m['test_r2']:.2f}"
        else:
            m1 = one_se_metrics.loc[model]
            txt = f"1-SE selected={int(m1['selected_total'])}\n1-SE test R2={m1['test_r2']:.2f}\nmin-MSE selected=11\nmin-MSE test R2={min_metrics.loc['test','r2']:.2f}"
        ax.text(0.05, 0.95, txt, transform=ax.transAxes, ha='left', va='top', fontsize=5.8,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
        ax.legend(frameon=False, loc='lower right')
        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
    fig.suptitle('q2-classo R3 five-fold CV, top-300 elevation models', y=0.965)
    base = OUT / 'atacama-top300-elevation-r3-cv-curves-one-se-and-minmse'
    fig.savefig(base.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(base.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)


def style_pred(ax):
    ax.plot(lims, lims, linestyle='--', color=ref, linewidth=0.85, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color=grid, linewidth=0.55)
    ax.set_xlabel('Observed elevation (m)')
    ax.set_ylabel('Predicted elevation (m)')


def metric_label(split):
    row = min_metrics.loc[split]
    return f"n={int(row['n'])}\nr={row['r']:.2f}\nR2={row['r2']:.2f}\nRMSE={row['rmse']:.0f} m"


def display_predictor(row):
    pred = str(row['predictor'])
    if pred.startswith('asvscaled-'):
        raw = pred.replace('asvscaled-', '')
        return {
            'log10-extract-concen': 'log10 extract conc.',
            'ph': 'pH',
            'ec': 'EC',
            'log10-amplicon-concentration': 'log10 amplicon conc.',
            'toc': 'TOC',
            'depth': 'depth',
        }.get(raw, raw.replace('-', ' '))
    if pred in label_map.index:
        return label_map.loc[pred, 'display_label']
    return row['display_label']


def plot_minmse_three_panel():
    dcoef = min_coefs.copy()
    dcoef['display_label'] = dcoef.apply(display_predictor, axis=1)
    dcoef = dcoef.sort_values('coefficient')
    fig = plt.figure(figsize=(8.7, 4.75))
    fig.subplots_adjust(left=0.31, right=0.985, bottom=0.15, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    y = np.arange(len(dcoef))
    vals = dcoef['coefficient'].to_numpy(float)
    colors = [blue if v < 0 else red for v in vals]
    ax0.barh(y, vals, color=colors, alpha=0.98)
    ax0.axvline(0, color=ref, linewidth=0.8)
    ax0.set_yticks(y)
    ax0.set_yticklabels(dcoef['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.2)
    ax0.tick_params(axis='y', pad=3)
    ax0.grid(True, axis='x', color=grid, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)
    for ax, split, color, title in [(ax1, 'train', blue, 'In-sample refit'), (ax2, 'test', red, 'Held-out prediction')]:
        dp = min_pred[min_pred['split'].eq(split)]
        for transect, dd in dp.groupby('transect-name'):
            ax.scatter(dd['observed_elevation'], dd['predicted_elevation'], s=22,
                       marker=markers.get(transect, 'o'), facecolor=color,
                       edgecolor='white', linewidth=0.45, alpha=0.95)
        style_pred(ax)
        ax.set_title(title, pad=4)
        ax.text(0.05, 0.95, metric_label(split), transform=ax.transAxes,
                ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    handles = [
        Line2D([0], [0], marker=markers['Baquedano'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Baquedano'),
        Line2D([0], [0], marker=markers['Yungay'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Yungay'),
    ]
    ax2.legend(handles=handles, title='Transect', title_fontsize=6.4, loc='lower right', frameon=True)
    fig.suptitle('q2-classo R3 min-MSE refit, top-300 covariate-augmented elevation model', y=0.965)
    base = OUT / 'atacama-top300-elevation-r3-covariate-augmented-minmse-three-panel'
    fig.savefig(base.with_suffix('.png'), dpi=300, bbox_inches='tight')
    fig.savefig(base.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)


plot_cv_curves()
plot_minmse_three_panel()
print(OUT / 'atacama-top300-elevation-r3-cv-curves-one-se-and-minmse.png')
print(OUT / 'atacama-top300-elevation-r3-covariate-augmented-minmse-three-panel.png')
