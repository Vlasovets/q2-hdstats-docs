from pathlib import Path
import io
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import zarr

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
TOP15 = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'
TAXONOMY = ROOT / 'work/atacama-full/taxonomy-silva138-export/taxonomy.tsv'
SCRIPT = Path(__file__)
sys.path.insert(0, str(ROOT / 'outputs/fresh-atacama-workflow/tutorial-figure-generation'))
from provenance_utils import write_provenance

pred = pd.read_csv(OUT / 'top300-elevation-asv-vs-covariate-augmented-predictions.tsv', sep='\t')
coefs = pd.read_csv(OUT / 'top300-elevation-asv-vs-covariate-augmented-selected-coefficients.tsv', sep='\t')
if 'predictor' not in coefs.columns and 'feature' in coefs.columns:
    coefs = coefs.rename(columns={'feature': 'predictor'})
metrics = pd.read_csv(OUT / 'top300-elevation-asv-vs-covariate-augmented-metrics.tsv', sep='\t').set_index('model')
ranking = pd.read_csv(OUT / 'top300-prevalence-relative-abundance-ranking-from-53-samples.tsv', sep='\t')
taxonomy = pd.read_csv(TAXONOMY, sep='\t', index_col=0)
metadata = pd.read_csv(TOP15 / 'atacama-53-sample-metadata.tsv', sep='\t', comment='#').set_index('sample-id')
top15_labels = pd.read_csv(TOP15 / 'top15-asv-label-map.tsv', sep='\t').set_index('feature_id')

if 'transect-name' not in pred.columns or 'site-name' not in pred.columns:
    pred = pred.merge(metadata[['transect-name', 'site-name']], left_on='sample-id', right_index=True, how='left')
if 'set' not in pred.columns and 'split' in pred.columns:
    pred['set'] = pred['split']


def parse_tax(tax):
    outd = {}
    rank_names = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    for part in str(tax).split(';'):
        part = part.strip()
        if '__' not in part:
            continue
        k, v = part.split('__', 1)
        v = v.replace('_', ' ').strip()
        if v:
            outd[rank_names.get(k, k)] = v
    return outd


def informative_label_from_tax(tax):
    ranks = parse_tax(tax)
    species = ranks.get('species', '')
    genus = ranks.get('genus', '')
    generic_species = species.lower() in {'uncultured bacterium', 'uncultured archaeon', 'uncultured organism', 'metagenome'}
    if species and not generic_species:
        return f'species: {species}'
    if genus and genus.lower() != 'uncultured':
        return f'genus: {genus}'
    for r in ['family', 'order', 'class', 'phylum', 'domain']:
        if ranks.get(r) and ranks[r].lower() != 'uncultured':
            return f'{r}: {ranks[r]}'
    return 'taxonomy: unassigned'

label_rows = []
for _, row in ranking.iterrows():
    fid = row['feature_id']
    asv_label = f"ASV-{int(row['prevalence_rank'])}"
    if fid in top15_labels.index:
        tax_label = top15_labels.loc[fid, 'taxonomic_label']
    else:
        tax_label = informative_label_from_tax(taxonomy.loc[fid, 'Taxon']) if fid in taxonomy.index else 'taxonomy: unassigned'
    label_rows.append({
        'feature_id': fid,
        'asv_label': asv_label,
        'prevalence_rank': int(row['prevalence_rank']),
        'taxonomic_label': tax_label,
        'display_label': f'{asv_label} ({tax_label})',
    })
label_map = pd.DataFrame(label_rows).set_index('feature_id')
label_map.reset_index().to_csv(OUT / 'top300-asv-label-map.tsv', sep='\t', index=False)

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
grid = '#E6E6E6'
ref = '#404040'
markers = {'Baquedano': 'o', 'Yungay': 's'}
lims = [850, 4850]


def read_problem(path):
    with zipfile.ZipFile(path) as artifact:
        name = next(n for n in artifact.namelist() if n.endswith('/data/problem.zip'))
        data = artifact.read(name)
    d = tempfile.mkdtemp()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(d)
    return zarr.open(d, mode='r')


def metric_label(model, subset):
    m = metrics.loc[model]
    prefix = 'train' if subset == 'train' else 'test'
    return f"n={int(m[f'{prefix}_n'])}\nr={m[f'{prefix}_pearson_r']:.2f}\nR2={m[f'{prefix}_r2']:.2f}\nRMSE={m[f'{prefix}_rmse']:.0f} m"


def style_pred(ax):
    ax.plot(lims, lims, linestyle='--', color=ref, linewidth=0.85, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color=grid, linewidth=0.55)
    ax.set_xlabel('Observed elevation (m)')
    ax.set_ylabel('Predicted elevation (m)')


def display_predictor(row):
    predictor = row['predictor']
    cov_labels = {
        'asvscaled-log10-extract-concen': 'log10 extract conc.',
        'asvscaled-ph': 'pH',
        'asvscaled-ec': 'EC',
        'asvscaled-log10-amplicon-concentration': 'log10 amplicon conc.',
        'asvscaled-toc': 'TOC',
        'asvscaled-depth': 'depth',
    }
    if str(predictor) in cov_labels:
        return cov_labels[str(predictor)]
    return label_map.loc[predictor, 'display_label']


def plot_three_panel(model, filename, title):
    dcoef = coefs[coefs['model'].eq(model)].copy()
    dcoef['display_label'] = dcoef.apply(display_predictor, axis=1)
    dcoef = dcoef.sort_values('coefficient')
    height = 3.15 if len(dcoef) <= 6 else 4.55
    left = 0.205 if len(dcoef) <= 6 else 0.29
    fig = plt.figure(figsize=(8.7, height))
    fig.subplots_adjust(left=left, right=0.985, bottom=0.16, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)

    y = np.arange(len(dcoef))
    colors = [blue if v < 0 else red for v in dcoef['coefficient'].to_numpy(float)]
    ax0.barh(y, dcoef['coefficient'].to_numpy(float), color=colors, alpha=0.98)
    ax0.axvline(0, color=ref, linewidth=0.8)
    ax0.set_yticks(y)
    ax0.set_yticklabels(dcoef['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.5)
    ax0.tick_params(axis='y', pad=3)
    ax0.grid(True, axis='x', color=grid, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    for ax, subset, color, panel_title in [(ax1, 'train', blue, 'In-sample refit'), (ax2, 'test', red, 'Held-out prediction')]:
        dp = pred[pred['model'].eq(model) & pred['set'].eq(subset)]
        for transect, dd in dp.groupby('transect-name'):
            ax.scatter(dd['observed_elevation'], dd['predicted_elevation'], s=22,
                       marker=markers.get(transect, 'o'), facecolor=color,
                       edgecolor='white', linewidth=0.45, alpha=0.95)
        style_pred(ax)
        ax.set_title(panel_title, pad=4)
        ax.text(0.05, 0.95, metric_label(model, subset), transform=ax.transAxes,
                ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    handles = [
        Line2D([0], [0], marker=markers['Baquedano'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Baquedano'),
        Line2D([0], [0], marker=markers['Yungay'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Yungay'),
    ]
    ax2.legend(handles=handles, title='Transect', title_fontsize=6.4, loc='lower right', frameon=True)
    fig.suptitle(title, y=0.965)
    base = OUT / filename
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(
        base,
        description=title,
        project_root=ROOT,
        script=SCRIPT,
        inputs=[
            OUT / 'top300-elevation-asv-vs-covariate-augmented-predictions.tsv',
            OUT / 'top300-elevation-asv-vs-covariate-augmented-selected-coefficients.tsv',
            OUT / 'top300-elevation-asv-vs-covariate-augmented-metrics.tsv',
            OUT / 'top300-asv-label-map.tsv',
            OUT / ('atacama-top300-elevation-asv-only-r3-cv.qza' if model == 'ASV-only' else 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza'),
        ],
        outputs=outputs,
        parameters={'model': model, 'figure_type': 'q2-classo R3 elevation three-panel'},
        notes=['Input TSVs are extracted from regenerated QIIME artifacts by extract_top300_classo_qiime_results.py.'],
    )


def cv_curve_data(model):
    path = OUT / ('atacama-top300-elevation-asv-only-r3-cv.qza' if model == 'ASV-only' else 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza')
    p = read_problem(path)
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
    for i, model in enumerate(['ASV-only', 'ASV + covariates']):
        ax = fig.add_subplot(gs[0, i])
        d = cv_curve_data(model)
        order = np.argsort(d['x'])
        x, y, se = d['x'][order], d['y'][order], d['se'][order]
        color = blue if model == 'ASV-only' else red
        ax.plot(x, y, color=color, linewidth=1.05)
        ax.fill_between(x, y - se, y + se, color=color, alpha=0.16, linewidth=0)
        ax.axvline(d['lambda_min'], color=ref, linewidth=0.75, linestyle=':', label='min')
        ax.axvline(d['lambda_1se'], color=ref, linewidth=0.85, linestyle='--', label='1-SE')
        ax.set_xscale('log')
        ax.invert_xaxis()
        ax.grid(True, color=grid, linewidth=0.55)
        ax.set_title(model, pad=4)
        ax.set_xlabel('lambda / lambda max')
        ax.set_ylabel('CV prediction error')
        m = metrics.loc[model]
        ax.text(0.05, 0.95, f"selected={int(m['selected_total'])}\ntest R2={m['test_r2']:.2f}",
                transform=ax.transAxes, ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
        ax.legend(frameon=False, loc='lower right')
        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
    fig.suptitle('q2-classo R3 five-fold CV, top-300 elevation models', y=0.965)
    base = OUT / 'atacama-top300-elevation-asv-vs-covariate-r3-cv-curves-matched-style'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(
        base,
        description='q2-classo R3 five-fold CV curves, top-300 elevation models',
        project_root=ROOT,
        script=SCRIPT,
        inputs=[
            OUT / 'top300-elevation-asv-vs-covariate-augmented-metrics.tsv',
            OUT / 'atacama-top300-elevation-asv-only-r3-cv.qza',
            OUT / 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv.qza',
        ],
        outputs=outputs,
        parameters={'models': ['ASV-only R3', 'ASV + covariates R3'], 'figure_type': 'CV curves'},
        notes=['CV curves are read from regenerated q2-classo QIIME artifacts.'],
    )

plot_cv_curves()
plot_three_panel('ASV-only', 'atacama-top300-elevation-asv-only-r3-cv-three-panel-coefficients-in-oos-matched-style', 'q2-classo R3 CV refit, top-300 ASV-only elevation model')
plot_three_panel('ASV + covariates', 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv-three-panel-coefficients-in-oos-matched-style', 'q2-classo R3 CV refit, top-300 covariate-augmented elevation model')
print(OUT / 'atacama-top300-elevation-asv-vs-covariate-r3-cv-curves-matched-style.png')
print(OUT / 'atacama-top300-elevation-asv-only-r3-cv-three-panel-coefficients-in-oos-matched-style.png')
print(OUT / 'atacama-top300-elevation-covariate-augmented-equalweights-r3-cv-three-panel-coefficients-in-oos-matched-style.png')
