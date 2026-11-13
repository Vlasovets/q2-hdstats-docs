"""Replot supplementary Figures S1-S3 from the archived QIIME-derived TSVs."""

from pathlib import Path
import argparse
import os

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(os.environ['ATACAMA_PROJECT_ROOT']) if 'ATACAMA_PROJECT_ROOT' in os.environ else Path(__file__).resolve().parents[3]
DATA = ROOT / 'outputs/fresh-atacama-workflow'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5,
    'legend.fontsize': 8.5,
    'pdf.fonttype': 42,
    'axes.linewidth': 0.7,
})


def read_tsv(path, index='sample-id'):
    return pd.read_csv(path, sep='\t').set_index(index)


def style_heatmap(ax, arr, *, cmap, vmin=None, vmax=None):
    im = ax.imshow(arr, aspect='auto', interpolation='nearest', cmap=cmap,
                   vmin=vmin, vmax=vmax, rasterized=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(length=0)
    return im


def figure_s1(out):
    d = DATA / 'metadata'
    meta = read_tsv(d / 'atacama-54-transformed-covariate-zscores.tsv')
    missing = read_tsv(d / 'atacama-54-covariate-missingness-mask.tsv')
    corr = pd.read_csv(d / 'atacama-54-transformed-covariate-pearson-correlation-clustered.tsv',
                       sep='\t', index_col=0)
    labels = pd.read_csv(d / 'atacama-54-aligned-covariate-transform-and-display-scaling.tsv', sep='\t')
    summary = pd.read_csv(d / 'atacama-54-numeric-covariate-summary.tsv', sep='\t').set_index('column')
    cols = corr.columns.tolist()
    names = labels.set_index('covariate').loc[cols, 'label'].tolist()
    meta = meta.loc[:, cols]
    missing = missing.loc[meta.index, cols]
    baq = np.array([str(x).startswith('BAQ') for x in meta.index])
    assert baq.sum() + (~baq).sum() == 54

    fig = plt.figure(figsize=(14.2, 6.4))
    gs = fig.add_gridspec(1, 3, left=0.085, right=0.95, bottom=0.37,
                          top=0.82, width_ratios=[6.8, 1.8, 4.1], wspace=0.34)
    ax = fig.add_subplot(gs[0, 0])
    im = style_heatmap(ax, meta.T.to_numpy(), cmap='RdBu_r', vmin=-2.8, vmax=2.8)
    ax.set_yticks(np.arange(len(cols)), names)
    ticks = np.arange(len(meta))
    ax.set_xticks(ticks, meta.index, rotation=86, ha='right',
                  rotation_mode='anchor', fontsize=8.8)
    ax.set_title('Numeric metadata by sample', pad=16)
    ax.set_xlabel('Samples ordered by transect and elevation')
    ax.axvline(baq.sum() - 0.5, color='black', linewidth=1.1)
    for j, col in enumerate(cols):
        for i in np.flatnonzero(missing[col].to_numpy(dtype=bool)):
            ax.add_patch(Rectangle((i - .5, j - .5), 1, 1, fill=False,
                                   edgecolor='black', linewidth=.65))
    strip = ax.inset_axes([0, 1.025, 1, .055])
    strip.imshow(baq[None, :], cmap=ListedColormap(['#4F79A4', '#DC8A3A']),
                 aspect='auto', interpolation='nearest')
    strip.set_axis_off()
    ax.text(.25, 1.105, 'Baquedano', ha='center', transform=ax.transAxes,
            color='#A75313', fontsize=9)
    ax.text(.75, 1.105, 'Yungay', ha='center', transform=ax.transAxes,
            color='#326596', fontsize=9)

    stat = fig.add_subplot(gs[0, 1])
    stat.set_xlim(0, 3)
    stat.set_ylim(ax.get_ylim())
    stat.set_xticks([.5, 1.5, 2.5], ['min', 'median', 'max'])
    stat.xaxis.tick_top()
    stat.tick_params(axis='x', labelsize=8.2, pad=5, length=0)
    stat.set_yticks([])
    for j, col in enumerate(cols):
        row = summary.loc[col]
        for k, field in enumerate(('min', 'median', 'max')):
            value = row[field]
            string = f'{value:.2f}' if abs(value) < 100 else f'{value:.0f}'
            stat.text(k + .5, j, string, ha='center', va='center', fontsize=7.5)
        stat.axhline(j + .5, color='#DDDDDD', linewidth=.35)
    for spine in stat.spines.values():
        spine.set_visible(False)

    axr = fig.add_subplot(gs[0, 2])
    imr = style_heatmap(axr, corr.loc[cols, cols].to_numpy(), cmap='RdBu_r', vmin=-1, vmax=1)
    axr.set_aspect('equal', adjustable='box')
    axr.set_xticks(np.arange(len(cols)), names, rotation=65, ha='right', rotation_mode='anchor')
    axr.set_yticks(np.arange(len(cols)), names)
    axr.set_title('Pearson correlation', pad=12)
    cb = fig.colorbar(imr, ax=axr, fraction=.045, pad=.02)
    cb.set_label('Pearson r')
    cax = fig.add_axes([.085, .08, .32, .018])
    cb2 = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb2.set_label('Within-covariate z-score')
    fig.savefig(out / 'atacama_covariates_overview.pdf', bbox_inches='tight')
    plt.close(fig)


def figure_s2(out):
    d = DATA / 'q2-classo-top15'
    clr = read_tsv(d / 'atacama-top15-53-samples-q2-classo-clr-coef05-matrix.tsv')
    mclr = read_tsv(d / 'atacama-top15-53-samples-q2-gglasso-mclr-matrix.tsv')
    assert clr.shape == mclr.shape == (53, 15)
    mclr = mclr.loc[clr.index, clr.columns]
    fig, axes = plt.subplots(1, 2, figsize=(8.7, 8.8), sharey=True)
    fig.subplots_adjust(left=.16, right=.96, bottom=.14, top=.94, wspace=.10)
    for ax, table, cmap, limits, title in [
        (axes[0], clr, 'RdBu_r', (-6, 6), 'q2-classo CLR (coef = 0.5)'),
        (axes[1], mclr, 'viridis', (0, 2.7), 'q2-gglasso mCLR'),
    ]:
        im = style_heatmap(ax, table.to_numpy(), cmap=cmap, vmin=limits[0], vmax=limits[1])
        ax.set_title(title, pad=10)
        ax.set_xticks(np.arange(15), [f'ASV-{i}' for i in range(1, 16)],
                      rotation=55, ha='right', rotation_mode='anchor')
        ax.set_xlabel('Prevalence-ranked ASV')
        ax.axhline(sum(str(x).startswith('BAQ') for x in table.index) - .5,
                   color='white' if ax is axes[1] else '#222222', linewidth=1)
        cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.05, pad=.10)
        cb.set_label('CLR value' if ax is axes[0] else 'mCLR value')
    ticks = np.arange(len(clr))
    axes[0].set_yticks(ticks, clr.index[ticks])
    axes[0].tick_params(axis='y', labelsize=7.8, pad=3)
    axes[0].set_ylabel('Samples ordered by transect and elevation')
    axes[1].tick_params(axis='y', labelleft=False)
    fig.savefig(out / 'atacama_top15_clr_mclr_heatmaps.pdf', bbox_inches='tight')
    plt.close(fig)


def figure_s3(out, selection_threshold=0.4):
    d = DATA / 'q2-classo-top15'
    cv = pd.read_csv(d / 'atacama-top15-elevation-r1-r3-cv-curves.tsv', sep='\t')
    prof = pd.read_csv(d / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv', sep='\t')
    fig = plt.figure(figsize=(10.4, 3.6))
    gs = fig.add_gridspec(1, 3, left=.075, right=.985, top=.82, bottom=.37,
                          width_ratios=[1, 1, 1.18], wspace=.37)
    for k, model in enumerate(('R1 CV', 'R3 CV')):
        ax = fig.add_subplot(gs[0, k])
        dcv = cv.loc[cv.model.eq(model)].sort_values('lambda')
        x = dcv['lambda'].to_numpy(float)
        y = dcv['cv_error'].to_numpy(float) / 1e6
        se = dcv['standard_error'].to_numpy(float) / 1e6
        color = '#356F9F' if k == 0 else '#C43B4D'
        ax.plot(x, y, color=color, linewidth=1.5)
        ax.fill_between(x, y - se, y + se, color=color, alpha=.16)
        for flag, line, label in [('is_lambda_min', ':', 'Min'), ('is_lambda_1se', '--', '1SE')]:
            selected = dcv.loc[dcv[flag].astype(bool), 'lambda']
            if not selected.empty:
                ax.axvline(selected.iloc[0], color='#444444', linestyle=line,
                           linewidth=1, label=label)
        ax.set_xscale('log')
        ax.invert_xaxis()
        ax.set_title(model.replace(' CV', ''))
        ax.set_xlabel(r'$\lambda/\lambda_{\max}$')
        if k == 0:
            ax.set_ylabel('CV error (million m$^2$)')
        ax.grid(color='#E6E6E6', linewidth=.5)
        ax.legend(frameon=False, loc='lower right', ncol=2)
        ax.spines[['top', 'right']].set_visible(False)
    ax = fig.add_subplot(gs[0, 2])
    prof = prof.sort_values('stability', ascending=False)
    short = prof['display_label'].str.extract(r'^(ASV-\d+)')[0].tolist()
    vals = prof['stability'].to_numpy(float)
    ax.bar(np.arange(len(vals)), vals, color=['#C43B4D' if v > selection_threshold else '#356F9F' for v in vals])
    threshold_label = '1/3' if np.isclose(selection_threshold, 1 / 3) else f'{selection_threshold:g}'
    ax.axhline(selection_threshold, color='#777777', linestyle=':', linewidth=1, label=f'{threshold_label} cutoff')
    ax.axhline(.7, color='#C43B4D', linestyle='--', linewidth=1, label='0.7 reference')
    ax.set_ylim(0, 1)
    ax.set_ylabel('Selection probability')
    ax.set_xlabel('ASV')
    ax.set_xticks(np.arange(len(vals)), short, rotation=70, ha='right')
    ax.set_title('R3 stability selection')
    ax.legend(frameon=False, fontsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(out / 'atacama_top15_classo_selection_large_type.pdf', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_s1(args.output_dir)
    figure_s2(args.output_dir)
    figure_s3(args.output_dir)
