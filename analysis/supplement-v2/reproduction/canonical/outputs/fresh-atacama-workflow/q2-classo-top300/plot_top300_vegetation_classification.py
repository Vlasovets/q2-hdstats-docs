from __future__ import annotations

from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
SCRIPT = Path(__file__)
sys.path.insert(0, str(ROOT / 'outputs/fresh-atacama-workflow/tutorial-figure-generation'))
from provenance_utils import write_provenance

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7.0,
    'axes.titlesize': 7.8,
    'axes.labelsize': 7.3,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 5.7,
    'legend.fontsize': 6.4,
    'figure.titlesize': 9.0,
    'axes.linewidth': 0.75,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
BLUE = '#356F9F'
RED = '#C43B4D'
GREEN = '#3E8A5B'
GRID = '#E6E6E6'
REF = '#404040'


def plot_cv():
    cv_path = OUT / 'atacama-top300-vegetation-classification-cv-curve.tsv'
    cv = pd.read_csv(cv_path, sep='\t')
    fig, ax = plt.subplots(figsize=(4.3, 3.05))
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.18, top=0.84)
    order = np.argsort(cv['lambda'].to_numpy(float))
    x = cv['lambda'].to_numpy(float)[order]
    y = cv['cv_error'].to_numpy(float)[order]
    se = cv['standard_error'].to_numpy(float)[order]
    ax.plot(x, y, color=GREEN, linewidth=1.1)
    ax.fill_between(x, y - se, y + se, color=GREEN, alpha=0.17, linewidth=0)
    for flag, ls, lab in [('is_lambda_min', ':', 'min'), ('is_lambda_1se', '--', '1-SE')]:
        if flag in cv.columns and cv[flag].astype(bool).any():
            ax.axvline(float(cv.loc[cv[flag].astype(bool), 'lambda'].iloc[0]), color=REF, linewidth=0.85, linestyle=ls, label=lab)
    ax.set_xscale('log')
    ax.invert_xaxis()
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.set_xlabel('lambda / lambda max')
    ax.set_ylabel('CV classification error')
    ax.set_title('q2-classo vegetation classification CV', pad=4)
    ax.legend(frameon=False, loc='upper right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    base = OUT / 'atacama-top300-vegetation-classification-cv-curve'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(base, description='q2-classo vegetation classification CV curve, top-300 Atacama ASVs', project_root=ROOT, script=SCRIPT,
                     inputs=[cv_path, OUT / 'atacama-top300-vegetation-classification-cv.qza'], outputs=outputs,
                     parameters={'figure_type': 'q2-classo classification CV curve', 'feature_set': 'top-300 prevalence/relative-abundance ASVs'},
                     notes=['Uses a TSV extracted from the q2-classo CV artifact; plotting does not refit the model.'])


def class_metric_label(metrics: pd.DataFrame, subset: str) -> str:
    row = metrics.loc[metrics['set'].eq(subset)].iloc[0]
    return f"n={int(row['n'])}\nacc={row['accuracy']:.2f}\nbal. acc={row['balanced_accuracy']:.2f}"


def plot_three_panel():
    pred_path = OUT / 'atacama-top300-vegetation-classification-cv-train-test-predictions.tsv'
    metrics_path = OUT / 'atacama-top300-vegetation-classification-cv-metrics.tsv'
    coef_path = OUT / 'atacama-top300-vegetation-classification-cv-coefficients.tsv'
    pred = pd.read_csv(pred_path, sep='\t')
    metrics = pd.read_csv(metrics_path, sep='\t')
    coefs = pd.read_csv(coef_path, sep='\t')
    coefs = coefs[~coefs['feature'].eq('intercept') & coefs['selected'].astype(bool)].copy()
    coefs = coefs.sort_values('refit_beta')

    fig = plt.figure(figsize=(9.6, 5.35))
    fig.subplots_adjust(left=0.34, right=0.985, bottom=0.13, top=0.83, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.34, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)

    y = np.arange(len(coefs))
    vals = coefs['refit_beta'].to_numpy(float)
    ax0.barh(y, vals, color=[BLUE if v < 0 else RED for v in vals], alpha=0.98)
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(y)
    ax0.set_yticklabels(coefs['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.05)
    ax0.tick_params(axis='y', pad=3)
    ax0.grid(True, axis='x', color=GRID, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    jitter = {'no': -0.09, 'yes': 0.09}
    colors = {'no': BLUE, 'yes': RED}
    labels = {'no': 'no vegetation', 'yes': 'vegetation'}
    for ax, subset, title in [(ax1, 'train', 'In-sample refit'), (ax2, 'test', 'Held-out prediction')]:
        dp = pred[pred['set'].eq(subset)].copy()
        for cls, dd in dp.groupby('observed_vegetation'):
            ax.scatter(dd['elevation'], dd['score'] + jitter.get(cls, 0.0), s=24, marker='o',
                       facecolor=colors.get(cls, '#777777'), edgecolor='white', linewidth=0.45,
                       alpha=0.95, label=labels.get(cls, cls))
        wrong = dp[~dp['correct'].astype(bool)]
        if len(wrong):
            ax.scatter(wrong['elevation'], wrong['score'], s=44, marker='x', color=REF, linewidth=0.9, label='misclassified')
        ax.axhline(0, color=REF, linewidth=0.85, linestyle='--')
        ax.grid(True, color=GRID, linewidth=0.55)
        ax.set_xlabel('Elevation (m)')
        ax.set_ylabel('Classification score')
        ax.set_title(title, pad=4)
        ax.text(0.95, 0.95, class_metric_label(metrics, subset), transform=ax.transAxes,
                ha='right', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    observed_handles = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=BLUE, markeredgecolor='white', markersize=4.7, label='no vegetation'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor=RED, markeredgecolor='white', markersize=4.7, label='vegetation'),
    ]
    wrong_handle = [Line2D([0], [0], marker='x', color=REF, linestyle='none', markersize=4.7, label='misclassified')]
    ax1.legend(handles=observed_handles, loc='lower left', frameon=True, title='Observed', title_fontsize=6.4)
    ax2.legend(handles=wrong_handle, loc='upper right', frameon=True)
    fig.suptitle('q2-classo CV refit, top-300 vegetation classification', y=0.965)

    base = OUT / 'atacama-top300-vegetation-classification-cv-three-panel-coefficients-in-oos'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(base, description='q2-classo CV refit, top-300 vegetation classification', project_root=ROOT, script=SCRIPT,
                     inputs=[pred_path, metrics_path, coef_path, OUT / 'top300-asv-label-map.tsv', OUT / 'atacama-top300-vegetation-classification-cv.qza', OUT / 'atacama-top300-vegetation-classification-cv-train-predictions.qza', OUT / 'atacama-top300-vegetation-classification-cv-test-predictions.qza'], outputs=outputs,
                     parameters={'figure_type': 'q2-classo vegetation classification three-panel', 'class_score_threshold': 0.0, 'positive_class': 'no'},
                     notes=['Uses saved q2-classo-derived prediction, metric, and coefficient tables; no model refit is performed by this plotting script.'])


def main():
    plot_cv()
    plot_three_panel()


if __name__ == '__main__':
    main()
