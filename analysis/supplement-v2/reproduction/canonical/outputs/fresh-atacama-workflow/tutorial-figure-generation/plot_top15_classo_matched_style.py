from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

from provenance_utils import write_provenance

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
SCRIPT = Path(__file__)
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'

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
BLUE = '#356F9F'
RED = '#C43B4D'
GREEN = '#3E8A5B'
GRID = '#E6E6E6'
REF = '#404040'
MARKERS = {'Baquedano': 'o', 'Yungay': 's'}
LIMS = [850, 4850]


def style_pred(ax):
    ax.plot(LIMS, LIMS, linestyle='--', color=REF, linewidth=0.85, zorder=0)
    ax.set_xlim(LIMS)
    ax.set_ylim(LIMS)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.set_xlabel('Observed elevation (m)')
    ax.set_ylabel('Predicted elevation (m)')


def metric_label(metrics: pd.DataFrame, subset: str) -> str:
    row = metrics.loc[metrics['set'].eq(subset)].iloc[0]
    return f"n={int(row['n'])}\nr={row['pearson_r']:.2f}\nR2={row['r2']:.2f}\nRMSE={row['rmse_m']:.0f} m"


def _feature_column(df: pd.DataFrame) -> str:
    if 'feature' in df.columns:
        return 'feature'
    return 'feature_id'


def plot_elevation_three_panel(*, pred_path: Path, metrics_path: Path, coef_path: Path, output_stem: str, title: str, coefficient_column: str, coefficient_filter: str | None = None):
    pred = pd.read_csv(pred_path, sep='\t')
    metrics = pd.read_csv(metrics_path, sep='\t')
    coefs = pd.read_csv(coef_path, sep='\t')
    if coefficient_filter is not None and 'model' in coefs.columns:
        coefs = coefs[coefs['model'].eq(coefficient_filter)]
    feature_col = _feature_column(coefs)
    coefs = coefs[~coefs[feature_col].eq('intercept')].copy()
    if 'selected' in coefs.columns:
        coefs = coefs[coefs['selected'].astype(bool)]
    if 'display_label' not in coefs.columns:
        labels = pd.read_csv(OUT / 'top15-asv-label-map.tsv', sep='\t').set_index('feature_id')
        coefs['display_label'] = [labels.loc[f, 'display_label'] for f in coefs[feature_col]]
    coefs = coefs.sort_values(coefficient_column)

    height = 3.15 if len(coefs) <= 6 else 4.1
    left = 0.22 if len(coefs) <= 6 else 0.28
    fig = plt.figure(figsize=(8.7, height))
    fig.subplots_adjust(left=left, right=0.985, bottom=0.16, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)

    y = np.arange(len(coefs))
    vals = coefs[coefficient_column].to_numpy(float)
    ax0.barh(y, vals, color=[BLUE if v < 0 else RED for v in vals], alpha=0.98)
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(y)
    ax0.set_yticklabels(coefs['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.5)
    ax0.tick_params(axis='y', pad=3)
    ax0.grid(True, axis='x', color=GRID, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    for ax, subset, color, panel_title in [(ax1, 'train', BLUE, 'In-sample refit'), (ax2, 'test', RED, 'Held-out prediction')]:
        dp = pred[pred['set'].eq(subset)]
        for transect, dd in dp.groupby('transect-name'):
            ax.scatter(dd['observed_elevation'], dd['predicted_elevation'], s=22,
                       marker=MARKERS.get(transect, 'o'), facecolor=color,
                       edgecolor='white', linewidth=0.45, alpha=0.95)
        style_pred(ax)
        ax.set_title(panel_title, pad=4)
        ax.text(0.05, 0.95, metric_label(metrics, subset), transform=ax.transAxes,
                ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    handles = [
        Line2D([0], [0], marker=MARKERS['Baquedano'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Baquedano'),
        Line2D([0], [0], marker=MARKERS['Yungay'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Yungay'),
    ]
    ax2.legend(handles=handles, title='Transect', title_fontsize=6.4, loc='lower right', frameon=True)
    fig.suptitle(title, y=0.965)
    base = OUT / output_stem
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(base, description=title, project_root=ROOT, script=SCRIPT,
                     inputs=[pred_path, metrics_path, coef_path, OUT / 'top15-asv-label-map.tsv'], outputs=outputs,
                     parameters={'figure_type': 'q2-classo elevation regression three-panel', 'coefficient_column': coefficient_column, 'coefficient_filter': coefficient_filter},
                     notes=['Uses saved tables exported from q2-classo artifacts; no model refit is performed by this plotting script.'])


def plot_vegetation_cv():
    cv_path = OUT / 'atacama-top15-vegetation-classification-cv-curve.tsv'
    cv = pd.read_csv(cv_path, sep='\t')
    fig, ax = plt.subplots(figsize=(4.3, 3.05))
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.18, top=0.84)
    x = cv['lambda'].to_numpy(float)
    y = cv['cv_error'].to_numpy(float)
    se = cv['standard_error'].to_numpy(float)
    order = np.argsort(x)
    x, y, se = x[order], y[order], se[order]
    ax.plot(x, y, color=GREEN, linewidth=1.1)
    ax.fill_between(x, y - se, y + se, color=GREEN, alpha=0.17, linewidth=0)
    for flag, ls, lab in [('is_lambda_min', ':', 'min'), ('is_lambda_1se', '--', '1-SE')]:
        if flag in cv.columns and cv[flag].any():
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
    base = OUT / 'atacama-top15-vegetation-classification-cv-curve'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(base, description='q2-classo vegetation classification CV curve, top-15 Atacama ASVs', project_root=ROOT, script=SCRIPT,
                     inputs=[cv_path], outputs=outputs,
                     parameters={'figure_type': 'q2-classo classification CV curve'},
                     notes=['Uses saved CV curve table derived from q2-classo output.'])


def class_metric_label(metrics: pd.DataFrame, subset: str) -> str:
    row = metrics.loc[metrics['set'].eq(subset)].iloc[0]
    return f"n={int(row['n'])}\nacc={row['accuracy']:.2f}\nbal. acc={row['balanced_accuracy']:.2f}"


def plot_vegetation_three_panel():
    pred_path = OUT / 'atacama-top15-vegetation-classification-cv-train-test-predictions.tsv'
    metrics_path = OUT / 'atacama-top15-vegetation-classification-cv-metrics.tsv'
    coef_path = OUT / 'atacama-top15-vegetation-classification-cv-coefficients.tsv'
    pred = pd.read_csv(pred_path, sep='\t')
    metrics = pd.read_csv(metrics_path, sep='\t')
    coefs = pd.read_csv(coef_path, sep='\t')
    coefs = coefs[~coefs['feature'].eq('intercept') & coefs['selected'].astype(bool)].copy()
    coefs = coefs.sort_values('refit_beta')
    fig = plt.figure(figsize=(8.7, 4.35))
    fig.subplots_adjust(left=0.285, right=0.985, bottom=0.15, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    y = np.arange(len(coefs))
    vals = coefs['refit_beta'].to_numpy(float)
    ax0.barh(y, vals, color=[BLUE if v < 0 else RED for v in vals], alpha=0.98)
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(y)
    ax0.set_yticklabels(coefs['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.3)
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
    handles, labels_ = ax1.get_legend_handles_labels()
    ax1.legend(handles[:2], labels_[:2], loc='lower left', frameon=True, title='Observed', title_fontsize=6.4)
    ax2.legend(handles[2:], labels_[2:], loc='upper right', frameon=True)
    fig.suptitle('q2-classo CV refit, top-15 vegetation classification', y=0.965)
    base = OUT / 'atacama-top15-vegetation-classification-cv-three-panel-coefficients-in-oos'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(base, description='q2-classo CV refit, top-15 vegetation classification', project_root=ROOT, script=SCRIPT,
                     inputs=[pred_path, metrics_path, coef_path], outputs=outputs,
                     parameters={'figure_type': 'q2-classo vegetation classification three-panel'},
                     notes=['Uses saved q2-classo-derived prediction, metric, and coefficient tables.'])


def main():
    plot_elevation_three_panel(
        pred_path=OUT / 'atacama-top15-elevation-r3-cv-two-asv-train-test-predictions.tsv',
        metrics_path=OUT / 'atacama-top15-elevation-r3-cv-two-asv-train-test-metrics.tsv',
        coef_path=OUT / 'atacama-top15-elevation-r1-r3-cv-coefficients.tsv',
        output_stem='atacama-top15-elevation-r3-cv-two-asv-three-panel-coefficients-in-oos',
        title='q2-classo R3 CV refit, top-15 two-ASV elevation model',
        coefficient_column='refit_beta',
        coefficient_filter='R3',
    )
    plot_elevation_three_panel(
        pred_path=OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-predictions.tsv',
        metrics_path=OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-metrics.tsv',
        coef_path=OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-coefficients.tsv',
        output_stem='atacama-top15-r3-stabsel-threshold04-three-panel-coefficients-in-oos',
        title='q2-classo R3 stability-selection refit, top-15 four-ASV elevation model',
        coefficient_column='stabsel_refit_beta',
    )
    plot_vegetation_cv()
    plot_vegetation_three_panel()


if __name__ == '__main__':
    main()
