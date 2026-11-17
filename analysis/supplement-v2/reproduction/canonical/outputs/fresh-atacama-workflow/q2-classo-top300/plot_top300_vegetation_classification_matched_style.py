from __future__ import annotations

from pathlib import Path
import json
import sys

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
from tutorial_plotting import draw_stability_profile

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7.0,
    'axes.titlesize': 7.8,
    'axes.labelsize': 7.3,
    'xtick.labelsize': 6.4,
    'ytick.labelsize': 6.3,
    'legend.fontsize': 6.3,
    'figure.titlesize': 9.0,
    'axes.linewidth': 0.75,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

BLUE = '#356F9F'
RED = '#C43B4D'
GREEN = '#3E8A5B'
GRAY = '#9CA8B3'
GRID = '#E6E6E6'
REF = '#404040'

CV_CURVE = OUT / 'atacama-top300-vegetation-classification-cv-curve.tsv'
CV_MODEL = OUT / 'atacama-top300-vegetation-classification-cv.qza'
STAB_PROFILE = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-profile.tsv'
STAB_MODEL = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04.qza'
STAB_PRED = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-train-test-predictions.tsv'
STAB_METRICS = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-metrics.tsv'
STAB_COEF = OUT / 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-coefficients.tsv'
CV_PRED = OUT / 'atacama-top300-vegetation-classification-cv-train-test-predictions.tsv'
CV_METRICS = OUT / 'atacama-top300-vegetation-classification-cv-metrics.tsv'
CV_COEF = OUT / 'atacama-top300-vegetation-classification-cv-coefficients.tsv'


def short_asv_label(label: str) -> str:
    text = str(label)
    if ' (' in text and text.endswith(')'):
        asv, tax = text[:-1].split(' (', 1)
        tax = tax.replace(': ', ' ', 1)
        return f'{asv}\n{tax}'
    return text.replace(': ', ' ')


def class_metric_label(metrics: pd.DataFrame, subset: str) -> str:
    row = metrics.loc[metrics['set'].eq(subset)].iloc[0]
    return f"n={int(row['n'])}\nacc={row['accuracy']:.2f}\nbal. acc={row['balanced_accuracy']:.2f}"


def plot_cv_and_stabsel() -> None:
    cv = pd.read_csv(CV_CURVE, sep='\t')
    profile = pd.read_csv(STAB_PROFILE, sep='\t')
    profile = profile.sort_values('stability', ascending=False).head(20).copy()
    profile['plot_label'] = profile['display_label'].map(short_asv_label)

    fig = plt.figure(figsize=(10.2, 4.9))
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.38, top=0.83, wspace=0.24)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[0.82, 1.38])

    ax0 = fig.add_subplot(gs[0, 0])
    order = np.argsort(cv['lambda'].to_numpy(float))
    x = cv['lambda'].to_numpy(float)[order]
    y = cv['cv_error'].to_numpy(float)[order]
    se = cv['standard_error'].to_numpy(float)[order]
    ax0.plot(x, y, color=GREEN, linewidth=1.1)
    ax0.fill_between(x, y - se, y + se, color=GREEN, alpha=0.17, linewidth=0)
    for flag, ls, lab in [('is_lambda_min', ':', 'min'), ('is_lambda_1se', '--', '1-SE')]:
        if flag in cv.columns and cv[flag].astype(bool).any():
            ax0.axvline(float(cv.loc[cv[flag].astype(bool), 'lambda'].iloc[0]), color=REF, linewidth=0.85, linestyle=ls, label=lab)
    ax0.set_xscale('log')
    ax0.invert_xaxis()
    ax0.grid(True, color=GRID, linewidth=0.55)
    ax0.set_xlabel('lambda / lambda max')
    ax0.set_ylabel('CV classification error')
    ax0.set_title('Five-fold CV', pad=4)
    ax0.legend(frameon=False, loc='upper right')
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    ax1 = fig.add_subplot(gs[0, 1])
    draw_stability_profile(
        ax1,
        profile['plot_label'],
        profile['stability'],
        title='C1 stability selection',
        xlabel='Prevalence-ranked ASV label',
        label_size=5.25,
        highlight_threshold=0.2,
        annotate_thresholds=True,
    )

    fig.suptitle('q2-classo vegetation classification model selection, top-300 ASVs', y=0.965)
    base = OUT / 'atacama-top300-vegetation-classification-cv-and-stabsel-selection'
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(
        base,
        description='q2-classo top-300 vegetation classification CV and stability-selection model-selection figure',
        project_root=ROOT,
        script=SCRIPT,
        inputs=[CV_CURVE, STAB_PROFILE, CV_MODEL, STAB_MODEL],
        outputs=outputs,
        parameters={'stabsel_B': 100, 'stabsel_q': 10, 'stabsel_threshold': 0.4, 'stabsel_label_threshold': 0.2},
        notes=['Stability-selection panel uses the vertical bar style of the top-15 stability-selection profile.'],
    )


def plot_class_axis_three_panel(kind: str) -> None:
    if kind == 'stabsel':
        pred_path = STAB_PRED
        metrics_path = STAB_METRICS
        coef_path = STAB_COEF
        model_path = STAB_MODEL
        coef_col = 'stabsel_refit_beta'
        title = 'q2-classo C1 stability-selection refit, top-300 vegetation classification'
        stem = 'atacama-top300-vegetation-classification-c1-stabsel-q10-b100-threshold04-three-panel-class-axis'
        selected_title = 'Stable predictors'
    elif kind == 'cv':
        pred_path = CV_PRED
        metrics_path = CV_METRICS
        coef_path = CV_COEF
        model_path = CV_MODEL
        coef_col = 'refit_beta'
        title = 'q2-classo CV refit, top-300 vegetation classification'
        stem = 'atacama-top300-vegetation-classification-cv-three-panel-class-axis'
        selected_title = 'Selected predictors'
    else:
        raise ValueError(kind)

    pred = pd.read_csv(pred_path, sep='\t')
    metrics = pd.read_csv(metrics_path, sep='\t')
    coefs = pd.read_csv(coef_path, sep='\t')
    coefs = coefs[~coefs['feature'].eq('intercept') & coefs['selected'].astype(bool)].copy()
    coefs = coefs.sort_values(coef_col)

    height = max(3.8, 2.3 + 0.22 * len(coefs))
    fig = plt.figure(figsize=(9.1, height))
    fig.subplots_adjust(left=0.335, right=0.985, bottom=0.17, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.30, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharey=ax1)

    ypos = np.arange(len(coefs))
    vals = coefs[coef_col].to_numpy(float)
    ax0.barh(ypos, vals, color=[BLUE if v < 0 else RED for v in vals], alpha=0.98)
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(ypos)
    ax0.set_yticklabels(coefs['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.05)
    ax0.tick_params(axis='y', pad=3)
    ax0.grid(True, axis='x', color=GRID, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title(selected_title, pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    class_order = ['no', 'yes']
    class_labels = {'no': 'no vegetation', 'yes': 'vegetation'}
    class_colors = {'no': BLUE, 'yes': RED}
    rng = np.random.default_rng(7)

    class_to_x = {c: i for i, c in enumerate(class_order)}
    for ax, subset, panel_title in [(ax1, 'train', 'In-sample refit'), (ax2, 'test', 'Held-out prediction')]:
        dp = pred[pred['set'].eq(subset)].copy()
        dp['_class_x'] = dp['observed_vegetation'].map(class_to_x).astype(float)
        dp['_plot_x'] = dp['_class_x'] + rng.uniform(-0.085, 0.085, len(dp))
        for cls in class_order:
            dd = dp[dp['observed_vegetation'].eq(cls)].copy()
            if dd.empty:
                continue
            ax.scatter(dd['_plot_x'], dd['score'], s=24, marker='o', facecolor=class_colors[cls], edgecolor='white', linewidth=0.45, alpha=0.95, label=class_labels[cls])
        wrong = dp[~dp['correct'].astype(bool)]
        if len(wrong):
            ax.scatter(wrong['_plot_x'], wrong['score'], s=46, marker='x', color=REF, linewidth=0.95, label='misclassified', zorder=4)
        ax.axhline(0, color=REF, linewidth=0.85, linestyle='--')
        ax.grid(True, axis='y', color=GRID, linewidth=0.55)
        ax.set_xlim(-0.45, len(class_order) - 0.55)
        ax.set_xticks(range(len(class_order)))
        ax.set_xticklabels([class_labels[c] for c in class_order], rotation=0)
        ax.set_xlabel('Observed vegetation label')
        ax.set_ylabel('Classification score')
        ax.set_title(panel_title, pad=4)
        ax.text(0.95, 0.95, class_metric_label(metrics, subset), transform=ax.transAxes, ha='right', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')

    observed_handles = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=BLUE, markeredgecolor='white', markersize=4.7, label='no vegetation'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor=RED, markeredgecolor='white', markersize=4.7, label='vegetation'),
    ]
    wrong_handle = [Line2D([0], [0], marker='x', color=REF, linestyle='none', markersize=4.7, label='misclassified')]
    ax1.legend(handles=observed_handles, loc='lower left', frameon=True, title='Observed', title_fontsize=6.4)
    ax2.legend(handles=wrong_handle, loc='upper left', frameon=True)
    fig.suptitle(title, y=0.965)

    base = OUT / stem
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(
        base,
        description=title + ' with observed vegetation class on the x-axis',
        project_root=ROOT,
        script=SCRIPT,
        inputs=[pred_path, metrics_path, coef_path, model_path, OUT / 'top300-asv-label-map.tsv'],
        outputs=outputs,
        parameters={'figure_type': 'q2-classo classification three-panel', 'x_axis': 'observed vegetation label', 'score_threshold': 0.0, 'model_kind': kind},
        notes=['This supersedes earlier classification score plots that used elevation on the x-axis.'],
    )


def main() -> None:
    plot_cv_and_stabsel()
    plot_class_axis_three_panel('cv')
    plot_class_axis_three_panel('stabsel')


if __name__ == '__main__':
    main()
