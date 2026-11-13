from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
OUT = ROOT
SCRIPT = Path(__file__)
sys.path.insert(0, str(ROOT))
from provenance_utils import write_provenance
from tutorial_plotting import draw_stability_profile

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7.0,
    'axes.titlesize': 7.8,
    'axes.labelsize': 7.3,
    'xtick.labelsize': 6.4,
    'ytick.labelsize': 6.1,
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
MARKERS = {'Baquedano': 'o', 'Yungay': 's'}
LIMS = [850, 4850]


def save(fig, base: Path, *, description: str, inputs: list[Path], parameters: dict | None = None, notes: list[str] | None = None):
    outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)
    write_provenance(
        base,
        description=description,
        project_root=ROOT,
        script=SCRIPT,
        inputs=inputs,
        outputs=outputs,
        parameters=parameters or {},
        notes=notes or ['Input TSVs are extracted from regenerated QIIME artifacts by extract_top15_classo_qiime_results.py.'],
    )


def style_elevation(ax):
    ax.plot(LIMS, LIMS, linestyle='--', color=REF, linewidth=0.85, zorder=0)
    ax.set_xlim(LIMS)
    ax.set_ylim(LIMS)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.set_xlabel('Observed elevation (m)')
    ax.set_ylabel('Predicted elevation (m)')


def reg_metric_label(metrics: pd.DataFrame, model: str, subset: str) -> str:
    row = metrics.set_index('model').loc[model]
    return (
        f"n={int(row[f'{subset}_n'])}\n"
        f"r={row[f'{subset}_pearson_r']:.2f}\n"
        f"R2={row[f'{subset}_r2']:.2f}\n"
        f"RMSE={row[f'{subset}_rmse']:.0f} m"
    )


def short_label(x: str) -> str:
    text = str(x)
    if ' (' in text and text.endswith(')'):
        asv, tax = text[:-1].split(' (', 1)
        tax = tax.replace(': ', ' ', 1)
        return f'{asv}\n{tax}'
    return text.replace(': ', ' ')


def plot_cv_curves():
    cv = pd.read_csv(OUT / 'atacama-top15-elevation-r1-r3-cv-curves.tsv', sep='\t')
    metrics = pd.read_csv(OUT / 'atacama-top15-elevation-r1-r3-train-test-metrics.tsv', sep='\t').set_index('model')
    fig = plt.figure(figsize=(8.7, 3.1))
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.18, top=0.80, wspace=0.30)
    gs = GridSpec(1, 2, figure=fig)
    for i, model in enumerate(['R1 CV', 'R3 CV']):
        ax = fig.add_subplot(gs[0, i])
        d = cv[cv['model'].eq(model)]
        order = np.argsort(d['lambda'].to_numpy(float))
        x = d['lambda'].to_numpy(float)[order]
        y = d['cv_error'].to_numpy(float)[order]
        se = d['standard_error'].to_numpy(float)[order]
        color = BLUE if model == 'R1 CV' else RED
        ax.plot(x, y, color=color, linewidth=1.05)
        ax.fill_between(x, y - se, y + se, color=color, alpha=0.16, linewidth=0)
        for flag, ls, lab in [('is_lambda_min', ':', 'min'), ('is_lambda_1se', '--', '1-SE')]:
            if d[flag].astype(bool).any():
                ax.axvline(float(d.loc[d[flag].astype(bool), 'lambda'].iloc[0]), color=REF, linewidth=0.85, linestyle=ls, label=lab)
        ax.set_xscale('log')
        ax.invert_xaxis()
        ax.grid(True, color=GRID, linewidth=0.55)
        ax.set_title(model.replace(' CV', ''), pad=4)
        ax.set_xlabel('lambda / lambda max')
        ax.set_ylabel('CV prediction error')
        ax.text(0.05, 0.95, f"selected={int(metrics.loc[model, 'selected_total'])}\ntest R2={metrics.loc[model, 'test_r2']:.2f}",
                transform=ax.transAxes, ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
        ax.legend(frameon=False, loc='lower right')
        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
    fig.suptitle('q2-classo five-fold CV, top-15 elevation models', y=0.965)
    save(
        fig,
        OUT / 'atacama-top15-elevation-r1-r3-cv-curves',
        description='q2-classo R1/R3 five-fold CV curves for top-15 elevation regression',
        inputs=[OUT / 'atacama-top15-elevation-r1-r3-cv-curves.tsv', OUT / 'atacama-top15-elevation-r1-cv.qza', OUT / 'atacama-top15-elevation-r3-cv.qza'],
        parameters={'models': ['R1 CV', 'R3 CV']},
    )


def plot_regression_three_panel(model: str, pred_path: Path, metrics_path: Path, coef_path: Path, coef_col: str, stem: str, title: str):
    pred = pd.read_csv(pred_path, sep='\t')
    metrics = pd.read_csv(metrics_path, sep='\t')
    coefs = pd.read_csv(coef_path, sep='\t')
    coefs = coefs[coefs['model'].eq(model) & ~coefs['feature'].eq('intercept') & coefs['selected'].astype(bool)].copy()
    coefs = coefs.sort_values(coef_col)
    height = max(3.2, 2.5 + 0.22 * len(coefs))
    fig = plt.figure(figsize=(8.7, height))
    fig.subplots_adjust(left=0.25, right=0.985, bottom=0.16, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.20, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)

    ypos = np.arange(len(coefs))
    vals = coefs[coef_col].to_numpy(float)
    ax0.barh(ypos, vals, color=[BLUE if v < 0 else RED for v in vals], alpha=0.98)
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(ypos)
    ax0.set_yticklabels(coefs['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.6)
    ax0.grid(True, axis='x', color=GRID, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    for ax, subset, color, panel_title in [(ax1, 'train', BLUE, 'In-sample refit'), (ax2, 'test', RED, 'Held-out prediction')]:
        d = pred[pred['model'].eq(model) & pred['set'].eq(subset)]
        for transect, dd in d.groupby('transect-name'):
            ax.scatter(dd['observed_elevation'], dd['predicted_elevation'], s=23,
                       marker=MARKERS.get(transect, 'o'), facecolor=color,
                       edgecolor='white', linewidth=0.45, alpha=0.95)
        style_elevation(ax)
        ax.set_title(panel_title, pad=4)
        ax.text(0.05, 0.95, reg_metric_label(metrics, model, subset), transform=ax.transAxes,
                ha='left', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    handles = [
        Line2D([0], [0], marker=MARKERS['Baquedano'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Baquedano'),
        Line2D([0], [0], marker=MARKERS['Yungay'], color='none', markerfacecolor='#777777', markeredgecolor='white', markersize=4.3, label='Yungay'),
    ]
    ax2.legend(handles=handles, title='Transect', title_fontsize=6.4, loc='lower right', frameon=True)
    fig.suptitle(title, y=0.965)
    save(
        fig,
        OUT / stem,
        description=title,
        inputs=[pred_path, metrics_path, coef_path, OUT / 'top15-asv-label-map.tsv'],
        parameters={'model': model},
    )


def plot_stability_profile():
    prof = pd.read_csv(OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv', sep='\t')
    prof = prof.sort_values('stability', ascending=False)
    fig, ax = plt.subplots(figsize=(6.4, 3.25))
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.40, top=0.82)
    draw_stability_profile(
        ax,
        prof['display_label'].map(short_label),
        prof['stability'],
        title='q2-classo R3 stability selection, top-15 elevation model',
        xlabel='Prevalence-ranked ASV label',
    )
    save(
        fig,
        OUT / 'atacama-top15-elevation-r3-stabsel-profile-asv-taxlabels',
        description='q2-classo R3 stability-selection profile for top-15 elevation regression',
        inputs=[OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv', OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-threshold04.qza'],
        parameters={'threshold': 0.4, 'label_threshold': 0.3},
    )


def class_metric_label(metrics: pd.DataFrame, subset: str) -> str:
    row = metrics.loc[metrics['set'].eq(subset)].iloc[0]
    return f"n={int(row['n'])}\nacc={row['accuracy']:.2f}\nbal. acc={row['balanced_accuracy']:.2f}"


def plot_vegetation():
    cv = pd.read_csv(OUT / 'atacama-top15-vegetation-classification-cv-curve.tsv', sep='\t')
    pred = pd.read_csv(OUT / 'atacama-top15-vegetation-classification-cv-train-test-predictions.tsv', sep='\t')
    metrics = pd.read_csv(OUT / 'atacama-top15-vegetation-classification-cv-metrics.tsv', sep='\t')
    coef = pd.read_csv(OUT / 'atacama-top15-vegetation-classification-cv-coefficients.tsv', sep='\t')

    fig, ax = plt.subplots(figsize=(4.3, 3.05))
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.18, top=0.84)
    order = np.argsort(cv['lambda'].to_numpy(float))
    x, y, se = cv['lambda'].to_numpy(float)[order], cv['cv_error'].to_numpy(float)[order], cv['standard_error'].to_numpy(float)[order]
    ax.plot(x, y, color=GREEN, linewidth=1.1)
    ax.fill_between(x, y - se, y + se, color=GREEN, alpha=0.17, linewidth=0)
    for flag, ls, lab in [('is_lambda_min', ':', 'min'), ('is_lambda_1se', '--', '1-SE')]:
        if cv[flag].astype(bool).any():
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
    save(fig, OUT / 'atacama-top15-vegetation-classification-cv-curve',
         description='q2-classo top-15 vegetation classification CV curve',
         inputs=[OUT / 'atacama-top15-vegetation-classification-cv-curve.tsv', OUT / 'atacama-top15-vegetation-classification-cv.qza'])

    sel = coef.query("feature != 'intercept' and selected").sort_values('refit_beta')
    fig = plt.figure(figsize=(8.7, max(3.2, 2.5 + 0.22 * len(sel))))
    fig.subplots_adjust(left=0.28, right=0.985, bottom=0.17, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.25, 1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2], sharey=ax1)
    vals = sel['refit_beta'].to_numpy(float)
    ypos = np.arange(len(sel))
    ax0.barh(ypos, vals, color=[BLUE if v < 0 else RED for v in vals])
    ax0.axvline(0, color=REF, linewidth=0.8)
    ax0.set_yticks(ypos)
    ax0.set_yticklabels(sel['display_label'], rotation=18, ha='right', rotation_mode='anchor', fontsize=5.5)
    ax0.grid(True, axis='x', color=GRID, linewidth=0.55)
    ax0.set_xlabel('Refit coefficient')
    ax0.set_title('Selected predictors', pad=4)
    for sp in ['top', 'right']:
        ax0.spines[sp].set_visible(False)

    class_order = ['no', 'yes']
    class_labels = {'no': 'no vegetation', 'yes': 'vegetation'}
    class_colors = {'no': BLUE, 'yes': RED}
    rng = np.random.default_rng(7)
    class_to_x = {c: i for i, c in enumerate(class_order)}
    for ax, subset, panel_title in [(ax1, 'train', 'In-sample refit'), (ax2, 'test', 'Held-out prediction')]:
        d = pred[pred['set'].eq(subset)].copy()
        d['_x'] = d['observed_vegetation'].map(class_to_x).astype(float) + rng.uniform(-0.085, 0.085, len(d))
        for cls in class_order:
            dd = d[d['observed_vegetation'].eq(cls)]
            ax.scatter(dd['_x'], dd['score'], s=24, marker='o', facecolor=class_colors[cls], edgecolor='white', linewidth=0.45, alpha=0.95)
        wrong = d[~d['correct'].astype(bool)]
        if len(wrong):
            ax.scatter(wrong['_x'], wrong['score'], s=46, marker='x', color=REF, linewidth=0.95, zorder=4)
        ax.axhline(0, color=REF, linewidth=0.85, linestyle='--')
        ax.grid(True, axis='y', color=GRID, linewidth=0.55)
        ax.set_xlim(-0.45, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([class_labels[c] for c in class_order])
        ax.set_xlabel('Observed vegetation label')
        ax.set_ylabel('Classification score')
        ax.set_title(panel_title, pad=4)
        ax.text(0.95, 0.95, class_metric_label(metrics, subset), transform=ax.transAxes, ha='right', va='top', fontsize=5.9,
                bbox=dict(facecolor='white', edgecolor='#CCCCCC', linewidth=0.45, boxstyle='round,pad=0.18'))
    ax2.set_ylabel('')
    ax1.legend(handles=[
        Line2D([0], [0], marker='o', color='none', markerfacecolor=BLUE, markeredgecolor='white', markersize=4.7, label='no vegetation'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor=RED, markeredgecolor='white', markersize=4.7, label='vegetation'),
    ], loc='lower left', frameon=True, title='Observed', title_fontsize=6.4)
    ax2.legend(handles=[Line2D([0], [0], marker='x', color=REF, linestyle='none', markersize=4.7, label='misclassified')],
               loc='upper right', frameon=True)
    fig.suptitle('q2-classo CV refit, top-15 vegetation classification', y=0.965)
    save(fig, OUT / 'atacama-top15-vegetation-classification-cv-three-panel-coefficients-in-oos',
         description='q2-classo top-15 vegetation classification three-panel figure',
         inputs=[OUT / 'atacama-top15-vegetation-classification-cv-train-test-predictions.tsv',
                 OUT / 'atacama-top15-vegetation-classification-cv-metrics.tsv',
                 OUT / 'atacama-top15-vegetation-classification-cv-coefficients.tsv',
                 OUT / 'atacama-top15-vegetation-classification-cv.qza'])


def main():
    plot_cv_curves()
    plot_regression_three_panel(
        'R3 CV',
        OUT / 'atacama-top15-elevation-r1-r3-train-test-predictions.tsv',
        OUT / 'atacama-top15-elevation-r1-r3-train-test-metrics.tsv',
        OUT / 'atacama-top15-elevation-r1-r3-cv-coefficients.tsv',
        'coefficient',
        'atacama-top15-elevation-r3-cv-two-asv-three-panel-coefficients-in-oos',
        'q2-classo R3 CV refit, top-15 elevation model',
    )
    plot_stability_profile()
    plot_regression_three_panel(
        'R3 stabsel q4 threshold 0.4',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-predictions.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-train-test-metrics.tsv',
        OUT / 'atacama-top15-elevation-r3-stabsel-threshold04-coefficients.tsv',
        'stabsel_refit_beta',
        'atacama-top15-r3-stabsel-threshold04-three-panel-coefficients-in-oos',
        'q2-classo R3 stability-selection refit, top-15 elevation model',
    )
    plot_vegetation()


if __name__ == '__main__':
    main()
