"""QIIME-only R3 stability selection at threshold 1/3 and matched figures.

Run with the established QIIME environment's Python. Previous threshold-0.4
artifacts and plotting scripts are retained.
"""

import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15'
FIGROOT = ROOT / 'outputs/fresh-atacama-workflow/tutorial-figure-generation'
APPENDIX = ROOT / 'outputs/_Application_Note__High_dimensional_regression_and_network_inference_in_QIIME_2_with_q2_classo_and_q2_gglasso/03_supplement'
THRESHOLD = 1 / 3
MODEL = 'R3 stabsel q4 threshold 1/3'
STEM = 'atacama-top15-elevation-r3-stabsel-threshold-third'
ARTIFACT = OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-threshold-third.qza'
sys.path[:0] = [str(OUT), str(FIGROOT)]
import extract_top15_classo_qiime_results as extract
import plot_top15_classo_qiime_figures as plots
from provenance_utils import write_provenance
from tutorial_plotting import draw_stability_profile


def main():
    commands = [[str(Path(sys.executable).with_name('qiime')), 'classo', 'regress',
        '--i-features', str(OUT / 'atacama-top15-train-53-clr.qza'),
        '--m-y-file', str(OUT / 'elevation-train-53-metadata-q2.tsv'),
        '--m-y-column', 'elevation', '--p-no-path', '--p-no-cv', '--p-stabsel',
        '--p-stabsel-seed', '1', '--p-stabsel-b', '100', '--p-stabsel-q', '4',
        '--p-stabsel-percent-ns', '0.5', '--p-stabsel-threshold', repr(THRESHOLD),
        '--p-stabsel-threshold-label', '0.3', '--p-no-lamfixed',
        '--p-concomitant', '--p-no-huber', '--p-intercept', '--o-result', str(ARTIFACT)]]
    for split in ('train', 'test'):
        commands.append([commands[0][0], 'classo', 'predict', '--i-features',
            str(OUT / f'atacama-top15-{split}-53-clr.qza'), '--i-problem', str(ARTIFACT),
            '--o-predictions', str(OUT / f'{STEM}-{split}-predictions.qza')])
    command_file = OUT / f'{STEM}-commands.sh'
    command_file.write_text('#!/usr/bin/env bash\nset -euo pipefail\n' + '\n'.join(shlex.join(c) for c in commands) + '\n')
    for command in commands:
        subprocess.run(command, check=True)

    problem = extract.open_problem(ARTIFACT)
    coefs = extract.solution_table(problem, 'StabSel', MODEL, 'stabsel_refit_beta')
    profile = pd.DataFrame({'feature': np.asarray(problem['solution/StabSel/label']).astype(str),
        'stability': np.asarray(problem['solution/StabSel/distribution']).astype(float)})
    profile['display_label'] = profile.feature.map(extract.display_label)
    profile = profile.query("feature != 'intercept'").sort_values('stability', ascending=False)
    # q and seed are fixed, so the threshold must not change subsampling scores.
    previous = pd.read_csv(OUT / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv', sep='\t')
    np.testing.assert_allclose(profile.set_index('feature').stability.sort_index(), previous.set_index('feature').stability.sort_index())
    expected = set(profile.loc[profile.stability > THRESHOLD, 'feature'])
    actual = set(coefs.loc[coefs.selected & coefs.feature.ne('intercept'), 'feature'])
    assert expected == actual, (expected, actual)
    parts, metrics = [], {'model': MODEL, 'selected_total': len(actual)}
    for split in ('train', 'test'):
        part = extract.extract_prediction(OUT / f'{STEM}-{split}-predictions.qza', 'elevation', split)
        part['model'] = MODEL
        part['observed_elevation'] = part.observed_elevation.astype(float)
        part['predicted_elevation'] = part.pop('prediction')
        for col in ('transect-name', 'site-name'):
            part[col] = extract.META.loc[part['sample-id'], col].to_numpy()
        metrics.update({f'{split}_{k}': v for k, v in extract.regression_metrics(part.observed_elevation.to_numpy(), part.predicted_elevation.to_numpy()).items()})
        parts.append(part)
    paths = [OUT / f'{STEM}-{suffix}.tsv' for suffix in ('coefficients', 'profile', 'train-test-predictions', 'train-test-metrics')]
    for frame, path in zip((coefs, profile, pd.concat(parts), pd.DataFrame([metrics])), paths):
        frame.to_csv(path, sep='\t', index=False)
    inputs = [ARTIFACT, command_file, OUT / f'{STEM}-train-predictions.qza', OUT / f'{STEM}-test-predictions.qza']
    write_provenance(OUT / f'{STEM}-extraction', description='QIIME-derived stability-selection refit at threshold 1/3', project_root=ROOT,
        script=Path(__file__), inputs=inputs, outputs=paths, parameters={'threshold': THRESHOLD, 'q': 4, 'B': 100, 'seed': 1})

    refit_stem = 'atacama-top15-r3-stabsel-threshold-third-three-panel-coefficients-in-oos'
    plots.plot_regression_three_panel(MODEL, paths[2], paths[3], paths[0], 'stabsel_refit_beta', refit_stem,
        'q2-classo R3 stability-selection refit, top-15 elevation model')
    fig, ax = plots.plt.subplots(figsize=(6.4, 3.25))
    fig.subplots_adjust(left=.10, right=.985, bottom=.40, top=.82)
    draw_stability_profile(ax, profile.display_label.map(plots.short_label), profile.stability,
        title='q2-classo R3 stability selection, top-15 elevation model', xlabel='Prevalence-ranked ASV label', selection_threshold=THRESHOLD)
    plots.save(fig, OUT / f'{STEM}-profile-figure', description='R3 stability profile with 1/3 cutoff', inputs=[paths[1], ARTIFACT], parameters={'threshold': THRESHOLD})
    spec = importlib.util.spec_from_file_location('replot_first_three', FIGROOT / 'replot_first_three.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.figure_s3(APPENDIX / 'figs', selection_threshold=THRESHOLD)
    shutil.copy2(OUT / f'{refit_stem}.pdf', APPENDIX / 'figs/atacama_top15_classo_stabsel_refit.pdf')
    for name in ('atacama_top15_classo_selection_large_type', 'atacama_top15_classo_stabsel_refit'):
        target = APPENDIX / 'figs' / name
        write_provenance(target, description='Top-15 R3 stability selection at threshold 1/3 using established plotters', project_root=ROOT,
            script=Path(__file__), inputs=inputs + paths + [FIGROOT / 'replot_first_three.py', OUT / 'plot_top15_classo_qiime_figures.py'],
            outputs=[target.with_suffix('.pdf')], parameters={'threshold': THRESHOLD, 'q': 4, 'B': 100, 'backend_commit': 'd22da8d5ea54eb499de2fce82fa6c0e1f0c2dcc0'})
    print(json.dumps(metrics, indent=2))
    print(coefs.loc[coefs.selected, ['display_label', 'stabsel_refit_beta']].to_string(index=False))


if __name__ == '__main__':
    main()
