"""Reproduce tutorial plots from QIIME artifacts using the established plotters.

This script extracts stored solutions and predictions. It does not fit models.
Run in the tutorial QIIME environment; zarr must be version 2.
"""
import argparse
from contextlib import ExitStack
import io
from pathlib import Path
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import zarr

import extract_atacama as extract
import style_atacama as plots
from provenance_utils import write_provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', choices=['regression', 'classification'])
    parser.add_argument('--input-dir', type=Path, default=Path('.'))
    parser.add_argument('--output-dir', type=Path, default=Path('figures'))
    args = parser.parse_args()
    source = args.input_dir.resolve()
    out = args.output_dir.resolve()
    if source == out:
        raise ValueError('Choose a separate output directory to preserve the input artifacts.')
    out.mkdir(parents=True, exist_ok=True)
    if args.analysis == 'regression':
        names = ['atacama-top15-elevation-r1-cv.qza', 'atacama-top15-elevation-r3-cv.qza',
                 'atacama-top15-elevation-r3-cv-train-predictions.qza',
                 'atacama-top15-elevation-r3-test-predictions.qza',
                 'atacama-top15-elevation-r3-stabsel-q4-b100-threshold-third.qza',
                 'atacama-top15-elevation-r3-stabsel-threshold-third-train-predictions.qza',
                 'atacama-top15-elevation-r3-stabsel-threshold-third-test-predictions.qza']
    else:
        names = ['atacama-top15-vegetation-classification-cv' + suffix + '.qza'
                 for suffix in ('', '-train-predictions', '-test-predictions')]
    names += ['atacama-53-sample-metadata.tsv', 'top15-asv-label-map.tsv']
    missing = [name for name in names if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError('Run the corresponding tutorial commands first. Missing: ' + ', '.join(missing))
    for name in names:
        shutil.copy2(source / name, out / name)
    extract.ROOT = plots.ROOT = out
    extract.OUT = plots.OUT = out
    extract.META = pd.read_csv(out / 'atacama-53-sample-metadata.tsv', sep='\t', comment='#').set_index('sample-id')
    extract.LABELS = pd.read_csv(out / 'top15-asv-label-map.tsv', sep='\t').set_index('feature_id')
    if len(extract.META) != 53 or len(extract.LABELS) != 15:
        raise ValueError('Expected the 53-sample, 15-ASV tutorial annotations.')
    with ExitStack() as stack:
        def open_problem(path):
            directory = Path(stack.enter_context(tempfile.TemporaryDirectory())).resolve()
            with zipfile.ZipFile(path) as artifact:
                members = [n for n in artifact.namelist() if n.endswith('/data/problem.zip')]
                if len(members) != 1:
                    raise ValueError('Not a supported q2-classo problem artifact: ' + str(path))
                payload = artifact.read(members[0])
            with zipfile.ZipFile(io.BytesIO(payload)) as problem:
                for member in problem.infolist():
                    if not (directory / member.filename).resolve().is_relative_to(directory):
                        raise ValueError('Unsafe path in problem archive')
                problem.extractall(directory)
            return zarr.open(str(directory), mode='r')
        extract.open_problem = open_problem
        if args.analysis == 'classification':
            extract.extract_classification()
            plots.plot_vegetation()
        else:
            regression(out, open_problem)
    write_provenance(out / ('reproduce-' + args.analysis),
                     description='Portable plotting from QIIME artifacts; no model fitting',
                     project_root=out, script=Path(__file__),
                     inputs=[source / n for n in names] + list(Path(__file__).parent.glob('*.py')),
                     outputs=sorted(out.glob('*.pdf')) + sorted(out.glob('*.png')),
                     parameters={'analysis': args.analysis}, notes=['Uses established Atacama plotting functions.'])
    print('Figures and provenance written to:', out)


def regression(out, open_problem):
    curves = []
    for model, stem in [('R1 CV', 'r1-cv'), ('R3 CV', 'r3-cv')]:
        curves.append(extract.cv_curve(open_problem(out / f'atacama-top15-elevation-{stem}.qza'), model))
    pd.concat(curves).to_csv(out / 'atacama-top15-elevation-r1-r3-cv-curves.tsv', sep='\t', index=False)
    specifications = [('R3 CV', 'CV', 'atacama-top15-elevation-r3-cv.qza',
                       'atacama-top15-elevation-r3-cv-train-predictions.qza',
                       'atacama-top15-elevation-r3-test-predictions.qza', 'coefficient',
                       'atacama-top15-elevation-r3-cv-two-asv-three-panel-coefficients-in-oos',
                       'q2-classo R3 CV refit, top-15 elevation model'),
                      ('R3 stabsel q4 threshold 1/3', 'StabSel',
                       'atacama-top15-elevation-r3-stabsel-q4-b100-threshold-third.qza',
                       'atacama-top15-elevation-r3-stabsel-threshold-third-train-predictions.qza',
                       'atacama-top15-elevation-r3-stabsel-threshold-third-test-predictions.qza',
                       'stabsel_refit_beta',
                       'atacama-top15-r3-stabsel-threshold-third-three-panel-coefficients-in-oos',
                       'q2-classo R3 stability-selection refit, top-15 elevation model')]
    for model, key, artifact, train, test, column, stem, title in specifications:
        problem = open_problem(out / artifact)
        coefficients = extract.solution_table(problem, key, model, column)
        selected = coefficients.loc[coefficients.selected & coefficients.feature.ne('intercept')]
        if key == 'StabSel':
            profile = pd.DataFrame({'feature': np.asarray(problem['solution/StabSel/label']).astype(str),
                                    'stability': np.asarray(problem['solution/StabSel/distribution']).astype(float)})
            profile['display_label'] = profile.feature.map(extract.display_label)
            profile = profile.query("feature != 'intercept'").sort_values('stability', ascending=False)
            assert set(profile.loc[profile.stability > 1 / 3, 'feature']) == set(selected.feature)
            profile.to_csv(out / 'atacama-top15-elevation-r3-stabsel-q4-b100-profile.tsv', sep='\t', index=False)
        predictions = []
        metric = {'model': model, 'selected_total': len(selected)}
        for split, name in [('train', train), ('test', test)]:
            part = extract.extract_prediction(out / name, 'elevation', split)
            part['model'] = model
            part['observed_elevation'] = part.observed_elevation.astype(float)
            part['predicted_elevation'] = part.pop('prediction')
            for field in ('transect-name', 'site-name'):
                part[field] = extract.META.loc[part['sample-id'], field].to_numpy()
            metric.update({f'{split}_{k}': v for k, v in extract.regression_metrics(
                part.observed_elevation.to_numpy(), part.predicted_elevation.to_numpy()).items()})
            predictions.append(part)
        paths = [out / f'{stem}-{kind}.tsv' for kind in ('predictions', 'metrics', 'coefficients')]
        for frame, path in zip((pd.concat(predictions), pd.DataFrame([metric]), coefficients), paths):
            frame.to_csv(path, sep='\t', index=False)
        plots.plot_regression_three_panel(model, paths[0], paths[1], paths[2], column, stem, title)
    import selection_atacama as selection
    selection.figure_s3(out, selection_threshold=1 / 3, data_dir=out)


if __name__ == '__main__':
    main()
