"""Check packaged source inputs and unchanged preparation in a clean directory."""
from pathlib import Path
import io
import json
import subprocess
import sys
import tempfile
import zipfile
import biom
import numpy as np
import pandas as pd
import zarr

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'outputs/standalone-tutorial'
checks = []


def table(path):
    with zipfile.ZipFile(path) as archive:
        name = next(n for n in archive.namelist() if n.endswith('/data/feature-table.biom'))
        with tempfile.NamedTemporaryFile(suffix='.biom') as file:
            file.write(archive.read(name))
            file.flush()
            value = biom.load_table(file.name)
    return pd.DataFrame(value.matrix_data.toarray(), index=value.ids(axis='observation'), columns=value.ids(axis='sample'))


with tempfile.TemporaryDirectory(prefix='q2-standalone-') as temp:
    folder = Path(temp)
    for bundle in ('atacama-top300-inputs', 'motus-inputs'):
        dest = folder / bundle
        with zipfile.ZipFile(ROOT / (bundle + '.zip')) as archive:
            archive.extractall(dest)
        if bundle.startswith('atacama'):
            before = {name: pd.read_csv(dest / name, sep='\t') for name in (
                'top300-training-supported-feature-ids.tsv',
                'top300-elevation-augmented-covariates-normalized-from-train.tsv',
                'top300-covariate-normalization-summary.tsv',
                'top300-prevalence-relative-abundance-ranking-from-53-samples.tsv')}
            for script in ('prepare_top300_inputs.py', 'create_classo_constraint_weights.py'):
                subprocess.run([sys.executable, script], cwd=dest, check=True, stdout=subprocess.DEVNULL)
            for name, frame in before.items():
                pd.testing.assert_frame_equal(frame, pd.read_csv(dest / name, sep='\t'), check_exact=False, rtol=1e-12, atol=1e-12)
                checks.append({'check': 'portable preparation matches audited ' + name, 'passed': True})
            assert np.asarray(zarr.load(str(dest / 'top300-zero-sum-constraint.zip'))).shape == (1, 285)
            assert np.asarray(zarr.load(str(dest / 'top300-unit-weights.zip'))).shape == (285,)
            matrix = table(dest / 'atacama-top300-prevalence-53-table.qza')
            assert matrix.shape == (300, 53) and (matrix.sum(axis=1) > 0).all()
            checks.append({'check': '300 network features, 285 prediction weights and constraint; no globally empty features', 'passed': True})
        else:
            train, test, full = [table(dest / name) for name in ('motus-train.qza', 'motus-test.qza', 'motus-prev10pct.qza')]
            assert train.shape == (60, 26) and test.shape == (60, 8) and full.shape == (60, 34)
            assert list(train.index) == list(test.index)
            assert (train.sum(axis=1) > 0).all()
            pd.testing.assert_frame_equal(full[train.columns], train)
            pd.testing.assert_frame_equal(full[test.columns], test)
            metadata = pd.read_csv(dest / 'outcomes.tsv', sep='\t', index_col=0)
            assert not set(metadata.loc[train.columns, 'host_subject_id']) & set(metadata.loc[test.columns, 'host_subject_id'])
            checks.append({'check': 'mOTU 34x60; 26/8 subject-disjoint split; training-supported features and full-table agreement', 'passed': True})
(ROOT / 'input-validation.json').write_text(json.dumps(checks, indent=2) + '\n')
print(f'PASS: {len(checks)} packaged-input and clean-directory preparation checks; no models refitted.')
