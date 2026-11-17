"""Regression tests for shortening QIIME examples without changing settings."""
from pathlib import Path
import json
import shlex
import unittest

from minimize_qiime_commands import resolved, simplify


SNAPSHOT = json.loads((Path(__file__).resolve().parents[1] /
                      'supplement-v2/command-harmonization/qiime-interface-defaults.json').read_text())


class MinimalCommands(unittest.TestCase):
    def check_command(self, text, removed, retained):
        new, report = simplify(text, SNAPSHOT)
        self.assertTrue(report['effective_settings_equal'])
        tokens = shlex.split(new.replace('\\\n', ' '))
        self.assertEqual(resolved(shlex.split(text), SNAPSHOT), resolved(tokens, SNAPSHOT))
        for flag in removed:
            self.assertNotIn(flag, tokens)
        for flag in retained:
            self.assertIn(flag, tokens)

    def test_r3_cv(self):
        self.check_command('qiime classo regress --p-concomitant --p-no-huber --p-intercept --p-cv --p-cv-one-se --p-cv-subsets 5 --p-cv-seed 1 --p-no-path --p-no-stabsel --p-no-lamfixed',
                           ['--p-concomitant', '--p-no-huber', '--p-intercept', '--p-cv', '--p-cv-one-se', '--p-cv-subsets', '--p-cv-seed'],
                           ['--p-no-path', '--p-no-stabsel', '--p-no-lamfixed'])

    def test_r1(self):
        self.check_command('qiime classo regress --p-no-concomitant --p-no-huber', ['--p-no-huber'], ['--p-no-concomitant'])

    def test_c2_no_refit(self):
        self.check_command('qiime classo classify --p-huber --p-rho 0 --p-no-refit', ['--p-rho'], ['--p-huber', '--p-no-refit'])

    def test_nondefault_grid(self):
        self.check_command('qiime classo regress --p-cv--nlam 120 --p-cv-lamin 0.0001 --p-cv-seed 1', ['--p-cv-seed'], ['--p-cv--nlam', '--p-cv-lamin'])

    def test_stability_seed_not_default(self):
        self.check_command('qiime classo regress --p-stabsel --p-stabsel-seed 1 --p-stabsel-q 4 --p-stabsel-percent-ns 0.5',
                           ['--p-stabsel', '--p-stabsel-percent-ns'], ['--p-stabsel-seed', '--p-stabsel-q'])

    def test_transform_and_covariance(self):
        self.check_command('qiime classo transform-features --p-transformation clr --p-coef 0.5', ['--p-transformation', '--p-coef'], [])
        self.check_command('qiime gglasso transform-features --p-transformation mclr --p-pseudo-count 1', ['--p-pseudo-count'], ['--p-transformation'])
        self.check_command('qiime gglasso calculate-covariance --p-method scaled --p-bias', ['--p-method', '--p-bias'], [])

    def test_latent_auto_is_not_false(self):
        self.check_command('qiime gglasso solve-problem --p-no-latent --p-gamma 0.01', ['--p-gamma'], ['--p-no-latent'])
        self.check_command('qiime gglasso solve-problem --p-latent --p-gamma 0.25', [], ['--p-latent', '--p-gamma'])

    def test_metadata_quoting(self):
        self.check_command('qiime feature-table filter-samples --p-where "[split]=\'train\'" --p-no-filter-empty-features', [], ['--p-where', '--p-no-filter-empty-features'])


if __name__ == '__main__':
    unittest.main()
