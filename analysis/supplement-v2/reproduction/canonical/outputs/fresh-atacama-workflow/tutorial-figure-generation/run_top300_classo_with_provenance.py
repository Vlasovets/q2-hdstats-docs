from __future__ import annotations

from pathlib import Path
import runpy

from provenance_utils import write_provenance

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
SCRIPT = Path(__file__)
TOP300 = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'
PLOT_SCRIPT = TOP300 / 'plot_top300_classo_elevation_matched_style.py'


def main():
    runpy.run_path(str(PLOT_SCRIPT), run_name='__main__')
    common_inputs = [
        TOP300 / 'top300-elevation-asv-vs-covariate-augmented-predictions.tsv',
        TOP300 / 'top300-elevation-asv-vs-covariate-augmented-selected-coefficients.tsv',
        TOP300 / 'top300-elevation-asv-vs-covariate-augmented-metrics.tsv',
        TOP300 / 'top300-prevalence-relative-abundance-ranking-from-54-samples.tsv',
        TOP300 / 'atacama-top300-elevation-asv-only-r1-cv.qza',
        TOP300 / 'atacama-top300-elevation-covariate-augmented-r1-cv.qza',
        ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top15/top15-asv-label-map.tsv',
        ROOT / 'work/atacama-full/taxonomy-silva138-export/taxonomy.tsv',
    ]
    figures = [
        ('atacama-top300-elevation-asv-vs-covariate-cv-curves-matched-style', 'q2-classo R1 five-fold CV curves for top-300 ASV-only and covariate-augmented elevation models'),
        ('atacama-top300-elevation-asv-only-r1-cv-three-panel-coefficients-in-oos-matched-style', 'q2-classo R1 CV refit, top-300 ASV-only elevation model'),
        ('atacama-top300-elevation-covariate-augmented-r1-cv-three-panel-coefficients-in-oos-matched-style', 'q2-classo R1 CV refit, top-300 covariate-augmented elevation model'),
    ]
    for stem, desc in figures:
        base = TOP300 / stem
        outputs = [base.with_suffix('.png'), base.with_suffix('.pdf')]
        write_provenance(base, description=desc, project_root=ROOT, script=SCRIPT,
                         inputs=[PLOT_SCRIPT, *common_inputs], outputs=outputs,
                         parameters={'figure_type': 'q2-classo top-300 matched-style figure'},
                         notes=['Wrapper records provenance for the plotting script that regenerates the figure from saved q2-classo artifacts and derived TSVs.'])


if __name__ == '__main__':
    main()
