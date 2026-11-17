# High-Dimensional Statistics with QIIME 2

[Published documentation](https://bio-datascience.github.io/q2-hdstats-docs/intro.html)

This revision rebuilds the q2-classo/q2-gglasso tutorial from the audited
supplement. Its mathematical models, commands, results, and figure captions
are the source of truth. The website puts the applications before the detailed
statistical background; the supplement retains its original S1-S5 order.
The original light/dark logos, group
attribution, and Jupyter Book theme are retained.

## Build and Validate

Use Python 3.11 and the repository's Jupyter Book 1.x requirements:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analysis/scripts/check_supplement_migration.py .
jupyter-book build docs --warningiserror
```

Open `docs/_build/html/intro.html` to preview the book. The initial conversion
uses Pandoc and Poppler. Some pages now have reviewed tutorial adaptations;
their original conversion baselines and edited hashes are retained separately.
Do not run the converter over those edits. Convert a new supplement into a clean
staging copy, compare against the baselines, and review changes before updating
the tutorial. Validate the edited book with:

```bash
python analysis/scripts/check_supplement_migration.py .
```

## Source and Provenance

| Path | Contents |
| --- | --- |
| `docs/tutorial/` | Thirteen pages derived from the supplement |
| `docs/images/tutorial-v2/` | Original figure PDFs, browser PNGs, provenance JSONs |
| `analysis/supplement-v2/source/` | Authoritative supplement snapshot |
| `analysis/supplement-v2/conversion-manifest.json` | Section mapping, command blocks, and figure hashes |
| `analysis/supplement-v2/reproduction/` | Preserved plotting and workflow scripts |
| `analysis/supplement-v2/script-path-map.json` | Original-to-bundled script paths and hashes |
| `docs/reproducibility.md` | Exact source pins and publication limitations |

Legacy analysis and tutorial files are retained for history, but are not the
source of the revised examples. Their old URLs redirect to the new subject areas.
The style-only invariants in `STYLE.md` concern prose passes on existing content;
this migration replaces superseded examples with the audited supplement and
checks that source using `check_supplement_migration.py`.

## Release Status

This branch is a development review version, not a published release. The top-15
Atacama walkthrough supplies prepared inputs and portable regression/classification
plotting scripts. Other input downloads, the full installation recipe, and
remaining plotting adapters still need packaging. Public
distribution of the supplied Bracken artifacts must be confirmed first. Older
Zenodo placeholder URLs must not be used as release inputs.

Branch pushes and pull requests build the documentation; only a push to `main`
deploys GitHub Pages. Do not merge until the input-distribution requirements
have been resolved.
