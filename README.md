[![GitHub Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://Vlasovets.github.io/q2-hdstats-docs/)
[![deploy-book](https://github.com/Vlasovets/q2-hdstats-docs/actions/workflows/ci.yml/badge.svg)](https://github.com/Vlasovets/q2-hdstats-docs/actions/workflows/ci.yml)


# High-dimensional Statistics with QIIME 2

[![Github Pages](https://img.shields.io/badge/github%20pages-121013?style=for-the-badge&logo=github&logoColor=white)](https://Vlasovets.github.io/q2-hdstats-docs/)

Welcome to the **High-dimensional Statistics with QIIME 2** documentation!  
This site features tutorials and usage guides for the [`q2-gglasso`](https://github.com/Vlasovets/q2-gglasso) and [`q2-classo`](https://github.com/Vlasovets/q2-classo) plugins, which bring advanced high-dimensional statistical modeling to microbiome analysis within the QIIME 2 framework.

---

## Repository layout

Despite the name, this repository holds **both** the book and the analysis that
produces the numbers in it.

| path | what it is |
|---|---|
| `docs/` | the Jupyter Book. `jupyter-book` reads only this directory (`path_to_book: docs`) |
| `docs/_data/` | tables the chapters render, **generated** by the analysis — never edited by hand |
| `analysis/` | the recompute pipeline: SLURM stages, scripts, pinned environment, reports |
| `main.tex` | the F1000Research Application Note |

The two used to be separate repositories, which meant the analysis wrote its
generated tables and the data manifest into a *sibling checkout by absolute
path*. The pipeline consequently ran on one machine only. Keeping them together
makes those paths relative and gives the manuscript's "Underlying data" section a
single artifact to cite.

### Running the analysis

Every SLURM stage resolves paths from one variable, so a different checkout needs
no edits:

```shell
export Q2_HDSTATS_REPO=/path/to/this/repo    # defaults to the author's checkout
sbatch analysis/slurm/01_lambda_path.sh
```

`analysis/data/`, `analysis/results/` and `analysis/publish/` are gitignored:
the repository tracks the **recipe**, not the output. Everything regenerates from
the committed scripts against the pinned environment in `analysis/envs/`.

Start with `analysis/README.md` for the stage-by-stage description, and
`analysis/reports/DECISIONS_NEEDED.md` for the decisions taken and the evidence
behind them.

---

## 🚀 Building the Documentation

**Create and activate conda environment:**
```shell
mamba create -n jupyter-book -c conda-forge jupyter-book

conda activate jupyter-book
```

**Build the book:**

```shell
jupyter-book build --all docs
```