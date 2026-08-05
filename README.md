[![GitHub Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://Vlasovets.github.io/q2-hdstats-docs/)
[![deploy-book](https://github.com/Vlasovets/q2-hdstats-docs/actions/workflows/ci.yml/badge.svg)](https://github.com/Vlasovets/q2-hdstats-docs/actions/workflows/ci.yml)


# High-dimensional Statistics with QIIME 2

[![Github Pages](https://img.shields.io/badge/github%20pages-121013?style=for-the-badge&logo=github&logoColor=white)](https://Vlasovets.github.io/q2-hdstats-docs/)

Welcome to the **High-dimensional Statistics with QIIME 2** documentation!  
This site features tutorials and usage guides for the [`q2-gglasso`](https://github.com/Vlasovets/q2-gglasso) and [`q2-classo`](https://github.com/Vlasovets/q2-classo) plugins, which bring advanced high-dimensional statistical modeling to microbiome analysis within the QIIME 2 framework.

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