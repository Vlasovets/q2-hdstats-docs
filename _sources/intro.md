# High-Dimensional Statistics with QIIME 2

We extend the [QIIME 2](https://qiime2.org/) microbiome multi-omics data science platform {cite}`bolyen2019reproducible` with advanced methods for high-dimensional statistical modeling. This suite currently features two powerful plugins:

- [**q2-gglasso**](https://github.com/bio-datascience/q2-gglasso): Solves General Graphical Lasso problems via sparse inverse covariance estimation {cite}`Schaipp2021, friedman2008sparse`.
- [**q2-classo**](https://github.com/bio-datascience/q2-classo): Implements sparse log-contrast regression and classification models {cite}`Simpson2021, combettes2021regression`.

<!-- ```{figure} images/png/overview.png
---
name: q2-overview
alt: High-dimensional statistics with QIIME2 overview
width: 600px
align: center
---
High-dimensional statistics with QIIME2.
``` -->

With these plugins, you can:

- **Analyze high-dimensional data** using regularization to control overfitting and improve interpretability.
- **Infer microbial association networks** with graphical lasso and group graphical lasso.
- **Perform sparse, interpretable classification and regression** with log-contrast models.
- **Build reproducible workflows** leveraging QIIME 2 artifact system.

Both plugins are fully integrated with the QIIME 2 framework, enabling reproducible workflows and seamless interoperability. They aim to help researchers model, visualize, and interpret complex microbiome data using state-of-the-art high-dimensional statistical methods.

---

This documentation will guide you through setup, installation, and example analyses using q2-gglasso and q2-classo, demonstrating how high-dimensional statistics can advance your microbiome research.