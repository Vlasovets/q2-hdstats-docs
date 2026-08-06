# Installation

## Prerequisites

Before installing q2-classo and q2-gglasso, ensure you have:

- QIIME 2 version **2026.7** or later
- Python **3.10** or later (the 2026.7 distribution ships Python 3.12)

Both plugins were migrated to the 2026.7 stack (Python 3.12, NumPy 2.x,
pandas 2.3). They will not install into a NumPy 1.x environment.

```{note}
**Two upstream renames you will run into.** In QIIME 2 **2026.4** the `amplicon`
distribution was renamed to **`qiime2`**, and the environment files were renamed
from `qiime2-*` to **`rachis-*`** — the framework package itself was rebranded
from `qiime2` to `rachis`. A compatibility shim keeps `import qiime2` working, so
existing analysis scripts do not need to change, but every install URL does.
```

## Installing QIIME 2

If you don't have QIIME 2 installed, follow the official installation guide at
[library.qiime2.org](https://library.qiime2.org/quickstart/qiime2).

Available distributions are `qiime2` (what this tutorial uses), `moshpit`,
`pathogenome` and `tiny`. Platform support for the `qiime2` distribution is
`linux-64` and `osx-64`; there is no `osx-arm64` build.

---

## Getting Help

If you encounter issues:

1. Check the [QIIME 2 Forum](https://forum.qiime2.org/)
2. Visit the plugin repositories:
   - [q2-gglasso GitHub](https://github.com/Vlasovets/q2-gglasso)
   - [q2-classo GitHub](https://github.com/Vlasovets/q2-classo)
3. Review the documentation:
   - [gglasso documentation](https://gglasso.readthedocs.io/en/latest/#)
   - [classo documentation](https://c-lasso.readthedocs.io/en/latest/index.html#)

## Next Steps

Once installation is complete, proceed to the installation of q2-gglasso and q2-classo. Both plugins are independent from each other, so you can install each one of them separately if you only want to do network analysis (q2-gglasso) or classification/regression tasks (q2-classo).
