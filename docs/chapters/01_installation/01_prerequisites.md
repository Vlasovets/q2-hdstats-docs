# Installation

## Prerequisites

Before installing q2-classo and q2-gglasso, you need:

- QIIME 2 version 2026.7 or later
- Python 3.10 or later (the 2026.7 distribution ships Python 3.12)

Both plugins were migrated to the 2026.7 stack (Python 3.12, NumPy 2.x,
pandas 2.3). They will not install into a NumPy 1.x environment.

```{note}
**Two upstream renames.** In QIIME 2 2026.4 the `amplicon` distribution was
renamed to `qiime2`, and the environment files were renamed from `qiime2-*` to
`rachis-*` — the framework package itself was rebranded from `qiime2` to
`rachis`. A compatibility shim keeps `import qiime2` working, so existing
analysis scripts need no change, but every install URL does.
```

## Installing QIIME 2

If you do not have QIIME 2 installed, follow the official installation guide at
[library.qiime2.org](https://library.qiime2.org/quickstart/qiime2).

The distributions are `qiime2` (the one this book uses), `moshpit`,
`pathogenome` and `tiny`. The `qiime2` distribution builds for `linux-64` and
`osx-64`. There is no `osx-arm64` build.

## Support

When an install fails or an action misbehaves:

1. Check the [QIIME 2 Forum](https://forum.qiime2.org/)
2. Visit the plugin repositories:
   - [q2-gglasso GitHub](https://github.com/Vlasovets/q2-gglasso)
   - [q2-classo GitHub](https://github.com/Vlasovets/q2-classo)
3. Review the documentation:
   - [gglasso documentation](https://gglasso.readthedocs.io/en/latest/#)
   - [classo documentation](https://c-lasso.readthedocs.io/en/latest/index.html#)

## Installing the plugins

Install q2-gglasso and q2-classo next. The plugins are independent of each
other: install q2-gglasso alone for network analysis, or q2-classo alone for
classification and regression.
