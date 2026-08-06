# Verifying Your Installation

Two minutes of checking here saves a lot of confusion later. Plugin registration
failures do not appear until you actually invoke an action, and several of them
produce errors that point somewhere other than the real cause.

## 1. The framework

```bash
qiime info
```

Expect QIIME 2 **2026.7**. Since the 2026.1 rebrand the framework package is
called `rachis`, so `qiime info` reports that name — this is expected, and
`import qiime2` still works through a compatibility shim.

## 2. Both plugins are registered

```bash
qiime dev refresh-cache
qiime gglasso --help
qiime classo --help
```

`qiime gglasso --help` must list **six** actions:

```
  build-groups          build-groups
  calculate-covariance  calculate_covariance
  pca                   Principal component analysis (PCA)
  solve-problem         solve_problem
  summarize             Summary table
  transform-features    transform-features
```

The one-line descriptions are not all English sentences: q2cli uses each action's
registered `name=` as its short help, and four of the six register that field as
the action name itself. A description column that repeats the action name is
therefore expected, not a sign of a broken registration.

`qiime classo --help` must list **eight**: `add-covariates`, `add-taxa`,
`classify`, `generate-data`, `predict`, `regress`, `summarize`,
`transform-features`.

```{warning}
If two actions are both shown as **regress**, your q2-classo predates the fix for
`classify` having been registered under the wrong name. The action works; it is
only mislabelled. Update to a current checkout.
```

## 3. The scientific stack is the one you expect

The most common failure mode is not a missing package but a *wrong* one: `pip`
cannot see conda's pins and will happily install a wheel over the distribution's
NumPy.

```bash
python -c "import numpy, pandas, scipy, numba, bokeh, zarr; \
print('numpy', numpy.__version__); print('pandas', pandas.__version__); \
print('scipy', scipy.__version__); print('numba', numba.__version__); \
print('bokeh', bokeh.__version__); print('zarr', zarr.__version__)"
```

Expected on a clean 2026.7 environment:

| Package | Expected | Why it matters |
|---|---|---|
| numpy | 2.4.x | q2-gglasso previously pinned `<=1.27`; 2.x is required now |
| pandas | 2.3.x | 3.x changes Copy-on-Write semantics further |
| scipy | 1.17.x | distribution pin |
| numba | 0.66.x | compiles GGLasso's JIT solver kernels |
| bokeh | 3.x | 2.4.3 cannot render the visualizations |
| zarr | 2.18.x | **must be < 3** — zarr 3 removed `zarr.hierarchy.Group` |

And the two solver libraries:

```bash
python -c "import gglasso, classo; print('gglasso', gglasso.__version__)"
```

`gglasso` should be **0.3.0** or later.

## 4. A solver actually runs

Registration succeeding does not mean the numerics work — the JIT kernels are
compiled on first call, and that is where a numba/NumPy mismatch surfaces.

```bash
python - <<'PY'
import numpy as np
from gglasso.problem import glasso_problem
from gglasso.helper.data_generation import (
    generate_precision_matrix, sample_covariance_matrix)

Sigma, Theta = generate_precision_matrix(p=20, M=2, style="erdos", prob=0.1, seed=1)
S, _ = sample_covariance_matrix(Sigma, 100)
P = glasso_problem(S, N=100, reg_params={"lambda1": 0.05}, latent=False)
P.solve()
print("solver OK")
PY
```

The first call is slow — that is numba compiling, not a hang.

```{note}
If this raises a `TypingError` or `LoweringError`, you have hit a numba/NumPy
incompatibility rather than a q2-gglasso bug. Re-run with `NUMBA_DISABLE_JIT=1`
to confirm: the kernels are valid pure Python and will work, just slowly.
```

## 5. Read an artifact

```bash
qiime tools peek data/atacama-counts.qza
```

This confirms the artifact API and the type system agree with what the tutorial
expects.

## Known rough edges

Neither plugin declares `Choices()` on its string parameters, so a misspelled
enum value is accepted by the CLI and only fails inside the function at runtime:

```
ValueError: Unknown transformation name, use clr and not 'clrr'
```

Affected: `--p-transformation`, `--p-method`, `--p-reg`, `--p-path-scale` and the
`--p-*-numerical-method` family. See
[Troubleshooting](../90_reference/04_troubleshooting.md) for the full list.
