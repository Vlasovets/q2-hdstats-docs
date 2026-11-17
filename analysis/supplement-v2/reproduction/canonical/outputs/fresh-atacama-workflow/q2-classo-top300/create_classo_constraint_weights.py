from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from zipfile import ZipFile

ROOT = Path('/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial')
OUT = ROOT / 'outputs/fresh-atacama-workflow/q2-classo-top300'

n_features = len(pd.read_csv(OUT / 'top300-training-supported-feature-ids.tsv', sep='\t'))
C = np.ones((1, n_features), dtype=float)
w = np.ones(n_features, dtype=float)

for name, arr in [('top300-zero-sum-constraint.zip', C), ('top300-unit-weights.zip', w)]:
    path = OUT / name
    if path.exists():
        path.unlink()
    zarr.save(str(path), arr)
    print(path)
