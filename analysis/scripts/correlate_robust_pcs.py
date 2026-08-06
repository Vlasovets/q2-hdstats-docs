import argparse
import numpy as np
import pandas as pd
import zarr
from biom import load_table
from q2_gglasso.utils import PCA
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


parser = argparse.ArgumentParser()
parser.add_argument(
    "--solution-zarr",
    default="work/atacama-full/top-300-single-zarr",
)
parser.add_argument(
    "--output-prefix",
    default="outputs/atacama-q2-gglasso/atacama-top-300-robust-pc5",
)
args = parser.parse_args()

table = load_table("work/atacama-full/top-300-clr-export/feature-table.biom")
counts = table.to_dataframe(dense=True)
low_rank = np.asarray(
    zarr.open(args.solution_zarr, mode="r")["solution/lowrank_"]
)
scores, _, eigenvalues = PCA(counts, low_rank, inverse=True)

pcs = pd.DataFrame(
    scores[:, :5],
    index=counts.index,
    columns=[f"robust-PC{i}" for i in range(1, 6)],
)
pcs.index.name = "sample-id"
pcs.to_csv(f"{args.output_prefix}-scores.tsv", sep="\t")

metadata = pd.read_csv(
    "work/atacama-full/sample-metadata.tsv",
    sep="\t",
    index_col=0,
    skiprows=[1],
)
numeric = metadata.apply(pd.to_numeric, errors="coerce")

rows = []
for pc in pcs.columns:
    for covariate in numeric.columns:
        paired = pd.concat([pcs[pc], numeric[covariate]], axis=1).dropna()
        if len(paired) < 3 or paired[covariate].nunique() < 2:
            continue
        rho, p_value = spearmanr(paired[pc], paired[covariate])
        rows.append(
            {
                "pc": pc,
                "covariate": covariate,
                "n": len(paired),
                "spearman-rho": rho,
                "p-value": p_value,
            }
        )

results = pd.DataFrame(rows)
results["q-value"] = multipletests(results["p-value"], method="fdr_bh")[1]
results["abs-rho"] = results["spearman-rho"].abs()
results = results.sort_values(
    ["q-value", "abs-rho"], ascending=[True, False]
).drop(columns="abs-rho")
results.to_csv(f"{args.output_prefix}-correlations.tsv", sep="\t", index=False)

variance = eigenvalues / eigenvalues.sum()
pd.DataFrame(
    {
        "pc": pcs.columns,
        "eigenvalue": eigenvalues[:5],
        "variance-fraction": variance[:5],
    }
).to_csv(f"{args.output_prefix}-variance.tsv", sep="\t", index=False)
