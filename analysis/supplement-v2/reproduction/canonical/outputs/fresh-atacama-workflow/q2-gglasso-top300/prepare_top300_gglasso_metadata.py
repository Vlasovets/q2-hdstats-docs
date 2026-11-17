from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial")
OUT = ROOT / "outputs/fresh-atacama-workflow/q2-gglasso-top300"
CLASSO300 = ROOT / "outputs/fresh-atacama-workflow/q2-classo-top300"
META = ROOT / "outputs/fresh-atacama-workflow/metadata/atacama-54-sample-metadata.tsv"

COVARIATES = [
    "elevation",
    "ph",
    "average-soil-relative-humidity",
    "average-soil-temperature",
    "percentcover",
    "log10-extract-concen",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    split = pd.read_csv(CLASSO300 / "train-test-split-53-q2.tsv", sep="\t")
    samples = split.loc[split["sample-id"] != "#q2:types", "sample-id"].tolist()

    md = pd.read_csv(META, sep="\t")
    md = md.loc[md["sample-id"] != "#q2:types"].set_index("sample-id")
    raw = pd.DataFrame(index=samples)
    raw["elevation"] = pd.to_numeric(md.loc[samples, "elevation"], errors="coerce")
    raw["ph"] = pd.to_numeric(md.loc[samples, "ph"], errors="coerce")
    raw["average-soil-relative-humidity"] = pd.to_numeric(
        md.loc[samples, "average-soil-relative-humidity"], errors="coerce"
    )
    raw["average-soil-temperature"] = pd.to_numeric(
        md.loc[samples, "average-soil-temperature"], errors="coerce"
    )
    raw["percentcover"] = pd.to_numeric(md.loc[samples, "percentcover"], errors="coerce")
    extract = pd.to_numeric(md.loc[samples, "extract-concen"], errors="coerce")
    raw["log10-extract-concen"] = np.log10(extract)

    means = raw.mean(axis=0, skipna=True)
    imputed = raw.fillna(means)

    q2 = imputed.copy()
    q2.index.name = "sample-id"
    q2 = pd.concat(
        [pd.DataFrame({col: ["numeric"] for col in q2.columns}, index=["#q2:types"]), q2]
    )
    q2.index.name = "sample-id"
    q2.to_csv(OUT / "atacama-top300-53-curated-metadata-for-gglasso.tsv", sep="\t")

    report = pd.DataFrame(
        {
            "covariate": COVARIATES,
            "mean_used_for_imputation": [means[c] for c in COVARIATES],
            "missing_values": [int(raw[c].isna().sum()) for c in COVARIATES],
            "minimum_after_imputation": [float(imputed[c].min()) for c in COVARIATES],
            "maximum_after_imputation": [float(imputed[c].max()) for c in COVARIATES],
        }
    )
    report.to_csv(
        OUT / "atacama-top300-53-curated-metadata-imputation-summary.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "samples": len(samples),
        "covariates": COVARIATES,
        "source_metadata": str(META),
        "missing_values": dict(zip(COVARIATES, report["missing_values"].tolist())),
        "note": "Missing numeric covariate values are mean-imputed before passing metadata to q2-gglasso.",
    }
    (OUT / "atacama-top300-gglasso-metadata-manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
