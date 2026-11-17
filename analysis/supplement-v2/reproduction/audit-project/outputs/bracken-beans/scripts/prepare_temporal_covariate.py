#!/usr/bin/env python3
"""Prepare a leakage-free time covariate for the Bracken temporal split."""

from __future__ import annotations

import json
import tempfile
import zipfile
from pathlib import Path

import biom
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"
METADATA = ROOT / "metadata.tsv"
TRAIN_CLR = ART / "bracken-species-prev14-total50k-train-clr05.qza"


def clr_table(path: Path) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as temp:
        with zipfile.ZipFile(path) as archive:
            member = next(name for name in archive.namelist() if name.endswith("data/feature-table.biom"))
            archive.extract(member, temp)
        table = biom.load_table(str(Path(temp, member)))
        return pd.DataFrame(table.matrix_data.toarray().T, index=table.ids(axis="sample"))


def main() -> None:
    metadata = pd.read_csv(METADATA, sep="\t", skiprows=[1]).set_index("sample-id")
    train_ids = metadata.index[metadata["analysis_set"].eq("train")]
    test_ids = metadata.index[metadata["analysis_set"].eq("test")]
    x = clr_table(TRAIN_CLR).loc[train_ids]
    time = metadata.loc[train_ids, "timepoint_hours"].astype(float)
    center = float(time.mean())
    centered = time - center
    target_norm = float(np.median(np.linalg.norm(x.to_numpy(float), axis=0)))
    scale = target_norm / float(np.linalg.norm(centered.to_numpy(float)))

    for name, ids in [("train", train_ids), ("test", test_ids)]:
        result = metadata.loc[ids, ["bean_type", "timepoint_hours", "analysis_set"]].copy()
        result["time_scaled"] = (result["timepoint_hours"].astype(float) - center) * scale
        result.index.name = "sample-id"
        result.to_csv(ROOT / f"metadata-{name}-time-scaled.tsv", sep="\t")

    (ROOT / "temporal-covariate-scaling.json").write_text(json.dumps({
        "training_time_mean_hours": center,
        "training_centered_time_l2_norm": float(np.linalg.norm(centered.to_numpy(float))),
        "training_median_clr_column_l2_norm": target_norm,
        "multiplicative_scale": scale,
        "rule": "Time was centered and scaled from training data only to the median CLR microbial-column L2 norm.",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
