#!/usr/bin/env python3
"""Plot the rank-one QIIME low-rank component and its sample-score associations."""

from __future__ import annotations

import hashlib
import json
import tempfile
import zipfile
from pathlib import Path

import biom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr


ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"
EXPORTS = ROOT / "exports"
FIGURES = ROOT / "figures"
PROVENANCE = ROOT / "provenance"
SOLUTION = ART / "motus-60-clr-slr-ebic-gamma010.qza"
TRANSFORMED_TABLE = ART / "motus-60-clr-pseudocount1.qza"
METADATA = EXPORTS / "motus-outcomes-all36.tsv"


def extract_member(path: Path, suffix: str, destination: Path) -> Path:
    with zipfile.ZipFile(path) as archive:
        member = next(name for name in archive.namelist() if name.endswith(suffix))
        archive.extract(member, destination)
    return destination / member


def open_solution(path: Path):
    temporary = tempfile.TemporaryDirectory()
    nested = extract_member(path, "problem.zip", Path(temporary.name))
    with zipfile.ZipFile(nested) as archive:
        archive.extractall(Path(temporary.name, "problem"))
    return zarr.open_group(Path(temporary.name, "problem"), mode="r"), temporary


def solution_labels(root) -> list[str]:
    return [str(root["labels"][str(index)][()]) for index in range(len(root["labels"]))]


def table_dataframe(path: Path) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as temporary:
        biom_path = extract_member(path, "data/feature-table.biom", Path(temporary))
        # BIOM is features x samples; q2-gglasso stores transformed tables as samples x features.
        return biom.load_table(str(biom_path)).to_dataframe(dense=True).T


def display_label(feature: str) -> str:
    return "mOTU " + feature.rsplit("_", 1)[-1]


def main():
    for directory in (EXPORTS, FIGURES, PROVENANCE):
        directory.mkdir(parents=True, exist_ok=True)
    solution, temporary = open_solution(SOLUTION)
    lowrank = np.asarray(solution["solution/lowrank_"], dtype=float)
    features = solution_labels(solution)
    eigenvalues, eigenvectors = np.linalg.eigh(lowrank)
    component = eigenvectors[:, np.argmax(np.abs(eigenvalues))]
    # Fix arbitrary eigenvector sign deterministically for a stable figure.
    component *= np.sign(component[np.argmax(np.abs(component))])
    if np.linalg.matrix_rank(lowrank, tol=1e-7) != 1:
        raise ValueError("This diagnostic is defined only for a rank-one selected low-rank component.")

    transformed = table_dataframe(TRANSFORMED_TABLE).reindex(columns=features)
    metadata = pd.read_csv(METADATA, sep="\t").set_index("sample-id")
    samples = transformed.index.intersection(metadata.index)
    if len(samples) != 34:
        raise ValueError(f"Expected 34 selected samples, found {len(samples)} shared samples.")
    score = transformed.loc[samples].to_numpy() @ component
    data = metadata.loc[samples, ["host_age_days", "diagnosis"]].copy()
    data["LV1_score"] = score
    data.index.name = "sample-id"
    data.reset_index().to_csv(EXPORTS / "motus-60-slr-selected-lv1-scores.tsv", sep="\t", index=False)
    pd.DataFrame(lowrank, index=features, columns=features).rename_axis("mOTU").to_csv(
        EXPORTS / "motus-60-slr-selected-lowrank-matrix.tsv", sep="\t"
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), constrained_layout=True)

    age = data["host_age_days"].to_numpy(float)
    axes[0].scatter(age, score, s=38, color="#2166ac", edgecolor="white", linewidth=0.5)
    slope, intercept = np.polyfit(age, score, 1)
    xline = np.linspace(age.min(), age.max(), 100)
    axes[0].plot(xline, slope * xline + intercept, color="#555555", linewidth=1)
    rho = float(pd.Series(age).corr(pd.Series(score), method="spearman"))
    axes[0].text(0.04, 0.95, f"Spearman rho = {rho:.2f}", transform=axes[0].transAxes, va="top", fontsize=8)
    axes[0].set(xlabel="Infant age (days)", ylabel="LV1 score", title="LV1 and infant age")
    axes[0].grid(color="#dddddd", linewidth=0.6)

    groups = ["control", "asthma"]
    values = [data.loc[data["diagnosis"] == group, "LV1_score"].to_numpy() for group in groups]
    axes[1].boxplot(values, tick_labels=["Control", "Asthma"], widths=0.5, patch_artist=True, boxprops={"facecolor": "#d9d9d9", "edgecolor": "#555555"}, medianprops={"color": "#222222"})
    rng = np.random.default_rng(20261021)
    for position, values_for_group, color in zip((1, 2), values, ("#2166ac", "#b2182b")):
        axes[1].scatter(position + rng.uniform(-0.09, 0.09, len(values_for_group)), values_for_group, s=31, color=color, edgecolor="white", linewidth=0.45, zorder=3)
    axes[1].set(ylabel="LV1 score", title="LV1 and maternal asthma")
    axes[1].grid(axis="y", color="#dddddd", linewidth=0.6)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)

    output = FIGURES / "motus-60-slr-selected-lowrank-and-covariates"
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    provenance = {
        "figure": str(output.relative_to(ROOT)),
        "selected_solution": str(SOLUTION.relative_to(ROOT)),
        "transformed_table": str(TRANSFORMED_TABLE.relative_to(ROOT)),
        "metadata": str(METADATA.relative_to(ROOT)),
        "solution_sha256": hashlib.sha256(SOLUTION.read_bytes()).hexdigest(),
        "source_of_truth": "q2-gglasso selected sparse-plus-low-rank solution and the matching QIIME CLR-transformed table.",
        "score_definition": "CLR-transformed sample vector multiplied by the unit eigenvector associated with the nonzero eigenvalue of L.",
        "rank": 1,
        "sign_convention": "The loading of largest absolute magnitude is positive.",
        "n_samples": int(len(samples)),
        "spearman_rho_age": rho,
    }
    (PROVENANCE / "motus-60-slr-selected-lowrank-and-covariates.json").write_text(json.dumps(provenance, indent=2) + "\n")
    temporary.cleanup()


if __name__ == "__main__":
    main()
