#!/usr/bin/env python3
"""Export and plot selected q2-gglasso mOTU solutions without refitting them."""

from __future__ import annotations

import hashlib
import json
import tempfile
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr


ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"
EXPORTS = ROOT / "exports"
FIGURES = ROOT / "figures"
PROVENANCE = ROOT / "provenance"
SGL = ART / "motus-60-mclr-sgl-ebic-gamma010.qza"
SLR = ART / "motus-60-clr-slr-ebic-gamma010.qza"
CLR_COVARIANCE = ART / "motus-60-clr-pseudocount1-correlation.qza"


def open_problem(path: Path):
    temporary = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(path) as outer:
        member = next(name for name in outer.namelist() if name.endswith("problem.zip"))
        outer.extract(member, temporary.name)
    with zipfile.ZipFile(Path(temporary.name, member)) as inner:
        inner.extractall(Path(temporary.name, "problem"))
    return zarr.open_group(Path(temporary.name, "problem"), mode="r"), temporary


def labels(root) -> list[str]:
    label_group = root["labels"]
    return [str(label_group[str(index)][()]) for index in range(len(label_group))]


def path_values(root, name: str) -> np.ndarray:
    values = np.asarray(root[f"modelselect_params/{name}"]).ravel()
    if values.size == 0 or (values.size == 1 and str(values[0]) == ""):
        return np.array([], dtype=float)
    return values.astype(float)


def scalar(root, path: str):
    value = root[path][()]
    return value.decode() if isinstance(value, bytes) else value


def partial_correlation(precision: np.ndarray) -> np.ndarray:
    diagonal = np.sqrt(np.outer(np.diag(precision), np.diag(precision)))
    matrix = -precision / diagonal
    np.fill_diagonal(matrix, 0.0)
    return matrix


def covariance_matrix(path: Path, feature_ids: list[str]) -> np.ndarray:
    with tempfile.TemporaryDirectory() as temporary:
        with zipfile.ZipFile(path) as archive:
            member = next(name for name in archive.namelist() if name.endswith("data/pairwise_comparisons.tsv"))
            archive.extract(member, temporary)
        frame = pd.read_csv(Path(temporary, member), sep="\t", index_col=0)
    return frame.reindex(index=feature_ids, columns=feature_ids).to_numpy(float)


def model_record(path: Path, name: str):
    root, temporary = open_problem(path)
    precision = np.asarray(root["solution/precision_"], dtype=float)
    lowrank = np.asarray(root["solution/lowrank_"], dtype=float) if "lowrank_" in root["solution"] else None
    lambda_grid = path_values(root, "lambda1_range")
    mu_grid = path_values(root, "mu1_range")
    selected_lambda = float(scalar(root, "modelselect_stats/BEST/lambda1"))
    selected_mu_raw = scalar(root, "modelselect_stats/BEST/mu1")
    selected_mu = None if lowrank is None or selected_mu_raw in ("", None) else float(selected_mu_raw)
    edges = int(np.count_nonzero(np.triu(np.abs(precision) > 1e-8, k=1)))
    rank = int(np.linalg.matrix_rank(lowrank, tol=1e-7)) if lowrank is not None else 0
    result = {
        "model": name,
        "artifact": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "labels": labels(root),
        "precision": precision,
        "lowrank": lowrank,
        "partial": partial_correlation(precision),
        "lambda1": selected_lambda,
        "mu1": selected_mu,
        "lambda1_path": lambda_grid.tolist(),
        "mu1_path": mu_grid.tolist(),
        "lambda1_on_boundary": bool(np.isclose(selected_lambda, lambda_grid[0]) or np.isclose(selected_lambda, lambda_grid[-1])),
        "mu1_on_boundary": None if selected_mu is None or mu_grid.size == 0 else bool(np.isclose(selected_mu, mu_grid[0]) or np.isclose(selected_mu, mu_grid[-1])),
        "rank": rank,
        "edge_count": edges,
    }
    temporary.cleanup()
    return result


def display_label(feature: str) -> str:
    return "mOTU " + feature.rsplit("_", 1)[-1]


def export_matrix(matrix: np.ndarray, feature_ids: list[str], path: Path):
    pd.DataFrame(matrix, index=feature_ids, columns=feature_ids).rename_axis("mOTU").to_csv(path, sep="\t")


def plot_partial_correlations(sgl: dict, slr: dict):
    if sgl["labels"] != slr["labels"]:
        raise ValueError("SGL and SLR feature orders differ; matched figure would be invalid.")
    feature_ids = sgl["labels"]
    readable = [display_label(label) for label in feature_ids]
    vmax = max(np.abs(sgl["partial"]).max(), np.abs(slr["partial"]).max())
    vmax = max(0.05, float(np.ceil(vmax * 20) / 20))

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2), constrained_layout=True)
    images = []
    for axis, record, title in zip(
        axes,
        (sgl, slr),
        ("SGL partial correlations", "Sparse component partial correlations"),
    ):
        image = axis.imshow(record["partial"], cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
        images.append(image)
        ticks = np.arange(0, len(readable), 5)
        axis.set_xticks(ticks, [readable[index] for index in ticks], rotation=90, fontsize=6)
        axis.set_yticks(ticks, [readable[index] for index in ticks], fontsize=6)
        axis.set_title(title, fontsize=10)
        axis.set_aspect("equal")
    colorbar = fig.colorbar(images[-1], ax=axes, shrink=0.79, pad=0.02)
    colorbar.set_label("Partial correlation", fontsize=9)
    colorbar.ax.tick_params(labelsize=7)
    output = FIGURES / "motus-60-sgl-vs-slr-selected-partial-correlation"
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return output, vmax


def plot_matrix_overview(covariance: np.ndarray, sgl: dict, slr: dict):
    lowrank = slr["lowrank"]
    partial_vmax = max(np.abs(sgl["partial"]).max(), np.abs(slr["partial"]).max())
    partial_vmax = max(0.05, float(np.ceil(partial_vmax * 20) / 20))
    covariance_vmax = 1.0
    lowrank_vmax = max(0.001, float(np.abs(lowrank).max()))

    # Mirror the compact fixed-cell S8 tutorial matrix style.  At p=60, the
    # feature IDs belong in the exported order table rather than on the axes.
    def add_grid(axis, size):
        axis.set_xticks(np.arange(-0.5, size, 1), minor=True)
        axis.set_yticks(np.arange(-0.5, size, 1), minor=True)
        axis.grid(which="minor", color="#f0f0f0", linewidth=0.22)
        axis.tick_params(which="minor", bottom=False, left=False)

    def draw(axis, matrix, title, vmax):
        image = axis.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
        axis.set_title(title, fontsize=8.5, pad=6)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        add_grid(axis, matrix.shape[0])
        for spine in axis.spines.values():
            spine.set_color("#333333")
            spine.set_linewidth(0.75)
        return image

    fig_width, fig_height = 8.27, 11.69  # A4 portrait, matching the supplement page.
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)
    cell_width = 0.053
    square = 60 * cell_width / fig_width
    left, right = 0.080, 0.535
    bottom, top = 0.190, 0.560
    ax_cov = fig.add_axes([left, top, square, square * fig_width / fig_height])
    ax_sgl = fig.add_axes([right, top, square, square * fig_width / fig_height])
    ax_lowrank = fig.add_axes([left, bottom, square, square * fig_width / fig_height])
    ax_sparse = fig.add_axes([right, bottom, square, square * fig_width / fig_height])
    cbar_height = square * fig_width / fig_height
    ax_cbar_cov = fig.add_axes([0.475, top, 0.013, cbar_height])
    ax_cbar_graph = fig.add_axes([0.932, top, 0.013, cbar_height])
    ax_cbar_lowrank = fig.add_axes([0.475, bottom, 0.013, cbar_height])
    ax_cbar_sparse = fig.add_axes([0.932, bottom, 0.013, cbar_height])

    image_cov = draw(ax_cov, covariance, "Empirical CLR correlation", covariance_vmax)
    draw(ax_sgl, sgl["partial"], "SGL partial correlations", partial_vmax)
    image_lowrank = draw(ax_lowrank, lowrank, "Low-rank component (L)", lowrank_vmax)
    image_sparse = draw(ax_sparse, slr["partial"], "Sparse partial correlations", partial_vmax)

    cbar_cov = fig.colorbar(image_cov, cax=ax_cbar_cov)
    cbar_cov.ax.tick_params(labelsize=6.0, length=2)
    cbar_graph = fig.colorbar(image_sparse, cax=ax_cbar_graph)
    cbar_graph.ax.tick_params(labelsize=6.0, length=2)
    cbar_lowrank = fig.colorbar(image_lowrank, cax=ax_cbar_lowrank)
    cbar_lowrank.ax.tick_params(labelsize=6.0, length=2)
    cbar_sparse = fig.colorbar(image_sparse, cax=ax_cbar_sparse)
    cbar_sparse.ax.tick_params(labelsize=6.0, length=2)
    output = FIGURES / "motus-60-clr-sgl-lowrank-sparse-matrices"
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    return output, {"covariance": covariance_vmax, "partial_correlation": partial_vmax, "lowrank": lowrank_vmax}


def main():
    for directory in (EXPORTS, FIGURES, PROVENANCE):
        directory.mkdir(parents=True, exist_ok=True)
    sgl = model_record(SGL, "SGL, mCLR")
    slr = model_record(SLR, "Sparse + low-rank, CLR(pseudocount=1)")
    if sgl["labels"] != slr["labels"]:
        raise ValueError("The QIIME artifacts do not have matched mOTU labels.")

    pd.DataFrame(
        [{key: value for key, value in record.items() if key not in {"labels", "precision", "lowrank", "partial", "lambda1_path", "mu1_path"}} for record in (sgl, slr)]
    ).to_csv(EXPORTS / "motus-60-gglasso-selection-summary.tsv", sep="\t", index=False)
    pd.DataFrame({"order": np.arange(1, len(sgl["labels"]) + 1), "feature_id": sgl["labels"], "display_label": [display_label(label) for label in sgl["labels"]]}).to_csv(
        EXPORTS / "motus-60-gglasso-motu-order.tsv", sep="\t", index=False
    )
    export_matrix(sgl["partial"], sgl["labels"], EXPORTS / "motus-60-sgl-selected-partial-correlation.tsv")
    export_matrix(slr["partial"], slr["labels"], EXPORTS / "motus-60-slr-selected-sparse-partial-correlation.tsv")
    figure, vmax = plot_partial_correlations(sgl, slr)
    covariance = covariance_matrix(CLR_COVARIANCE, sgl["labels"])
    overview, overview_scales = plot_matrix_overview(covariance, sgl, slr)

    provenance = {
        "figure": str(figure.relative_to(ROOT)),
        "source_of_truth": "Selected q2-gglasso GGLassoProblem artifacts; this script only extracts and converts precision matrices to partial correlations.",
        "gamma": 0.10,
        "covariance_normalization": "scaled covariance with bias=True (1/n)",
        "models": [{key: value for key, value in record.items() if key not in {"precision", "lowrank", "partial"}} for record in (sgl, slr)],
        "matched_motu_order_export": "exports/motus-60-gglasso-motu-order.tsv",
        "color_scale": {"vmin": -vmax, "vmax": vmax, "cmap": "RdBu_r"},
    }
    (PROVENANCE / "motus-60-gglasso-matrix-figure.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (PROVENANCE / "motus-60-clr-sgl-lowrank-sparse-matrices.json").write_text(json.dumps({
        "figure": str(overview.relative_to(ROOT)),
        "source_artifacts": [str(CLR_COVARIANCE.relative_to(ROOT)), str(SGL.relative_to(ROOT)), str(SLR.relative_to(ROOT))],
        "source_of_truth": "Selected QIIME covariance and GGLassoProblem artifacts; this script only converts precision matrices to partial correlations.",
        "mOTU_order": "exports/motus-60-gglasso-motu-order.tsv",
        "color_scales": overview_scales,
    }, indent=2) + "\n")
    print(pd.read_csv(EXPORTS / "motus-60-gglasso-selection-summary.tsv", sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
