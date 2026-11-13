from __future__ import annotations

import io
import json
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import zarr


PROJECT_ROOT = Path("/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial")
WORKFLOW_DIR = PROJECT_ROOT / "outputs/fresh-atacama-workflow"
OUT = WORKFLOW_DIR / "q2-gglasso-top300"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402

SGL_QZA = OUT / "atacama-top300-53-sgl-ebic-gamma030.qza"
CORR_QZA = OUT / "atacama-top300-53-mclr-correlation.qza"


def open_problem(path: Path):
    with zipfile.ZipFile(path) as artifact:
        member = next(name for name in artifact.namelist() if name.endswith("/data/problem.zip"))
        payload = artifact.read(member)
    tmp = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(io.BytesIO(payload)) as problem_zip:
        problem_zip.extractall(tmp.name)
    return zarr.open(tmp.name, mode="r"), tmp


def read_labels(root) -> list[str]:
    labels = root["labels"]
    return [str(np.asarray(labels[str(i)])) for i in range(len(labels))]


def read_pairwise_qza(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as artifact:
        member = next(
            name for name in artifact.namelist() if name.endswith("/data/pairwise_comparisons.tsv")
        )
        df = pd.read_csv(artifact.open(member), sep="\t")
    return df.set_index(df.columns[0])


def edge_count(theta: np.ndarray, tol: float = 1e-8) -> int:
    off = np.abs(theta) > tol
    np.fill_diagonal(off, False)
    return int(off.sum() // 2)


def connected_nodes(theta: np.ndarray, tol: float = 1e-8) -> int:
    off = np.abs(theta) > tol
    np.fill_diagonal(off, False)
    return int(off.any(axis=0).sum())


def partial_correlation(theta: pd.DataFrame) -> pd.DataFrame:
    values = theta.to_numpy(dtype=float)
    denom = np.sqrt(np.outer(np.diag(values), np.diag(values)))
    with np.errstate(divide="ignore", invalid="ignore"):
        partial = -values / denom
    partial[~np.isfinite(partial)] = 0.0
    np.fill_diagonal(partial, 0.0)
    return pd.DataFrame(partial, index=theta.index, columns=theta.columns)


def main() -> None:
    outputs: list[Path] = []
    root, tmp = open_problem(SGL_QZA)
    try:
        labels = read_labels(root)
        theta = pd.DataFrame(
            np.asarray(root["solution/precision_"], dtype=float),
            index=labels,
            columns=labels,
        )
        adjacency = np.asarray(root["solution/adjacency_"], dtype=float)
        best_lambda = float(np.asarray(root["modelselect_stats/BEST/lambda1"]))
        lambda_path = np.asarray(root["modelselect_stats/LAMBDA"], dtype=float).ravel()
        ebic = np.asarray(root["modelselect_stats/BIC/0.3"], dtype=float).ravel()
    finally:
        tmp.cleanup()

    theta_path = OUT / "atacama-top300-sgl-gamma030-selected-precision.tsv"
    partial_path = OUT / "atacama-top300-sgl-gamma030-selected-partial-correlation.tsv"
    theta.to_csv(theta_path, sep="\t")
    partial_correlation(theta).to_csv(partial_path, sep="\t")
    outputs.extend([theta_path, partial_path])

    path_df = pd.DataFrame({"lambda1": lambda_path, "ebic": ebic})
    path_path = OUT / "atacama-top300-sgl-gamma030-lambda-ebic-path.tsv"
    path_df.to_csv(path_path, sep="\t", index=False)
    outputs.append(path_path)

    corr = read_pairwise_qza(CORR_QZA)
    corr_path = OUT / "atacama-top300-mclr-correlation.tsv"
    corr.to_csv(corr_path, sep="\t")
    outputs.append(corr_path)

    summary = pd.DataFrame(
        [
            {
                "model": "SGL",
                "gamma": 0.30,
                "selected_lambda1": best_lambda,
                "n_features": theta.shape[0],
                "possible_off_diagonal_pairs": theta.shape[0] * (theta.shape[0] - 1) // 2,
                "edge_count_precision_nonzero": edge_count(theta.to_numpy()),
                "edge_count_adjacency_nonzero": edge_count(adjacency),
                "connected_nodes": connected_nodes(theta.to_numpy()),
                "min_offdiagonal_precision": float(
                    np.min(theta.to_numpy()[~np.eye(theta.shape[0], dtype=bool)])
                ),
                "max_offdiagonal_precision": float(
                    np.max(theta.to_numpy()[~np.eye(theta.shape[0], dtype=bool)])
                ),
            }
        ]
    )
    summary_path = OUT / "atacama-top300-sgl-gamma030-summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    outputs.append(summary_path)

    prov = write_provenance(
        OUT / "atacama-top300-sgl-gamma030-qza-extraction",
        description="Extract selected top-300 SGL baseline results from q2-gglasso QIIME artifacts.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=[SGL_QZA, CORR_QZA],
        outputs=outputs,
        parameters={
            "gamma": 0.30,
            "partial_correlation": "-Theta_ij / sqrt(Theta_ii * Theta_jj), diagonal set to zero",
            "edge_threshold": 1e-8,
        },
        notes=[
            "Inputs are QIIME 2 artifacts produced by run_top300_sgl_qiime_workflow.sh.",
            "No direct GGLasso solver call is made in this extraction.",
        ],
    )

    print(json.dumps({"summary": summary.to_dict(orient="records")[0], "provenance": str(prov)}, indent=2))


if __name__ == "__main__":
    main()
