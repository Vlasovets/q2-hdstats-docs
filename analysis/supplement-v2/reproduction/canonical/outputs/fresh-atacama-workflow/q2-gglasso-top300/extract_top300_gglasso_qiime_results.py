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
LABEL_MAP = WORKFLOW_DIR / "q2-classo-top300/top300-asv-label-map.tsv"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402


QZA = {
    "correlation": OUT / "atacama-top300-53-clr-pseudocount1-correlation.qza",
    "sgl": OUT / "atacama-top300-53-clr-pseudocount1-sgl-ebic-gamma025.qza",
    "mixed": OUT / "atacama-top300-53-clr-pseudocount1-mixed-curated-covariates-weight005-ebic-gamma025.qza",
    "slr": OUT / "atacama-top300-53-clr-pseudocount1-slr-lambda065-mu-sweep-gamma025.qza",
}

COVARIATES = [
    "elevation",
    "ph",
    "average-soil-relative-humidity",
    "average-soil-temperature",
    "percentcover",
    "log10-extract-concen",
]

COVARIATE_LABELS = {
    "elevation": "elevation",
    "ph": "pH",
    "average-soil-relative-humidity": "soil RH",
    "average-soil-temperature": "soil temp.",
    "percentcover": "cover",
    "log10-extract-concen": "log10 extract",
}


def read_asv_map() -> pd.DataFrame:
    mapping = pd.read_csv(LABEL_MAP, sep="\t")
    mapping["rank"] = mapping["asv_label"].str.extract(r"(\d+)").astype(int)
    return mapping.sort_values("rank").reset_index(drop=True)


def read_pairwise_qza(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as artifact:
        member = next(
            name for name in artifact.namelist() if name.endswith("/data/pairwise_comparisons.tsv")
        )
        df = pd.read_csv(artifact.open(member), sep="\t")
    return df.set_index(df.columns[0])


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


def matrix_from_solution(path: Path, key: str) -> tuple[pd.DataFrame, dict[str, float]]:
    root, tmp = open_problem(path)
    try:
        labels = read_labels(root)
        values = np.asarray(root[f"solution/{key}"], dtype=float)
        stats = {
            "lambda1": float(np.asarray(root["modelselect_stats/BEST/lambda1"]))
            if "modelselect_stats/BEST/lambda1" in root
            else np.nan,
            "mu1": float(np.asarray(root["modelselect_stats/BEST/mu1"]))
            if "modelselect_stats/BEST/mu1" in root
            else np.nan,
        }
        if key == "lowrank_":
            stats["rank"] = int(np.linalg.matrix_rank(values))
        return pd.DataFrame(values, index=labels, columns=labels), stats
    finally:
        tmp.cleanup()


def precision_to_partial(theta: pd.DataFrame) -> pd.DataFrame:
    values = theta.to_numpy(dtype=float)
    denom = np.sqrt(np.outer(np.diag(values), np.diag(values)))
    with np.errstate(divide="ignore", invalid="ignore"):
        partial = -values / denom
    partial[~np.isfinite(partial)] = 0.0
    np.fill_diagonal(partial, 0.0)
    return pd.DataFrame(partial, index=theta.index, columns=theta.columns)


def relabel_asv_square(matrix: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    features = mapping["feature_id"].tolist()
    labels = mapping["asv_label"].tolist()
    out = matrix.loc[features, features].copy()
    out.index = labels
    out.columns = labels
    return out


def relabel_asv_covariate_square(matrix: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    features = mapping["feature_id"].tolist()
    labels = mapping["asv_label"].tolist()
    present_covariates = [cov for cov in COVARIATES if cov in matrix.index]
    ordered = features + present_covariates
    ordered_labels = labels + [COVARIATE_LABELS[cov] for cov in present_covariates]
    out = matrix.loc[ordered, ordered].copy()
    out.index = ordered_labels
    out.columns = ordered_labels
    return out


def lowrank_loadings(lowrank: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    values = (lowrank.to_numpy(dtype=float) + lowrank.to_numpy(dtype=float).T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(values)
    order = np.argsort(np.abs(eigvals))[::-1][:n_components]
    loadings = {}
    for i, eig_idx in enumerate(order, start=1):
        vec = eigvecs[:, eig_idx] * np.sqrt(abs(eigvals[eig_idx]))
        if np.nansum(vec) < 0:
            vec = -vec
        loadings[f"LV{i}"] = vec
    return pd.DataFrame(loadings, index=lowrank.index)


def sparse_plus_latent_matrix(sparse_partial: pd.DataFrame, loadings: pd.DataFrame) -> pd.DataFrame:
    asvs = sparse_partial.index.tolist()
    latent = loadings.columns.tolist()
    out = pd.DataFrame(0.0, index=asvs + latent, columns=asvs + latent)
    out.loc[asvs, asvs] = sparse_partial.to_numpy(dtype=float)
    for lv in latent:
        out.loc[asvs, lv] = loadings.loc[asvs, lv].to_numpy(dtype=float)
        out.loc[lv, asvs] = loadings.loc[asvs, lv].to_numpy(dtype=float)
    return out


def edge_nodes(matrix: pd.DataFrame, tol: float = 1e-8) -> set[str]:
    arr = np.abs(matrix.to_numpy(dtype=float)) > tol
    np.fill_diagonal(arr, False)
    active = arr.any(axis=0) | arr.any(axis=1)
    return set(matrix.index[active])


def edge_count(matrix: pd.DataFrame, tol: float = 1e-8) -> int:
    arr = np.abs(matrix.to_numpy(dtype=float)) > tol
    np.fill_diagonal(arr, False)
    return int(np.triu(arr, 1).sum())


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    df.to_csv(path, sep="\t")
    return path


def main() -> None:
    mapping = read_asv_map()
    outputs: list[Path] = []

    corr = relabel_asv_square(read_pairwise_qza(QZA["correlation"]), mapping)
    sgl_theta, sgl_stats = matrix_from_solution(QZA["sgl"], "precision_")
    mixed_theta, mixed_stats = matrix_from_solution(QZA["mixed"], "precision_")
    slr_theta, slr_stats = matrix_from_solution(QZA["slr"], "precision_")
    slr_lowrank, lowrank_stats = matrix_from_solution(QZA["slr"], "lowrank_")

    sgl = relabel_asv_square(precision_to_partial(sgl_theta), mapping)
    mixed_full = relabel_asv_covariate_square(precision_to_partial(mixed_theta), mapping)
    mixed = mixed_full.loc[corr.index, corr.columns].copy()
    slr = relabel_asv_square(precision_to_partial(slr_theta), mapping)
    lowrank = relabel_asv_square(slr_lowrank, mapping)
    lv_loadings = lowrank_loadings(lowrank, n_components=2)
    slr_augmented = sparse_plus_latent_matrix(slr, lv_loadings)

    matrices = {
        "atacama-top300-clr-correlation-asv-order.tsv": corr,
        "atacama-top300-sgl-gamma025-partial-correlation-asv-order.tsv": sgl,
        "atacama-top300-mixed-gamma025-asv-block-partial-correlation-asv-order.tsv": mixed,
        "atacama-top300-mixed-gamma025-full-asv-covariate-partial-correlation.tsv": mixed_full,
        "atacama-top300-slr-gamma025-sparse-partial-correlation-asv-order.tsv": slr,
        "atacama-top300-slr-gamma025-lowrank-asv-order.tsv": lowrank,
        "atacama-top300-slr-gamma025-lowrank-top2-loadings-asv-order.tsv": lv_loadings,
        "atacama-top300-slr-gamma025-sparse-plus-lv12-augmented-matrix.tsv": slr_augmented,
    }
    for name, matrix in matrices.items():
        outputs.append(write_tsv(matrix, OUT / name))

    active_union = sorted(
        edge_nodes(sgl) | edge_nodes(mixed) | edge_nodes(slr),
        key=lambda label: int(label.split("-")[1]),
    )
    active_path = OUT / "atacama-top300-union-active-asvs-sgl-mixed-slr.tsv"
    pd.DataFrame(
        {
            "asv_label": active_union,
            "active_in_sgl": [label in edge_nodes(sgl) for label in active_union],
            "active_in_mixed": [label in edge_nodes(mixed) for label in active_union],
            "active_in_slr": [label in edge_nodes(slr) for label in active_union],
        }
    ).to_csv(active_path, sep="\t", index=False)
    outputs.append(active_path)

    summary = pd.DataFrame(
        [
            {
                "model": "SGL",
                "gamma": 0.25,
                "lambda1": sgl_stats["lambda1"],
                "mu1": np.nan,
                "rank": 0,
                "asv_asv_edges": edge_count(sgl),
                "active_asvs": len(edge_nodes(sgl)),
            },
            {
                "model": "mixed",
                "gamma": 0.25,
                "lambda1": mixed_stats["lambda1"],
                "mu1": np.nan,
                "rank": 0,
                "asv_asv_edges": edge_count(mixed),
                "active_asvs": len(edge_nodes(mixed)),
            },
            {
                "model": "SLR",
                "gamma": 0.25,
                "lambda1": slr_stats["lambda1"],
                "mu1": slr_stats["mu1"],
                "rank": lowrank_stats["rank"],
                "asv_asv_edges": edge_count(slr),
                "active_asvs": len(edge_nodes(slr)),
            },
            {
                "model": "union",
                "gamma": 0.25,
                "lambda1": np.nan,
                "mu1": np.nan,
                "rank": np.nan,
                "asv_asv_edges": np.nan,
                "active_asvs": len(active_union),
            },
        ]
    )
    summary_path = OUT / "atacama-top300-gglasso-sgl-mixed-slr-summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    outputs.append(summary_path)

    prov = write_provenance(
        OUT / "atacama-top300-gglasso-qza-extraction",
        description="Extract top-300 q2-gglasso SGL, mixed, and SLR matrices from QIIME artifacts.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=[LABEL_MAP, *QZA.values()],
        outputs=outputs,
        parameters={
            "gamma": 0.25,
            "sgl_lambda1": sgl_stats["lambda1"],
            "mixed_lambda1": mixed_stats["lambda1"],
            "slr_lambda1": slr_stats["lambda1"],
            "slr_mu1": slr_stats["mu1"],
            "partial_correlation": "-Theta_ij / sqrt(Theta_ii * Theta_jj), diagonal set to zero",
            "mixed_covariate_labels": COVARIATE_LABELS,
            "slr_latent_loading_extraction": "Top two eigencomponents of the QIIME-extracted low-rank matrix, scaled by sqrt(abs(eigenvalue)) and sign-oriented to positive sum.",
            "active_asv_union_definition": "ASVs with at least one non-zero ASV-ASV partial correlation in SGL, mixed, or SLR.",
        },
        notes=[
            "All fitted matrices are extracted from q2-gglasso QIIME 2 artifacts.",
            "No direct GGLasso solver call is made in this extraction.",
            "The reduced plotting set is the union of edge-active ASVs across the three fitted models.",
            "The mixed full matrix appends the six tutorial covariates as rows and columns after ASV-300.",
            "The SLR augmented matrix appends LV1 and LV2 as rows and columns after ASV-300.",
        ],
    )

    print(
        json.dumps(
            {
                "summary": summary.to_dict(orient="records"),
                "active_union_n": len(active_union),
                "provenance": str(prov),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
