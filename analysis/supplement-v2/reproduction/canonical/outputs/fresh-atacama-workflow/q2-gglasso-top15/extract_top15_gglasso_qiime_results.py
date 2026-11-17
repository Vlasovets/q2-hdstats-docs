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
OUT_DIR = WORKFLOW_DIR / "q2-gglasso-top15"
LABEL_MAP = WORKFLOW_DIR / "q2-classo-top15/top15-asv-label-map.tsv"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402


QZA = {
    "mclr_correlation": OUT_DIR / "atacama-top15-53-mclr-correlation.qza",
    "sgl": OUT_DIR / "atacama-top15-53-sgl-ebic-gamma001.qza",
    "mixed": OUT_DIR / "atacama-top15-53-mixed-curated-covariates-weight005-ebic-gamma001.qza",
    "slr": OUT_DIR / "atacama-top15-53-slr-ebic-gamma001-mu-grid.qza",
}

COVARIATE_LABELS = {
    "elevation": "elevation",
    "ph": "pH",
    "average-soil-relative-humidity": "soil humidity",
    "average-soil-temperature": "soil temperature",
    "percentcover": "percent cover",
    "log10-extract-concen": "log10 extract conc.",
}


def read_asv_map() -> pd.DataFrame:
    mapping = pd.read_csv(LABEL_MAP, sep="\t")
    mapping["rank"] = mapping["asv_label"].str.extract(r"(\d+)").astype(int)
    return mapping.sort_values("rank").reset_index(drop=True)


def read_pairwise_qza(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as artifact:
        member = next(
            name
            for name in artifact.namelist()
            if name.endswith("/data/pairwise_comparisons.tsv")
        )
        df = pd.read_csv(artifact.open(member), sep="\t")
    return df.set_index(df.columns[0])


def open_problem(path: Path):
    with zipfile.ZipFile(path) as artifact:
        member = next(
            name for name in artifact.namelist() if name.endswith("/data/problem.zip")
        )
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
        stats = {}
        for stat in ("lambda1", "mu1"):
            stat_key = f"modelselect_stats/BEST/{stat}"
            if stat_key in root:
                stats[stat] = float(np.asarray(root[stat_key]))
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


def relabel_asvs(matrix: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    features = mapping["feature_id"].tolist()
    labels = mapping["asv_label"].tolist()
    out = matrix.loc[features, features].copy()
    out.index = labels
    out.columns = labels
    return out


def relabel_full_mixed(matrix: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    features = mapping["feature_id"].tolist()
    feature_labels = mapping["asv_label"].tolist()
    covariates = [name for name in matrix.index if name not in set(features)]
    ordered = features + covariates
    labels = feature_labels + [COVARIATE_LABELS.get(name, name) for name in covariates]
    out = matrix.loc[ordered, ordered].copy()
    out.index = labels
    out.columns = labels
    return out


def leading_lowrank_loading(lowrank: pd.DataFrame) -> pd.Series:
    sym = (lowrank.to_numpy(dtype=float) + lowrank.to_numpy(dtype=float).T) / 2
    eigvals, eigvecs = np.linalg.eigh(sym)
    idx = int(np.argmax(np.abs(eigvals)))
    loading = eigvecs[:, idx] * np.sqrt(abs(eigvals[idx]))
    if np.nansum(loading) < 0:
        loading = -loading
    return pd.Series(loading, index=lowrank.index, name="LV1")


def edge_count(partial: pd.DataFrame, tol: float = 1e-9) -> int:
    values = partial.to_numpy(dtype=float)
    upper = np.triu(np.abs(values) > tol, k=1)
    return int(upper.sum())


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    df.to_csv(path, sep="\t", index=True)
    return path


def main() -> None:
    mapping = read_asv_map()
    outputs: list[Path] = []
    summaries: list[dict[str, object]] = []

    corr = read_pairwise_qza(QZA["mclr_correlation"])
    corr_asv = relabel_asvs(corr, mapping)
    outputs.append(write_tsv(corr_asv, OUT_DIR / "atacama-top15-mclr-correlation-asv-order.tsv"))

    sgl_theta, sgl_stats = matrix_from_solution(QZA["sgl"], "precision_")
    sgl_partial = relabel_asvs(precision_to_partial(sgl_theta), mapping)
    outputs.append(write_tsv(sgl_partial, OUT_DIR / "atacama-top15-sgl-gamma001-partial-correlation-asv-order.tsv"))
    summaries.append({"model": "SGL", **sgl_stats, "edges": edge_count(sgl_partial)})

    mixed_theta, mixed_stats = matrix_from_solution(QZA["mixed"], "precision_")
    mixed_partial = precision_to_partial(mixed_theta)
    mixed_asv = relabel_asvs(mixed_partial, mapping)
    mixed_full = relabel_full_mixed(mixed_partial, mapping)
    outputs.append(write_tsv(mixed_asv, OUT_DIR / "atacama-top15-mixed-weight005-gamma001-asv-block-partial-correlation-asv-order.tsv"))
    outputs.append(write_tsv(mixed_full, OUT_DIR / "atacama-top15-mixed-weight005-gamma001-full-asv-covariate-partial-correlation.tsv"))
    summaries.append({"model": "mixed covariate-adjusted", **mixed_stats, "asv_edges": edge_count(mixed_asv)})

    slr_theta, slr_stats = matrix_from_solution(QZA["slr"], "precision_")
    slr_lowrank, _ = matrix_from_solution(QZA["slr"], "lowrank_")
    slr_partial = relabel_asvs(precision_to_partial(slr_theta), mapping)
    slr_lowrank_asv = relabel_asvs(slr_lowrank, mapping)
    loading = leading_lowrank_loading(slr_lowrank_asv)

    augmented = slr_partial.copy()
    augmented["LV1"] = loading
    augmented.loc["LV1"] = list(loading) + [0.0]

    outputs.append(write_tsv(slr_partial, OUT_DIR / "atacama-top15-slr-gamma001-sparse-partial-correlation-asv-order.tsv"))
    outputs.append(write_tsv(slr_lowrank_asv, OUT_DIR / "atacama-top15-slr-gamma001-lowrank-rank1-matrix-asv-order.tsv"))
    outputs.append(write_tsv(loading.to_frame(), OUT_DIR / "atacama-top15-slr-gamma001-lowrank-rank1-loading-asv-order.tsv"))
    outputs.append(write_tsv(augmented, OUT_DIR / "atacama-top15-slr-gamma001-sparse-plus-lv1-augmented-matrix.tsv"))
    summaries.append({"model": "SLR latent rank-1", **slr_stats, "sparse_edges": edge_count(slr_partial)})

    summary_path = OUT_DIR / "atacama-top15-q2-gglasso-model-selection-summary.tsv"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
    outputs.append(summary_path)

    extraction_provenance = OUT_DIR / "atacama-top15-q2-gglasso-qza-extraction"
    prov_path = write_provenance(
        extraction_provenance,
        description="Extract top-15 q2-gglasso tutorial matrices from QIIME 2 artifacts.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=[LABEL_MAP, *QZA.values()],
        outputs=outputs,
        parameters={
            "asv_order": "ASV-1 through ASV-15 from top15-asv-label-map.tsv",
            "partial_correlation": "-Theta_ij / sqrt(Theta_ii * Theta_jj), diagonal set to zero",
            "lowrank_loading": "leading eigenvector of q2-gglasso lowrank_ matrix, scaled by sqrt(abs(eigenvalue))",
        },
        notes=[
            "No direct GGLasso solver call is made here; q2-gglasso model outputs are read from QIIME 2 .qza artifacts.",
            "The mCLR correlation input is the q2-gglasso calculate-covariance scaled output.",
        ],
    )

    print(json.dumps({"outputs": [str(p) for p in outputs], "provenance": str(prov_path)}, indent=2))


if __name__ == "__main__":
    main()
