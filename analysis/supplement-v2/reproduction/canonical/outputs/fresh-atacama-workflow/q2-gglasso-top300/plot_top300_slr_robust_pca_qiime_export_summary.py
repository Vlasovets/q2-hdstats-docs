from __future__ import annotations

import ast
import base64
import gzip
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path("/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial")
WORKFLOW_DIR = PROJECT_ROOT / "outputs/fresh-atacama-workflow"
OUT = WORKFLOW_DIR / "q2-gglasso-top300"
EXPORT_DIR = OUT / "exported-robust-pca-qzv"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402


COVARIATES = [
    ("elevation", "Elevation (m)"),
    ("ph", "pH"),
    ("average-soil-relative-humidity", "Soil relative humidity (%)"),
    ("average-soil-temperature", "Soil temperature (deg C)"),
    ("percentcover", "Vegetation cover (%)"),
    ("log10-extract-concen", "log10 extract concentration"),
]

HTML_INPUTS = {
    covariate: EXPORT_DIR
    / f"atacama-top300-slr-rank2-robust-pca-color-by-{covariate}"
    / "index.html"
    for covariate, _ in COVARIATES
}
QZV_INPUTS = {
    covariate: OUT / f"atacama-top300-slr-rank2-robust-pca-color-by-{covariate}.qzv"
    for covariate, _ in COVARIATES
}

PNG_OUT = OUT / "atacama-top300-slr-rank2-robust-pca-q2-export-2x3-covariates.png"
PDF_OUT = OUT / "atacama-top300-slr-rank2-robust-pca-q2-export-2x3-covariates.pdf"
SCORES_OUT = OUT / "atacama-top300-slr-rank2-robust-pca-q2-export-scores-and-covariates.tsv"
CORR_OUT = OUT / "atacama-top300-slr-rank2-robust-pca-q2-export-spearman-summary.tsv"


def parse_docs_json(html_path: Path) -> dict:
    text = html_path.read_text(errors="replace")
    start = text.find("const docs_json = ")
    if start < 0:
        raise ValueError(f"Could not find Bokeh docs_json in {html_path}")
    start += len("const docs_json = ")
    end = text.find(";\n        const render_items", start)
    if end < 0:
        raise ValueError(f"Could not find end of Bokeh docs_json in {html_path}")
    return json.loads(ast.literal_eval(text[start:end]))


def decode_ndarray(value):
    if not (isinstance(value, dict) and value.get("type") == "ndarray"):
        return value
    payload = value["array"]
    if isinstance(payload, list):
        return payload
    raw = gzip.decompress(base64.b64decode(payload["data"]))
    return np.frombuffer(raw, dtype=np.dtype(value["dtype"])).tolist()


def find_sample_source(html_path: Path) -> pd.DataFrame:
    docs = parse_docs_json(html_path)
    document = next(iter(docs.values()))
    sources: list[dict] = []

    def walk(obj):
        if isinstance(obj, dict):
            if obj.get("name") == "ColumnDataSource":
                sources.append(obj)
            for item in obj.values():
                walk(item)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(document)
    for source in sources:
        data = source.get("attributes", {}).get("data")
        if not (isinstance(data, dict) and data.get("type") == "map"):
            continue
        entries = data["entries"]
        keys = {key for key, _ in entries}
        if {"index", "x", "y"}.issubset(keys) and {name for name, _ in COVARIATES}.issubset(keys):
            decoded = {key: decode_ndarray(value) for key, value in entries}
            return pd.DataFrame(decoded)
    raise ValueError(f"Could not find q2-gglasso PCA sample ColumnDataSource in {html_path}")


def load_q2_exported_scores() -> pd.DataFrame:
    frames = []
    for covariate, _ in COVARIATES:
        frame = find_sample_source(HTML_INPUTS[covariate])
        frame = frame[["index", "x", "y", *[name for name, _ in COVARIATES]]].copy()
        frame["q2_export_color_by"] = covariate
        frames.append(frame)

    reference = frames[0].drop(columns=["q2_export_color_by"]).copy()
    for frame in frames[1:]:
        comparable = frame.drop(columns=["q2_export_color_by"]).copy()
        if not reference.equals(comparable):
            raise ValueError("The q2 exported PCA sources are not identical across color-by visualizations.")
    reference = reference.rename(columns={"index": "sample-id", "x": "PC1", "y": "PC2"})
    return reference


def correlation_summary(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for covariate, label in COVARIATES:
        values = pd.to_numeric(scores[covariate], errors="coerce")
        for component in ["PC1", "PC2"]:
            result = spearmanr(scores[component], values, nan_policy="omit")
            rows.append(
                {
                    "covariate": covariate,
                    "label": label,
                    "component": component,
                    "spearman_rho": result.statistic,
                    "p_value": result.pvalue,
                    "n": int((scores[component].notna() & values.notna()).sum()),
                }
            )
    return pd.DataFrame(rows)


def draw_panel(ax, scores: pd.DataFrame, covariate: str, label: str, corr: pd.DataFrame):
    values = pd.to_numeric(scores[covariate], errors="coerce")
    scatter = ax.scatter(
        scores["PC1"],
        scores["PC2"],
        c=values,
        cmap="viridis",
        s=34,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.92,
    )
    ax.axhline(0, color="#d7d7d7", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#d7d7d7", linewidth=0.7, zorder=0)
    ax.set_title(label, fontsize=9.5, pad=4)
    ax.set_xlabel("Latent PC1", fontsize=8)
    ax.set_ylabel("Latent PC2", fontsize=8)
    ax.tick_params(labelsize=6.8, length=2)
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.7)

    rho_pc1 = corr.query("covariate == @covariate and component == 'PC1'")["spearman_rho"].iloc[0]
    rho_pc2 = corr.query("covariate == @covariate and component == 'PC2'")["spearman_rho"].iloc[0]
    ax.text(
        0.03,
        0.97,
        f"rho(PC1)={rho_pc1:+.2f}\nrho(PC2)={rho_pc2:+.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.86},
    )
    return scatter


def main() -> None:
    scores = load_q2_exported_scores()
    scores.to_csv(SCORES_OUT, sep="\t", index=False)
    corr = correlation_summary(scores)
    corr.to_csv(CORR_OUT, sep="\t", index=False)

    xlim = (scores["PC1"].min() - 0.4, scores["PC1"].max() + 0.4)
    ylim = (scores["PC2"].min() - 0.4, scores["PC2"].max() + 0.4)

    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.1), dpi=350)
    plt.subplots_adjust(left=0.060, right=0.965, top=0.885, bottom=0.080, wspace=0.34, hspace=0.36)
    fig.suptitle(
        "Atacama top-300 sparse + low-rank model: q2-gglasso PCA visualization summary",
        fontsize=11,
        y=0.965,
    )
    for ax, (covariate, label) in zip(axes.ravel(), COVARIATES):
        scatter = draw_panel(ax, scores, covariate, label, corr)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.025)
        cbar.ax.tick_params(labelsize=6.3, length=2)

    fig.savefig(PNG_OUT, dpi=350, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)

    prov = write_provenance(
        PNG_OUT,
        description="Render a static 2x3 tutorial summary from exported q2-gglasso PCA visualizations for the top-300 sparse+low-rank Atacama model.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=[*HTML_INPUTS.values(), *QZV_INPUTS.values()],
        outputs=[PNG_OUT, PDF_OUT, SCORES_OUT, CORR_OUT],
        parameters={
            "source_of_coordinates": "Bokeh ColumnDataSource embedded in exported q2-gglasso PCA .qzv HTML files.",
            "q2_action": "qiime gglasso pca",
            "covariates": [name for name, _ in COVARIATES],
            "layout": "2x3",
            "axis_labels_from_q2_visualization": ["PC1 (55.0%)", "PC2 (44.9%)"],
            "axis_labels_in_static_figure": ["Latent PC1", "Latent PC2"],
            "q2_percentage_display": "Suppressed in the static tutorial figure because the q2 visualizer percentages are normalized within the rank-2 low-rank component and are not total CLR feature-table variance explained.",
            "correlation_annotation": "Spearman correlation of each covariate with PC1 and PC2, computed from q2-exported coordinates.",
        },
        notes=[
            "This script does not recompute the latent projection; it decodes the sample coordinates already embedded by q2-gglasso in the exported visualizations.",
            "The q2-gglasso .qzv files remain the workflow source of truth.",
        ],
    )

    print(f"wrote {PNG_OUT}")
    print(f"wrote {PDF_OUT}")
    print(f"wrote {SCORES_OUT}")
    print(f"wrote {CORR_OUT}")
    print(f"wrote {prov}")
    print(corr.sort_values("p_value").head(8).to_string(index=False))


if __name__ == "__main__":
    main()
