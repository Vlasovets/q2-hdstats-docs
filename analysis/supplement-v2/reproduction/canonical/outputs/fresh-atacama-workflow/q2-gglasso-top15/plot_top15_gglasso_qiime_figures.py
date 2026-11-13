from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path("/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial")
WORKFLOW_DIR = PROJECT_ROOT / "outputs/fresh-atacama-workflow"
OUT_DIR = WORKFLOW_DIR / "q2-gglasso-top15"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402


INPUTS = {
    "mclr": OUT_DIR / "atacama-top15-mclr-correlation-asv-order.tsv",
    "sgl": OUT_DIR / "atacama-top15-sgl-gamma001-partial-correlation-asv-order.tsv",
    "mixed": OUT_DIR / "atacama-top15-mixed-weight005-gamma001-full-asv-covariate-partial-correlation.tsv",
    "slr": OUT_DIR / "atacama-top15-slr-gamma001-sparse-plus-lv1-augmented-matrix.tsv",
    "summary": OUT_DIR / "atacama-top15-q2-gglasso-model-selection-summary.tsv",
}

PNG_OUT = OUT_DIR / "atacama-top15-q2-gglasso-four-matrix-comparison-final-tight.png"
PDF_OUT = OUT_DIR / "atacama-top15-q2-gglasso-four-matrix-comparison-final-tight.pdf"
PNG_SPACING_ALIAS = OUT_DIR / "atacama-top15-q2-gglasso-four-matrix-comparison-final-spacing.png"
PDF_SPACING_ALIAS = OUT_DIR / "atacama-top15-q2-gglasso-four-matrix-comparison-final-spacing.pdf"

COVARIATE_TO_CODE = {
    "elevation": "elevation",
    "pH": "pH",
    "soil humidity": "soil RH",
    "soil temperature": "soil temp.",
    "percent cover": "cover",
    "log10 extract conc.": "log10\nextract",
}


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def coded_covariates(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix.rename(index=COVARIATE_TO_CODE, columns=COVARIATE_TO_CODE)


def slr_with_lv_rowcol(slr_augmented: pd.DataFrame) -> pd.DataFrame:
    asvs = [idx for idx in slr_augmented.index if idx != "LV1"]
    return slr_augmented.loc[asvs + ["LV1"], asvs + ["LV1"]].copy()


def add_grid(ax, nrows: int, ncols: int) -> None:
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="#f0f0f0", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)


def draw_matrix(
    ax,
    matrix: pd.DataFrame,
    *,
    title: str,
    vmin: float,
    vmax: float,
    title_size: float = 8.5,
    tick_size: float = 5.6,
    boundary: int | None = None,
):
    arr = matrix.to_numpy(dtype=float)
    image = ax.imshow(arr, cmap="RdBu_r", vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title, fontsize=title_size, pad=7)
    ax.set_aspect("equal")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=tick_size)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=tick_size)
    ax.tick_params(length=0, pad=1.2)
    add_grid(ax, matrix.shape[0], matrix.shape[1])
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.75)
    if boundary is not None:
        line = boundary - 0.5
        ax.axhline(line, color="#333333", linewidth=1.1)
        ax.axvline(line, color="#333333", linewidth=1.1)
    return image


def axes_from_cells(
    *,
    left: float,
    bottom: float,
    ncols: int,
    nrows: int,
    cell_width: float,
    fig_width: float,
    fig_height: float,
) -> list[float]:
    """Return axes rectangle with one shared physical cell size."""
    width = ncols * cell_width / fig_width
    height = nrows * cell_width / fig_height
    return [left, bottom, width, height]


def main() -> None:
    mclr = read_matrix(INPUTS["mclr"])
    sgl = read_matrix(INPUTS["sgl"])
    mixed = coded_covariates(read_matrix(INPUTS["mixed"]))
    slr = slr_with_lv_rowcol(read_matrix(INPUTS["slr"]))

    fig_width, fig_height = 16.2, 5.95
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)
    fig.suptitle(
        "Atacama top-15 ASVs: empirical mCLR correlation and selected q2-gglasso graph matrices",
        fontsize=9.2,
        y=0.895,
    )

    cell_width = 0.176
    asv_bottom = 0.30
    mixed_bottom = asv_bottom - 6 * cell_width / fig_height
    slr_bottom = asv_bottom - cell_width / fig_height

    ax_mclr = fig.add_axes(
        axes_from_cells(
            left=0.048,
            bottom=asv_bottom,
            ncols=15,
            nrows=15,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_cbar_mclr = fig.add_axes([0.228, asv_bottom, 0.010, 15 * cell_width / fig_height])
    ax_sgl = fig.add_axes(
        axes_from_cells(
            left=0.318,
            bottom=asv_bottom,
            ncols=15,
            nrows=15,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_mixed = fig.add_axes(
        axes_from_cells(
            left=0.522,
            bottom=mixed_bottom,
            ncols=21,
            nrows=21,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_slr = fig.add_axes(
        axes_from_cells(
            left=0.785,
            bottom=slr_bottom,
            ncols=16,
            nrows=16,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_cbar_graph = fig.add_axes([0.965, slr_bottom, 0.010, 16 * cell_width / fig_height])

    img_mclr = draw_matrix(
        ax_mclr,
        mclr,
        title="Empirical mCLR correlation",
        vmin=-1,
        vmax=1,
    )
    img_sgl = draw_matrix(
        ax_sgl,
        sgl,
        title="SGL partial correlations",
        vmin=-0.35,
        vmax=0.35,
    )
    draw_matrix(
        ax_mixed,
        mixed,
        title="Mixed graph partial correlations",
        vmin=-0.35,
        vmax=0.35,
        tick_size=4.9,
        boundary=15,
    )
    img_slr = draw_matrix(
        ax_slr,
        slr,
        title="Sparse partial correlation + latent variable (LV1)",
        vmin=-0.35,
        vmax=0.35,
    )

    cbar_mclr = fig.colorbar(img_mclr, cax=ax_cbar_mclr)
    cbar_mclr.set_label("mCLR correlation", rotation=90, labelpad=8, fontsize=7.5)
    cbar_mclr.ax.tick_params(labelsize=6.2, length=2)

    cbar_graph = fig.colorbar(img_slr, cax=ax_cbar_graph)
    cbar_graph.set_label("signed partial correlation / LV1 loading", rotation=90, labelpad=14, fontsize=7.5)
    cbar_graph.ax.tick_params(labelsize=6.2, length=2)

    fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.035)
    shutil.copyfile(PNG_OUT, PNG_SPACING_ALIAS)
    shutil.copyfile(PDF_OUT, PDF_SPACING_ALIAS)
    plt.close(fig)

    write_provenance(
        PNG_OUT,
        description="Render the S8 top-15 q2-gglasso four-panel matrix comparison from QIIME-derived TSV extracts using fixed ASV-ASV cell geometry.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=INPUTS.values(),
        outputs=[PNG_OUT, PDF_OUT, PNG_SPACING_ALIAS, PDF_SPACING_ALIAS],
        parameters={
            "figure_size_inches": [16.2, 5.95],
            "dpi": 300,
            "savefig": {"bbox_inches": "tight", "pad_inches": 0.035},
            "shared_asv_cell_width_inches": cell_width,
            "mclr_correlation_color_scale": [-1, 1],
            "graph_color_scale": [-0.35, 0.35],
            "mixed_model_covariate_labels": COVARIATE_TO_CODE,
            "lv1_display": "LV1 is appended as both a row and a column; the ASV-ASV block keeps the same cell size and vertical alignment as the 15 x 15 panels.",
        },
        notes=[
            "The final-tight filename is kept because the tutorial figure registry points to it.",
            "The final-spacing aliases are also written for compatibility with the accepted S8 naming history.",
            "The mCLR, SGL, mixed-model ASV block, and sparse+latent ASV block use the same physical cell size.",
            "All matrix values originate from QIIME 2 .qza extraction, not from a direct solver run.",
        ],
    )

    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {PNG_SPACING_ALIAS}")
    print(f"Wrote {PDF_SPACING_ALIAS}")


if __name__ == "__main__":
    main()
