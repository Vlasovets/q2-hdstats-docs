from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


PROJECT_ROOT = Path("/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial")
WORKFLOW_DIR = PROJECT_ROOT / "outputs/fresh-atacama-workflow"
OUT = WORKFLOW_DIR / "q2-gglasso-top300"
PROVENANCE_UTILS = WORKFLOW_DIR / "tutorial-figure-generation"

sys.path.insert(0, str(PROVENANCE_UTILS))
from provenance_utils import write_provenance  # noqa: E402


INPUTS = {
    "active": OUT / "atacama-top300-union-active-asvs-sgl-mixed-slr.tsv",
    "correlation": OUT / "atacama-top300-clr-correlation-asv-order.tsv",
    "sgl": OUT / "atacama-top300-sgl-gamma025-partial-correlation-asv-order.tsv",
    "mixed_full": OUT / "atacama-top300-mixed-gamma025-full-asv-covariate-partial-correlation.tsv",
    "slr_augmented": OUT / "atacama-top300-slr-gamma025-sparse-plus-lv12-augmented-matrix.tsv",
    "summary": OUT / "atacama-top300-gglasso-sgl-mixed-slr-summary.tsv",
}

PNG_OUT = OUT / "atacama-top300-gglasso-reduced-active-asv-augmented-matrix-comparison.png"
PDF_OUT = OUT / "atacama-top300-gglasso-reduced-active-asv-augmented-matrix-comparison.pdf"
ORDER_OUT = OUT / "atacama-top300-gglasso-reduced-active-asv-cluster-order.tsv"

COVARIATE_LABELS = ["elevation", "pH", "soil RH", "soil temp.", "cover", "log10 extract"]
LV_LABELS = ["LV1", "LV2"]


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def cluster_order(correlation: pd.DataFrame) -> list[str]:
    values = correlation.to_numpy(dtype=float)
    clipped = np.clip(values, -1.0, 1.0)
    distance = 1.0 - np.abs(clipped)
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    linked = linkage(condensed, method="average")
    return correlation.index[leaves_list(linked)].tolist()


def axes_from_cells(
    *,
    left: float,
    top: float,
    ncols: int,
    nrows: int,
    cell_width: float,
    fig_width: float,
    fig_height: float,
) -> list[float]:
    width = ncols * cell_width / fig_width
    height = nrows * cell_width / fig_height
    bottom = top - height
    return [left, bottom, width, height]


def draw_heatmap(
    ax,
    matrix: pd.DataFrame,
    *,
    title: str,
    vmin: float,
    vmax: float,
    tick_labels: list[str] | None = None,
    boundary: int | None = None,
):
    image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title, fontsize=7.8, pad=4)
    ax.set_aspect("equal")

    if tick_labels:
        tick_positions = [matrix.index.get_loc(label) for label in tick_labels if label in matrix.index]
        visible_labels = [label for label in tick_labels if label in matrix.index]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(visible_labels, rotation=90, fontsize=5.2)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(visible_labels, fontsize=5.2)
        ax.tick_params(length=0, pad=1.0)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(0.55)
        spine.set_color("#333333")

    if boundary is not None:
        line = boundary - 0.5
        ax.axhline(line, color="#333333", linewidth=0.65)
        ax.axvline(line, color="#333333", linewidth=0.65)
    return image


def main() -> None:
    active = pd.read_csv(INPUTS["active"], sep="\t")["asv_label"].tolist()
    corr_all = read_matrix(INPUTS["correlation"])
    corr_active = corr_all.loc[active, active]
    order = cluster_order(corr_active)
    n_asv = len(order)

    corr = corr_all.loc[order, order]
    sgl = read_matrix(INPUTS["sgl"]).loc[order, order]
    mixed_full = read_matrix(INPUTS["mixed_full"])
    mixed = mixed_full.loc[order + COVARIATE_LABELS, order + COVARIATE_LABELS]
    slr_augmented_full = read_matrix(INPUTS["slr_augmented"])
    slr = slr_augmented_full.loc[order + LV_LABELS, order + LV_LABELS]

    pd.DataFrame({"cluster_order": np.arange(1, n_asv + 1), "asv_label": order}).to_csv(
        ORDER_OUT, sep="\t", index=False
    )

    sgl_vmax = max(0.05, float(np.nanmax(np.abs(sgl.to_numpy(dtype=float)))))
    other_vmax = max(
        0.05,
        float(np.nanmax(np.abs(mixed.to_numpy(dtype=float)))),
        float(np.nanmax(np.abs(slr.to_numpy(dtype=float)))),
    )

    fig_width, fig_height = 8.27, 8.8
    cell_width = 0.01335
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=450)
    fig.suptitle(
        f"Atacama top-300 q2-gglasso networks on the {n_asv} edge-active ASVs",
        fontsize=9.4,
        y=0.965,
    )

    left_col = 0.070
    right_col = 0.555
    top_row_top = 0.890
    bottom_row_top = 0.500
    asv_height = n_asv * cell_width / fig_height

    ax_corr = fig.add_axes(
        axes_from_cells(
            left=left_col,
            top=top_row_top,
            ncols=n_asv,
            nrows=n_asv,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_sgl = fig.add_axes(
        axes_from_cells(
            left=right_col,
            top=top_row_top,
            ncols=n_asv,
            nrows=n_asv,
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_mixed = fig.add_axes(
        axes_from_cells(
            left=left_col,
            top=bottom_row_top,
            ncols=n_asv + len(COVARIATE_LABELS),
            nrows=n_asv + len(COVARIATE_LABELS),
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )
    ax_slr = fig.add_axes(
        axes_from_cells(
            left=right_col,
            top=bottom_row_top,
            ncols=n_asv + len(LV_LABELS),
            nrows=n_asv + len(LV_LABELS),
            cell_width=cell_width,
            fig_width=fig_width,
            fig_height=fig_height,
        )
    )

    ax_cbar_corr = fig.add_axes([0.448, top_row_top - asv_height, 0.012, asv_height])
    ax_cbar_sgl = fig.add_axes([0.940, top_row_top - asv_height, 0.012, asv_height])
    ax_cbar_other = fig.add_axes(
        [
            0.940,
            bottom_row_top - ((n_asv + len(LV_LABELS)) * cell_width / fig_height),
            0.012,
            (n_asv + len(LV_LABELS)) * cell_width / fig_height,
        ]
    )

    img_corr = draw_heatmap(
        ax_corr,
        corr,
        title="Empirical CLR correlation",
        vmin=-1,
        vmax=1,
    )
    img_sgl = draw_heatmap(
        ax_sgl,
        sgl,
        title="SGL partial correlations",
        vmin=-sgl_vmax,
        vmax=sgl_vmax,
    )
    draw_heatmap(
        ax_mixed,
        mixed,
        title="Mixed graph partial correlations",
        vmin=-other_vmax,
        vmax=other_vmax,
        boundary=n_asv,
    )
    img_graph = draw_heatmap(
        ax_slr,
        slr,
        title="Sparse partial correlation + latent variables (LV1, LV2)",
        vmin=-other_vmax,
        vmax=other_vmax,
        boundary=n_asv,
    )

    fig.text(
        left_col,
        bottom_row_top - ((n_asv + len(COVARIATE_LABELS)) * cell_width / fig_height) - 0.028,
        "covariates appended: " + ", ".join(COVARIATE_LABELS),
        fontsize=5.8,
        ha="left",
        va="center",
    )
    fig.text(
        right_col,
        bottom_row_top - ((n_asv + len(LV_LABELS)) * cell_width / fig_height) - 0.028,
        "latent variables appended: " + ", ".join(LV_LABELS),
        fontsize=5.8,
        ha="left",
        va="center",
    )

    cbar_corr = fig.colorbar(img_corr, cax=ax_cbar_corr)
    cbar_corr.set_label("CLR correlation", rotation=90, labelpad=7, fontsize=7.0)
    cbar_corr.ax.tick_params(labelsize=5.8, length=2)

    cbar_sgl = fig.colorbar(img_sgl, cax=ax_cbar_sgl)
    cbar_sgl.set_label("SGL partial correlation", rotation=90, labelpad=10, fontsize=7.0)
    cbar_sgl.ax.tick_params(labelsize=5.8, length=2)

    cbar_graph = fig.colorbar(img_graph, cax=ax_cbar_other)
    cbar_graph.set_label("mixed partial correlation / LV loading", rotation=90, labelpad=10, fontsize=7.0)
    cbar_graph.ax.tick_params(labelsize=5.8, length=2)

    fig.savefig(PNG_OUT, dpi=450, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)

    write_provenance(
        PNG_OUT,
        description="Plot clustered top-300 q2-gglasso matrix comparison with fixed ASV-ASV geometry and appended covariate/LV rows and columns.",
        project_root=PROJECT_ROOT,
        script=Path(__file__),
        inputs=INPUTS.values(),
        outputs=[PNG_OUT, PDF_OUT, ORDER_OUT],
        parameters={
            "active_asv_definition": "Union of ASVs with at least one ASV-ASV edge in SGL, mixed, or SLR.",
            "n_active_asvs": n_asv,
            "appended_mixed_covariates": COVARIATE_LABELS,
            "appended_slr_latent_variables": LV_LABELS,
            "appended_label_display": "Exact appended labels are written as panel annotations because individual tick labels overlap at the 228-ASV cell scale.",
            "cluster_distance": "1 - abs(empirical CLR correlation), average linkage",
            "same_order_all_panels": True,
            "layout": "A4 portrait 2x2 with tightened vertical panel spacing",
            "same_physical_asv_cell_size_all_panels": True,
            "shared_asv_cell_width_inches": cell_width,
            "sgl_partial_correlation_color_scale": [-sgl_vmax, sgl_vmax],
            "mixed_slr_partial_correlation_lv_color_scale": [-other_vmax, other_vmax],
            "correlation_color_scale": [-1, 1],
        },
        notes=[
            "The plotted matrices are QIIME-derived TSV extracts from q2-gglasso artifacts.",
            "Rows and columns are reduced to the active ASV union and ordered by clustering the empirical CLR correlation matrix.",
            "The ASV-ASV block has the same physical cell size in all four panels; covariates and latent variables are appended below and to the right.",
        ],
    )

    print(f"active ASVs: {n_asv}")
    print(f"SGL color scale: {-sgl_vmax:.4g} to {sgl_vmax:.4g}")
    print(f"mixed/SLR color scale: {-other_vmax:.4g} to {other_vmax:.4g}")
    print(f"wrote {PNG_OUT}")
    print(f"wrote {PDF_OUT}")


if __name__ == "__main__":
    main()
