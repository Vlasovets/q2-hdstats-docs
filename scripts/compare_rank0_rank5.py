import colorsys
import html

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from bokeh.io import output_file, save
from bokeh.layouts import gridplot, row
from bokeh.models import (
    BasicTicker,
    ColorBar,
    ColumnDataSource,
    Div,
    HoverTool,
    LinearColorMapper,
    MultiLine,
    TabPanel,
    Tabs,
)
from bokeh.palettes import RdBu11
from bokeh.plotting import figure
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform


# Retargeted from the original hardcoded lambda=0.95 / rank-5 stems.
#
# lambda=0.95 was the EARLIER exploratory fit. Gate C1 reproduced the whole
# linear lambda path through the CLI and confirmed lambda=0.8 as the eBIC
# selection (216 edges, min eBIC 16130.0988 at gamma=0.3), and the mu1 scout
# confirmed mu1=15 -> rank 2 with 202 edges. So the canonical comparison is
# rank 0 vs rank 2 at lambda=0.8, not rank 0 vs rank 5 at lambda=0.95.
#
# Overridable from the command line so this is not hardcoded a second time:
#   compare_rank0_rank5.py [OUTPUT_STEM] [SLR_STEM] [SGL_STEM]
import sys as _sys

_ROOT = "results/gglasso/"

OUTPUT_STEM = (
    _sys.argv[1] if len(_sys.argv) > 1
    else "results/figures/atacama-top-300-rank0-rank2-comparison"
)
SLR_STEM = (
    _sys.argv[2] if len(_sys.argv) > 2
    else _ROOT + "atacama-top-300-network-slr-lambda0.8-rank2"
)
SGL_STEM = (
    _sys.argv[3] if len(_sys.argv) > 3
    else _ROOT + "atacama-top-300-network-sgl-lambda0.8"
)


def edge_key(source, target):
    return tuple(sorted((source, target)))


def collision_layout(graph, node_sizes):
    nodes = list(graph.nodes)
    raw = nx.spring_layout(
        graph,
        weight=None,
        seed=42,
        k=2.2 / np.sqrt(max(len(nodes), 1)),
        iterations=600,
    )
    raw_x = np.array([raw[node][0] for node in nodes])
    raw_y = np.array([raw[node][1] for node in nodes])
    x_span = max(raw_x.max() - raw_x.min(), 1e-12)
    y_span = max(raw_y.max() - raw_y.min(), 1e-12)
    values = np.column_stack(
        (
            70 + 860 * (raw_x - raw_x.min()) / x_span,
            70 + 710 * (raw_y - raw_y.min()) / y_span,
        )
    )
    radii = np.array([node_sizes[node] / 2 for node in nodes])
    for _ in range(800):
        max_overlap = 0.0
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                delta = values[j] - values[i]
                distance = np.linalg.norm(delta)
                minimum = radii[i] + radii[j] + 7
                if distance >= minimum:
                    continue
                if distance < 1e-12:
                    angle = (i * 37 + j * 17) * np.pi / 180
                    direction = np.array([np.cos(angle), np.sin(angle)])
                else:
                    direction = delta / distance
                overlap = minimum - distance
                values[i] -= direction * overlap / 2
                values[j] += direction * overlap / 2
                max_overlap = max(max_overlap, overlap)
        values[:, 0] = np.clip(values[:, 0], radii + 15, 985 - radii)
        values[:, 1] = np.clip(values[:, 1], radii + 15, 835 - radii)
        if max_overlap < 0.05:
            break
    return {node: values[i] for i, node in enumerate(nodes)}


def edge_source(edge_frame, positions):
    return ColumnDataSource(
        {
            "xs": [
                [positions[row.source][0], positions[row.target][0]]
                for row in edge_frame.itertuples()
            ],
            "ys": [
                [positions[row.source][1], positions[row.target][1]]
                for row in edge_frame.itertuples()
            ],
            "pair": [
                f"{row.source} - {row.target}" for row in edge_frame.itertuples()
            ],
            "partial": edge_frame["partial_correlation"].tolist(),
            "width": (
                0.8
                + 3.2
                * edge_frame["partial_correlation"].abs()
                / max(edge_frame["partial_correlation"].abs().max(), 1e-12)
            ).tolist(),
        }
    )


def network_plot(
    title,
    nodes,
    positions,
    shared_edges,
    extra_edges,
    active_nodes,
    show_extra,
    core_nodes=None,
):
    plot = figure(
        width=790,
        height=720,
        title=title,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
        x_range=(0, 1000),
        y_range=(0, 850),
        x_axis_type=None,
        y_axis_type=None,
    )
    plot.grid.visible = False

    shared_source = edge_source(shared_edges, positions)
    shared_renderer = plot.multi_line(
        xs="xs",
        ys="ys",
        source=shared_source,
        line_color="#3b78a8",
        line_width="width",
        line_alpha=0.65,
    )
    renderers = [shared_renderer]
    if show_extra:
        extra_source = edge_source(extra_edges, positions)
        extra_renderer = plot.multi_line(
            xs="xs",
            ys="ys",
            source=extra_source,
            line_color="#d05b34",
            line_width="width",
            line_dash="dashed",
            line_alpha=0.85,
        )
        renderers.append(extra_renderer)

    node_order = list(positions)
    if core_nodes is None:
        core_nodes = set()
    node_source = ColumnDataSource(
        {
            "x": [positions[node][0] for node in node_order],
            "y": [positions[node][1] for node in node_order],
            "asv": node_order,
            "genus": nodes.loc[node_order, "genus"].tolist(),
            "mean_clr": nodes.loc[node_order, "mean_clr"].tolist(),
            "degree": nodes.loc[node_order, "degree"].fillna(0).tolist(),
            "size": [
                nodes.loc[node, "node_size"] if node in active_nodes else 7
                for node in node_order
            ],
            "color": [
                nodes.loc[node, "color"] if node in active_nodes else "#d5d8dc"
                for node in node_order
            ],
            "alpha": [0.9 if node in active_nodes else 0.25 for node in node_order],
            "line_color": [
                "#111827" if node in core_nodes else "#ffffff"
                for node in node_order
            ],
            "line_width": [
                2.2 if node in core_nodes else 0.7 for node in node_order
            ],
            "membership": [
                "rank 5" if node in core_nodes else "rank 0 only"
                for node in node_order
            ],
        }
    )
    node_renderer = plot.scatter(
        x="x",
        y="y",
        source=node_source,
        size="size",
        fill_color="color",
        fill_alpha="alpha",
        line_color="line_color",
        line_width="line_width",
    )
    plot.add_tools(
        HoverTool(
            renderers=[node_renderer],
            tooltips=[
                ("ASV", "@asv"),
                ("Genus", "@genus"),
                ("Mean CLR", "@mean_clr{0.000}"),
                ("Degree", "@degree{0}"),
                ("Network membership", "@membership"),
            ],
        ),
        HoverTool(
            renderers=renderers,
            line_policy="interp",
            tooltips=[
                ("Pair", "@pair"),
                ("Partial correlation", "@partial{0.0000}"),
            ],
        ),
    )
    return plot


def off_diagonal_values(matrix):
    mask = ~np.eye(matrix.shape[0], dtype=bool)
    return np.abs(matrix[mask])


def matrix_plot(matrix, title, color_limit):
    clipped = np.clip(matrix, -color_limit, color_limit)
    mapper = LinearColorMapper(
        palette=RdBu11,
        low=-color_limit,
        high=color_limit,
    )
    plot = figure(
        width=540,
        height=500,
        title=title,
        x_range=(0, matrix.shape[0]),
        y_range=(0, matrix.shape[0]),
        x_axis_type=None,
        y_axis_type=None,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )
    plot.image(
        image=[clipped[::-1]],
        x=0,
        y=0,
        dw=matrix.shape[0],
        dh=matrix.shape[0],
        color_mapper=mapper,
    )
    plot.grid.visible = False
    plot.add_layout(
        ColorBar(
            color_mapper=mapper,
            ticker=BasicTicker(),
            label_standoff=8,
            width=12,
        ),
        "right",
    )
    return plot


def dendrogram_plot(linkage_matrix):
    data = dendrogram(linkage_matrix, no_plot=True)
    plot = figure(
        width=540,
        height=500,
        title="Hierarchical clustering of empirical covariance",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
        x_axis_type=None,
        y_axis_label="Linkage distance",
    )
    plot.multi_line(
        xs=data["icoord"],
        ys=data["dcoord"],
        line_color="#30343b",
        line_width=1.1,
        line_alpha=0.8,
    )
    plot.grid.visible = False
    return plot


slr_edges = pd.read_csv(f"{SLR_STEM}-edges.tsv", sep="\t")
sgl_edges = pd.read_csv(f"{SGL_STEM}-edges.tsv", sep="\t")
nodes = pd.read_csv(f"{SGL_STEM}-nodes.tsv", sep="\t").set_index("asv")

slr_keys = {
    edge_key(row.source, row.target) for row in slr_edges.itertuples()
}
sgl_keys = {
    edge_key(row.source, row.target) for row in sgl_edges.itertuples()
}
extra_keys = sgl_keys - slr_keys
extra_edges = sgl_edges[
    [
        edge_key(row.source, row.target) in extra_keys
        for row in sgl_edges.itertuples()
    ]
].copy()
shared_edges = sgl_edges[
    [
        edge_key(row.source, row.target) in slr_keys
        for row in sgl_edges.itertuples()
    ]
].copy()

union_graph = nx.Graph()
union_graph.add_edges_from(
    [(row.source, row.target) for row in sgl_edges.itertuples()]
)
active_rank0 = set(union_graph.nodes)
active_rank5 = set(slr_edges["source"]) | set(slr_edges["target"])

genera = sorted(nodes.loc[list(active_rank0), "genus"].unique())
colors = {
    genus: "#{:02x}{:02x}{:02x}".format(
        *[
            int(255 * channel)
            for channel in colorsys.hsv_to_rgb(i / max(len(genera), 1), 0.62, 0.82)
        ]
    )
    for i, genus in enumerate(genera)
}
nodes["color"] = nodes["genus"].map(colors)
node_sizes = nodes["node_size"].fillna(10).to_dict()
positions = collision_layout(union_graph, node_sizes)

overlay_plot = network_plot(
    "Joint overlay: 90 shared edges + 55 rank-0-only edges",
    nodes,
    positions,
    shared_edges,
    extra_edges,
    active_rank0,
    True,
    active_rank5,
)

legend_items = "".join(
    (
        "<div style='display:flex;align-items:center;gap:6px;margin:3px 0'>"
        f"<span style='width:11px;height:11px;background:{colors[genus]};"
        "display:inline-block;border-radius:50%'></span>"
        f"<span>{html.escape(genus)}</span></div>"
    )
    for genus in genera
)
legend = Div(
    width=240,
    height=720,
    text=(
        "<h3 style='margin-top:0'>Comparison</h3>"
        "<p style='font-size:12px'>Blue solid: shared edge<br>"
        "<b style='color:#d05b34'>Orange dashed: rank-0-only edge</b><br>"
        "Black node outline: connected at rank 5<br>"
        "Node size: mean CLR abundance</p>"
        "<h3>Genus</h3>"
        f"{legend_items}"
    ),
)

# Retargeted: the originals pointed at unpacked lambda=0.95 zarr directories that
# only existed in the author's working tree. Read straight out of the .qza instead,
# so this runs against the artifacts the recompute actually produces.
def _open_qza(path):
    """Open the zarr payload inside a QIIME 2 GGLassoProblem artifact."""
    import os
    import tempfile
    import zipfile

    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")]
        if not inner:
            raise SystemExit(f"{path}: no problem.zip inside the artifact")
        z.extract(inner[0], tmp)
    return zarr.open(zarr.ZipStore(os.path.join(tmp, inner[0]), mode="r"))


SGL_QZA = (
    _sys.argv[4] if len(_sys.argv) > 4
    else _ROOT + "atacama-top-300-sgl-linear-path.qza"
)
SLR_QZA = (
    _sys.argv[5] if len(_sys.argv) > 5
    else _ROOT + "atacama-top-300-slr-lambda0.8-mu15.qza"
)

sgl_solution = _open_qza(SGL_QZA)
slr_solution = _open_qza(SLR_QZA)
rank0_precision = np.asarray(sgl_solution["solution/precision_"])
rank5_sparse = np.asarray(slr_solution["solution/precision_"])
rank5_lowrank = np.asarray(slr_solution["solution/lowrank_"])
rank5_total = rank5_sparse - rank5_lowrank
empirical_covariance = np.asarray(
    slr_solution["solution/sample_covariance_"]
)

# Convert the scaled covariance (correlation) to a distance matrix, then apply
# one optimal average-linkage ordering to every matrix panel.
covariance_symmetric = (empirical_covariance + empirical_covariance.T) / 2
covariance_distance = np.sqrt(
    np.maximum(0, 2 * (1 - np.clip(covariance_symmetric, -1, 1)))
)
np.fill_diagonal(covariance_distance, 0)
linkage_matrix = linkage(
    squareform(covariance_distance, checks=False),
    method="average",
    optimal_ordering=True,
)
permutation = leaves_list(linkage_matrix)
empirical_covariance = empirical_covariance[np.ix_(permutation, permutation)]
rank0_precision = rank0_precision[np.ix_(permutation, permutation)]
rank5_total = rank5_total[np.ix_(permutation, permutation)]
rank5_sparse = rank5_sparse[np.ix_(permutation, permutation)]
rank5_lowrank = rank5_lowrank[np.ix_(permutation, permutation)]

precision_limit = np.quantile(
    np.concatenate(
        [
            off_diagonal_values(rank0_precision),
            off_diagonal_values(rank5_total),
            off_diagonal_values(rank5_sparse),
        ]
    ),
    0.99,
)
lowrank_limit = np.quantile(off_diagonal_values(rank5_lowrank), 0.99)

covariance_matrix_plot = matrix_plot(
    empirical_covariance,
    "Empirical covariance (correlation)",
    1.0,
)
rank0_matrix_plot = matrix_plot(
    rank0_precision,
    "Rank 0 inverse covariance",
    precision_limit,
)
rank5_total_plot = matrix_plot(
    rank5_total,
    "Rank 5 total inverse covariance = sparse - low-rank",
    precision_limit,
)
rank5_sparse_plot = matrix_plot(
    rank5_sparse,
    "Rank 5 sparse component",
    precision_limit,
)
rank5_lowrank_plot = matrix_plot(
    rank5_lowrank,
    "Rank 5 low-rank component",
    lowrank_limit,
)
clustering_plot = dendrogram_plot(linkage_matrix)
matrix_note = Div(
    width=1080,
    text=(
        "<p><b>Matrix ordering:</b> average-linkage hierarchical clustering "
        "of the empirical covariance, using correlation distance "
        "<code>sqrt(2(1-C))</code>. "
        "The three precision panels share one diverging color scale, clipped "
        "at the 99th percentile of absolute off-diagonal values. The low-rank "
        "panel uses its own scale.</p>"
    ),
)

network_panel = TabPanel(
    child=row(overlay_plot, legend),
    title="Joint network",
)
matrix_panel = TabPanel(
    child=gridplot(
        [
            [matrix_note, None],
            [covariance_matrix_plot, clustering_plot],
            [rank0_matrix_plot, rank5_total_plot],
            [rank5_sparse_plot, rank5_lowrank_plot],
        ]
    ),
    title="Inverse covariance",
)

output_file(f"{OUTPUT_STEM}.html", title="Rank-0 versus rank-5 comparison")
save(Tabs(tabs=[network_panel, matrix_panel]))

labels = [f"ASV-{i + 1}" for i in permutation]
for name, matrix in {
    "empirical-covariance": empirical_covariance,
    "rank0-precision": rank0_precision,
    "slr-total-precision": rank5_total,
    "slr-sparse": rank5_sparse,
    "slr-lowrank": rank5_lowrank,
}.items():
    pd.DataFrame(matrix, index=labels, columns=labels).to_csv(
        f"{OUTPUT_STEM}-{name}.tsv",
        sep="\t",
    )

figure, axes = plt.subplots(3, 2, figsize=(12, 15), constrained_layout=True)
matrix_specs = [
    (empirical_covariance, "Empirical covariance (correlation)", 1.0),
    (rank0_precision, "Rank 0 inverse covariance", precision_limit),
    (
        rank5_total,
        "Rank 5 total inverse covariance\n(sparse - low-rank)",
        precision_limit,
    ),
    (rank5_sparse, "Rank 5 sparse component", precision_limit),
    (rank5_lowrank, "Rank 5 low-rank component", lowrank_limit),
]
for axis, (matrix, title, limit) in zip(axes.flat[:5], matrix_specs):
    image = axis.imshow(
        np.clip(matrix, -limit, limit),
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        interpolation="nearest",
        aspect="equal",
    )
    axis.set_title(title)
    axis.set_xticks([])
    axis.set_yticks([])
    figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
axes.flat[5].set_title("Hierarchical clustering of empirical covariance")
dendrogram(
    linkage_matrix,
    ax=axes.flat[5],
    no_labels=True,
    color_threshold=0,
    above_threshold_color="#30343b",
)
axes.flat[5].set_xlabel("Clustered ASVs")
axes.flat[5].set_ylabel("Linkage distance")
figure.suptitle(
    "Atacama covariance and precision comparison\n"
    "(hierarchical covariance ordering)",
    fontsize=14,
)
figure.savefig(f"{OUTPUT_STEM}-matrices.png", dpi=200)
plt.close(figure)

network_figure, network_axis = plt.subplots(figsize=(14, 11))
network_axis.set_title(
    "Rank-0 versus rank-5 network overlay\n"
    "solid blue = shared; dashed orange = rank-0 only",
    fontsize=15,
)
network_axis.axis("off")
nx.draw_networkx_edges(
    union_graph,
    positions,
    edgelist=[
        (row.source, row.target) for row in shared_edges.itertuples()
    ],
    edge_color="#3b78a8",
    width=1.4,
    alpha=0.55,
    ax=network_axis,
)
nx.draw_networkx_edges(
    union_graph,
    positions,
    edgelist=[
        (row.source, row.target) for row in extra_edges.itertuples()
    ],
    edge_color="#d05b34",
    width=2.2,
    alpha=0.9,
    style="dashed",
    ax=network_axis,
)
rank0_only_nodes = sorted(active_rank0 - active_rank5)
nx.draw_networkx_nodes(
    union_graph,
    positions,
    nodelist=rank0_only_nodes,
    node_color=[nodes.loc[node, "color"] for node in rank0_only_nodes],
    node_size=[
        nodes.loc[node, "node_size"] * 12 for node in rank0_only_nodes
    ],
    edgecolors="#ffffff",
    linewidths=0.8,
    alpha=0.9,
    ax=network_axis,
)
rank5_nodes = sorted(active_rank5)
nx.draw_networkx_nodes(
    union_graph,
    positions,
    nodelist=rank5_nodes,
    node_color=[nodes.loc[node, "color"] for node in rank5_nodes],
    node_size=[
        nodes.loc[node, "node_size"] * 12 for node in rank5_nodes
    ],
    edgecolors="#111827",
    linewidths=1.8,
    alpha=0.95,
    ax=network_axis,
)
network_figure.savefig(
    f"{OUTPUT_STEM}-network.png",
    dpi=180,
    bbox_inches="tight",
)
plt.close(network_figure)

pd.DataFrame(
    {
        "metric": [
            "slr_edges",
            "rank0_edges",
            "shared_edges",
            "rank0_only_edges",
            "slr_connected_nodes",
            "rank0_connected_nodes",
        ],
        "value": [
            len(slr_keys),
            len(sgl_keys),
            len(slr_keys & sgl_keys),
            len(extra_keys),
            len(active_rank5),
            len(active_rank0),
        ],
    }
).to_csv(f"{OUTPUT_STEM}-summary.tsv", sep="\t", index=False)
