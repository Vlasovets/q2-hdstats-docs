import argparse
import colorsys
import html

import networkx as nx
import numpy as np
import pandas as pd
import zarr
from biom import load_table
from bokeh.io import output_file, save
from bokeh.layouts import row
from bokeh.models import ColumnDataSource, Div, HoverTool, MultiLine
from bokeh.plotting import figure


parser = argparse.ArgumentParser()
parser.add_argument(
    "--solution-zarr",
    default="work/atacama-full/lambda0.95-mu10.5-zarr",
)
parser.add_argument(
    "--output-stem",
    default=(
        "outputs/atacama-q2-gglasso/"
        "atacama-top-300-network-lambda0.95-mu10.5"
    ),
)
parser.add_argument(
    "--title",
    default="Sparse + low-rank network: lambda=0.95, mu=10.5",
)
args = parser.parse_args()
OUTPUT = f"{args.output_stem}.html"


def genus_from_taxonomy(value):
    levels = [part.strip() for part in str(value).split(";")]
    genera = [part[3:] for part in levels if part.startswith("g__") and len(part) > 3]
    return genera[-1] if genera else "Unassigned"


frequency = load_table(
    "work/atacama-full/top-300-frequency-export/feature-table.biom"
).to_dataframe(dense=True)
clr = load_table(
    "work/atacama-full/top-300-clr-export/feature-table.biom"
).to_dataframe(dense=True)
taxonomy = pd.read_csv(
    "work/atacama-full/taxonomy-silva138-export/taxonomy.tsv",
    sep="\t",
    index_col=0,
)

# q2-gglasso renames features ASV-1...ASV-p after sorting by total count.
feature_order = frequency.sum(axis=1).sort_values(kind="stable").index
mapping = pd.DataFrame(
    {
        "asv": [f"ASV-{i + 1}" for i in range(len(feature_order))],
        "feature_id": feature_order,
        "total_abundance": frequency.loc[feature_order].sum(axis=1).to_numpy(),
    }
).set_index("asv")
mapping["genus"] = [
    genus_from_taxonomy(taxonomy.loc[feature_id, "Taxon"])
    if feature_id in taxonomy.index
    else "Unassigned"
    for feature_id in mapping["feature_id"]
]
mapping["mean_clr"] = clr.mean(axis=0).reindex(mapping.index)

solution = zarr.open(
    args.solution_zarr, mode="r"
)
precision = np.asarray(solution["solution/precision_"])
adjacency = np.asarray(solution["solution/adjacency_"]).astype(bool)
diagonal = np.diag(precision)
partial = -precision / np.sqrt(np.outer(diagonal, diagonal))
np.fill_diagonal(partial, 0)

graph = nx.Graph()
for i, asv in enumerate(mapping.index):
    graph.add_node(asv, matrix_index=i)
for i, j in zip(*np.where(np.triu(adjacency, 1))):
    graph.add_edge(
        mapping.index[i],
        mapping.index[j],
        partial_correlation=float(partial[i, j]),
    )

connected = [node for node, degree in graph.degree() if degree > 0]
network = graph.subgraph(connected).copy()

clr_values = mapping.loc[connected, "mean_clr"]
clr_range = clr_values.max() - clr_values.min()
if clr_range == 0:
    sizes = pd.Series(18.0, index=connected)
else:
    sizes = 10 + 28 * (clr_values - clr_values.min()) / clr_range

raw_positions = nx.spring_layout(
    network,
    weight=None,
    seed=42,
    k=2.2 / np.sqrt(max(network.number_of_nodes(), 1)),
    iterations=500,
)

# Convert to pixel-like plot coordinates, then enforce separation based on
# rendered node diameters. This keeps the force-directed topology readable
# without allowing large abundance-scaled nodes to overlap.
raw_x = np.array([raw_positions[node][0] for node in connected])
raw_y = np.array([raw_positions[node][1] for node in connected])
x_span = max(raw_x.max() - raw_x.min(), 1e-12)
y_span = max(raw_y.max() - raw_y.min(), 1e-12)
position_array = np.column_stack(
    (
        70 + 860 * (raw_x - raw_x.min()) / x_span,
        70 + 710 * (raw_y - raw_y.min()) / y_span,
    )
)
radii = sizes.loc[connected].to_numpy() / 2

for _ in range(600):
    max_overlap = 0.0
    for i in range(len(connected)):
        for j in range(i + 1, len(connected)):
            delta = position_array[j] - position_array[i]
            distance = np.linalg.norm(delta)
            minimum = radii[i] + radii[j] + 8
            if distance >= minimum:
                continue
            if distance < 1e-12:
                angle = (i * 37 + j * 17) * np.pi / 180
                direction = np.array([np.cos(angle), np.sin(angle)])
            else:
                direction = delta / distance
            overlap = minimum - distance
            position_array[i] -= direction * overlap / 2
            position_array[j] += direction * overlap / 2
            max_overlap = max(max_overlap, overlap)
    position_array[:, 0] = np.clip(position_array[:, 0], radii + 15, 985 - radii)
    position_array[:, 1] = np.clip(position_array[:, 1], radii + 15, 835 - radii)
    if max_overlap < 0.05:
        break

positions = {
    node: position_array[i] for i, node in enumerate(connected)
}

genera = sorted(mapping.loc[connected, "genus"].unique())
colors = {
    genus: "#{:02x}{:02x}{:02x}".format(
        *[
            int(255 * channel)
            for channel in colorsys.hsv_to_rgb(i / max(len(genera), 1), 0.62, 0.82)
        ]
    )
    for i, genus in enumerate(genera)
}

node_source = ColumnDataSource(
    {
        "x": [positions[node][0] for node in connected],
        "y": [positions[node][1] for node in connected],
        "asv": connected,
        "feature_id": mapping.loc[connected, "feature_id"].tolist(),
        "genus": mapping.loc[connected, "genus"].tolist(),
        "mean_clr": mapping.loc[connected, "mean_clr"].tolist(),
        "total_abundance": mapping.loc[connected, "total_abundance"].tolist(),
        "degree": [network.degree(node) for node in connected],
        "size": sizes.loc[connected].tolist(),
        "color": [colors[mapping.loc[node, "genus"]] for node in connected],
    }
)

edge_xs, edge_ys, edge_values, edge_colors, edge_widths, edge_names = (
    [],
    [],
    [],
    [],
    [],
    [],
)
edge_magnitudes = [
    abs(data["partial_correlation"]) for _, _, data in network.edges(data=True)
]
edge_max = max(edge_magnitudes, default=1)
for source, target, data in network.edges(data=True):
    value = data["partial_correlation"]
    edge_xs.append([positions[source][0], positions[target][0]])
    edge_ys.append([positions[source][1], positions[target][1]])
    edge_values.append(value)
    edge_colors.append("#2474b5" if value > 0 else "#c43c39")
    edge_widths.append(0.8 + 4.2 * abs(value) / edge_max)
    edge_names.append(f"{source} - {target}")

edge_source = ColumnDataSource(
    {
        "xs": edge_xs,
        "ys": edge_ys,
        "partial": edge_values,
        "color": edge_colors,
        "width": edge_widths,
        "pair": edge_names,
    }
)

plot = figure(
    width=1000,
    height=850,
    title=(
        f"{args.title} | "
        f"{network.number_of_nodes()} connected ASVs, "
        f"{network.number_of_edges()} edges"
    ),
    tools="pan,wheel_zoom,box_zoom,reset,save",
    active_scroll="wheel_zoom",
    x_range=(0, 1000),
    y_range=(0, 850),
    x_axis_type=None,
    y_axis_type=None,
)
plot.grid.visible = False
edge_renderer = plot.multi_line(
    xs="xs",
    ys="ys",
    source=edge_source,
    line_color="color",
    line_width="width",
    line_alpha=0.65,
)
node_renderer = plot.scatter(
    x="x",
    y="y",
    source=node_source,
    size="size",
    fill_color="color",
    line_color="#ffffff",
    line_width=0.8,
    fill_alpha=0.9,
)
plot.add_tools(
    HoverTool(
        renderers=[node_renderer],
        tooltips=[
            ("ASV", "@asv"),
            ("Genus", "@genus"),
            ("Mean CLR", "@mean_clr{0.000}"),
            ("Total count", "@total_abundance{0}"),
            ("Degree", "@degree"),
            ("Feature ID", "@feature_id"),
        ],
    ),
    HoverTool(
        renderers=[edge_renderer],
        line_policy="interp",
        tooltips=[
            ("Pair", "@pair"),
            ("Partial correlation", "@partial{0.0000}"),
        ],
    ),
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
    height=850,
    text=(
        "<h3 style='margin-top:0'>Genus</h3>"
        "<p style='font-size:12px'>Blue edges: positive partial correlation<br>"
        "Red edges: negative partial correlation<br>"
        "Node size: mean CLR abundance</p>"
        f"{legend_items}"
    ),
)

output_file(OUTPUT, title="Atacama sparse genus network")
save(row(plot, legend))

pd.DataFrame(
    {
        "source": [source for source, _, _ in network.edges(data=True)],
        "target": [target for _, target, _ in network.edges(data=True)],
        "partial_correlation": [
            data["partial_correlation"] for _, _, data in network.edges(data=True)
        ],
    }
).to_csv(
    f"{args.output_stem}-edges.tsv",
    sep="\t",
    index=False,
)
mapping["layout_x"] = pd.Series(
    {node: positions[node][0] for node in connected}
)
mapping["layout_y"] = pd.Series(
    {node: positions[node][1] for node in connected}
)
mapping["node_size"] = sizes
mapping["degree"] = pd.Series(
    {node: network.degree(node) for node in connected}
)
mapping.to_csv(
    f"{args.output_stem}-nodes.tsv",
    sep="\t",
)
