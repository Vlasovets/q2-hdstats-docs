# Figure audit — 2026-08-14

Reviewed every figure in the book against two questions:

1. **Can a reader produce this by running a documented command?** A figure that
   illustrates this repository's own plotting code, rather than plugin output,
   teaches the reader something they cannot repeat.
2. **Does it follow the [cnsplots](https://github.com/faridrashidi/cnsplots)
   conventions?** Specifically: no legend over the plotting area, no label
   colliding with data, no hardcoded style that contradicts the house defaults.

The book went from **27 figures to 22**.

## What the plugins can actually draw

This is the finding that drove the removals, and it is checkable from source.

`q2-gglasso` registers **six actions, of which two are visualizers**:

| action | kind |
|---|---|
| `transform-features` | method |
| `build-groups` | method |
| `calculate-covariance` | method |
| `solve-problem` | method |
| `pca` | **visualizer** |
| `summarize` | **visualizer** |

`q2-classo` registers **eight actions, of which one is a visualizer**
(`summarize`; the rest are methods).

Reading `q2_gglasso/_summarize/_visualizer.py` and `q2_gglasso/_pca/_visualizer.py`:
both build **bokeh** output — heatmaps, scatter plots, data tables and hover
tooltips. Neither imports `networkx` or any graph-layout library, and neither
emits a node-link diagram.

**So no `qiime` command in either plugin draws a network.** The book's two network
figures came from `analysis/scripts/visualize_sparse_network.py`, a script in this
repository. They showed real results, but they showed them through a rendering
path a reader has no access to.

## Removed

| Was | Figure | Why |
|---|---|---|
| Fig 7 | `atacama-network-rank0-vs-rank2.png` in the tutorial | No plugin command produces a network drawing. |
| Fig 23 | the same figure in `04_highdim_atacama/06_interpretation` | As above. |
| Fig 25 | `pipeline_overview_placeholder.png` | Self-declared placeholder; caption began "**Placeholder.**". |
| Fig 26 | `compositional_regression_placeholder.png` | Placeholder, and duplicated the simplex figure that already exists. |
| Fig 27 | `toy_vs_hd_network_placeholder.png` | Placeholder, and would have been a third network drawing. |

Nothing quantitative was lost. Both network captions carried edge-set facts —
202 shared, 14 removed by the rank-2 latent block, 0 added; the removed edges
falling into three components of 6, 4 and 2 nodes — and those facts are now stated
as prose in the chapters that used to hold the figures. `qiime gglasso summarize`
reports the per-solution edge sets, so a reader can verify them.

`fig_network_rank0_vs_rank2()` was deleted from
`analysis/scripts/make_docs_figures.py` as well; leaving it would have regenerated
a figure the book no longer references. The script's docstring records why.

## Kept, and why

Everything remaining is either **plugin output** or a **concept diagram**:

- `atacama-ebic-lambda-selection`, `atacama-rank-tradeoff`,
  `toy-lambda-path-gamma` — read directly from `solve-problem` results committed
  under `analysis/results/tables/`.
- `atacama-top-300-*` heatmaps, `scatter_pc` — the `summarize` and `pca`
  visualizers' own output.
- `classo_reg`, `classo_class`, `reg`, `slc_fig`, `reg_tree`, `slr_example`,
  `example_gglasso`, `overview`, `ph` — schematics and `.qzv` captures.
- `compositional-simplex-zero-sum` — the one concept every chapter depends on.

## Style fixes

Three defects found by rendering the figures and measuring the output, not by
reading the code.

**Two in-plot legends.** `fig_ebic_path` and `fig_rank_tradeoff` both called
`ax.legend(loc=...)` with no `bbox_to_anchor`, putting the legend inside the axes
on top of the data. The cnsplots convention — read from `cns.settings`, not
assumed — is `loc="upper left"`, `bbox_to_anchor=(1, 1.02)`, `frameon=False`,
7 pt. `fig_rank_tradeoff` now uses it.

**A legend that had to go rather than move.** In `fig_ebic_path` the anchored
position landed on the twin right-hand axis's tick labels and axis title. Both
series are already colour-matched to their own axis, so the legend only restated
the axes; its one piece of information, the value of $\gamma$, moved into the
y-axis label. The figure now carries no legend.

**Point labels sitting on the data.** In `fig_rank_tradeoff`, `where="post"` draws
a vertical step segment to the right of every dropping point, so labels centred
above their markers landed on those segments; the first label also collided with
the top y-tick, and the third with the pale "connected nodes" series. Offsets are
now per-point and documented inline with the obstruction each one avoids.

`cnsplots` is deliberately **not** added as a dependency. The generator must stay
installable from `analysis/requirements-figures.txt` alone, so the three settings
values are hardcoded with a comment recording their source.

## Verification

- `jupyter-book build docs --warningiserror` — clean.
- `python analysis/scripts/make_docs_figures.py` — runs, and reproduces the same
  numbers it did before the edits: eBIC 16130.0988 at λ=0.8 with 216 edges; ranks
  0/2/5/10 → 216/202/158/110 edges.
- Legend and label placement checked by rendering each figure and measuring ink
  extents in the saved raster, because `bbox_inches="tight"` re-lays out at save
  time and pre-save bounding boxes do not reflect the shipped image.

## Note on figure numbers

Removing Fig 7 shifted everything after it down by one. The two `classo_cv`
heatmaps that were **Fig 21 and Fig 22** are now **Fig 20 and Fig 21**. They are
the same two images:

- `atacama-top-300-r1-cv5-selected-taxa-heatmap.png`
- `atacama-filtered-r1-cv5-selected-predictors-heatmap.png`
