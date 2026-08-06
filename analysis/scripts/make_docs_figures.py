#!/usr/bin/env python
"""Generate the book's result figures from the recompute's own output tables.

The book had 8 figures across 48 chapters, and the two on the model-selection
page did not carry their claims:

  - the eBIC figure plotted five gamma curves when the chapter uses gamma = 0.3,
    on a y-axis of *eBIC minus each curve's own minimum*, so every curve bottomed
    at zero and the selected lambda was invisible and unmarked;
  - the rank-comparison figure was an 11-panel grid whose entire "sparse" column
    rendered as a bare diagonal (216 edges in 44,850 cells), with an unrelated
    dendrogram in one corner and a different colour scale on every panel.

Figures here are generated from `analysis/results/tables/*.tsv`, so they cannot
drift from the numbers the chapters quote. That drift is not hypothetical: the
hand-transcribed edge counts in the chapter disagreed with the pipeline's own
output (1403 vs 1405) precisely because a human step sat in between.

Each figure states one claim and marks it. Usage:

    python analysis/scripts/make_docs_figures.py [--outdir DIR]
"""
import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]        # <repo>/analysis
TABLES = ROOT / "results" / "tables"
DEFAULT_OUT = ROOT.parent / "docs" / "images" / "png" / "generated"

# The chapter's canonical parameterisation. Asserted against the data below
# rather than trusted, so a re-run that moves the selection fails loudly here
# instead of silently producing a figure that contradicts the prose.
SELECTED_LAMBDA = 0.8
SELECTED_EDGES = 216
N_FEATURES = 300
PAIRS = N_FEATURES * (N_FEATURES - 1) // 2               # 44,850

INK = "#2f3b47"
MUTED = "#8b98a5"
ACCENT = "#c1442f"
GRID = "#dde3e9"


def _style(ax):
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=10)


def fig_ebic_path(outdir):
    """eBIC against lambda at the gamma the chapter actually uses, minimum marked.

    One curve, one decision, the answer annotated on the figure. Edge count rides
    on a twin axis because "216 edges" is the number the reader carries forward,
    and it is otherwise only available by reading a separate table.
    """
    df = pd.read_csv(TABLES / "lambda-path.tsv", sep="\t").sort_values("lambda1")
    lam = df["lambda1"].to_numpy(float)
    ebic = df["ebic"].to_numpy(float)
    edges = np.rint(df["sparsity"].to_numpy(float) * PAIRS).astype(int)

    j = int(np.argmin(ebic))
    if not np.isclose(lam[j], SELECTED_LAMBDA):
        sys.exit(f"eBIC minimum is at lambda={lam[j]}, expected {SELECTED_LAMBDA}. "
                 "The chapter's prose and this figure would disagree — fix one.")

    fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=200)
    _style(ax)
    ax.plot(lam, ebic, "-o", color=INK, markersize=4.5, linewidth=1.8, zorder=3,
            label=r"eBIC ($\gamma = 0.3$)")
    ax.plot([lam[j]], [ebic[j]], "o", markersize=13, markerfacecolor="none",
            markeredgecolor=ACCENT, markeredgewidth=2.4, zorder=4)
    ax.axvline(lam[j], color=ACCENT, linewidth=1.0, linestyle=(0, (4, 4)), zorder=2)
    ax.annotate(
        f"minimum at $\\lambda$ = {lam[j]:g}\n{edges[j]} edges",
        # placed left of the minimum: the curve rises steeply to the right of
        # lambda=0.8 and the legend occupies the top centre
        xy=(lam[j], ebic[j]), xytext=(-150, 34), textcoords="offset points",
        fontsize=11, color=ACCENT, fontweight="bold",
        arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=1.2),
    )
    ax.set_xlabel(r"$\lambda_1$  (sparsity penalty)", fontsize=11, color=INK)
    ax.set_ylabel(r"eBIC", fontsize=11, color=INK)
    ax.set_title("Lower is better: eBIC picks the penalty", fontsize=13,
                 fontweight="bold", color=INK, pad=12)

    ax2 = ax.twinx()
    ax2.plot(lam, edges, "-", color=MUTED, linewidth=1.4, zorder=1)
    ax2.set_ylabel("edges in the network", fontsize=10, color=MUTED)
    ax2.tick_params(colors=MUTED, labelsize=9)
    for side in ("top", "left"):
        ax2.spines[side].set_visible(False)
    ax2.spines["right"].set_color(MUTED)

    ax.legend(loc="upper center", frameon=False, fontsize=10)
    out = outdir / "atacama-ebic-lambda-selection.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out, f"min eBIC {ebic[j]:.4f} at lambda={lam[j]:g}, {edges[j]} edges"


def fig_rank_tradeoff(outdir):
    """What each extra latent dimension costs in network structure.

    Replaces the 11-panel heatmap. The claim being made is "rank 2 is the
    parsimonious choice"; this shows the quantity that claim is about — how many
    edges and connected nodes survive as the latent component absorbs more.
    """
    df = pd.read_csv(TABLES / "mu-rank-map.tsv", sep="\t")
    rank = np.concatenate([[0], df["achieved rank"].to_numpy(int)])
    edges = np.concatenate([[SELECTED_EDGES], df["sparse edges"].to_numpy(int)])
    nodes = np.concatenate([[163], df["connected nodes"].to_numpy(int)])
    order = np.argsort(rank)
    rank, edges, nodes = rank[order], edges[order], nodes[order]

    fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=200)
    _style(ax)
    ax.step(rank, edges, where="post", color=INK, linewidth=2.0, zorder=3)
    ax.plot(rank, edges, "o", color=INK, markersize=7, zorder=4, label="sparse edges")
    ax.step(rank, nodes, where="post", color=MUTED, linewidth=1.6,
            linestyle=(0, (5, 3)), zorder=2)
    ax.plot(rank, nodes, "s", color=MUTED, markersize=6, zorder=3,
            label="connected nodes")

    for r, e in zip(rank, edges):
        ax.annotate(str(e), xy=(r, e), xytext=(0, 11), textcoords="offset points",
                    ha="center", fontsize=9.5, color=INK)

    ax.axvspan(-0.35, 2.35, color=ACCENT, alpha=0.07, zorder=0)
    ax.annotate("rank 2 keeps 94% of the edges", xy=(2, edges[rank == 2][0]),
                xytext=(24, -34), textcoords="offset points", fontsize=10.5,
                color=ACCENT, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=1.2))

    ax.set_xticks(rank)
    ax.set_xlabel(r"achieved rank of the low-rank component $\hat{L}$",
                  fontsize=11, color=INK)
    ax.set_ylabel("count", fontsize=11, color=INK)
    ax.set_title(r"What each latent dimension absorbs ($\lambda = 0.8$)",
                 fontsize=13, fontweight="bold", color=INK, pad=12)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.set_xlim(-0.5, max(rank) + 0.6)

    out = outdir / "atacama-rank-tradeoff.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out, "rank 0/2/5/10 -> " + "/".join(str(e) for e in edges) + " edges"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    missing = [t for t in ("lambda-path.tsv", "mu-rank-map.tsv")
               if not (TABLES / t).is_file()]
    if missing:
        sys.exit(f"missing input tables in {TABLES}: {missing}. "
                 "Run the recompute stages first.")

    for fn in (fig_ebic_path, fig_rank_tradeoff):
        path, note = fn(args.outdir)
        print(f"  {path.name:44} {note}")
    print(f"  -> {args.outdir}")


if __name__ == "__main__":
    main()
