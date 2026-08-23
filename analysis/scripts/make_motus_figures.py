#!/usr/bin/env python
"""Figures for the shotgun (mOTUs) section of the high-dimensional chapter.

One figure per plugin action, and no more. The figure audit
(`analysis/reports/figure-audit.md`) removed every image the plugins cannot actually
produce, including two node-link network drawings -- neither q2-gglasso nor q2-classo
imports networkx or emits a graph layout, so a network picture in this book would come
from a path the reader has no access to. Nothing here draws a network.

Each figure is generated from a COMMITTED TSV under analysis/results/tables/, so it
cannot drift from the numbers the chapter quotes, and it rebuilds from a clean clone
with no cluster and no QIIME 2:

    pip install -r analysis/requirements-figures.txt
    python analysis/scripts/make_motus_figures.py

Style follows the cnsplots conventions already adopted for the tier-2 figures: legends
anchored outside the axes at (1, 1.02), frameless, 7pt. cnsplots is deliberately not a
dependency -- the generator must stay installable from requirements-figures.txt alone --
so those three values are hardcoded here with their source recorded.
"""
import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
DEFAULT_OUT = ROOT.parent / "docs" / "images" / "png" / "generated"

INK, MUTED, ACCENT, GRID = "#2f3b47", "#8b98a5", "#c1442f", "#dde3e9"

# Asserted rather than trusted: a re-run that moves the selection should fail here
# instead of quietly producing a figure that contradicts the prose beside it.
SELECTED_LAMBDA = 0.30
SELECTED_EDGES = 481
N_FEATURES = 100
PAIRS = N_FEATURES * (N_FEATURES - 1) // 2      # 4,950


def _style(ax):
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=9)


def fig_lambda_path(outdir):
    """`solve-problem` -- where eBIC puts the penalty, and what it costs in edges."""
    d = pd.read_csv(TABLES / "motus-lambda-path.tsv", sep="\t").sort_values("lambda1")
    j = int(d["ebic_gamma0.15"].idxmin())
    lam, edges = float(d.loc[j, "lambda1"]), int(d.loc[j, "edges"])
    if abs(lam - SELECTED_LAMBDA) > 1e-9 or edges != SELECTED_EDGES:
        sys.exit(f"selection moved: lambda={lam} edges={edges}, expected "
                 f"{SELECTED_LAMBDA}/{SELECTED_EDGES}. Update the chapter, not this check.")

    fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=200)
    _style(ax)
    ax.plot(d.lambda1, d["ebic_gamma0.15"], "-", color=INK, linewidth=1.9, zorder=3)
    ax.plot([lam], [d.loc[j, "ebic_gamma0.15"]], "o", markersize=9, zorder=4,
            markerfacecolor="white", markeredgecolor=ACCENT, markeredgewidth=2.2)
    ax.set_xlabel(r"$\lambda_1$", fontsize=10, color=INK)
    ax.set_ylabel(r"eBIC at $\gamma = 0.15$", fontsize=10, color=INK)

    ax2 = ax.twinx()
    ax2.plot(d.lambda1, d.edges, "--", color=ACCENT, linewidth=1.5, zorder=2)
    ax2.set_ylabel("edges", fontsize=10, color=ACCENT)
    ax2.tick_params(axis="y", colors=ACCENT, labelsize=9)
    for s in ("top", "left"):
        ax2.spines[s].set_visible(False)
    ax2.spines["right"].set_color(ACCENT)

    ax.set_title(r"eBIC selects $\lambda_1 = 0.30$: %d edges of %s"
                 % (edges, f"{PAIRS:,}"),
                 fontsize=11.5, fontweight="bold", color=INK, pad=10)
    ax.annotate(f"{edges} edges\n({100*edges/PAIRS:.1f}% of pairs)",
                xy=(lam, d.loc[j, "ebic_gamma0.15"]),
                xytext=(lam - 0.19, d["ebic_gamma0.15"].max() * 0.90),
                fontsize=9, color=INK, ha="left",
                arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=1,
                                shrinkB=8, connectionstyle="angle,angleA=0,angleB=90"))
    out = outdir / "motus-lambda-path.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return out, f"eBIC min at lambda={lam}, {edges} edges"


def fig_mu_rank(outdir):
    """`solve-problem --p-latent` -- what each latent dimension costs the network."""
    d = pd.read_csv(TABLES / "motus-top100-mu-rank-map.tsv", sep="\t").sort_values("rank")
    fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=200)
    _style(ax)
    ax.step(d["rank"], d.edges, where="post", color=INK, linewidth=1.9, zorder=3)
    ax.plot(d["rank"], d.edges, "o", color=INK, markersize=5, zorder=4)
    ax.axhline(SELECTED_EDGES, color=ACCENT, linestyle=":", linewidth=1.4, zorder=2)
    ax.text(d["rank"].max() * 0.50, SELECTED_EDGES * 0.90,
            f"rank 0 = {SELECTED_EDGES} edges", fontsize=9, color=ACCENT, va="top")
    ax.set_xlabel("achieved rank of the latent block", fontsize=10, color=INK)
    ax.set_ylabel("edges in the sparse component", fontsize=10, color=INK)
    ax.set_title("Three latent dimensions absorb more than half the network",
                 fontsize=11.5, fontweight="bold", color=INK, pad=10)
    for _, r in d.iterrows():
        if int(r["rank"]) in (3, 7, 16):
            # the rank-16 point sits on the right spine, so its label goes leftwards
            dx, ha = (-8, "right") if int(r["rank"]) == 16 else (6, "left")
            ax.annotate(rf"$\mu_1$={r.mu1:g}", (r["rank"], r.edges),
                        textcoords="offset points", xytext=(dx, 9),
                        fontsize=8.5, color=MUTED, ha=ha)
    out = outdir / "motus-mu-rank.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return out, f"{len(d)} mu values, ranks {d['rank'].min()}-{d['rank'].max()}"


def fig_trac_coefficients(outdir):
    """`classo regress` -- which clades the log-contrast model keeps."""
    d = pd.read_csv(TABLES / "motus-trac-selected.tsv", sep="\t")
    d = d.sort_values("coefficient")
    fig, ax = plt.subplots(figsize=(7.4, 2.9), dpi=200)
    _style(ax)
    colors = [ACCENT if v > 0 else "#3b6ea5" for v in d.coefficient]
    ax.barh(range(len(d)), d.coefficient, color=colors, height=0.55, zorder=3)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.clade, fontsize=9.5, color=INK)
    ax.axvline(0, color=MUTED, linewidth=1)
    ax.set_xlabel("log-contrast coefficient", fontsize=10, color=INK)
    ax.set_title("trac keeps three clades for postnatal age",
                 fontsize=11.5, fontweight="bold", color=INK, pad=10)
    out = outdir / "motus-trac-coefficients.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return out, f"{len(d)} clades selected"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    need = ["motus-lambda-path.tsv", "motus-top100-mu-rank-map.tsv",
            "motus-trac-selected.tsv"]
    missing = [t for t in need if not (TABLES / t).is_file()]
    if missing:
        sys.exit(f"missing input tables in {TABLES}: {missing}. "
                 "Run analysis/slurm/44_motus_tutorial_run.sh first.")

    for fn in (fig_lambda_path, fig_mu_rank, fig_trac_coefficients):
        path, note = fn(args.outdir)
        print(f"  {path.name:36} {note}")
    print(f"  -> {args.outdir}")


if __name__ == "__main__":
    main()
