"""Shared, publication-style plotting primitives for tutorial result figures."""

from __future__ import annotations

import numpy as np


BLUE = "#356F9F"
RED = "#C43B4D"
GRAY = "#9CA8B3"
GRID = "#E6E6E6"
REF = "#404040"


def draw_stability_profile(
    ax,
    labels,
    probabilities,
    *,
    title: str,
    xlabel: str,
    label_size: float = 5.3,
    label_rotation: float = 50,
    selection_threshold: float = 0.4,
    highlight_threshold: float = 0.3,
    annotate_thresholds: bool = False,
):
    """Draw a reusable stability-selection bar profile on an existing axis."""
    probabilities = np.asarray(probabilities, dtype=float)
    positions = np.arange(len(probabilities))
    colors = np.where(
        probabilities > selection_threshold,
        RED,
        np.where(probabilities >= highlight_threshold, BLUE, GRAY),
    )
    ax.bar(positions, probabilities, color=colors, width=0.78)
    ax.axhline(selection_threshold, color=RED, linestyle="--", linewidth=0.85)
    ax.axhline(highlight_threshold, color=REF, linestyle=":", linewidth=0.85)
    if annotate_thresholds:
        label_x = len(probabilities) - 0.25
        ax.text(label_x, selection_threshold + 0.015, f"{selection_threshold:g}", ha="right", va="bottom", fontsize=6.2, color=RED)
        ax.text(label_x, highlight_threshold + 0.015, f"{highlight_threshold:g}", ha="right", va="bottom", fontsize=6.2, color=REF)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        labels,
        rotation=label_rotation,
        ha="right",
        rotation_mode="anchor",
        fontsize=label_size,
    )
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Selection probability")
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=4)
    ax.grid(True, axis="y", color=GRID, linewidth=0.55)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
