"""Shared plotting conventions for the q2-hdstats tutorial result figures."""

from __future__ import annotations

import numpy as np


BLUE = "#356F9F"
RED = "#C43B4D"
GRAY = "#9CA8B3"
GRID = "#E6E6E6"
REF = "#404040"


def draw_cv_curve(ax, x, error, standard_error, *, ylabel):
    """Draw a CV curve without estimator-specific wording in the panel."""
    order = np.argsort(x)
    x, error, standard_error = (np.asarray(value)[order] for value in (x, error, standard_error))
    ax.plot(x, error, color=BLUE, linewidth=1.1)
    ax.fill_between(x, error - standard_error, error + standard_error, color=BLUE, alpha=0.17, linewidth=0)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set(xlabel="lambda / lambda max", ylabel=ylabel)
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.spines[["top", "right"]].set_visible(False)


def add_cv_rules(ax, lambda_min, lambda_one_se):
    """Add consistent minimum-error and 1SE rules to a CV axis."""
    ax.axvline(lambda_min, color=REF, linestyle=":", linewidth=0.85, label="minimum")
    ax.axvline(lambda_one_se, color=REF, linestyle="--", linewidth=0.85, label="1SE")
    ax.legend(frameon=False, loc="lower right", fontsize=6.3)


def draw_stability_profile(
    ax,
    labels,
    probabilities,
    *,
    xlabel,
    label_size=5.3,
    label_rotation=50,
    selection_threshold=0.4,
    highlight_threshold=0.3,
):
    """Draw a reusable stability-selection profile without a panel subtitle."""
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
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=label_rotation, ha="right", rotation_mode="anchor", fontsize=label_size)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Selection probability")
    ax.set_xlabel(xlabel)
    ax.grid(True, axis="y", color=GRID, linewidth=0.55)
    ax.spines[["top", "right"]].set_visible(False)


def signed_descending(items):
    """Sort (label, coefficient) pairs from largest positive to largest negative."""
    return sorted(items, key=lambda item: item[1], reverse=True)


def matched_limits(*arrays, pad_fraction=0.05):
    """Return one square plotting range shared by paired prediction panels."""
    values = np.concatenate([np.asarray(array, dtype=float).ravel() for array in arrays])
    lower, upper = float(values.min()), float(values.max())
    pad = max((upper - lower) * pad_fraction, 1e-8)
    return lower - pad, upper + pad
