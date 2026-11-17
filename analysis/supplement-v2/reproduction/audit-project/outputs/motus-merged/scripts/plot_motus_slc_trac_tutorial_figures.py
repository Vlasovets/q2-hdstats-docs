#!/usr/bin/env python3
"""Render reproducible model-selection and refit figures from QIIME 2 artifacts."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[1]
SHARED_PLOTTERS = ROOT.parent / "tutorial-figure-generation"
sys.path.insert(0, str(SHARED_PLOTTERS))
from tutorial_plotting import (
    BLUE,
    GRID,
    RED,
    REF,
    add_cv_rules,
    draw_cv_curve,
    draw_stability_profile,
    matched_limits,
    signed_descending,
)

ART = ROOT / "artifacts"
OUT = ROOT / "figures"
PROVENANCE = ROOT / "provenance"
TAXONOMY = Path(
    "/Users/cmueller/Documents/GitHub/q2-hdstats-docs/analysis/results/"
    "motus-merged/taxonomy-378.tsv"
)
COLORS = {"negative": BLUE, "positive": RED, "neutral": "#666666"}


def open_problem(path: Path):
    temp = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(path) as outer:
        member = next(name for name in outer.namelist() if name.endswith("problem.zip"))
        outer.extract(member, temp.name)
    with zipfile.ZipFile(Path(temp.name, member)) as inner:
        inner.extractall(Path(temp.name, "problem"))
    return zarr.open_group(Path(temp.name, "problem"), mode="r"), temp


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_input(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def predictions(path: Path, selection: str):
    root, temp = open_problem(path)
    prediction_key = {"cv": "YhatCV", "stabsel": "YhatStabSel"}[selection]
    return np.asarray(root["sample_labels"]), np.asarray(root[prediction_key]), temp


def taxonomy_labels() -> dict[str, str]:
    taxonomy = pd.read_csv(TAXONOMY, sep="\t").set_index("Feature ID")
    result = {}
    for feature_id, lineage in taxonomy["Taxon"].items():
        labels = [part.strip() for part in lineage.split(";")]
        terminal = next((part[3:] for part in reversed(labels) if part.startswith("g__") and part != "g__"), None)
        result[feature_id] = terminal or feature_id.replace("ref_mOTU_v31_", "mOTU ")
    return result


def readable_label(label: str, method: str, labels: dict[str, str]) -> str:
    if method == "slc":
        return labels.get(label, label.replace("ref_mOTU_v31_", "mOTU "))
    if label.startswith("m__ref_mOTU_v31_"):
        return f"mOTU {label.rsplit('_', 1)[-1]}"
    ranks = {"p__": "Phylum: ", "c__": "Class: ", "o__": "Order: ", "f__": "Family: ", "g__": "Genus: "}
    for prefix, display in ranks.items():
        if label.startswith(prefix):
            return f"{display}{label[len(prefix):]}"
    return label


def metrics(observed: np.ndarray, predicted: np.ndarray) -> tuple[float, float, float]:
    correlation = float(np.corrcoef(observed, predicted)[0, 1])
    rmse = float(np.sqrt(np.mean((observed - predicted) ** 2)))
    r_squared = float(1 - np.sum((observed - predicted) ** 2) / np.sum((observed - observed.mean()) ** 2))
    return correlation, rmse, r_squared


def scatter(ax, observed: np.ndarray, predicted: np.ndarray, title: str, limits: tuple[float, float]):
    lo, hi = limits
    correlation, rmse, r_squared = metrics(observed, predicted)
    ax.scatter(observed, predicted, s=42, color=COLORS["negative"], edgecolor="white", linewidth=0.6)
    ax.plot([lo, hi], [lo, hi], color="#555555", linewidth=1)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel="Observed infant age (days)", ylabel="Predicted infant age (days)", title=title)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="#dddddd", linewidth=0.6)
    ax.text(0.04, 0.95, f"n = {len(observed)}\nr = {correlation:.2f}\nR$^2$ = {r_squared:.2f}\nRMSE = {rmse:.1f} days", transform=ax.transAxes, va="top", fontsize=8)
    return {"r": correlation, "rmse_days": rmse, "r_squared": r_squared}


def selected_terms(solution, labels: dict[str, str], method: str):
    beta = np.asarray(solution["refit"])
    names = np.asarray(solution["label"])
    keep = (names != "intercept") & (np.abs(beta) > 1e-10)
    selected = [(readable_label(name, method, labels), value) for name, value in zip(names[keep], beta[keep])]
    return signed_descending(selected)


def stability_terms(stability, method: str, label_map: dict[str, str], label_threshold: float = 0.2):
    frequencies = np.asarray(stability["distribution"])[1:]
    names = np.asarray(stability["label"])[1:]
    keep = frequencies >= label_threshold
    return sorted(
        [(readable_label(name, method, label_map), value) for name, value in zip(names[keep], frequencies[keep])],
        key=lambda item: item[1], reverse=True,
    )


def plot_model_selection(root, method: str, label_map: dict[str, str], path: Path):
    cv = root["solution/CV"]
    stability = root["solution/StabSel"]
    x = np.asarray(cv["xGraph"])
    y = np.asarray(cv["yGraph"])
    standard_error = np.asarray(cv["standard_error"])
    lambda_min = float(cv.attrs["lambda_min"])
    lambda_one_se = float(cv.attrs["lambda_1SE"])
    selected = selected_terms(cv, label_map, method)
    stable = stability_terms(stability, method, label_map)

    fig, (ax_cv, ax_stability) = plt.subplots(
        1, 2, figsize=(11.0, 3.7), gridspec_kw={"width_ratios": [1.05, 1.25]}
    )
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.43, top=0.83, wspace=0.30)
    draw_cv_curve(ax_cv, x, y, standard_error, ylabel="CV prediction error")
    add_cv_rules(ax_cv, lambda_min, lambda_one_se)

    labels, probabilities = zip(*stable)
    draw_stability_profile(
        ax_stability,
        labels,
        probabilities,
        xlabel="mOTU / taxonomic term",
        label_size=4.6,
        label_rotation=55,
    )
    fig.savefig(path.with_suffix(".png"), dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return {
        "lambda_min": lambda_min,
        "lambda_1SE": lambda_one_se,
        "cv_selected_terms": [{"term": term, "coefficient": float(value)} for term, value in selected],
        "stable_terms": [{"term": term, "probability": float(value)} for term, value in stable],
    }


def plot_refit(root, selection: str, method: str, label_map: dict[str, str], train_y, train_hat, test_y, test_hat, path: Path):
    selected = selected_terms(root[f"solution/{selection}"], label_map, method)
    fig = plt.figure(figsize=(11.2, 3.8), constrained_layout=True)
    grid = GridSpec(1, 3, figure=fig, width_ratios=[0.9, 1, 1])
    ax_coef = fig.add_subplot(grid[0, 0])
    names, values = zip(*selected)
    colors = [COLORS["positive"] if value > 0 else COLORS["negative"] for value in values]
    ax_coef.barh(np.arange(len(values)), values, color=colors)
    ax_coef.axvline(0, color="#555555", linewidth=0.8)
    panel_title = "CV-selected terms" if selection == "CV" else "Stability-selected terms"
    ax_coef.set(yticks=np.arange(len(values)), yticklabels=names, xlabel="Refit coefficient", title=panel_title)
    ax_coef.invert_yaxis()
    ax_coef.tick_params(axis="y", labelsize=8)
    ax_coef.grid(axis="x", color="#dddddd", linewidth=0.6)

    limits = matched_limits(train_y, train_hat, test_y, test_hat)
    train_metrics = scatter(fig.add_subplot(grid[0, 1]), train_y, train_hat, "Training set", limits)
    test_metrics = scatter(fig.add_subplot(grid[0, 2]), test_y, test_hat, "Held-out set", limits)
    fig.savefig(path.with_suffix(".png"), dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return {"train": train_metrics, "test": test_metrics}


def main(method: str):
    if method not in {"slc", "trac"}:
        raise SystemExit("Usage: plot_slc_results.py {slc|trac}")
    OUT.mkdir(exist_ok=True)
    PROVENANCE.mkdir(exist_ok=True)
    stem = f"motus-age-{method}-r3"
    model = ART / f"{stem}.qza"
    train_file = ART / f"{stem}-train-predictions.qza"
    test_file = ART / f"{stem}-test-predictions.qza"
    outcome = ROOT / "exports/motus-outcomes-all36.tsv"
    label_map = taxonomy_labels()
    root, model_temp = open_problem(model)
    metadata = pd.read_csv(outcome, sep="\t").set_index("sample-id")
    train_ids, train_cv_hat, train_cv_temp = predictions(train_file, "cv")
    train_stability_ids, train_stability_hat, train_stability_temp = predictions(train_file, "stabsel")
    test_ids, test_cv_hat, test_temp = predictions(test_file, "cv")
    test_stability_ids, test_stability_hat, test_stability_temp = predictions(test_file, "stabsel")
    train_y = metadata.loc[train_ids, "host_age_days"].astype(float).to_numpy()
    test_y = metadata.loc[test_ids, "host_age_days"].astype(float).to_numpy()
    if not (np.array_equal(train_ids, train_stability_ids) and np.array_equal(test_ids, test_stability_ids)):
        raise ValueError("CV and stability-selection prediction artifacts have inconsistent sample ordering.")

    model_selection = plot_model_selection(root, method, label_map, OUT / f"{stem}-model-selection")
    cv_refit = plot_refit(root, "CV", method, label_map, train_y, train_cv_hat, test_y, test_cv_hat, OUT / f"{stem}-cv-refit")
    stability_refit = None
    if method == "trac":
        stability_refit = plot_refit(
            root,
            "StabSel",
            method,
            label_map,
            train_y,
            train_stability_hat,
            test_y,
            test_stability_hat,
            OUT / f"{stem}-stabsel-refit",
        )
    provenance = {
        "figures": [f"{stem}-model-selection", f"{stem}-cv-refit"] + ([f"{stem}-stabsel-refit"] if stability_refit else []),
        "generator": str(Path(__file__).resolve()),
        "inputs": {relative_input(path): sha256(path) for path in [model, train_file, test_file, outcome, TAXONOMY]},
        "qiime_model": f"q2-classo {method.upper()}: CLR with zeros replaced by 0.5, R3 concomitant, train-only 5-fold CV and stability selection",
        "model_selection": model_selection,
        "cv_refit_metrics": cv_refit,
        "stability_refit_metrics": stability_refit,
    }
    (PROVENANCE / f"{stem}-figures.json").write_text(json.dumps(provenance, indent=2) + "\n")
    model_temp.cleanup()
    train_cv_temp.cleanup()
    train_stability_temp.cleanup()
    test_temp.cleanup()
    test_stability_temp.cleanup()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) == 2 else "slc")
