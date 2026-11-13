#!/usr/bin/env python3
"""Render tutorial-style plots from the QIIME-native Bracken bean classifier."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import zipfile
from pathlib import Path

import biom
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import zarr


ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"
OUT = ROOT / "figures"
EXPORT = ROOT / "exports"
PROVENANCE = ROOT / "provenance"
SCRIPT = Path(__file__).resolve()
SHARED = Path(
    "/Users/cmueller/Documents/ChatGPT/Atacama-q2-gglasso-q2-classo-tutorial/"
    "outputs/fresh-atacama-workflow/tutorial-figure-generation"
)
sys.path.insert(0, str(SHARED))
from tutorial_plotting import draw_stability_profile


plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.0, "axes.titlesize": 7.8,
    "axes.labelsize": 7.3, "xtick.labelsize": 6.1, "ytick.labelsize": 5.6,
    "legend.fontsize": 6.3, "figure.titlesize": 9.0, "axes.linewidth": 0.75,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
BLUE, RED, GREEN, GRID, REF = "#356F9F", "#C43B4D", "#3E8A5B", "#E6E6E6", "#404040"
MODEL = ART / "bracken-bean-type-clr05-classo.qza"
PREDICTIONS = ART / "bracken-bean-type-clr05-classo-predictions.qza"
CLR = ART / "bracken-species-prev14-total50k-clr05.qza"
TAXONOMY = ART / "bracken-species-taxonomy.qza"
METADATA = ROOT / "metadata.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def open_problem(path: Path):
    temp = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(path) as outer:
        member = next(name for name in outer.namelist() if name.endswith("problem.zip"))
        outer.extract(member, temp.name)
    with zipfile.ZipFile(Path(temp.name, member)) as inner:
        inner.extractall(Path(temp.name, "problem"))
    return zarr.open_group(Path(temp.name, "problem"), mode="r"), temp


def qza_member(path: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(path) as archive:
        member = next(name for name in archive.namelist() if name.endswith(suffix))
        return archive.read(member)


def taxonomy_labels() -> dict[str, str]:
    lines = qza_member(TAXONOMY, "data/taxonomy.tsv").decode().splitlines()
    taxonomy = pd.read_csv(pd.io.common.StringIO("\n".join(lines)), sep="\t")
    labels = {}
    for feature, lineage in taxonomy[["Feature ID", "Taxon"]].itertuples(index=False):
        ranks = [part.strip() for part in lineage.split(";")]
        species = next((part[3:] for part in reversed(ranks) if part.startswith("s__") and part != "s__"), None)
        genus = next((part[3:] for part in reversed(ranks) if part.startswith("g__") and part != "g__"), None)
        labels[str(feature)] = species or genus or str(feature)
    return labels


def clr_table() -> pd.DataFrame:
    temp = tempfile.TemporaryDirectory()
    qza_path = Path(temp.name, "table.qza")
    qza_path.write_bytes(CLR.read_bytes())
    with zipfile.ZipFile(qza_path) as archive:
        member = next(name for name in archive.namelist() if name.endswith("data/feature-table.biom"))
        biom_path = Path(temp.name, "feature-table.biom")
        biom_path.write_bytes(archive.read(member))
    table = biom.load_table(str(biom_path))
    result = pd.DataFrame(table.matrix_data.toarray().T, index=table.ids(axis="sample"), columns=table.ids(axis="observation"))
    temp.cleanup()
    return result


def save(fig, stem: str, description: str, inputs: list[Path], parameters: dict | None = None) -> None:
    OUT.mkdir(exist_ok=True)
    PROVENANCE.mkdir(exist_ok=True)
    outputs = [OUT / f"{stem}.png", OUT / f"{stem}.pdf"]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    plt.close(fig)
    payload = {
        "description": description, "generator": str(SCRIPT),
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
        "outputs": [str(path.relative_to(ROOT)) for path in outputs],
        "parameters": parameters or {},
    }
    (PROVENANCE / f"{stem}.json").write_text(json.dumps(payload, indent=2) + "\n")


def export_tables(root, pred_root, labels: dict[str, str], metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    EXPORT.mkdir(exist_ok=True)
    cv, ss = root["solution/CV"], root["solution/StabSel"]
    cv_table = pd.DataFrame({
        "lambda_over_lambda_max": np.asarray(cv["xGraph"]),
        "cv_misclassification_error": np.asarray(cv["yGraph"]),
        "standard_error": np.asarray(cv["standard_error"]),
    })
    cv_table["is_lambda_min"] = np.isclose(cv_table["lambda_over_lambda_max"], float(cv.attrs["lambda_min"]))
    cv_table["is_lambda_1se"] = np.isclose(cv_table["lambda_over_lambda_max"], float(cv.attrs["lambda_1SE"]))
    cv_table.to_csv(EXPORT / "bracken-bean-classification-cv-curve.tsv", sep="\t", index=False)

    feature_labels = np.asarray(ss["label"])
    stability = np.asarray(ss["distribution"])
    stability_table = pd.DataFrame({
        "feature_id": feature_labels, "stability": stability,
        "display_label": ["intercept" if x == "intercept" else labels.get(str(x), str(x)) for x in feature_labels],
    })
    stability_table = stability_table.query("feature_id != 'intercept'").sort_values("stability", ascending=False)
    stability_table.to_csv(EXPORT / "bracken-bean-classification-stability-profile.tsv", sep="\t", index=False)

    sample_ids = np.asarray(pred_root["sample_labels"]).astype(str)
    score_table = metadata.loc[sample_ids, ["bean_type", "timepoint_hours"]].copy()
    score_table.index.name = "sample-id"
    score_table["cv_score"] = np.asarray(pred_root["YhatCV"])
    score_table["stabsel_score"] = np.asarray(pred_root["YhatStabSel"])
    score_table.to_csv(EXPORT / "bracken-bean-classification-in-sample-scores.tsv", sep="\t")
    return cv_table, stability_table, score_table.reset_index()


def coefficient_table(solution, labels: dict[str, str]) -> pd.DataFrame:
    names = np.asarray(solution["label"]).astype(str)
    beta = np.asarray(solution["refit"])
    selection = np.asarray(solution["selected_param"]).astype(bool)
    result = pd.DataFrame({"feature_id": names, "refit_coefficient": beta, "selected": selection})
    result = result.query("feature_id != 'intercept' and selected").copy()
    result["display_label"] = result["feature_id"].map(labels).fillna(result["feature_id"])
    return result.sort_values("refit_coefficient")


def plot_clr_heatmap(clr: pd.DataFrame, metadata: pd.DataFrame, labels: dict[str, str]) -> None:
    order = metadata.sort_values(["bean_type", "timepoint_hours"]).index
    data = clr.loc[order].T
    species = pd.DataFrame({"feature_id": data.index, "display_label": [labels.get(x, x) for x in data.index]})
    species.to_csv(EXPORT / "bracken-filtered-species-labels.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 8.0))
    fig.subplots_adjust(left=0.30, right=0.87, bottom=0.18, top=0.86)
    image = ax.imshow(data.to_numpy(float), aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-5, vmax=5)
    for boundary in np.flatnonzero(metadata.loc[order, "bean_type"].to_numpy()[1:] != metadata.loc[order, "bean_type"].to_numpy()[:-1]):
        ax.axvline(boundary + 0.5, color="white", linewidth=1.4)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=55, ha="right", fontsize=5.4)
    ax.set_yticks(np.arange(len(data.index)))
    ax.set_yticklabels([labels.get(x, x) for x in data.index], fontsize=3.6)
    ax.set_xlabel("Samples, grouped by bean type and ordered by fermentation time")
    ax.set_ylabel("61 prevalence- and abundance-filtered species")
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("CLR abundance, pseudocount 0.5", fontsize=6.6)
    top = ax.secondary_xaxis("top")
    top.set_xticks(np.arange(len(order)))
    top.set_xticklabels([f"{int(metadata.loc[sample, 'timepoint_hours'])} h" for sample in order], fontsize=5.1)
    top.set_xlabel("Fermentation time", fontsize=6.7, labelpad=5)
    top.tick_params(length=0, pad=2)
    for name, start, end, color in [("Forasteiro", 0, 7, BLUE), ("Hybrid", 7, 14, RED)]:
        x = (start + end - 1) / 2
        ax.text(x, -4.8, name, color=color, ha="center", va="bottom", fontsize=7.0, fontweight="bold", clip_on=False)
    save(fig, "bracken-bean-clr-filtered-data-sorted-by-label", "CLR-transformed filtered Bracken species table, sorted by bean type and fermentation time", [CLR, METADATA, TAXONOMY], {"filter": "prevalence = 14/14; total frequency >= 50000", "clr_pseudocount": 0.5})


def plot_selection(cv: pd.DataFrame, stability: pd.DataFrame) -> None:
    display = stability[stability["stability"] >= 0.20].copy()
    fig = plt.figure(figsize=(8.7, 3.25))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.43, top=0.82, wspace=0.28)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.32])
    ax = fig.add_subplot(gs[0, 0])
    d = cv.sort_values("lambda_over_lambda_max")
    ax.plot(d["lambda_over_lambda_max"], d["cv_misclassification_error"], color=GREEN, linewidth=1.1)
    ax.fill_between(d["lambda_over_lambda_max"], d["cv_misclassification_error"] - d["standard_error"], d["cv_misclassification_error"] + d["standard_error"], color=GREEN, alpha=0.17, linewidth=0)
    for flag, linestyle, label in [("is_lambda_min", ":", "min"), ("is_lambda_1se", "--", "1-SE")]:
        value = d.loc[d[flag], "lambda_over_lambda_max"]
        if len(value): ax.axvline(value.iloc[0], color=REF, linewidth=0.85, linestyle=linestyle, label=label)
    ax.set_xscale("log"); ax.invert_xaxis(); ax.grid(True, color=GRID, linewidth=0.55)
    ax.set(xlabel="lambda / lambda max", ylabel="CV classification error", title="q2-classo bean-type classification CV")
    ax.legend(frameon=False, loc="lower right")
    for spine in ["top", "right"]: ax.spines[spine].set_visible(False)
    ax = fig.add_subplot(gs[0, 1])
    draw_stability_profile(ax, display["display_label"], display["stability"], title="q2-classo stability selection", xlabel="Species", label_size=5.2, label_rotation=45)
    save(fig, "bracken-bean-classification-cv-and-stability", "Cross-validation and stability-selection diagnostics for QIIME-native Bracken bean-type classification", [MODEL, METADATA], {"cv_folds": 5, "stabsel_B": 100, "stabsel_threshold": 0.4, "label_threshold": 0.2})


def score_mapping(scores: pd.DataFrame, column: str) -> dict[str, str]:
    means = scores.groupby("bean_type")[column].mean().sort_values()
    return {means.index[0]: "negative", means.index[-1]: "positive"}


def classification_metrics(d: pd.DataFrame, column: str) -> tuple[float, pd.DataFrame]:
    signs = score_mapping(d, column)
    predicted = np.where(d[column] >= 0, next(k for k, v in signs.items() if v == "positive"), next(k for k, v in signs.items() if v == "negative"))
    result = d.copy(); result["predicted"] = predicted; result["correct"] = result["bean_type"] == result["predicted"]
    return float(result["correct"].mean()), result


def plot_refit(solution_name: str, solution, scores: pd.DataFrame, labels: dict[str, str]) -> None:
    coef = coefficient_table(solution, labels)
    if coef.empty:
        raise RuntimeError(f"{solution_name} refit selected no species; cannot render tutorial three-panel plot.")
    coef.to_csv(EXPORT / f"bracken-bean-classification-{solution_name.lower()}-coefficients.tsv", sep="\t", index=False)
    score_column = "cv_score" if solution_name == "CV" else "stabsel_score"
    accuracy, d = classification_metrics(scores, score_column)
    fig = plt.figure(figsize=(8.7, max(3.15, 2.55 + 0.24 * len(coef))))
    fig.subplots_adjust(left=0.26, right=0.985, bottom=0.18, top=0.80, wspace=0.34)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.28, 1.0, 1.0])
    ax0, ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
    pos = np.arange(len(coef)); values = coef["refit_coefficient"].to_numpy(float)
    ax0.barh(pos, values, color=[BLUE if x < 0 else RED for x in values])
    ax0.axvline(0, color=REF, linewidth=0.8); ax0.grid(True, axis="x", color=GRID, linewidth=0.55)
    ax0.set(yticks=pos, yticklabels=coef["display_label"], xlabel="Refit coefficient", title="Selected species")
    ax0.tick_params(axis="y", labelsize=5.5)
    for spine in ["top", "right"]: ax0.spines[spine].set_visible(False)
    order = ["Forasteiro", "Hybrid"]; colors = {"Forasteiro": BLUE, "Hybrid": RED}; rng = np.random.default_rng(20261021)
    for bean in order:
        dd = d[d["bean_type"] == bean]; xx = np.repeat(order.index(bean), len(dd)) + rng.uniform(-0.08, 0.08, len(dd))
        ax1.scatter(xx, dd[score_column], s=27, facecolor=colors[bean], edgecolor="white", linewidth=0.45)
    ax1.axhline(0, color=REF, linewidth=0.85, linestyle="--"); ax1.grid(True, axis="y", color=GRID, linewidth=0.55)
    ax1.set(xlim=(-0.45, 1.45), xticks=[0, 1], xticklabels=order, xlabel="Observed bean type", ylabel="Classification score", title="In-sample refit score")
    ax1.legend(handles=[Line2D([0],[0], marker='o', color='none', markerfacecolor=BLUE, markeredgecolor='white', markersize=4.7, label='Forasteiro'), Line2D([0],[0], marker='o', color='none', markerfacecolor=RED, markeredgecolor='white', markersize=4.7, label='Hybrid')], loc='lower left', frameon=True, title='Observed', title_fontsize=6.3)
    table = pd.crosstab(pd.Categorical(d["bean_type"], categories=order), pd.Categorical(d["predicted"], categories=order)).reindex(index=order, columns=order, fill_value=0)
    image = ax2.imshow(table.to_numpy(), cmap="Blues", vmin=0, vmax=max(1, table.to_numpy().max()))
    for (i, j), value in np.ndenumerate(table.to_numpy()): ax2.text(j, i, str(value), ha="center", va="center", fontsize=8, color=REF)
    ax2.set(xticks=[0, 1], xticklabels=order, yticks=[0, 1], yticklabels=order, xlabel="Predicted bean type", title="In-sample classifications")
    ax2.text(0.95, 0.95, f"n={len(d)}\naccuracy={accuracy:.2f}", transform=ax2.transAxes, ha="right", va="top", fontsize=5.9, bbox=dict(facecolor="white", edgecolor="#CCCCCC", linewidth=0.45, boxstyle="round,pad=0.18"))
    stem = f"bracken-bean-classification-{solution_name.lower()}-refit-three-panel"
    save(fig, stem, f"{solution_name} refit three-panel summary for QIIME-native Bracken bean-type classification", [MODEL, PREDICTIONS, METADATA, TAXONOMY], {"selection": solution_name, "in_sample_only": True})


def main() -> None:
    # Retain the QIIME metadata typing row for the action, but skip it for plotting.
    metadata = pd.read_csv(METADATA, sep="\t", skiprows=[1]).set_index("sample-id")
    labels = taxonomy_labels()
    model, model_temp = open_problem(MODEL)
    predicted, prediction_temp = open_problem(PREDICTIONS)
    clr = clr_table()
    plot_clr_heatmap(clr, metadata, labels)
    cv, stability, scores = export_tables(model, predicted, labels, metadata)
    plot_selection(cv, stability)
    plot_refit("CV", model["solution/CV"], scores, labels)
    plot_refit("StabSel", model["solution/StabSel"], scores, labels)
    model_temp.cleanup(); prediction_temp.cleanup()


if __name__ == "__main__":
    main()
