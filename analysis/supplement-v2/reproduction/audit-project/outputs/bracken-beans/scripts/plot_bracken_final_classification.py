#!/usr/bin/env python3
"""Render final tutorial figures from the CV-only MOSHPIT bean classifier."""

from __future__ import annotations

import hashlib
import json
import tempfile
import zipfile
from pathlib import Path
import sys

import biom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr


ROOT = Path(__file__).resolve().parents[1]
ART, FIG, EXP, PROV = (ROOT / name for name in ("artifacts", "figures", "exports", "provenance"))
sys.path.insert(0, str(ROOT.parent / "tutorial-figure-generation"))
from tutorial_plotting import BLUE, GRID, RED, REF, add_cv_rules, draw_cv_curve, signed_descending
MODEL = ART / "bracken-bean-type-train-clr05-c2-rho0-cv-no-refit-pr60.qza"
CLR = ART / "bracken-species-prev14-total50k-clr05.qza"
TAXONOMY = ART / "bracken-species-taxonomy.qza"
METADATA = ROOT / "metadata.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def open_problem(path: Path):
    temporary = tempfile.TemporaryDirectory()
    with zipfile.ZipFile(path) as archive:
        member = next(name for name in archive.namelist() if name.endswith("problem.zip"))
        archive.extract(member, temporary.name)
    with zipfile.ZipFile(Path(temporary.name, member)) as archive:
        archive.extractall(Path(temporary.name, "problem"))
    return zarr.open_group(Path(temporary.name, "problem"), mode="r"), temporary


def member(path: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(path) as archive:
        return archive.read(next(name for name in archive.namelist() if name.endswith(suffix)))


def taxonomy_labels() -> dict[str, str]:
    taxonomy = pd.read_csv(pd.io.common.StringIO(member(TAXONOMY, "data/taxonomy.tsv").decode()), sep="\t")
    result = {}
    for feature, lineage in taxonomy[["Feature ID", "Taxon"]].itertuples(index=False):
        ranks = [value.strip() for value in lineage.split(";")]
        species = next((value[3:] for value in reversed(ranks) if value.startswith("s__") and value != "s__"), None)
        genus = next((value[3:] for value in reversed(ranks) if value.startswith("g__") and value != "g__"), None)
        result[str(feature)] = species or genus or str(feature)
    return result


def clr_table() -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as temporary:
        biom_path = Path(temporary, "table.biom")
        biom_path.write_bytes(member(CLR, "data/feature-table.biom"))
        table = biom.load_table(str(biom_path))
        return pd.DataFrame(table.matrix_data.toarray().T, index=table.ids(axis="sample"), columns=table.ids(axis="observation"))


def save(fig, stem: str, description: str, inputs: list[Path], parameters: dict) -> None:
    FIG.mkdir(exist_ok=True); EXP.mkdir(exist_ok=True); PROV.mkdir(exist_ok=True)
    pdf, png = FIG / f"{stem}.pdf", FIG / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight"); fig.savefig(png, dpi=300, bbox_inches="tight"); plt.close(fig)
    (PROV / f"{stem}.json").write_text(json.dumps({
        "description": description,
        "generator": str(Path(__file__).resolve()),
        "inputs": {str(path.relative_to(ROOT)): digest(path) for path in inputs},
        "parameters": parameters,
    }, indent=2) + "\n")


def penalized_scores(model_path: Path, clr: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    root, temporary = open_problem(model_path)
    cv = root["solution/CV"]
    labels = np.asarray(cv["label"]).astype(str)
    coefficients = np.asarray(cv["beta"])
    intercept = coefficients[labels == "intercept"][0]
    microbial_labels = labels[labels != "intercept"]
    microbial_coefficients = coefficients[labels != "intercept"]
    result = metadata.loc[clr.index, ["bean_type", "timepoint_hours", "analysis_set"]].copy()
    result.index.name = "sample-id"
    result["score"] = intercept + clr.loc[:, microbial_labels].to_numpy() @ microbial_coefficients
    result["predicted"] = np.where(result["score"] >= 0, "Hybrid", "Forasteiro")
    result["correct"] = result["bean_type"] == result["predicted"]
    temporary.cleanup()
    return result.reset_index()


def plot_data(clr: pd.DataFrame, metadata: pd.DataFrame, labels: dict[str, str]) -> None:
    order = metadata.sort_values(["bean_type", "timepoint_hours"]).index
    fig, axis = plt.subplots(figsize=(7.2, 8.0))
    image = axis.imshow(clr.loc[order].T, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-5, vmax=5)
    for boundary in np.flatnonzero(metadata.loc[order, "bean_type"].to_numpy()[1:] != metadata.loc[order, "bean_type"].to_numpy()[:-1]):
        axis.axvline(boundary + 0.5, color="white", linewidth=1.4)
    axis.set(xticks=np.arange(len(order)), xticklabels=order, yticks=np.arange(clr.shape[1]), yticklabels=[labels.get(feature, feature) for feature in clr.columns], xlabel="Samples", ylabel="Species")
    axis.tick_params(axis="x", rotation=55, labelsize=5.4)
    axis.tick_params(axis="y", labelsize=3.6)
    bar = fig.colorbar(image, ax=axis, fraction=0.045, pad=0.02); bar.set_label("CLR abundance, pseudocount 0.5", fontsize=6.6)
    top = axis.secondary_xaxis("top")
    top.set_xticks(np.arange(len(order)))
    top.set_xticklabels([f"{int(metadata.loc[sample, 'timepoint_hours'])} h" for sample in order], fontsize=5.1)
    top.set_xlabel("Fermentation time", fontsize=6.7, labelpad=5)
    top.tick_params(length=0, pad=2)
    for name, start, end, color in (("Forasteiro", 0, 7, BLUE), ("Hybrid", 7, 14, RED)):
        axis.text((start + end - 1) / 2, -4.8, name, color=color, ha="center", va="bottom", fontweight="bold")
    fig.subplots_adjust(bottom=0.18, left=0.30, right=0.87, top=0.86)
    save(fig, "bracken-bean-clr-filtered-data-sorted-by-label", "Filtered CLR species table", [CLR, METADATA, TAXONOMY], {"filter": "present in 14/14 samples; total frequency >= 50000", "clr_pseudocount": 0.5})


def plot_cv(cv: pd.DataFrame) -> None:
    values = cv.sort_values("lambda_over_lambda_max")
    fig, axis = plt.subplots(figsize=(4.6, 3.25))
    draw_cv_curve(axis, values.lambda_over_lambda_max, values.error, values.se, ylabel="CV classification error")
    minimum = values.loc[values["minimum"], "lambda_over_lambda_max"].iloc[0]
    one_se = values.loc[values["one_se"], "lambda_over_lambda_max"].iloc[0]
    add_cv_rules(axis, minimum, one_se)
    save(fig, "bracken-bean-classification-cv", "Five-fold training cross-validation curve", [MODEL], {"cv_folds": 5, "selection": "1SE", "c_lasso_source": "bio-datascience/c-lasso@d22da8d5ea54eb499de2fce82fa6c0e1f0c2dcc0"})


def plot_penalized(coef: pd.DataFrame, train: pd.DataFrame, test: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.55), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    fig.subplots_adjust(left=.24, right=.98, bottom=.20, top=.82, wspace=.38)
    values = coef.penalized_coefficient.to_numpy(float); pos = np.arange(len(coef))
    axes[0].barh(pos, values, color=[BLUE if value < 0 else RED for value in values]); axes[0].axvline(0, color=REF, linewidth=.8); axes[0].grid(axis="x", color=GRID, linewidth=.55)
    axes[0].set(yticks=pos, yticklabels=coef.display_label, xlabel="Penalized coefficient", title="Selected species"); axes[0].tick_params(axis="y", labelsize=5.6)
    axes[0].invert_yaxis()
    score_limits = (min(train.score.min(), test.score.min()), max(train.score.max(), test.score.max()))
    score_pad = max((score_limits[1] - score_limits[0]) * 0.06, 0.02)
    score_limits = (score_limits[0] - score_pad, score_limits[1] + score_pad)
    for axis, frame, title in ((axes[1], train, "In-sample classification"), (axes[2], test, "Held-out prediction")):
        for bean, color, xpos in (("Forasteiro", BLUE, 0), ("Hybrid", RED, 1)):
            subset = frame.query("bean_type == @bean")
            jitter = np.linspace(-.06, .06, len(subset)) if len(subset) > 1 else np.array([0.])
            axis.scatter(np.repeat(xpos, len(subset)) + jitter, subset.score, s=30, facecolor=color, edgecolor="white", linewidth=.45)
        axis.axhline(0, color=REF, linewidth=.85, linestyle="--"); axis.grid(axis="y", color=GRID, linewidth=.55)
        accuracy = frame.correct.mean()
        axis.set(xlim=(-.45, 1.45), ylim=score_limits, xticks=[0, 1], xticklabels=["Forasteiro", "Hybrid"], xlabel="Observed bean type", ylabel="Classification score", title=title)
        axis.text(.96, .95, f"n={len(frame)}\naccuracy={accuracy:.2f}", transform=axis.transAxes, ha="right", va="top", fontsize=6.3)
        axis.spines[["top", "right"]].set_visible(False)
    save(fig, "bracken-bean-classification-cv-penalized-three-panel", "CV-selected penalized species coefficients and training/held-out classification scores", [MODEL, CLR, METADATA, TAXONOMY], {"selection": "1SE", "refit": False, "score_source": "intercept + CLR design matrix times the CV-selected penalized coefficient vector", "train_n": 10, "test_n": 4, "c_lasso_source": "bio-datascience/c-lasso@d22da8d5ea54eb499de2fce82fa6c0e1f0c2dcc0"})


def main() -> None:
    metadata = pd.read_csv(METADATA, sep="\t", skiprows=[1]).set_index("sample-id")
    labels = taxonomy_labels()
    model, temporary = open_problem(MODEL)
    cv = model["solution/CV"]
    curve = pd.DataFrame({"lambda_over_lambda_max": np.asarray(cv["xGraph"]), "error": np.asarray(cv["yGraph"]), "se": np.asarray(cv["standard_error"])})
    curve["minimum"] = np.isclose(curve.lambda_over_lambda_max, float(cv.attrs["lambda_min"])); curve["one_se"] = np.isclose(curve.lambda_over_lambda_max, float(cv.attrs["lambda_1SE"]))
    curve.to_csv(EXP / "bracken-bean-classification-cv-curve.tsv", sep="\t", index=False)
    selected = np.asarray(cv["selected_param"]).astype(bool); names = np.asarray(cv["label"]).astype(str); beta = np.asarray(cv["beta"])
    coef = pd.DataFrame({"feature_id": names[selected], "penalized_coefficient": beta[selected]}).query("feature_id != 'intercept'")
    coef["display_label"] = coef.feature_id.map(labels).fillna(coef.feature_id)
    coef = coef.set_index("display_label").loc[[label for label, _ in signed_descending(list(zip(coef.display_label, coef.penalized_coefficient)))]].reset_index()
    coef.to_csv(EXP / "bracken-bean-classification-cv-coefficients.tsv", sep="\t", index=False)
    scores = penalized_scores(MODEL, clr_table(), metadata)
    train, test = scores.query("analysis_set == 'train'").copy(), scores.query("analysis_set == 'test'").copy()
    pd.concat([train.assign(partition="train"), test.assign(partition="test")]).to_csv(EXP / "bracken-bean-classification-cv-scores.tsv", sep="\t", index=False)
    plot_data(clr_table(), metadata, labels); plot_cv(curve); plot_penalized(coef, train, test)
    temporary.cleanup()


if __name__ == "__main__":
    main()
