#!/usr/bin/env python3
"""Create QIIME 2 inputs from the complete merged mOTU result.

Run this only with the established q2-hdstats-qiime2-2026.7 interpreter.
The raw table is intentionally used here so removal of mOTUs' ``unassigned``
complement is explicit and recorded before QIIME 2 filtering begins.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import biom
import pandas as pd
import qiime2


SOURCE = Path("/Users/cmueller/Documents/GitHub/q2-hdstats-docs")
RAW_TABLE = SOURCE / "analysis/results/motus-merged/table-raw-with-unassigned.tsv"
TAXONOMY = SOURCE / "analysis/results/motus-merged/taxonomy-378.tsv"
METADATA = SOURCE / "analysis/config/map13241-fecal-wgs.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(output_dir: Path) -> None:
    artifacts = output_dir / "artifacts"
    exports = output_dir / "exports"
    artifacts.mkdir(parents=True, exist_ok=True)
    exports.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(RAW_TABLE, sep="\t", index_col=0)
    if "unassigned" not in table.columns:
        raise ValueError("The complete raw mOTU table must contain 'unassigned'.")
    table = table.drop(columns="unassigned")
    if table.shape != (36, 378):
        raise ValueError(f"Unexpected profiled table shape after removal: {table.shape}")
    if not all(float(value).is_integer() for value in table.to_numpy().ravel()):
        raise ValueError("mOTU counts must be integral for FeatureTable[Frequency].")

    biom_table = biom.Table(
        table.to_numpy(dtype=int).T,
        observation_ids=table.columns.astype(str),
        sample_ids=table.index.astype(str),
    )
    qiime2.Artifact.import_data("FeatureTable[Frequency]", biom_table).save(
        artifacts / "motus-complete-profiled.qza"
    )

    taxonomy = pd.read_csv(TAXONOMY, sep="\t", index_col=0)
    if set(taxonomy.index) != set(table.columns):
        raise ValueError("Taxonomy IDs do not match profiled mOTU IDs.")
    qiime2.Artifact.import_data("FeatureData[Taxonomy]", taxonomy).save(
        artifacts / "motus-complete-taxonomy.qza"
    )

    metadata = pd.read_csv(METADATA, sep="\t", dtype=str).set_index("run_accession")
    selected = metadata.loc[table.index, ["host_age_days", "host_subject_id", "diagnosis"]]
    # Fixed subject-level holdout: two asthma and two control infants, with
    # age coverage across the observed range. No infant is represented in both
    # training and test sets.
    test_subjects = {"J10", "J17", "J12", "J22"}
    selected["split"] = selected["host_subject_id"].map(
        lambda subject: "test" if subject in test_subjects else "train"
    )
    selected.index.name = "sample-id"
    selected.to_csv(exports / "motus-outcomes-all36.tsv", sep="\t")

    manifest = output_dir / "source-manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["source", "sha256", "purpose"])
        writer.writerow([RAW_TABLE, sha256(RAW_TABLE), "complete raw table; unassigned removed explicitly"])
        writer.writerow([TAXONOMY, sha256(TAXONOMY), "taxonomy for TRAC aggregation"])
        writer.writerow([METADATA, sha256(METADATA), "age outcome and subject metadata"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: prepare_complete_motus_inputs.py OUTPUT_DIRECTORY")
    main(Path(sys.argv[1]).resolve())
