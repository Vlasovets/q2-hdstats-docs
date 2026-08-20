#!/usr/bin/env python
"""Merge the 36 per-sample mOTUs profiles and answer Gate M3.

Gate M3 asks one question: does this dataset, profiled with mOTUs, support the
covariance estimate the network chapter would need? The smoke test on two samples
suggested no -- 705 assigned marker-gene inserts from a 2.0 M-read run -- but two
samples cannot settle it. This merges all 36 and measures.

Reads the BIOM payload straight out of each `.qza` with h5py rather than going through
`qiime2.Artifact.load`. That is deliberate: building a QIIME 2 PluginManager pulls in
q2_composition -> rpy2 -> R, which on this cluster's login node dies with
"9 arguments passed to .Internal(list.files) which requires 8". A `.qza` is a zip with a
BIOM inside; nothing here needs the framework.

Writes:
  analysis/results/tables/motus-profile-summary.tsv    per-sample depth and richness
  analysis/results/tables/motus-feature-prevalence.tsv per-feature prevalence and totals
  analysis/results/motus-merged/table.tsv              the merged samples x features table

Usage:
    python analysis/scripts/aggregate_motus_profiles.py
"""
import argparse
import csv
import io
import pathlib
import sys
import zipfile

import h5py
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_IN = ROOT / "results" / "motus-profiles"
DEFAULT_TABLES = ROOT / "results" / "tables"
DEFAULT_MERGED = ROOT / "results" / "motus-merged"
SHEET = ROOT / "config" / "map13241-fecal-wgs.tsv"

# The literal row mOTUs emits for reads it could not assign. It is the complement of
# everything profiled, so it must not enter a log-ratio: including it makes the
# "composition" a mixture of taxa and not-taxa.
UNASSIGNED = "unassigned"


def _member(qza: pathlib.Path, suffix: str) -> bytes:
    with zipfile.ZipFile(qza) as z:
        names = [n for n in z.namelist() if n.endswith(suffix)]
        if not names:
            sys.exit(f"{qza}: no member ending {suffix}")
        return z.read(names[0])


def read_biom(qza: pathlib.Path) -> pd.Series:
    """Return one sample's counts as a Series indexed by feature id."""
    with h5py.File(io.BytesIO(_member(qza, "/data/feature-table.biom")), "r") as f:
        ids = [i.decode() for i in f["observation/ids"][:]]
        data = f["observation/matrix/data"][:]
        indptr = f["observation/matrix/indptr"][:]
        vals = [float(data[indptr[i]:indptr[i + 1]].sum()) for i in range(len(ids))]
    return pd.Series(vals, index=ids, dtype=float)


def read_tax(qza: pathlib.Path) -> dict:
    txt = _member(qza, "/data/taxonomy.tsv").decode()
    rows = list(csv.reader(io.StringIO(txt), delimiter="\t"))
    return {r[0]: r[1] for r in rows if r and r[0] != "Feature ID"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=pathlib.Path, default=DEFAULT_IN)
    ap.add_argument("--min-prevalence", type=int, default=3,
                    help="feature must be nonzero in at least this many samples")
    args = ap.parse_args()

    tables = sorted(args.indir.glob("*-table.qza"))
    if not tables:
        sys.exit(f"no *-table.qza under {args.indir}")

    meta = {r["run_accession"]: r for r in
            csv.DictReader(open(SHEET), delimiter="\t")} if SHEET.is_file() else {}

    series, tax = {}, {}
    for p in tables:
        run = p.name[: -len("-table.qza")]
        series[run] = read_biom(p)
        tp = p.with_name(f"{run}-taxonomy.qza")
        if tp.is_file():
            tax.update(read_tax(tp))

    mat = pd.DataFrame(series).T.fillna(0.0)          # samples x features
    mat.index.name = "run_accession"
    print(f"  merged: {mat.shape[0]} samples x {mat.shape[1]} features "
          "(union, including 'unassigned')")

    unassigned = mat[UNASSIGNED] if UNASSIGNED in mat.columns else pd.Series(0.0, index=mat.index)
    bio = mat.drop(columns=[UNASSIGNED], errors="ignore")
    print(f"  after dropping '{UNASSIGNED}': {bio.shape[1]} features")

    depth = bio.sum(axis=1)
    rich = (bio > 0).sum(axis=1)
    prev = (bio > 0).sum(axis=0)
    totals = bio.sum(axis=0)

    DEFAULT_TABLES.mkdir(parents=True, exist_ok=True)
    DEFAULT_MERGED.mkdir(parents=True, exist_ok=True)

    per_sample = pd.DataFrame({
        "run_accession": bio.index,
        "sample_name": [meta.get(r, {}).get("sample_name", "") for r in bio.index],
        "host_subject_id": [meta.get(r, {}).get("host_subject_id", "") for r in bio.index],
        "diagnosis": [meta.get(r, {}).get("diagnosis", "") for r in bio.index],
        "host_age_days": [meta.get(r, {}).get("host_age_days", "") for r in bio.index],
        "read_count": [meta.get(r, {}).get("read_count", "") for r in bio.index],
        "assigned_counts": depth.values,
        "unassigned_counts": unassigned.reindex(bio.index).values,
        "richness": rich.values,
    })
    per_sample.to_csv(DEFAULT_TABLES / "motus-profile-summary.tsv", sep="\t", index=False)

    feat = pd.DataFrame({
        "feature_id": bio.columns,
        "prevalence": prev.values,
        "total_counts": totals.values,
        "max_in_one_sample": bio.max(axis=0).values,
        "taxon": [tax.get(f, "") for f in bio.columns],
    }).sort_values(["prevalence", "total_counts"], ascending=False)
    feat.to_csv(DEFAULT_TABLES / "motus-feature-prevalence.tsv", sep="\t", index=False)
    bio.to_csv(DEFAULT_MERGED / "table.tsv", sep="\t")

    print("\n  --- per-sample assigned depth ---")
    print(f"    min {depth.min():.0f}   median {depth.median():.0f}   max {depth.max():.0f}")
    print(f"    samples with < 100 assigned counts: {(depth < 100).sum()} of {len(depth)}")
    print("  --- per-sample richness (nonzero features) ---")
    print(f"    min {rich.min()}   median {rich.median():.0f}   max {rich.max()}")
    print("  --- feature prevalence ---")
    for k in (1, 2, 3, 5, 10, 18):
        print(f"    features present in >= {k:2d} of {bio.shape[0]} samples: {(prev >= k).sum()}")

    kept = bio.loc[:, prev >= args.min_prevalence]
    nz_frac = float((kept.to_numpy() > 0).mean()) if kept.shape[1] else 0.0
    print(f"\n  --- after prevalence filter (>= {args.min_prevalence} samples) ---")
    print(f"    p = {kept.shape[1]}   n = {kept.shape[0]}   nonzero cells = {100*nz_frac:.1f}%")
    if kept.shape[1]:
        km = kept.to_numpy()
        print(f"    counts: median {np.median(km):.0f}, "
              f"cells equal to 1: {(km == 1).sum()} of {km.size}")

    print("\n  --- GATE M3 ---")
    verdict_fail = kept.shape[1] < 20 or kept.shape[0] < 30
    if verdict_fail:
        print(f"    FAIL  p={kept.shape[1]} (need >= 20), n={kept.shape[0]} (need >= 30)")
    else:
        print(f"    p={kept.shape[1]}, n={kept.shape[0]} clear the stated thresholds")
    print(f"    -> {DEFAULT_TABLES / 'motus-profile-summary.tsv'}")
    print(f"    -> {DEFAULT_TABLES / 'motus-feature-prevalence.tsv'}")
    sys.exit(1 if verdict_fail else 0)


if __name__ == "__main__":
    main()
