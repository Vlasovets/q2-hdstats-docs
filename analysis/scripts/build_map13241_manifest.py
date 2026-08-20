#!/usr/bin/env python
"""Build the sample sheet for the MAP study (Qiita 13241) shotgun subset.

Joins two public sources and writes one committed TSV that every downstream stage reads,
so nothing has to re-derive which runs to fetch:

  * Qiita sample information (98 samples x 197 columns) -- the clinical metadata
  * ENA read_run report for PRJEB52147 (192 runs) -- the FASTQ URLs, sizes and MD5s

Three things here are not obvious and each one silently corrupts the result if missed.

1. THE QIITA ZIP LISTS THE SAME TEMPLATE MEMBER TWICE. Two ZipInfo entries, identical
   name, size and CRC. Iterating `infolist()` yields 196 rows instead of 98 and every
   count downstream doubles. This is not hypothetical -- it produced a wrong sample
   count during planning, and the wrong number survived into a written summary before
   anyone re-derived it. We deduplicate by name and assert the row count.

2. THE JOIN KEY IS THE STRIPPED ENA sample_alias. ENA carries
   `qiita_sid_13241:<sample_name>`; the Qiita template carries `<sample_name>`.
   `secondary_sample_accession` (ERS...) appears nowhere in the template, so joining on
   it silently yields zero rows.

3. THE STUDY IS HALF 16S. PRJEB52147 has 192 runs: 96 AMPLICON/SINGLE on MiSeq and 96
   WGS/PAIRED on NovaSeq. mOTUs can only profile the latter. The WGS discriminator used
   here is `experiment_alias` starting with `qiita_ptid_10582:`, which is set-equal to
   `library_strategy == "WGS"` but says explicitly which Qiita preparation it is.

Usage:
    python analysis/scripts/build_map13241_manifest.py [--out PATH] [--sample-type feces]
"""
import argparse
import csv
import io
import pathlib
import sys
import urllib.request
import zipfile

QIITA_SAMPLE_INFO = (
    "https://qiita.ucsd.edu/public_download/?data=sample_information&study_id=13241"
)
ENA_FIELDS = (
    "run_accession,experiment_alias,sample_alias,secondary_sample_accession,"
    "library_strategy,library_layout,instrument_model,read_count,base_count,"
    "fastq_bytes,fastq_md5,fastq_ftp"
)
ENA_REPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB52147"
    f"&result=read_run&format=tsv&limit=0&fields={ENA_FIELDS}"
)
WGS_PREP_PREFIX = "qiita_ptid_10582:"
QIITA_ALIAS_PREFIX = "qiita_sid_13241:"

EXPECTED_TEMPLATE_ROWS = 98
EXPECTED_ENA_RUNS = 192

# Columns kept in the sample sheet. Everything a downstream stage or chapter needs, and
# nothing else -- the template has 197 columns and most are irrelevant here.
KEEP = [
    "sample_name", "sample_type", "host_subject_id", "diagnosis", "timepoint",
    "host_age_days", "gestation_age_at_birth_weeks", "birth_weight_g", "infant_sex",
    "delivery_method",
]


def _get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=300) as fh:
        return fh.read()


def load_qiita_template() -> list[dict]:
    """Read the sample template, guarding the duplicate-member trap."""
    zf = zipfile.ZipFile(io.BytesIO(_get(QIITA_SAMPLE_INFO)))
    names = sorted(set(zf.namelist()))
    if len(names) != 1:
        sys.exit(f"expected exactly one template member, got {names}")
    if len(zf.namelist()) != len(names):
        print(f"  note: zip lists {len(zf.namelist())} members for "
              f"{len(names)} unique name(s) -- deduplicated", file=sys.stderr)
    with zf.open(names[0]) as fh:
        rows = list(csv.DictReader(io.TextIOWrapper(fh, "utf-8"), delimiter="\t"))
    if len(rows) != EXPECTED_TEMPLATE_ROWS:
        sys.exit(f"expected {EXPECTED_TEMPLATE_ROWS} template rows, got {len(rows)}. "
                 "If Qiita revised the study this is a real change, not a bug -- "
                 "re-check the counts quoted in the chapters before bumping this.")
    return rows


def load_ena_runs() -> list[dict]:
    text = _get(ENA_REPORT).decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(rows) != EXPECTED_ENA_RUNS:
        print(f"  warning: expected {EXPECTED_ENA_RUNS} ENA runs, got {len(rows)}",
              file=sys.stderr)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path(__file__).resolve().parents[1]
                    / "config" / "map13241-fecal-wgs.tsv")
    ap.add_argument("--sample-type", default="feces",
                    help="sample_type to keep (default: feces)")
    args = ap.parse_args()

    template = load_qiita_template()
    runs = load_ena_runs()
    by_name = {r["sample_name"]: r for r in template}
    print(f"  template rows: {len(template)}   ENA runs: {len(runs)}")

    wgs = [r for r in runs if r["experiment_alias"].startswith(WGS_PREP_PREFIX)]
    print(f"  WGS runs (prep {WGS_PREP_PREFIX.rstrip(':')}): {len(wgs)}")

    out_rows, unmatched = [], []
    for r in wgs:
        alias = r["sample_alias"]
        key = alias.split(":", 1)[1] if alias.startswith(QIITA_ALIAS_PREFIX) else alias
        meta = by_name.get(key)
        if meta is None:
            unmatched.append(alias)
            continue
        if meta.get("sample_type", "").strip() != args.sample_type:
            continue
        ftp = r["fastq_ftp"].split(";")
        md5 = r["fastq_md5"].split(";")
        size = [int(b) for b in r["fastq_bytes"].split(";") if b]
        if len(ftp) != 2:
            print(f"  warning: {r['run_accession']} has {len(ftp)} fastq files, "
                  "expected 2 (paired)", file=sys.stderr)
        row = {k: meta.get(k, "") for k in KEEP}
        row.update({
            "run_accession": r["run_accession"],
            "read_count": r["read_count"],
            "base_count": r["base_count"],
            "fastq_1_url": "https://" + ftp[0] if ftp else "",
            "fastq_2_url": "https://" + ftp[1] if len(ftp) > 1 else "",
            "fastq_1_md5": md5[0] if md5 else "",
            "fastq_2_md5": md5[1] if len(md5) > 1 else "",
            "bytes_total": sum(size),
        })
        out_rows.append(row)

    if unmatched:
        print(f"  warning: {len(unmatched)} WGS runs did not join to the template: "
              f"{unmatched[:3]}", file=sys.stderr)

    out_rows.sort(key=lambda r: r["run_accession"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(out_rows[0]) if out_rows else []
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    subjects = {r["host_subject_id"] for r in out_rows}
    total = sum(r["bytes_total"] for r in out_rows)
    print(f"  -> {args.out}")
    print(f"  {len(out_rows)} runs, {len(subjects)} subjects, "
          f"{total:,} bytes ({total / 2**30:.2f} GiB)")
    dx = {}
    for r in out_rows:
        dx[r["diagnosis"]] = dx.get(r["diagnosis"], 0) + 1
    print(f"  diagnosis: {dx}")
    per_subject = {}
    for r in out_rows:
        per_subject.setdefault(r["host_subject_id"], []).append(r["run_accession"])
    sizes = sorted({len(v) for v in per_subject.values()})
    print(f"  samples per subject: {sizes}")
    ages = [float(r["host_age_days"]) for r in out_rows
            if r["host_age_days"] not in ("", "na", "not applicable")]
    if ages:
        print(f"  host_age_days: n={len(ages)} min={min(ages):g} max={max(ages):g}")


if __name__ == "__main__":
    main()
