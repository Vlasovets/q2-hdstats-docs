#!/usr/bin/env python3
"""Assemble the tutorial data release and fill in the manifest.

Collects the tier 1 and tier 2 artifacts — which currently live scattered across
two plugin repos and the recompute tree — into one `publish/` directory, computes
real sizes and checksums, and rewrites `docs/_data/manifest.tsv` with them.

Only the **URL** column stays `ZENODO_DOI_PENDING`: bytes and sha256 are facts
about files we have, so there is no reason for a reader to wait on a DOI to
verify a download. Tier 3 rows already carry real values and are left untouched.

Run from anywhere:
    python scripts/package_release.py [--tar]
"""

import argparse
import hashlib
import pathlib
import shutil
import sys
import tarfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
GG = pathlib.Path("/home/itg/oleg.vlasovets/slr_example/q2-gglasso/data")
CL = pathlib.Path("/home/itg/oleg.vlasovets/slr_example/q2-classo/data")
DOCS = pathlib.Path("/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs")
REGEN = ROOT / "results" / "tier2-regen"

# tier -> [(published filename, source directory)] or
#         [(published filename, source directory, source filename)] when the
#         file is named differently at its source.
LAYOUT = {
    1: [
        ("atacama-counts.qza", GG),
        ("classification.qza", GG),
        ("selected-atacama-sample-metadata.tsv", GG),
        # NB: this one lives in q2-classo, not q2-gglasso.
        ("atacama-selected-covariates-veg.tsv", CL),
    ],
    2: [
        ("atacama-top-300-table.qza", ROOT / "data"),
        # The two derived tier-2 artifacts come from the REGENERATION, not from
        # data/. data/ still holds the June 2026 copies, whose features are
        # positional ASV-k labels rather than real feature IDs; they are kept as
        # the historical record and as the reference the permutation analysis was
        # run against, and must not be republished. Sourcing these two from data/
        # would silently revert the bundle and rewrite manifest.tsv back to the
        # old checksums while still exiting 0.
        ("atacama-top-300-clr.qza", REGEN, "clr.qza"),
        ("atacama-top-300-correlation.qza", REGEN, "correlation.qza"),
        ("atacama-taxonomy-silva138.qza", ROOT / "data"),
        ("sample-metadata.tsv", ROOT / "data"),
        ("top-300-asvs.tsv", ROOT / "data"),
        ("atacama-classo-outcomes-mean-imputed.tsv", ROOT / "data"),
    ],
}

PENDING = "ZENODO_DOI_PENDING"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tar", action="store_true", help="also write per-tier tarballs")
    args = ap.parse_args()

    pub = ROOT / "publish"
    collected = {}
    missing = []

    for tier, entries in LAYOUT.items():
        dest = pub / f"tier{tier}"
        dest.mkdir(parents=True, exist_ok=True)
        for entry in entries:
            name, src_dir = entry[0], entry[1]
            src_name = entry[2] if len(entry) > 2 else name
            src = pathlib.Path(src_dir) / src_name
            if not src.is_file():
                missing.append(f"tier{tier}/{name} (looked for {src_name} in {src_dir})")
                continue
            shutil.copy2(src, dest / name)
            collected[name] = (tier, (dest / name).stat().st_size, sha256(dest / name))

    if missing:
        print("MISSING — not packaged:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)

    # --- rewrite the manifest, preserving tier 3 rows verbatim ----------------
    man = DOCS / "docs" / "_data" / "manifest.tsv"
    rows, header = [], None
    for i, line in enumerate(man.read_text().splitlines()):
        parts = line.split("\t")
        if i == 0:
            header = parts
            continue
        if not parts or not parts[0]:
            continue
        name, tier = parts[0], parts[1]
        if tier == "3":
            rows.append(parts)                      # already real, leave alone
        elif name in collected:
            t, size, digest = collected[name]
            rows.append([name, str(t), str(size), digest, PENDING])
        else:
            rows.append(parts)

    # add anything collected that the manifest did not already list
    listed = {r[0] for r in rows}
    for name, (t, size, digest) in sorted(collected.items()):
        if name not in listed:
            rows.append([name, str(t), str(size), digest, PENDING])

    rows.sort(key=lambda r: (int(r[1]), r[0]))
    man.write_text("\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")

    n_real = sum(1 for r in rows if r[3] != PENDING)
    print(f"packaged {len(collected)} files into {pub}")
    print(f"manifest: {len(rows)} rows, {n_real} with real checksums, "
          f"{len(rows) - n_real} still pending a DOI")

    if args.tar:
        for tier in LAYOUT:
            tar_path = pub / f"q2-hdstats-tutorial-data-tier{tier}.tar.gz"
            with tarfile.open(tar_path, "w:gz") as tf:
                tf.add(pub / f"tier{tier}", arcname=f"tier{tier}")
            print(f"  {tar_path.name}  {tar_path.stat().st_size} bytes")

    if missing:
        sys.exit(1)


if __name__ == "__main__":
    main()
