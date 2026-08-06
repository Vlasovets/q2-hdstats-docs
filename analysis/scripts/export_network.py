#!/usr/bin/env python3
"""Export a GGLassoProblem solution as the edge/node TSVs the figure scripts read.

`compare_rank0_rank5.py` (from Christian's share bundle) consumes
`{stem}-edges.tsv` and `{stem}-nodes.tsv`, but the script that produced them was
not in the bundle. This reconstructs that step from the artifacts.

Contract, reverse-engineered from the consumer:
  {stem}-edges.tsv : source, target, partial_correlation
  {stem}-nodes.tsv : asv (index), genus, color, node_size

Partial correlations come from the precision matrix Theta as
    rho_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)
which is the standard conversion; the sign flip is why a positive partial
correlation corresponds to a negative off-diagonal precision entry.

Usage:
    export_network.py --solution S.qza --taxonomy T.qza --clr C.qza \
                      --stem results/gglasso/atacama-top-300-network-sgl-lambda0.8
"""

import argparse
import colorsys
import os
import tempfile
import zipfile

import numpy as np
import pandas as pd
import zarr


EDGE_TOL = 1e-8


def _open_problem(qza):
    """Return the zarr root of the problem.zip inside a .qza."""
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z:
        inner = [n for n in z.namelist() if n.endswith("problem.zip")]
        if not inner:
            raise SystemExit(f"{qza}: no problem.zip inside the artifact")
        z.extract(inner[0], tmp)
    return zarr.open(zarr.ZipStore(os.path.join(tmp, inner[0]), mode="r"))


def _labels(root, p):
    """Feature labels, in solver order.

    The solver stores them under `labels/<i>`; fall back to positional names so
    an artifact written before labels were injected still exports.
    """
    if "labels" not in root:
        return [f"F{i}" for i in range(p)]
    out = []
    for i in range(p):
        try:
            out.append(str(np.asarray(root[f"labels/{i}"]).item()))
        except (KeyError, IndexError):
            out.append(f"F{i}")
    return out


def _read_qza_tsv(qza, member_suffix):
    """Read a single TSV member out of a .qza into a DataFrame."""
    with zipfile.ZipFile(qza) as z:
        hits = [n for n in z.namelist() if n.endswith(member_suffix)]
        if not hits:
            raise SystemExit(f"{qza}: no member ending in {member_suffix}")
        with z.open(hits[0]) as fh:
            return pd.read_csv(fh, sep="\t", index_col=0)


def _genus_of(taxon):
    """Deepest informative rank at or above genus, from a SILVA-style string."""
    if not isinstance(taxon, str):
        return "Unassigned"
    parts = [p.strip() for p in taxon.split(";")]
    for prefix in ("g__", "f__", "o__", "c__", "p__"):
        for p in parts:
            if p.startswith(prefix) and len(p) > len(prefix):
                return p[len(prefix):]
    return "Unassigned"


def _display_name(label, genus, width=6):
    """Readable node label that still identifies the feature exactly.

    Real feature IDs are 32 hex characters, which are illegible as network node
    labels or heatmap axis ticks. Prefixing the deepest informative rank and
    keeping a short slice of the ID gives something like
    ``Rubrobacter (a7b877)`` -- readable, and still resolvable back to one
    feature.

    A 5-character prefix is already unique across the 300-ASV Atacama table; 6
    is used for margin. Uniqueness depends on the feature set, so the caller
    asserts it rather than trusting this default.
    """
    # Case-insensitive: QIIME 2 emits lowercase MD5, but an uppercase hash or a
    # hash from another pipeline would otherwise fall through and put the raw
    # 32-character ID on the figure with no warning.
    if len(label) == 32 and all(c in "0123456789abcdef" for c in label.lower()):
        return f"{genus} ({label[:width]})"
    if len(label) > 12:
        # Not a recognised hash, but still too long to use as a tick label.
        return f"{genus} ({label[:width]}…)"
    return label  # already short: ASV-k, or a human-readable name


def _palette(names):
    """Stable, evenly-spaced hues so the same genus keeps its colour across runs."""
    uniq = sorted(set(names))
    out = {}
    for i, name in enumerate(uniq):
        if name == "Unassigned":
            out[name] = "#bbbbbb"
            continue
        r, g, b = colorsys.hsv_to_rgb(i / max(len(uniq), 1), 0.55, 0.85)
        out[name] = "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solution", required=True)
    ap.add_argument("--taxonomy")
    ap.add_argument("--clr", help="clr/mclr table, for node sizes")
    ap.add_argument(
        "--id-map",
        help="top-300-asvs.tsv: maps abundance-rank -> feature-id. Needed when the "
             "solution carries ASV-n labels (transform-features with "
             "--p-keep-original-id False) but the taxonomy is keyed by feature ID. "
             "Without it the taxonomy join silently yields all-Unassigned.",
    )
    ap.add_argument("--stem", required=True)
    ap.add_argument("--tol", type=float, default=EDGE_TOL)
    args = ap.parse_args()

    root = _open_problem(args.solution)
    theta = np.asarray(root["solution/precision_"])
    p = theta.shape[0]
    labels = _labels(root, p)

    d = np.sqrt(np.diag(theta))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = -theta / np.outer(d, d)
    np.fill_diagonal(rho, 1.0)

    iu = np.triu_indices(p, 1)
    keep = np.abs(theta[iu]) > args.tol
    edges = pd.DataFrame({
        "source": [labels[i] for i in iu[0][keep]],
        "target": [labels[j] for j in iu[1][keep]],
        "partial_correlation": rho[iu][keep],
    })
    edges.to_csv(f"{args.stem}-edges.tsv", sep="\t", index=False)

    # --- nodes -----------------------------------------------------------------
    # Solver labels may be ASV-n rather than feature IDs; build the bridge first.
    lookup = {lab: lab for lab in labels}
    if args.id_map:
        m = pd.read_csv(args.id_map, sep="\t")
        m = m[m.iloc[:, 0].astype(str).str.startswith("#") == False]  # noqa: E712
        if "abundance-rank" in m.columns and "feature-id" in m.columns:
            for fid, rank in zip(m["feature-id"], m["abundance-rank"]):
                try:
                    lookup[f"ASV-{int(float(rank))}"] = str(fid)
                except (TypeError, ValueError):
                    continue

    genus = pd.Series("Unassigned", index=labels, name="genus")
    if args.taxonomy:
        tax = _read_qza_tsv(args.taxonomy, "taxonomy.tsv")
        col = "Taxon" if "Taxon" in tax.columns else tax.columns[0]
        mapped = tax[col].map(_genus_of)
        genus = pd.Series(
            [mapped.get(lookup.get(lab, lab), "Unassigned") for lab in labels],
            index=labels, name="genus",
        )
        n_hit = int((genus != "Unassigned").sum())
        if n_hit == 0:
            raise SystemExit(
                "taxonomy join matched NOTHING. The solution's labels "
                f"(e.g. {labels[0]!r}) do not appear in the taxonomy, whose keys "
                f"look like {list(tax.index[:1])!r}. Pass --id-map "
                "top-300-asvs.tsv, or rebuild the clr table with "
                "--p-keep-original-id. Refusing to emit an all-Unassigned table."
            )
        if n_hit < len(labels):
            print(f"  warning: taxonomy matched {n_hit}/{len(labels)} features")

    size = pd.Series(10.0, index=labels, name="node_size")
    mean_clr = pd.Series(np.nan, index=labels, name="mean_clr")
    if args.clr:
        clr = _read_qza_tsv(args.clr, "feature-table.tsv") if args.clr.endswith(".tsv") else None
        if clr is None:
            # biom-backed table: orient so features are rows, then take row means
            import biom
            with zipfile.ZipFile(args.clr) as z:
                hits = [n for n in z.namelist() if n.endswith("feature-table.biom")]
                if hits:
                    tmp = tempfile.mkdtemp()
                    z.extract(hits[0], tmp)
                    tab = biom.load_table(os.path.join(tmp, hits[0])).to_dataframe()
                    if set(labels) & set(tab.columns):
                        tab = tab.T          # stored samples x features
                    means = tab.reindex(labels).mean(axis=1)
                    mean_clr = means.rename("mean_clr")
                    rng = means.max() - means.min()
                    if rng > 0:
                        size = 6 + 18 * (means - means.min()) / rng
                    size = size.fillna(10.0)
                    size.name = "node_size"

    # Degree in THIS network, from the edges just written.
    deg = pd.Series(0, index=labels, name="degree", dtype=int)
    if len(edges):
        counts = pd.concat([edges["source"], edges["target"]]).value_counts()
        deg = counts.reindex(labels).fillna(0).astype(int).rename("degree")

    # Readable labels for figures. The index stays the exact feature ID -- this
    # column is for display only and must never be used as a join key.
    display = pd.Series(
        [_display_name(lab, genus[lab]) for lab in labels],
        index=labels, name="display",
    )
    if display.duplicated().any():
        dup = sorted(display[display.duplicated(keep=False)].unique())[:3]
        raise SystemExit(
            f"display names collide (e.g. {dup}). Widen the ID prefix in "
            "_display_name; two features would otherwise be indistinguishable "
            "on the figure."
        )

    colours = _palette(list(genus))
    nodes = pd.DataFrame({
        "display": display,
        "genus": genus,
        "mean_clr": mean_clr,
        "degree": deg,
        "color": [colours[g] for g in genus],
        "node_size": size,
    })
    nodes.index.name = "asv"
    nodes.to_csv(f"{args.stem}-nodes.tsv", sep="\t")

    print(f"{args.stem}-edges.tsv : {len(edges)} edges")
    print(f"{args.stem}-nodes.tsv : {len(nodes)} nodes, "
          f"{genus.nunique()} distinct genera")


if __name__ == "__main__":
    main()
