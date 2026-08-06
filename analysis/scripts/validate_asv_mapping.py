#!/usr/bin/env python
"""Does the documented `ASV-k` -> feature-ID recovery actually work?

`04_highdim_atacama/06_interpretation.md` tells readers that when a bundle was
built without `--p-keep-original-id`, they can recover feature identity from
`top-300-asvs.tsv` via

    ASV-n  <->  abundance-rank = 301 - n

That is only sound if `abundance-rank` breaks ties the same way the plugin's
`rename_index_with_sum` did when it assigned the labels. It does not obviously:
209 of the 300 features share a total-abundance value with another feature
(61 groups), so `abundance-rank` imposes an ARBITRARY total order on them.

This tests the claim end to end. The shipped correlation matrix is the
regenerated one reordered (verified exactly, residual 0.000e+00). So if route 2
describes the shipped ordering, permuting the regenerated matrix by route 2's
mapping must reproduce the shipped matrix exactly.

Usage:
    python scripts/validate_asv_mapping.py SHIPPED.qza REGEN.qza RAW.qza ASVS.tsv OUT.md
"""
import os
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd


def load_tsv(qza):
    with tempfile.TemporaryDirectory() as tmp, zipfile.ZipFile(qza) as z:
        name = [x for x in z.namelist() if x.endswith(".tsv") and "/data/" in x][0]
        z.extract(name, tmp)
        return pd.read_csv(os.path.join(tmp, name), sep="\t", index_col=0)


def main():
    ship_q, regen_q, raw_q, asvs_tsv, out = sys.argv[1:6]

    import biom
    import qiime2

    A = load_tsv(ship_q).to_numpy(float)
    B = load_tsv(regen_q).to_numpy(float)
    n = A.shape[0]

    raw = qiime2.Artifact.load(raw_q).view(biom.Table).to_dataframe()  # features x samples
    tot = raw.sum(axis=1)

    # Order the CURRENT plugin produces: ascending total abundance, stable sort,
    # so ties keep the input table's order. Position k <-> ASV-(k+1).
    regen_order = list(tot.sort_values(kind="stable").index)
    pos_in_regen = {fid: k for k, fid in enumerate(regen_order)}

    # Order route 2 claims the shipped bundle uses.
    tsv = pd.read_csv(asvs_tsv, sep="\t")
    tsv = tsv[tsv.iloc[:, 0] != "#q2:types"].copy()
    tsv.columns = [c.strip() for c in tsv.columns]
    fid_col = tsv.columns[0]
    tsv["abundance-rank"] = tsv["abundance-rank"].astype(int)
    tsv["total-abundance"] = tsv["total-abundance"].astype(float)
    route2_order = list(
        tsv.sort_values("abundance-rank", ascending=False)[fid_col]
    )  # rank 300 (least abundant) first == ASV-1

    L = ["# Does the documented `ASV-k` recovery actually work?", "",
         "`06_interpretation.md` route 2 maps `ASV-n` to `abundance-rank = 301 - n`",
         "using `top-300-asvs.tsv`. Testing that against the artifacts.", ""]

    missing = [f for f in route2_order if f not in pos_in_regen]
    if missing or len(route2_order) != n:
        L += [f"- route-2 order has {len(route2_order)} ids, "
              f"{len(missing)} not found in the table — cannot test.", ""]
        open(out, "w").write("\n".join(L) + "\n"); print("\n".join(L)); return 1

    # sanity: does the TSV's abundance agree with the table's?
    tsv_tot = tsv.set_index(fid_col)["total-abundance"]
    agree = np.allclose(tsv_tot.reindex(regen_order).to_numpy(float),
                        tot.reindex(regen_order).to_numpy(float))

    Q = np.array([pos_in_regen[f] for f in route2_order])
    resid = float(np.max(np.abs(A - B[np.ix_(Q, Q)])))
    works = resid < 1e-12

    # how far off is it, and is the damage confined to ties?
    n_disagree = int((Q != np.arange(n)).sum())
    uniq, counts = np.unique(tot.to_numpy(float), return_counts=True)
    tied_vals = set(uniq[counts > 1].tolist())
    unique_feats = [f for f in regen_order if float(tot[f]) not in tied_vals]
    # of the features with a UNIQUE abundance, does route 2 place them correctly?
    ship_pos_route2 = {f: k for k, f in enumerate(route2_order)}

    L += ["## Result", "",
          f"- TSV total-abundance matches the table : {agree}",
          f"- permuting the regenerated matrix by route 2 reproduces the shipped "
          f"matrix to : **{resid:.3e}**",
          f"- positions where route 2 differs from the current plugin order : "
          f"{n_disagree} / {n}",
          f"- features with a UNIQUE total abundance (route 2 is safe for these) : "
          f"**{len(unique_feats)} / {n}**",
          f"- features sharing an abundance value (route 2 is a guess for these) : "
          f"**{n - len(unique_feats)} / {n}**", ""]

    if works:
        L += ["**Route 2 is correct for the shipped bundle.** The `abundance-rank` "
              "column happens to break ties exactly the way the plugin did when it "
              "built these artifacts.", "",
              "That is a coincidence of provenance, not a guarantee. `abundance-rank` "
              "imposes an arbitrary total order on the "
              f"{n - len(unique_feats)} tied features, and nothing ties it to the "
              "sort inside `rename_index_with_sum`. It already diverged once: the "
              "stable-sort change in `9a4c08b` moved 158 features. Anyone applying "
              "route 2 to an artifact built with the current plugin would "
              "mis-assign those.", ""]
    else:
        L += ["**Route 2 does NOT recover the shipped ordering.** Readers following "
              "the documented procedure silently attach the wrong taxon to tied "
              f"features — up to {n - len(unique_feats)} of {n}. The join raises "
              "nothing, so the error is invisible.", ""]

    L += ["## Recommendation", "",
          "Route 1 (`--p-keep-original-id`, now the default) is the only reliable "
          "path and should be presented as such rather than as one of two options. "
          "Route 2 should carry an explicit warning that it is exact only for the "
          f"{len(unique_feats)} features whose total abundance is unique.", "",
          "The durable fix is in the plugin: break ties on the feature ID so that "
          "`ASV-k` becomes a pure function of the data rather than of the input "
          "row order — `df.sum(axis=1).sort_values(kind='stable')` replaced by a "
          "sort on `(row_sum, feature_id)`.", "",
          "## Verdict", "",
          f"**{'Route 2 works for THIS bundle but is not robust' if works else 'Route 2 is broken — silently mislabels tied features'}**"]

    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
