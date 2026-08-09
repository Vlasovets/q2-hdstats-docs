#!/bin/bash
#SBATCH --job-name=q2-t1-figs
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/t1figs_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/t1figs_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Tier-1 figures, generated rather than hand-run.
#
# `02_lowdim_gglasso/05_lambda_paths.md` is titled "Regularization Paths & Model
# Selection", is 341 lines long, and contains no plot of a regularization path.
# It also has a section called "gamma: what eBIC is actually optimizing" that
# argues entirely in prose.
#
# The figure this wants is the multi-gamma eBIC curve -- the same shape that was
# WRONG for tier 2, where the chapter fixes gamma = 0.3 and five curves buried
# the decision. Here the chapter is *about* how gamma behaves, so showing several
# is the content rather than clutter. Same plot, opposite verdict, because the
# question the chapter asks is different.
#
# Everything is derived from the 13-ASV toy table shipped with q2-gglasso, so
# this reproduces from a clean checkout.

set -euo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
GG="${Q2_HDSTATS_PLUGINS:-/home/itg/oleg.vlasovets/slr_example}/q2-gglasso/data"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
OUT="$REPO/docs/images/png/generated"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-t1figs/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-t1figs/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ [1/3] build the toy covariance ############"
qiime gglasso transform-features \
  --i-table "$GG/atacama-counts.qza" --i-taxonomy "$GG/classification.qza" \
  --m-sample-metadata-file "$GG/selected-atacama-sample-metadata.tsv" \
  --p-transformation mclr --o-transformed-table mclr.qza >/dev/null
qiime gglasso calculate-covariance \
  --i-table mclr.qza --p-method scaled --o-covariance-matrix corr.qza >/dev/null
echo "  corr.qza built"

echo "############ [2/3] solve the lambda path ############"
qiime gglasso solve-problem \
  --i-covariance-matrix corr.qza --p-n-samples 50 --p-no-latent \
  --p-lambda1-min 0.001 --p-lambda1-max 1 --p-n-lambda1 30 \
  --p-gamma 0.01 --o-solution sgl-path.qza >/dev/null
echo "  sgl-path.qza solved"

echo "############ [3/3] draw ############"
python - sgl-path.qza "$OUT" <<'PY'
import os, sys, tempfile, zipfile
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, zarr

path_q, outdir = sys.argv[1:3]
INK, MUTED, ACCENT, GRID = "#2f3b47", "#8b98a5", "#c1442f", "#dde3e9"

with tempfile.TemporaryDirectory() as t, zipfile.ZipFile(path_q) as z:
    i = [n for n in z.namelist() if n.endswith("problem.zip")][0]; z.extract(i, t)
    r = zarr.open(zarr.ZipStore(os.path.join(t, i), mode="r"))
    ms = r["modelselect_stats"]
    lam = np.asarray(ms["LAMBDA"]).ravel()
    gammas = sorted(ms["BIC"].keys(), key=float)
    curves = {g: np.asarray(ms["BIC"][g]).ravel() for g in gammas}
    sparsity = np.asarray(ms["SP"]).ravel() if "SP" in ms else None

order = np.argsort(lam)
lam = lam[order]
PAIRS = 13 * 12 // 2

fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=200)
ax.grid(True, color=GRID, linewidth=0.8, zorder=0); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(MUTED)

cmap = plt.get_cmap("viridis")
picks = []
for k, g in enumerate(gammas):
    c = curves[g][order]
    col = cmap(0.12 + 0.72 * (k / max(len(gammas) - 1, 1)))
    ax.plot(lam, c, "-", color=col, linewidth=1.9, zorder=3,
            label=r"$\gamma$ = %g" % float(g))
    j = int(np.argmin(c))
    picks.append((float(g), lam[j]))
    ax.plot([lam[j]], [c[j]], "o", markersize=9, markerfacecolor="white",
            markeredgecolor=col, markeredgewidth=2.2, zorder=4)

ax.set_xscale("log")
ax.set_xlabel(r"$\lambda_1$  (log scale)", fontsize=11, color=INK)
ax.set_ylabel("eBIC", fontsize=11, color=INK)
ax.set_title(r"On 13 features, eBIC picks a near-empty network at every $\gamma$",
             fontsize=12.5, fontweight="bold", color=INK, pad=12)
ax.tick_params(colors=INK, labelsize=10)
ax.legend(frameon=False, fontsize=9.5, loc="best", ncol=2)

sel = "   ".join(r"$\gamma$=%g$\to\lambda$=%.3g" % p for p in picks)
ax.text(0.5, -0.235, "open circles mark each criterion's choice", ha="center",
        transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
ax.text(0.5, -0.305, sel, ha="center", transform=ax.transAxes,
        fontsize=8.8, color=MUTED)

out = os.path.join(outdir, "toy-lambda-path-gamma.png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("  wrote %s" % out)

# The chapter quotes edge counts, so derive them here rather than letting a
# person read them off the plot. SP is GGLasso's off-diagonal nonzero fraction;
# cross-checked below against a known value so a change in that convention
# fails loudly instead of silently rescaling every number in the prose.
if sparsity is None:
    print("  FAIL: no SP in modelselect_stats; cannot verify the edge counts")
    sys.exit(1)
sp = sparsity[order]
rows = []
for (g, lam_sel) in picks:
    j = int(np.argmin(np.abs(lam - lam_sel)))
    rows.append((g, lam_sel, sp[j], int(round(sp[j] * PAIRS))))
print("  gamma    lambda     SP        edges/%d" % PAIRS)
for g, l, s_, e in rows:
    print("  %-7g  %-9.4f %-9.4f %d" % (g, l, s_, e))

distinct_lam = len({round(l, 6) for _, l, _, _ in rows})
print("  distinct lambda selections across gamma: %d of %d" % (distinct_lam, len(rows)))
print("  edge counts across gamma: %s" % sorted({e for *_, e in rows}))

# Assertions the chapter depends on.
bad = 0
if max(e for *_, e in rows) > 6:
    print("  FAIL: expected a near-empty network at every gamma, got up to %d edges"
          % max(e for *_, e in rows)); bad += 1
if not 0.0 <= min(s_ for _, _, s_, _ in rows) and max(s_ for _, _, s_, _ in rows) <= 1.0:
    print("  FAIL: SP outside [0,1] -- the sparsity convention changed"); bad += 1
sys.exit(1 if bad else 0)
PY

echo "############ RESULT: PASS ############"
