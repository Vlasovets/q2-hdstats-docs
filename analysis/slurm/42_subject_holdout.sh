#!/bin/bash
#SBATCH --job-name=q2-subj-holdout
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/subjholdout_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs/analysis/slurm/logs/subjholdout_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#
# Is the diagnosis classifier real, or is it memorising infants?
#
# Stage 41 found that classifying maternal-asthma `diagnosis` from the top-100 mOTUs
# table drops cross-validated error from 0.727 to 0.360 with five selected features --
# a 50% improvement, and by far the strongest signal anything in this study has produced.
#
# That number should not be believed, for a structural reason. `diagnosis` is a property
# of the INFANT, not of the sample, and each infant contributes two samples. q2-classo's
# --p-cv-subsets does a random split with no grouping parameter, so an infant's week-2
# sample lands in training while its week-5 sample lands in testing. The model does not
# have to learn anything about asthma; it only has to recognise that infant. With 18
# infants and a subject-constant label, that is the easiest possible shortcut.
#
# The honest estimate holds out WHOLE INFANTS. This stage does leave-one-subject-out:
# 18 folds, each removing one infant's two samples, fitting on the remaining 17 infants
# and predicting the held-out pair. Nothing about the held-out infant is ever seen.
#
# Both estimates are computed on the same data with the same settings, so the difference
# between them IS the leakage. That comparison is worth more to the book than either
# number alone -- it is a concrete, reproducible demonstration of why grouped splitting
# matters, on a real dataset, which is exactly the kind of thing a methods tutorial
# should show rather than assert.
#
# A baseline is included: predicting the majority class. Accuracy has to be read against
# it, not against 0.5, and with 9 v 9 infants a 95% interval on any accuracy estimate is
# roughly +/- 0.23 -- stated up front so the result is read with the right precision.

set -uo pipefail
REPO="${Q2_HDSTATS_REPO:-/home/itg/oleg.vlasovets/slr_example/q2-hdstats-docs}"
ROOT="$REPO/analysis"
PREFIX=/home/itg/oleg.vlasovets/.conda/envs/q2-2026.7-slr
CONDA=/home/itg/oleg.vlasovets/miniconda3/bin/conda
MERGED="$ROOT/results/motus-merged/table.tsv"
SHEET="$ROOT/config/map13241-fecal-wgs.tsv"
OUT="$ROOT/results/subject-holdout"
TABLES="$ROOT/results/tables"

for cand in "/localscratch/${USER}" "${LOCAL_SCRATCH:-}" "/var/tmp/${USER}"; do
  [[ -z "$cand" ]] && continue
  if mkdir -p "$cand" 2>/dev/null && [[ -w "$cand" ]]; then
    SCRATCH="$cand/q2-subjholdout/${SLURM_JOB_ID:-manual}"; break; fi
done
SCRATCH="${SCRATCH:-/lustre/scratch/users/oleg.vlasovets/q2-subjholdout/${SLURM_JOB_ID:-manual}}"
export TMPDIR="$SCRATCH/tmp"; mkdir -p "$TMPDIR" "$OUT" "$TABLES"
trap '[[ "${KEEP_SCRATCH:-0}" == "1" ]] || rm -rf "$SCRATCH"' EXIT

export PYTHONNOUSERSITE=1
# shellcheck source=/dev/null
source "$(dirname "$CONDA")/../etc/profile.d/conda.sh"; conda activate "$PREFIX"
cd "$SCRATCH"

echo "############ leave-one-subject-out, via the classo python API ############"
# The plugin is the same code path -- q2_classo._func.classify calls into classo -- but
# driving it directly is the only way to control the fold assignment, which the QIIME 2
# interface does not expose. The chapter should say exactly that.
python - "$MERGED" "$SHEET" "$OUT" "$TABLES" <<'PY'
import sys, csv, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
merged, sheet, out, tables = sys.argv[1:5]

# c-lasso 1.0.11 uses np.infty in five places (solve_R1:212, solve_R2:239/:293,
# solve_R3:205, solve_R4:211). NumPy 2.0 REMOVED that alias and this environment runs
# 2.4.2, so any code path computing lambdamax raises AttributeError -- which is what
# killed the regression half of the first run of this stage.
#
# Restoring the alias is a local workaround for a pinned dependency, not a fix. It is
# done here rather than by editing site-packages so the stage stays reproducible from a
# clean environment. The plugin itself does not reach the broken line (eight CLI
# combinations tested), so this only matters for direct API use like the loop below.
if not hasattr(np, "infty"):
    np.infty = np.inf

from classo import classo_problem

t = pd.read_csv(merged, sep="\t", index_col=0)
t = t.loc[t.sum(1) >= 100]
meta = {r["run_accession"]: r for r in csv.DictReader(open(sheet), delimiter="\t")}
base = t.loc[:, (t > 0).sum(0) >= 2]
X_raw = base[list(base.sum(0).sort_values(ascending=False).head(100).index)]
X_raw = X_raw.loc[:, X_raw.sum(0) > 0]

subj = np.array([meta[s]["host_subject_id"] for s in X_raw.index])
dx   = np.array([meta[s]["diagnosis"] for s in X_raw.index])
age  = np.array([float(meta[s]["host_age_days"]) for s in X_raw.index])
print(f"  X: {X_raw.shape[0]} samples x {X_raw.shape[1]} features, "
      f"{len(set(subj))} subjects")

# clr with a 0.5 zero-replacement, matching what q2-classo's transform-features does.
A = X_raw.to_numpy(dtype=float); A[A == 0] = 0.5
L = np.log(A); Xc = L - L.mean(axis=1, keepdims=True)
C = np.ones((1, Xc.shape[1]))                 # the zero-sum constraint

def fit_predict(Xtr, ytr, Xte, classification):
    p = classo_problem(Xtr, ytr, C=C)
    p.formulation.classification = classification
    p.formulation.concomitant = False
    p.model_selection.CV = False
    p.model_selection.StabSel = False
    p.model_selection.LAMfixed = True
    p.model_selection.LAMfixedparameters.rescaled_lam = True
    p.model_selection.LAMfixedparameters.lam = 0.1
    p.solve()
    b = np.asarray(p.solution.LAMfixed.beta).ravel()
    nsel = int((np.abs(b) > 1e-10).sum())
    return Xte @ b, nsel

# ---- classification of diagnosis -------------------------------------------------
y = np.where(dx == "asthma", 1.0, -1.0)
maj = max((y == 1).mean(), (y == -1).mean())
print(f"\n  DIAGNOSIS  majority-class baseline accuracy = {maj:.3f}")

rows = []
# (a) leave-one-SUBJECT-out -- honest
correct, nsel_all = [], []
for s in sorted(set(subj)):
    te = subj == s; tr = ~te
    pred, nsel = fit_predict(Xc[tr], y[tr], Xc[te], True)
    correct.extend(list((np.sign(pred) == y[te]).astype(int))); nsel_all.append(nsel)
acc_subj = float(np.mean(correct))
print(f"  leave-one-SUBJECT-out accuracy = {acc_subj:.3f} "
      f"({int(np.sum(correct))}/{len(correct)} samples), "
      f"median support {int(np.median(nsel_all))}")

# (b) leave-one-SAMPLE-out -- leaky, because the infant's other sample stays in training
correct2 = []
for i in range(len(y)):
    te = np.zeros(len(y), bool); te[i] = True; tr = ~te
    pred, _ = fit_predict(Xc[tr], y[tr], Xc[te], True)
    correct2.append(int(np.sign(pred)[0] == y[i]))
acc_samp = float(np.mean(correct2))
print(f"  leave-one-SAMPLE-out accuracy  = {acc_samp:.3f} "
      f"({int(np.sum(correct2))}/{len(correct2)} samples)   <-- LEAKY")
print(f"  leakage inflation = {acc_samp - acc_subj:+.3f}")
rows.append(dict(outcome="diagnosis", metric="accuracy", baseline=round(maj,3),
                 subject_holdout=round(acc_subj,3), sample_holdout=round(acc_samp,3),
                 leakage=round(acc_samp-acc_subj,3)))

# ---- regression of host_age_days --------------------------------------------------
print(f"\n  HOST_AGE_DAYS  variance baseline (predict the mean) R^2 = 0")
pred_s = np.zeros(len(age))
for s in sorted(set(subj)):
    te = subj == s; tr = ~te
    pr, _ = fit_predict(Xc[tr], age[tr] - age[tr].mean(), Xc[te], False)
    pred_s[te] = pr + age[tr].mean()
ss_res = float(((age - pred_s) ** 2).sum()); ss_tot = float(((age - age.mean()) ** 2).sum())
r2_subj = 1 - ss_res / ss_tot
pred_l = np.zeros(len(age))
for i in range(len(age)):
    te = np.zeros(len(age), bool); te[i] = True; tr = ~te
    pr, _ = fit_predict(Xc[tr], age[tr] - age[tr].mean(), Xc[te], False)
    pred_l[i] = pr[0] + age[tr].mean()
r2_samp = 1 - float(((age - pred_l) ** 2).sum()) / ss_tot
print(f"  leave-one-SUBJECT-out R^2 = {r2_subj:+.3f}")
print(f"  leave-one-SAMPLE-out  R^2 = {r2_samp:+.3f}   <-- LEAKY")
print(f"  leakage inflation = {r2_samp - r2_subj:+.3f}")
rows.append(dict(outcome="host_age_days", metric="R2", baseline=0.0,
                 subject_holdout=round(r2_subj,3), sample_holdout=round(r2_samp,3),
                 leakage=round(r2_samp-r2_subj,3)))

d = pd.DataFrame(rows)
d.to_csv(f"{tables}/subject-holdout-vs-leaky.tsv", sep="\t", index=False)
print("\n" + d.to_string(index=False))
print(f"\n  -> {tables}/subject-holdout-vs-leaky.tsv")
print("\n  --- READ THIS WITH THE RIGHT PRECISION ---")
print("    18 subjects, 9 v 9. A 95% interval on any accuracy here is roughly +/- 0.23,")
print("    so only a very large gap over the baseline would mean anything.")
PY
echo "############ RESULT: complete ############"
