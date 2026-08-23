# The Concomitant Formulation

The concomitant formulation estimates the model coefficients and the noise
level σ jointly, instead of holding σ fixed. Reach for it when the variance is
not homogeneous across your samples.

It applies to `qiime classo regress` only — see
[Robust Classification with the Huber Loss](#robust-classification-with-the-huber-loss)
below for the classification equivalent.

## The two formulations

Set `--p-concomitant True` on `qiime classo regress` and the solver switches
from the standard formulation to the **concomitant formulation**: R3 for least
squares, R4 when combined with the Huber loss.

### Mathematical formulation

**Standard Formulation (R1):**
```
min ||y - X·β||² + λ·||β||₁    s.t. C·β = 0
```
The noise level σ is held fixed rather than estimated.

**Concomitant Formulation (R3, least squares):**
```
min ||y - X·β||² / (2σ) + (n/2)·σ + λ·||β||₁    s.t. C·β = 0
```
Here σ is estimated jointly with the coefficients β, and there is no tuning parameter for σ. The coefficient `n/2`, with `n` the number of samples, is a fixed constant that follows from the perspective function of the squared loss — the plugin does not expose it.

**Concomitant Formulation with Huber loss (R4):**
```
min  h_{ρ,σ}(y - X·β) + n·σ + λ·||β||₁    s.t. C·β = 0
```
Combine `--p-concomitant True` with `--p-huber True` for the robust analogue. As in R3, σ is estimated with no tuning parameter; the σ-coefficient is `n` rather than `n/2`, again from the perspective function of the Huber loss.

### Key differences

| Aspect | Standard (R1) | Concomitant (R3) |
|--------|--------------|-----------------|
| Noise level | Fixed | Jointly estimated |
| Heteroscedasticity | Assumed equal variance | Adaptive to varying noise |
| Uncertainty estimates | May be inaccurate | More reliable |
| Computation | Faster | Slightly slower |
| Output label | "Formulation: R1" | "Formulation: R3 (concomitant)" |

## When to use the concomitant formulation

Use `--p-concomitant True` when:
- Your residual plots show variance increasing with the predicted values (heteroscedasticity)
- Your samples come from different environments with potentially different noise levels
- You need more reliable confidence intervals for your predictions
- You want feature selection that holds up across diverse data
- Computation time is not a binding constraint

## Regression with the concomitant formulation

### Log-contrast with concomitant

Replace the training step with:

```bash
qiime classo regress \
    --i-features data/regress-xtraining_lc.qza \
    --i-c data/ccovariates_lc.qza \
    --i-weights data/wcovariates_lc.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column average-soil-temperature \
    --p-concomitant True \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 1 \
    --p-no-cv-one-se \
    --o-result data/regresstaxa_lc_concomitant.qza
```

Then proceed with prediction and visualization using `regresstaxa_lc_concomitant.qza`.

### trac with concomitant

Replace the training step with:

```bash
qiime classo regress \
    --i-features data/regress-xtraining_trac.qza \
    --i-c data/ccovariates_trac.qza \
    --i-weights data/wcovariates_trac.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column average-soil-temperature \
    --p-concomitant True \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 1 \
    --p-no-cv-one-se \
    --o-result data/regresstaxa_trac_concomitant.qza
```

Then proceed with prediction and visualization using `regresstaxa_trac_concomitant.qza`.

## Robust Classification with the Huber Loss

```{important}
**The concomitant formulation is not available for classification.**
`qiime classo classify` has no `--p-concomitant` flag — the parameter is absent
from its signature, and the underlying `classo_problem` forces
`formulation.concomitant = False` for classification problems. Passing
`--p-concomitant` to `classify` fails with an unrecognised-parameter error.

The robust option for classification is the **Huber hinge loss**, enabled with
`--p-huber True` and tuned through `--p-rho` (default `0.0` for classification,
which is wired to `rho_classification`).
```

### Log-contrast classification with Huber loss

For both [log-contrast classification](../04_classification/01_logcontrast.md) and [trac classification](../04_classification/02_trac.md), replace the training step with:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_lc.qza \
    --i-c data/ccovariates_lc.qza \
    --i-weights data/wcovariates_lc.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column vegetation \
    --p-huber True \
    --p-rho 0.0 \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 42 \
    --p-no-cv-one-se \
    --o-result data/classifytaxa_lc_huber.qza
```

Then proceed with prediction and visualization using `classifytaxa_lc_huber.qza`.

### trac classification with Huber loss

Replace the training step with:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_trac.qza \
    --i-c data/ccovariates_trac.qza \
    --i-weights data/wcovariates_trac.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column vegetation \
    --p-huber True \
    --p-rho 0.0 \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 42 \
    --p-no-cv-one-se \
    --o-result data/classifytaxa_trac_huber.qza
```

Then proceed with prediction and visualization using `classifytaxa_trac_huber.qza`.

## Reading the formulation label

Run a regression with `--p-concomitant True` and the output reports:

```
Formulation: R3 (concomitant)
```

Add `--p-huber True` and the label becomes `Formulation: R4 (concomitant + Huber)`.
Either label tells you that:
- the fit estimates the noise level σ jointly with the coefficients
- σ comes out of the fit rather than out of a hyperparameter, so there is nothing to tune
- the coefficients and predictions are more robust to heteroscedasticity than under the standard formulation

## Practical notes

1. **Compare both formulations.** Run the analysis with `--p-concomitant False` and with `True`, and see which fits your data better

2. **Examine the residuals.** Plot residuals against fitted values before deciding:
   - Constant variance → the standard formulation is sufficient
   - Variance rising or falling with the fitted values → use the concomitant formulation

3. **Cross-validation performance.** Compare CV scores between the formulations; the concomitant fit tends to score better on heteroscedastic data

4. **Computational budget.** The concomitant formulation takes roughly 1.5-2x more computation time, which is worth weighing on large datasets

5. **Nothing to tune for σ.** The concomitant formulation estimates the noise level itself; there is no penalty parameter for σ to set
