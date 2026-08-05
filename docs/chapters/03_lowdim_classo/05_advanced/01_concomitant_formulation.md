# Advanced: Using Concomitant Formulation

The concomitant formulation provides a more robust approach to both regression and classification by jointly estimating the model coefficients and the noise level (σ). This is particularly useful when variance is not homogeneous across your samples.

## What is Concomitant Formulation?

By setting `--p-concomitant True`, you switch from the standard formulation to the **concomitant formulation** (R3 for least squares, R4 when combined with the Huber loss):

### Mathematical Formulation

**Standard Formulation (R1):**
```
min ||y - X·β||² + λ·||β||₁    s.t. C·β = 0
```
The noise level σ is treated as fixed and is not estimated.

**Concomitant Formulation (R3, least squares):**
```
min ||y - X·β||² / (2σ) + (n/2)·σ + λ·||β||₁    s.t. C·β = 0
```
Here σ is estimated jointly with the coefficients β. The concomitant formulation has **no tuning parameter for σ**: the coefficient `n/2` (with `n` the number of samples) is a fixed constant that follows from the perspective function of the squared loss — it is not a hyperparameter and is not exposed to the user.

**Concomitant Formulation with Huber loss (R4):**
```
min  h_{ρ,σ}(y - X·β) + n·σ + λ·||β||₁    s.t. C·β = 0
```
The robust (Huber) analogue, obtained by combining `--p-concomitant True` with `--p-huber True`. As in R3, σ is estimated with **no tuning parameter**; the σ-coefficient is `n` (rather than `n/2`), which again follows directly from the perspective function of the Huber loss.

### Key Differences

| Aspect | Standard (R1) | Concomitant (R3) |
|--------|--------------|-----------------|
| Noise level | Fixed | Jointly estimated |
| Heteroscedasticity | Assumed equal variance | Adaptive to varying noise |
| Uncertainty estimates | May be inaccurate | More reliable |
| Computation | Faster | Slightly slower |
| Output label | "Formulation: R1" | "Formulation: R3 (concomitant)" |

## When to Use Concomitant Formulation

Use `--p-concomitant True` when:
- Your **residual plots show increasing variance** with predicted values (heteroscedasticity)
- You have **samples from different environments** with potentially different noise levels
- You need **more reliable confidence intervals** for your predictions
- You want **robust feature selection** across diverse data
- Computation time is not a critical constraint

## Regression with Concomitant Formulation

### Log-Contrast with Concomitant

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

### trac with Concomitant

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

## Classification with Concomitant Formulation

### Log-Contrast with Huber Loss

For both [log-contrast classification](../03_classification/01_logcontrast.md) and [trac classification](../03_classification/02_trac.md), replace the training step with:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_lc.qza \
    --i-c data/ccovariates_lc.qza \
    --i-weights data/wcovariates_lc.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column vegetation \
    --p-huber False \
    --p-concomitant True \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 42 \
    --p-no-cv-one-se \
    --o-result data/classifytaxa_lc_concomitant.qza
```

Then proceed with prediction and visualization using `classifytaxa_lc_concomitant.qza`.

### trac with Concomitant

Replace the training step with:

```bash
qiime classo classify \
    --i-features data/classify-xtraining_trac.qza \
    --i-c data/ccovariates_trac.qza \
    --i-weights data/wcovariates_trac.qza \
    --m-y-file data/atacama-selected-covariates-veg.tsv \
    --m-y-column vegetation \
    --p-huber False \
    --p-concomitant True \
    --p-stabsel \
    --p-cv \
    --p-path \
    --p-lamfixed \
    --p-stabsel-threshold 0.5 \
    --p-cv-seed 42 \
    --p-no-cv-one-se \
    --o-result data/classifytaxa_trac_concomitant.qza
```

Then proceed with prediction and visualization using `classifytaxa_trac_concomitant.qza`.

## Interpreting Results with Concomitant Formulation

When you run your pipeline with `--p-concomitant True`, the output will display:

```
Formulation: R3 (concomitant)
```

(or `Formulation: R4 (concomitant + Huber)` when `--p-huber True` is also set). This indicates:
- The model is using the concomitant formulation, estimating the noise level σ jointly with the coefficients
- σ is determined by the model, not by a hyperparameter — there is nothing to tune
- Coefficients and predictions are more robust to heteroscedasticity

## Tips for Optimal Results

1. **Compare Both Formulations**: Run your analysis with both `--p-concomitant False` and `True` to see which produces better results for your data

2. **Examine Residuals**: Plot residuals vs. fitted values before deciding:
   - If variance is constant → Standard formulation is sufficient
   - If variance increases/decreases with fitted values → Use concomitant

3. **Cross-Validation Performance**: Compare CV scores between formulations; concomitant often performs better on heteroscedastic data

4. **Computational Budget**: Concomitant requires ~1.5-2x more computation time, so consider this for large datasets

5. **No σ Tuning**: The concomitant formulation estimates the noise level σ automatically — there is no penalty parameter for σ to set or tune
