# Advanced: Using Concomitant Formulation

The concomitant formulation provides a more robust approach to both regression and classification by jointly estimating the model coefficients and the noise level (σ). This is particularly useful when variance is not homogeneous across your samples.

## What is Concomitant Formulation?

By setting `--p-concomitant True`, you switch from standard formulation to the **concomitant formulation (R3)**:

### Mathematical Formulation

**Standard Formulation (R1):**
```
min ||y - X·β||² / σ + λ·||β||₁    s.t. C·β = 0
```
Assumes a fixed noise level across all samples.

**Concomitant Formulation (R3):**
```
min ||y - X·β||² / σ + e·σ + λ·||β||₁    s.t. C·β = 0
```
Where `e = 20.0` is the penalty parameter for the standard deviation (σ).

### Key Differences

| Aspect | Standard (R1) | Concomitant (R3) |
|--------|--------------|-----------------|
| Noise level | Fixed | Jointly estimated |
| Heteroscedasticity | Assumed equal variance | Adaptive to varying noise |
| Uncertainty estimates | May be inaccurate | More reliable |
| Computation | Faster | Slightly slower |
| Output label | "Formulation: R1" | "Formulation: R3 (concomitant with e = 20.0)" |

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

### TRAC with Concomitant

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

Replace the training step with:

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

### TRAC with Concomitant

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
Formulation: R3 (concomitant with e = 20.0)
```

This indicates:
- The model is using concomitant formulation
- The noise level parameter `e = 20.0` controls the penalty strength
- Coefficients and predictions are more robust to heteroscedasticity

## Tips for Optimal Results

1. **Compare Both Formulations**: Run your analysis with both `--p-concomitant False` and `True` to see which produces better results for your data

2. **Examine Residuals**: Plot residuals vs. fitted values before deciding:
   - If variance is constant → Standard formulation is sufficient
   - If variance increases/decreases with fitted values → Use concomitant

3. **Cross-Validation Performance**: Compare CV scores between formulations; concomitant often performs better on heteroscedastic data

4. **Computational Budget**: Concomitant requires ~1.5-2x more computation time, so consider this for large datasets

5. **Parameter Tuning**: The default `e = 20.0` works well for most cases, but can be adjusted if needed
