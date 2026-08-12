## Section 23: Baseline 5hmC as Predictor of DMR Susceptibility
**Key numbers:**
- Model A (baseline 5hmC) AUC = 0.762 [0.755, 0.769], McFadden R2 = 0.141 (baseline_hmc_model_comparison.tsv)
- Model B (K119ub signal) AUC = 0.573 [0.564, 0.581] — far weaker (baseline_hmc_model_comparison.tsv)
- Model C (combined) AUC = 0.790; adding K119ub barely improves on 5hmC alone (baseline_hmc_model_comparison.tsv)
- Standardized betas: wt_hmc = +1.070 (Model A); in Model C wt_hmc = +1.168, K119ub flips to -0.597 (baseline_hmc_model_coefficients.tsv)
- N = 18,402 genes, 10,909 DMR+ (59.3%) (baseline_hmc_model_comparison.tsv)

**What this shows:** Baseline wildtype 5hmC strongly predicts which genes become hmC DMRs (AUC 0.762) and vastly outperforms WT K119ub signal (AUC 0.573), supporting the substrate-availability model — genes with more 5hmC have more oxidized substrate for an impaired TET pathway to affect. K119ub enters the combined model with a negative sign and adds little discrimination.

**Figures:**
- 23a_roc_curves — ROC overlay of Models A/B/C.
- 23b_predicted_probability — Model A logistic probability vs baseline 5hmC.
- 23c_model_comparison — AIC bars + AUC dot/CI.
- 23d_dose_response_scatter — baseline 5hmC vs hmC mod_difference.
- 23e_dose_response_by_chromatin — dose-response faceted by chromatin state.
