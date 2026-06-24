## Section 25: Delta-Ratio Linear Model Refits (Continuous Response)
**Key numbers:**
- Section 24 linear refits: TET impediment R2 = 0.259 vs DNMT3A recruitment R2 = 0.044 (~6x more variance explained) (delta_ratio_linear_model_comparison.tsv)
- Section 23 linear refits: baseline 5hmC R2 = 0.267, combined R2 = 0.296, K119ub R2 = 0.023 (delta_ratio_binary_vs_continuous.tsv)
- Top-2 features unchanged (baseline_hmc rank 1, k119ub rank 2) and signs flip exactly as the response inversion predicts: baseline_hmc -0.0084, k119ub +0.0069 (delta_ratio_feature_importance.tsv)
- OLS vs QR(0.5) agree in rank for top features; baseline_hmc beta stable and negative across tau 0.25/0.5/0.75 (delta_ratio_ols_vs_qr_comparison.tsv; delta_ratio_qr_coefficients.tsv)

**What this shows:** Refitting the Section 23/24 models with continuous delta_ratio (instead of binary DMR status) reproduces the same feature hierarchy and exact sign-flips, and TET impediment again far out-explains DNMT3A recruitment. Quantile regression confirms OLS is robust (no sign changes across quantiles), so dichotomization did not create the result.

**Figures:**
- 25a_dose_response_delta_ratio — baseline 5hmC vs delta_ratio.
- 25b_binary_vs_continuous_s23 — Section 23 AUC vs R2.
- 25c_feature_importance_comparison — logistic vs linear |beta|.
- 25d_binary_vs_continuous_s24 — Section 24 AUC vs R2 (with delta-R2).
- 25e_residual_diagnostics — OLS residual 4-panel.
- 25f_predicted_vs_observed — predicted vs observed delta_ratio.
- 25g_ols_vs_qr_coefficients — OLS vs QR(0.5) betas.
- 25h_quantile_process — QR betas across tau.
- 25i_qr_residual_diagnostics — OLS QQ vs QR residuals.
- 25j_qr_validation_composite — QR robustness composite.
