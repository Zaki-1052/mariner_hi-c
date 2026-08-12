## Section 24: DNMT3A Recruitment vs TET Impediment (Hypermethylation Prediction)
**Key numbers:**
- TET impediment AUC = 0.793 [0.785, 0.801] >> DNMT3A recruitment AUC = 0.696 [0.687, 0.706]; DeLong p = 9.43e-49 (dnmt3a_model_comparison.tsv; dnmt3a_exclusive_model_comparison.tsv)
- Baseline 5hmC is the #1 predictor (standardized beta = +1.25); K119ub is NEGATIVE (beta = -1.05, rank 2) (dnmt3a_feature_importance.tsv)
- Exclusive models: TET (baseline_hmc only) AUC = 0.790 >> DNMT3A (k119ub+cpg) AUC = 0.658; DeLong p = 6.13e-75 (dnmt3a_exclusive_model_comparison.tsv)
- K119ub stays negative in every 5hmC tertile (OR 0.46/0.24/0.44); K119ubxhmC interaction beta = -0.151, p = 8.65e-05 (dnmt3a_interaction_results.tsv)
- 11,937 genes, 6,192 hyper-DMR (51.9%); 10-fold CV optimism ~0 (dnmt3a_model_comparison.tsv; dnmt3a_cv_results.tsv)

**What this shows:** mC hypermethylation in BAP1-KO is predicted by TET-impediment features (baseline 5hmC substrate), not DNMT3A-recruitment features. The TET model significantly out-discriminates the DNMT3A model, and H2AK119ub is a negative predictor — the opposite of the direct DNMT3A-UDR recruitment prediction (Chen et al. 2024). The negative K119ub direction is robust across non-promoter genes and all 5hmC tertiles, and the interaction is significantly negative, so there is no evidence for a dual mechanism.

**Figures:**
- 24a_feature_correlation — Spearman heatmap of 7 predictors.
- 24b_roc_curves — ROC of 5 logistic models.
- 24c_feature_importance — LR betas + RF importance.
- 24d_model_comparison — AIC + AUC with DeLong annotation.
- 24e_k119ub_predicted_probability — negative-slope P(hyper-DMR) vs K119ub.
- 24f_chromatin_state_predictions — predictions by chromatin state.
- 24g_cv_auc_comparison — 10-fold CV AUC vs in-sample.
- 24h_stratified_comparison — all genes vs non-promoter.
- 24i_interaction_k119ub_hmc — tertile sigmoids + OR forest.
- 24j_exclusive_model_comparison — shared vs exclusive ROC/AUC.
- 24k_dnmt3a_vs_tet_composite — composite summary figure.
