# Sections 22-26: Demethylation efficiency ratio (TET conversion), baseline-5hmC predictor, DNMT3A vs TET-impediment models, TET triple-KO comparison

## Summary
This group adjudicates the mechanism by which BAP1 loss restructures gene-body DNA methylation. It collapses the two marks into a single TET conversion-efficiency ratio, 5hmC/(5mC+5hmC), and asks: (i) does that ratio drop in BAP1-KO (Section 22); (ii) does baseline wildtype 5hmC predict which genes become hmC DMRs (Section 23); (iii) can mC hypermethylation be predicted better by DNMT3A-recruitment features or TET-impediment features (Section 24); (iv) do the same features predict the continuous delta-ratio, not just binary DMR status (Section 25); and (v) how does the BAP1-KO ratio shift compare quantitatively to a published TET triple-KO (Section 26). The answer across all five sections is consistent: [INTERPRETATION] BAP1 loss causes a partial, graded impairment of TET-mediated demethylation. Baseline 5hmC (the TET substrate) is the single strongest predictor of both DMR status (AUC 0.762) and hypermethylation (standardized beta = +1.25), the TET-impediment model significantly outperforms the DNMT3A-recruitment model (AUC 0.793 vs 0.696, DeLong p = 9.43e-49), H2AK119ub is a *negative* predictor of hypermethylation (argues against direct DNMT3A-UDR recruitment), and the BAP1-KO ratio shift reproduces only ~3-10% of the TET triple-KO effect with a graded (not binary) per-gene response.

## Section 22: section_22_demethylation_ratio
### Analysis question
Per gene, compute the TET conversion-efficiency ratio 5hmC/(5mC+5hmC) in control vs BAP1-KO and test whether it decreases in the mutant (impaired active demethylation).
### Key results
- Valid ratio data = 20,915 genes (source: demethylation_ratio_all_genes.tsv, row count).
- WT median ratio = 0.1275; KO median ratio = 0.1193; median delta_ratio = -0.0087 (source: demethylation_ratio_all_genes.tsv, computed from ratio_ctrl/ratio_mut/delta_ratio columns).
- Genes with decreased ratio (impaired TET) = 71.36%; increased = 28.45% (source: demethylation_ratio_all_genes.tsv, sign of delta_ratio column).
- Mean WT ratio = 0.1547, mean KO ratio = 0.1411, mean delta_ratio = -0.0136 (source: demethylation_ratio_all_genes.tsv).
- DMR-status breakdown of the 20,915 genes: Both mC+hmC DMR = 8,371; Not DMR = 7,027; hmC DMR only = 3,113; mC DMR only = 2,404 (source: demethylation_ratio_all_genes.tsv, dmr_status column).
- Per-chromatin-state median delta_ratio (Wilcoxon vs 0, BH-corrected): Active_Promoter = -0.0281 (n=6,356, ***); Active_Enhancer = -0.0248 (n=49, ***); Poised_Enhancer = -0.0073 (n=92, ***); Unmarked = -0.0051 (n=10,971, ***); Bivalent_Promoter = -0.0046 (n=150, ***); Polycomb = +0.0022 (n=56, ns); Repressed_Promoter = +0.0051 (n=3,241, ***) (source: demethylation_ratio_by_chromatin_state.tsv).
- Per-sample median ratios: ctrl-M = 0.1292, ctrl-F = 0.1278, mut-M = 0.1202, mut-F = 0.1168 (source: demethylation_ratio_per_sample.tsv, columns ratio_ctrl_M/F, ratio_mut_M/F).
- Top 5 most-decreased genes: Fndc10 (delta -0.543, WT 1.000 -> KO 0.457, Active_Promoter); H2ac4 (-0.500, 1.000 -> 0.500, Unmarked); Fam89b (-0.400, 0.885 -> 0.485, Active_Promoter); H1f0 (-0.363, 0.864 -> 0.502, Active_Promoter); Naa38 (-0.353, 0.740 -> 0.387, Active_Promoter) (source: demethylation_ratio_top_genes.tsv).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The dominant decrease in 5hmC/(5mC+5hmC) (71% of genes, both group-level and per-sample medians falling) indicates that BAP1 loss reduces the fraction of modified cytosines that have been oxidized by TET, i.e., a genome-wide attenuation of active demethylation. The effect is strongest at active chromatin (Active_Promoter, Active_Enhancer) and is reversed at Repressed_Promoter — consistent with TET impairment biting hardest where TET is normally most active, rather than where Polycomb already silences the locus. The partial magnitude (WT 0.128 -> KO 0.118, not -> 0) signals an indirect, not ablative, TET defect.
### Plot inventory
- 22a_wt_ratio_density — WT demethylation-ratio density distribution (all genes).
- 22b_wt_vs_ko_ratio — paired violin+box of control vs mutant ratio.
- 22c_delta_ratio_histogram — distribution of delta_ratio, colored by sign.
- 22d_delta_ratio_vs_diff — 2-panel scatter of delta_ratio vs mc_diff and vs hmc_diff.
- 22e_delta_ratio_by_dmr_status — delta_ratio across 4 DMR-status groups.
- 22f_delta_ratio_by_chromatin_state — delta_ratio across 7 chromatin states.
- 22g_top30_decreased_ratio — lollipop of 30 most negative delta_ratio genes.
- 22h_per_sample_ratio — per-sample (4 replicate) ratio distributions.

## Section 23: section_23_baseline_hmc_predictor
### Analysis question
Does higher wildtype 5hmC (substrate availability) predict which genes become hmC DMRs? Compare baseline 5hmC vs WT K119ub signal vs combined as logistic predictors of DMR status.
### Key results
- Model A (baseline 5hmC): AUC = 0.7617 [0.7545, 0.7689], McFadden R2 = 0.1406, AIC = 21,378.6 (source: baseline_hmc_model_comparison.tsv).
- Model B (K119ub signal): AUC = 0.5728 [0.5642, 0.5815], McFadden R2 = 0.0163, AIC = 24,471.6 (source: baseline_hmc_model_comparison.tsv).
- Model C (combined): AUC = 0.7900 [0.7833, 0.7968], McFadden R2 = 0.1693, AIC = 20,667.5 (source: baseline_hmc_model_comparison.tsv).
- N = 18,402 genes, 10,909 DMR+ (59.3%) used for all three models (source: baseline_hmc_model_comparison.tsv, n_genes/n_dmr).
- Model A standardized coefficient for wt_hmc: beta = +1.070 (raw OR = 1.107e6, p ~ 0); in combined Model C wt_hmc beta = +1.168 while K119ub beta flips to -0.597 (source: baseline_hmc_model_coefficients.tsv).
- Model B K119ub OR = 0.793 [0.773, 0.813], standardized beta = -0.335 — K119ub signal alone is a *negative*, weak predictor of hmC DMR status (source: baseline_hmc_model_coefficients.tsv).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Baseline 5hmC is a strong predictor of hmC-DMR susceptibility (AUC 0.762) and vastly outperforms K119ub signal (AUC 0.573), directly supporting the substrate-availability model: genes that carry more 5hmC in wildtype have more oxidized substrate for an impaired TET pathway to act on, so they show the largest hmC changes in BAP1-KO. Adding K119ub on top of 5hmC barely improves discrimination (0.762 -> 0.790) and enters with a negative sign, foreshadowing Section 24's conclusion that K119ub does not act as a positive DNMT3A-recruitment signal.
### Plot inventory
- 23_baseline_hmc_predictor — section output directory (created by script).
- 23a_roc_curves — ROC overlay of Models A/B/C.
- 23b_predicted_probability — Model A logistic probability curve vs baseline 5hmC.
- 23c_model_comparison — AIC bars + AUC dot/CI for the three models.
- 23d_dose_response_scatter — baseline 5hmC vs hmC mod_difference (dose-response).
- 23e_dose_response_by_chromatin — dose-response faceted by chromatin state.

## Section 24: section_24_dnmt3a_prediction
### Analysis question
Can mC hypermethylation (mC-up DMR) be predicted better by DNMT3A-recruitment features (K119ub + ATAC + CpG density) or by TET-impediment features (baseline 5hmC + ATAC)? Adjudicate the dual-mechanism hypothesis via logistic regression, random forest, DeLong test, interaction and exclusive-model comparisons.
### Key results
- Feature matrix = 11,937 genes, 6,192 hyper-DMR (51.9%) (source: dnmt3a_model_comparison.tsv, n_genes/n_hyper_dmr).
- Model AUCs: TET impediment = 0.7931 [0.7850, 0.8011]; DNMT3A recruitment = 0.6964 [0.6871, 0.7058]; Full = 0.8608; K119ub only = 0.6581; Stepwise = 0.8608 (source: dnmt3a_model_comparison.tsv).
- McFadden R2: TET = 0.1894 vs DNMT3A = 0.0799 vs K119ub-only = 0.0565 vs Full = 0.3130 (source: dnmt3a_model_comparison.tsv).
- DeLong test, DNMT3A recruitment vs TET impediment (shared models): p = 9.43e-49 (TET wins); exclusive models (DNMT3A = k119ub+cpg vs TET = baseline_hmc only): DeLong p = 6.13e-75 (source: dnmt3a_exclusive_model_comparison.tsv, delong_p column).
- Standardized full-model betas: baseline_hmc = +1.251 (rank 1, strongest); k119ub = -1.049 (rank 2, NEGATIVE); log_expression = +0.749; log_gene_length = -0.242; baseline_mc = +0.107; atac_count = -0.105; cpg_density = +0.082 (source: dnmt3a_feature_importance.tsv, lr_standardized_beta).
- Random forest importance (MeanDecreaseAccuracy) agrees: baseline_hmc = 158.2 (rank 1), k119ub = 134.9 (rank 2); LR-vs-RF rank correlation high (source: dnmt3a_feature_importance.tsv, rf_mean_decrease_accuracy/rf_rank).
- 10-fold CV AUCs near-identical to in-sample (negligible optimism): TET = 0.7928 +/- 0.0247 (optimism 0.00026); DNMT3A = 0.6961 +/- 0.0097; Full = 0.8603 +/- 0.0191; K119ub = 0.6581 +/- 0.0163 (source: dnmt3a_cv_results.tsv).
- Exclusive-model AUCs: TET exclusive (baseline_hmc only) = 0.7904 [0.7824, 0.7985] >> DNMT3A exclusive (k119ub+cpg) = 0.6580 [0.6483, 0.6678] (source: dnmt3a_exclusive_model_comparison.tsv).
- K119ub x baseline_hmc interaction: raw OR = 0.172 [0.071, 0.414], standardized beta = -0.151, p = 8.65e-05; interaction-model AUC = 0.8607 vs Full 0.8608, LR test vs Full p = 8.89e-05 (source: dnmt3a_interaction_results.tsv).
- K119ub is negative within every 5hmC tertile (univariate OR): Low = 0.458, Mid = 0.237, High = 0.436 (all p < 2.2e-16) — direction does NOT reverse at high 5hmC (source: dnmt3a_interaction_results.tsv).
- Non-promoter stratum (n=6,597, 2,283 hyper-DMR): TET = 0.7897 still beats DNMT3A = 0.7170; Full = 0.8712 (source: dnmt3a_stratified_comparison.tsv).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The mechanism that explains BAP1-KO gene-body hypermethylation is TET impediment, not DNMT3A-UDR recruitment. Baseline 5hmC (the TET substrate) is the single best predictor of which genes hypermethylate, and the TET model significantly out-discriminates the DNMT3A model (DeLong p = 9.43e-49). Critically, H2AK119ub is a *negative* predictor (beta = -1.05): genes with more K119ub are LESS likely to hypermethylate, the opposite of what direct DNMT3A recruitment by the K119ub-binding UDR (Chen et al. 2024) would predict. The negative direction is robust — it persists in the non-promoter subset (so it is not Active_Promoter confounding) and across all three 5hmC tertiles, and the K119ub x 5hmC interaction is significantly *negative* (the dual-mechanism prediction was a positive interaction), so there is no evidence that K119ub flips to a recruiting signal at high-substrate genes. Stepwise selection recovering the full model and near-zero CV optimism confirm the result is not overfit.
### Plot inventory
- 24_dnmt3a_prediction — section output directory.
- 24a_feature_correlation — Spearman correlation heatmap of the 7 predictors.
- 24b_roc_curves — ROC overlay of the 5 logistic models.
- 24c_feature_importance — LR standardized betas + RF MeanDecreaseAccuracy (2-panel).
- 24d_model_comparison — AIC bars + AUC dot/CI for 5 models with DeLong annotation.
- 24e_k119ub_predicted_probability — logistic curve P(hyper-DMR) vs K119ub (negative slope).
- 24f_chromatin_state_predictions — full-model predicted probability split by chromatin state.
- 24g_cv_auc_comparison — 10-fold CV AUC boxplots vs in-sample.
- 24h_stratified_comparison — all-genes vs non-promoter AUC and betas.
- 24i_interaction_k119ub_hmc — tertile-stratified sigmoid + forest of per-tertile K119ub OR.
- 24j_exclusive_model_comparison — shared vs exclusive ROC and AUC/CI.
- 24k_dnmt3a_vs_tet_composite — composite (24b + 24e + 24j-left) summary figure.

## Section 25: section_25_delta_ratio_models
### Analysis question
Refit the Section 23 and 24 models with the continuous delta_ratio as response (linear/OLS plus quantile regression) instead of binary DMR status — do the same features predict *how much* TET activity changes, with consistent sign, and is OLS robust?
### Key results
- Section 23 linear refits (R2 | logistic AUC for the same model): A baseline 5hmC R2 = 0.2666 (AUC 0.762); B K119ub R2 = 0.0231 (AUC 0.573); C combined R2 = 0.2956 (AUC 0.790) (source: delta_ratio_linear_model_comparison.tsv + delta_ratio_binary_vs_continuous.tsv).
- Section 24 linear refits: Full R2 = 0.3902 (AUC 0.861); TET impediment R2 = 0.2589 (AUC 0.793); DNMT3A recruitment R2 = 0.0440 (AUC 0.696); K119ub only R2 = 0.0311 (AUC 0.658) (source: delta_ratio_linear_model_comparison.tsv + delta_ratio_binary_vs_continuous.tsv).
- TET impediment explains ~6x more delta_ratio variance than DNMT3A recruitment in the linear frame (R2 0.259 vs 0.044), echoing the logistic AUC gap (source: delta_ratio_linear_model_comparison.tsv).
- Sign consistency (full linear model, standardized betas, delta_ratio is more-negative = more impairment): baseline_hmc = -0.0084 (logistic was +1.25 -> flipped, CONSISTENT); k119ub = +0.0069 (logistic -1.05 -> flipped, CONSISTENT); both top features flip sign exactly as the response inversion predicts (source: delta_ratio_feature_importance.tsv).
- Linear feature-importance ranking (|standardized beta|): baseline_hmc (rank 1), k119ub (rank 2), baseline_mc (3), log_expression (4), log_gene_length (5), atac_count (6), cpg_density (7) — top-2 match the logistic ranks 1 and 2 (source: delta_ratio_feature_importance.tsv, linear_rank vs lr_rank).
- OLS vs quantile-regression QR(0.5) coefficients agree in rank for the top features (baseline_hmc, k119ub rank 1 and 2 under both); QR(0.5) baseline_hmc beta = -0.0095 vs OLS -0.0084 (source: delta_ratio_ols_vs_qr_comparison.tsv).
- Quantile process (tau 0.25/0.5/0.75): baseline_hmc beta stable and negative across all quantiles (-0.0095, -0.0095, -0.0093); k119ub beta stays positive and grows with tau (+0.0062, +0.0086, +0.0110) (source: delta_ratio_qr_coefficients.tsv).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Recasting the analysis with a continuous response confirms the binary-DMR conclusions are not an artifact of dichotomization: the same two features (baseline 5hmC, then K119ub) dominate, their signs flip exactly as expected when the response is inverted (high-5hmC genes lose the most demethylation efficiency; high-K119ub genes lose less), and the TET model again far out-explains the DNMT3A model. Quantile regression shows the baseline-5hmC effect is homogeneous across the delta-ratio distribution while the K119ub effect strengthens in the upper tail, but neither changes sign — so OLS is robust and the mechanistic ranking is stable. The low absolute R2 values are expected (R2 and AUC are not on the same scale) and do not undercut the consistent feature hierarchy.
### Plot inventory
- 25a_dose_response_delta_ratio — baseline 5hmC vs delta_ratio scatter.
- 25b_binary_vs_continuous_s23 — Section 23 logistic AUC vs linear R2 bars.
- 25c_feature_importance_comparison — logistic |beta| vs linear |beta| scatter.
- 25d_binary_vs_continuous_s24 — Section 24 logistic AUC vs linear R2 (with delta-R2).
- 25e_residual_diagnostics — 4-panel OLS residual diagnostics.
- 25f_predicted_vs_observed — predicted vs observed delta_ratio, colored by chromatin state.
- 25g_ols_vs_qr_coefficients — OLS vs QR(0.5) standardized-beta scatter.
- 25h_quantile_process — QR betas across tau = 0.25/0.5/0.75 (faceted).
- 25i_qr_residual_diagnostics — OLS QQ vs QR(0.5) residual histogram.
- 25j_qr_validation_composite — composite QR robustness (25g + 25h).

## Section 26: section_26_tet_ko_comparison
### Analysis question
Quantitatively compare the BAP1-KO demethylation-ratio shift to a published TET triple-KO (GSE166423, Lopez-Moyado et al.): is BAP1-KO a partial/attenuated, graded version of the TET-KO signature, and are the same genes affected?
### Key results
- Matched genes = 20,851 (source: tet_ko_comparison_summary.tsv, "N matched genes").
- Median delta_ratio: BAP1-KO = -0.0085 vs TET-KO = -0.2597; absolute attenuation (BAP1/TET) = 0.0326 i.e. BAP1-KO reproduces ~3.3% of the TET-KO absolute effect (source: tet_ko_comparison_summary.tsv).
- Relative (baseline-normalized) attenuation = 0.0872 (~8.7% of TET-KO effect); BAP1-KO median relative delta = -0.0872 vs TET-KO = -1.000 (source: tet_ko_comparison_summary.tsv).
- Effect sizes: BAP1-KO Cliff's delta = -0.430, Cohen's d = -0.712; TET-KO Cliff's delta = -0.680, Cohen's d = -1.824 (source: tet_ko_comparison_summary.tsv).
- KS statistic = 0.6741 (p < 2.2e-16): the two delta distributions differ strongly (source: tet_ko_comparison_summary.tsv).
- QQ slope (absolute) = 0.0933; relative QQ slope = 0.3479 — BAP1 effect is a fraction of the TET-KO shift at every quantile (source: tet_ko_comparison_summary.tsv).
- Per-gene Spearman rho (absolute) = 0.2199 (p < 2.2e-16) BUT residualizing on WT baselines collapses it to 0.0975, and baseline accounts for 99.9% of explained variance (R2 baseline 0.8328 vs full 0.8337, gene-specific delta-R2 = 0.0008); relative rho = -0.0120 (p = 0.084, n.s.) (source: tet_ko_comparison_summary.tsv).
- Response decomposition: TET-KO is binary — 68.5% complete loss (relative delta = -1) + 27.0% no WT signal + 4.5% partial; BAP1-KO is graded — 0.8% strong (rel < -0.5), 45.1% moderate (-0.5 to -0.1), 54.1% weak/no change (source: tet_ko_comparison_summary.tsv).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] BAP1-KO is a partial, graded attenuation of the TET triple-KO demethylation signature, not a phenocopy. It reproduces only ~3% (absolute) to ~9% (baseline-normalized) of the TET-KO ratio shift, and where TET-KO drives an all-or-nothing collapse of 5hmC (68% of genes at complete loss), BAP1-KO produces a continuous, mostly moderate dimming of demethylation efficiency. The apparent positive per-gene correlation (rho 0.22) is almost entirely an artifact of shared WT baseline structure (residualized rho 0.10; 99.9% of explained variance is baseline-driven; relative rho ~ 0), so the two perturbations do not converge on the same individual genes through a shared mechanism. This pattern is exactly what an INDIRECT TET defect predicts: BAP1 loss (via PRC1/H2AK119ub) impairs TET *access/efficiency* rather than eliminating TET enzyme — a dimmer switch, not an off switch.
### Plot inventory
- 26a_wt_ratio_comparison — overlaid WT ratio densities, BAP1 WT vs TET-KO WT.
- 26b_delta_ratio_density — delta_ratio densities, BAP1-KO vs TET-KO (KS + attenuation).
- 26c_qq_plot — QQ of BAP1 vs TET-KO delta quantiles (absolute, slope annotation).
- 26d_effect_size_comparison — Cliff's delta / %decreased / Cohen's d / relative attenuation bars.
- 26e_per_gene_scatter — per-gene BAP1 delta vs TET-KO delta scatter.
- 26f_chromatin_stratified — BAP1 vs TET-KO delta boxplots per chromatin state.
- 26g_relative_delta_density — baseline-normalized relative-delta densities.
- 26h_relative_qq_plot — relative-delta QQ (step structure of TET-KO binary loss).
- 26i_response_decomposition — stacked bars: TET-KO binary vs BAP1-KO graded response.

## Cross-section synthesis
The five sections form a single mechanistic argument that converges on TET impediment as the driver of BAP1-KO methylation remodeling. Section 22 establishes the phenotype: a genome-wide, partial drop in the 5hmC/(5mC+5hmC) conversion ratio (71% of genes down, strongest at active chromatin). Section 23 shows that baseline 5hmC — the TET substrate — predicts hmC-DMR susceptibility far better than K119ub (AUC 0.762 vs 0.573), supporting substrate availability. Section 24 directly pits the two competing models against each other and finds TET impediment wins decisively (AUC 0.793 vs 0.696, DeLong p = 9.43e-49), with baseline 5hmC the #1 predictor and K119ub entering *negatively* — the opposite of the DNMT3A-UDR recruitment prediction — robustly across strata, tertiles, and an exclusive-feature comparison. Section 25 confirms the same feature hierarchy and exact sign-flips when delta_ratio is modeled continuously, and shows OLS is robust under quantile regression. Section 26 then bounds the magnitude: BAP1-KO recovers only ~3-10% of a TET triple-KO's effect and does so in a graded rather than binary fashion, with the per-gene "convergence" being a baseline artifact. Together these support the paper's thesis that [INTERPRETATION] BAP1 loss restructures DNA methylation by indirectly blocking TET-mediated active demethylation (a dimmer on the TET pathway), rather than by directly recruiting DNMT3A through H2AK119ub.

## Tables used
- demethylation_ratio_all_genes.tsv — per-gene mC/hmC means, ratio_ctrl/ratio_mut/delta_ratio, DMR status, chromatin state (20,915 genes); source of Section 22 medians, %-decreased, and DMR-status counts.
- demethylation_ratio_by_chromatin_state.tsv — per-state median delta_ratio, IQR, Wilcoxon-vs-0 p/p_adj, significance label (7 states).
- demethylation_ratio_per_sample.tsv — per-replicate ratios (ctrl-M/F, mut-M/F) used for per-sample medians.
- demethylation_ratio_top_genes.tsv — top 50 most-decreased delta_ratio genes (Section 22 lollipop / top-5 list).
- baseline_hmc_model_comparison.tsv — Section 23 Models A/B/C AIC, AUC+CI, McFadden R2, N/DMR counts.
- baseline_hmc_model_coefficients.tsv — Section 23 raw OR and standardized betas for wt_hmc, gb_ctrl_signal.
- baseline_hmc_predictor_all_genes.tsv — per-gene wt_hmc, K119ub signal, DMR q/mod_difference, predicted probability (dose-response source).
- dnmt3a_model_comparison.tsv — Section 24 five models: AIC, AUC+CI, McFadden R2, n_genes/n_hyper_dmr.
- dnmt3a_feature_importance.tsv — LR standardized betas + ranks and RF MeanDecreaseAccuracy/Gini + ranks for 7 features.
- dnmt3a_exclusive_model_comparison.tsv — shared vs exclusive DNMT3A/TET AUC+CI, AIC, R2, DeLong p.
- dnmt3a_cv_results.tsv — 10-fold CV mean/sd/min/max AUC and optimism per model.
- dnmt3a_interaction_results.tsv — K119ub x baseline_hmc interaction OR/beta/p, interaction-model AUC/LR-test, per-tertile K119ub OR.
- dnmt3a_stratified_comparison.tsv — all-genes vs non-promoter AUC+CI, AIC, R2 for all 5 models.
- dnmt3a_model_coefficients.tsv — full raw OR + standardized betas for every model term (k119ub, baseline_hmc, etc.).
- delta_ratio_linear_model_comparison.tsv — Section 23/24 linear refit R2, adj R2, F-stat, N, residual SE.
- delta_ratio_binary_vs_continuous.tsv — side-by-side logistic AUC/McFadden vs linear R2/adj R2.
- delta_ratio_feature_importance.tsv — logistic vs linear standardized betas and ranks (sign-consistency check).
- delta_ratio_ols_vs_qr_comparison.tsv — OLS vs QR(0.5) standardized betas and ranks.
- delta_ratio_qr_coefficients.tsv — QR betas+CIs at tau 0.25/0.5/0.75 per feature, with OLS reference.
- tet_ko_comparison_summary.tsv — 37-row Section 26 summary (attenuation, effect sizes, KS, rho, QQ slopes, decomposition, variance analysis).
- tet_ko_comparison_matched_genes.tsv — per-gene matched BAP1-KO and TET-KO ratios/deltas (20,851 genes).

## Data-quality flags
- **DeLong p-value precision (Section 24).** The TSV stores the exact DeLong p for DNMT3A-vs-TET as 9.43e-49 (shared) and 6.13e-75 (exclusive) in dnmt3a_exclusive_model_comparison.tsv. TODO.md row 44 and the script annotations report it only as the formatted display value "p < 2.2e-16" (the script's fmt_p() floor). Cite the exact TSV values; "p < 2.2e-16" is a display floor, not the true magnitude.
- **Section 24 hyper-DMR rate discrepancy.** The TSV-confirmed rate is 51.9% (6,192/11,937; dnmt3a_model_comparison.tsv). TODO.md rows 11a/11b state "46.2%" and "11,936 genes" — not confirmed in any table. Use 51.9% / 11,937 from the TSV.
- **Section 23 K119ub AUC discrepancy.** baseline_hmc_model_comparison.tsv gives Model B (K119ub) AUC = 0.5728. TODO.md row 23 states "K119ub AUC=0.592" — [UNVERIFIED: 0.592 per TODO.md, not confirmed in tables]. Use 0.573 from the TSV.
- **Section 23 dose-response rho not in any TSV.** The Spearman rho(wt_hmc vs hmC mod_difference) is computed in-script and printed to stdout only; it is NOT exported to a table. TODO.md row 23 cites "rho=-0.586" — [UNVERIFIED: -0.586 per TODO.md, not confirmed in tables]. A data spot-check (baseline_hmc_predictor_all_genes.tsv) confirms the negative direction: 18/20 of the highest-wt_hmc genes have negative hmC mod_difference, consistent with a negative dose-response, but the exact coefficient cannot be confirmed from a TSV.
- **Section 22 Cliff's delta / paired-Wilcoxon p not in a TSV.** Section 22's Cliff's delta and paired Wilcoxon p are printed to stdout / embedded in figure annotations only (not exported). TODO.md row 14b cites "Cliff's delta=0.455 (medium)" and "72.5% decreased"; the all-genes TSV yields 71.36% decreased (delta_ratio < 0). The 72.5% figure is [UNVERIFIED: 72.5% per TODO.md] — use 71.36% (sign of delta_ratio in demethylation_ratio_all_genes.tsv). Section 26 separately computes a BAP1-KO Cliff's delta of -0.430 vs 0 (tet_ko_comparison_summary.tsv) over the 20,851-gene matched set — a different baseline (vs 0, not paired ctrl-vs-mut), so it is not directly comparable to the row-14b value.
- **Row-count drift across sections (expected, not an error).** N differs by section due to different join cascades and complete-case filters: Section 22 = 20,915; Section 23 model set = 18,402; Section 24 = 11,937; Section 25 inherits these; Section 26 matched = 20,851. All are internally consistent with each section's documented merge logic.
- **Median delta_ratio rounding.** demethylation_ratio_all_genes.tsv gives median delta_ratio = -0.00867; tet_ko_comparison_summary.tsv reports BAP1-KO median delta = -0.008477 over its (slightly smaller, TET-matched) 20,851-gene set. The small difference is the different gene set, not a discrepancy.
- **Section 26 external dependency.** The TET-KO comparison depends on data/tet_ko_gene_signal.tsv (produced upstream by preprocess_tet_ko_bigwig.R on HPC from GSE166423). The output tables are present and populated (20,851 matched genes), so the dependency was satisfied for this run. TET-KO data is BS/OxBS-seq from a different tissue/protocol (Lopez-Moyado et al.); cross-platform comparison is qualitative for magnitude, which the script appropriately frames via baseline-normalized relative metrics.
