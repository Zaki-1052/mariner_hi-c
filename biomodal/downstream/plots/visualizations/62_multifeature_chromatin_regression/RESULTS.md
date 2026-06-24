## Section 62: Multi-Feature Chromatin Regression at MeCP2 Peaks
**Key numbers:**
- Binding-level R2: CG only = 0.017, Chromatin only = 0.246, Full = 0.260 (n = 202,574 peaks) (source: 62_model_comparison_summary.tsv)
- Differential R2: CG only = 0.013, Chromatin only = 0.170, Full = 0.173 (source: 62_model_comparison_summary.tsv)
- Variance partition (binding): Chromatin-unique = 24.3% vs CG-unique = 1.5% (source: 62_variance_partition.tsv)
- Top standardized betas: K119ub = +0.199, ATAC = +0.114, CG 5mC = +0.089, K27me3 = +0.061 (source: 62_standardized_coefficients.tsv)
- CG-unexplained residual correlates most with K119ub change (rho = 0.209) and K27me3 change (rho = 0.127) (source: 62_residual_chromatin_regression.tsv)

**What this shows:** Chromatin marks explain MeCP2 binding ~15x better than CG methylation alone (24.3% vs 1.5% unique variance), with H2AK119ub the dominant predictor in both standardized regression and LASSO and the strongest correlate of the CG-unexplained MeCP2 residual. The chromatin R2 advantage is largest at enhancers (Active_Enhancer gain +0.376), supporting H2AK119ub rather than CG methylation as the upstream signal positioning MeCP2 in BAP1-KO.

**Figures:**
- 62a_model_r2_comparison, 62b_variance_partition — model R2 + variance partition
- 62c_lasso_path, 62d_standardized_coefficients — feature selection + effect sizes
- 62e_chromatin_state_r2, 62f_residual_vs_k27me3, 62g_composite
