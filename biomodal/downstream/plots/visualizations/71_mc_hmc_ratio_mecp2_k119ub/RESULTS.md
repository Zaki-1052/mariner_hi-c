## Section 71: 5mC/5hmC Ratio vs MeCP2 at H2AK119ub Loci

**Key numbers:**
- n=12,149 genes with mc_diff, hmc_diff, MeCP2 fold, and K119ub data (71_statistics.tsv)
- log2(mC/hmC ratio) vs MeCP2: Spearman rho=0.014, p=0.113 — NOT significant (71_statistics.tsv)
- K119ub vs MeCP2: Spearman rho=0.187, p=3.53e-96; K119ub-only model R2=0.0738 (71_statistics.tsv; 71_model_comparison.tsv)
- Variance partition: K119ub unique=7.3%, ratio unique=0.0%, unexplained=92.6% (71_variance_partition.tsv)
- Standardized betas: K119ub +0.040 (p=1.6e-202) vs ratio +0.0020 (p=0.131, NS) (71_standardized_betas.tsv)

**What this shows:** Adjudicating between a 5mC-skew-gradient model and a K119ub-driven model, the per-gene 5mC/5hmC change ratio has effectively zero independent predictive power for MeCP2 binding (ratio-only R2=0.001, NS beta, flat quintile dose-response with Kruskal-Wallis p=5.2e-03 but no monotonic trend), whereas K119ub alone explains ~7.3%. MeCP2 does not titrate to how 5mC-skewed a gene is — it responds to the Polycomb mark. A small significant K119ub-by-ratio interaction (p=4.3e-14) indicates mild synergy, but K119ub is unambiguously the driver. The interaction model has the lowest AIC (-12,831.8).

**Figures:**
- 71a_k119ub_mecp2_ratio_scatter/ — K119ub vs MeCP2 colored by log2(mC/hmC ratio)
- 71b_mecp2_by_ratio_quintile/ — MeCP2 fold violins by ratio quintile (dose-response)
- 71c_ratio_vs_mecp2_direct/ — direct ratio vs MeCP2 with loess + OLS
- 71d_nested_model_comparison/ — R2 bars for 4 nested linear models
