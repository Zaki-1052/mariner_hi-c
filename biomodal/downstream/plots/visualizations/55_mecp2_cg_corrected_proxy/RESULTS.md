## Section 55: MeCP2 as CG-Corrected Proxy for Non-CG Methylation in TADs
**Key numbers:**
- Model A (MeCP2 ~ CG 5mC): R²=0.0211, slope=149.65 (p=3.6e-58), n=12,141 TADs (55_mecp2_regression_model_summary.tsv / _coefficients.tsv)
- Model B (+CG 5hmC): R²=0.0359; 5hmC term=359.88 (p=4.2e-42) steeper than 5mC term=74.70 (p=3.1e-12) (55_mecp2_regression_coefficients.tsv)
- Variance ratio (ctrl): MeCP2 raw=0.0680, residual Model A=0.0665, Model B=0.0656 — all far below CG 5mC=1.018 (55_mecp2_variance_ratio_summary.tsv)
- CG-stratified residual (Model A): all quartiles median ≤0, reject zero — Q2 median=-0.268 (p=6.6e-232), Q3=-0.245 (p=1.8e-183) (55_mecp2_cg_stratified_summary.tsv)
- DMR-stratified: Hyper-DMR TADs median resid=-0.080 (n=6,939); No-DMR=-0.269 (n=5,202) (55_mecp2_dmr_stratified_summary.tsv)

**What this shows:** CG methylation explains only ~2% of TAD-level MeCP2 enrichment (~3.6% with 5hmC). The CG-corrected residual is NOT TAD-organized (variance ratio ~0.066, like raw MeCP2 and ~15x below CG 5mC) and is median-negative, so MeCP2 sits below the CG prediction in most TADs — not a positive non-CG excess. The strong 5hmC regression coefficient again points to the TET axis.

**Figures:**
- `55a_mecp2_vs_cg_scatter/` — MeCP2 vs CG 5mC per TAD (Model A)
- `55b_model_diagnostics/` — Model A residuals/Q-Q
- `55c_residual_distributions/` — ctrl vs mut residuals
- `55d_variance_ratio_comparison/` — CG vs MeCP2-raw vs residuals (key panel)
- `55e_boundary_correlation/` — residual vs TAD boundary score
- `55f_cg_stratified_residuals/` — residual by CG-log2FC quartile (key test)
- `55g_dmr_stratified_residuals/` — residual by CG DMR class
- `55h_model_comparison/` — Model A vs B variance ratios
- `55i_composite/` — scatter + variance ratio + CG-stratified composite
