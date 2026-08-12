## Section 11: MeCP2 (CUT&RUN) Correlation at DMRs
**Key numbers:**
- DMR-overlap Fisher OR = 5.13, p = 1.27e-33 (hyper DMRs favor MeCP2-Up peaks) (source: mecp2_dmr_overlap_summary.tsv)
- Hyper DMRs 7.16% overlap MeCP2-Up vs 1.89% MeCP2-Down (source: mecp2_dmr_overlap_summary.tsv)
- Binary MeCP2-gain coordinated vs non-coordinated: OR = 1.04, p = 0.795 (NS) (source: mecp2_binary_gain_coordinated_fisher.tsv)
- Stratified Spearman (mC% vs MeCP2 fold): Coordinated rho = +0.024 (p=0.065), mC-only rho = +0.076 (p=5.5e-04) (source: mecp2_stratified_correlations.tsv)
- Delta-ratio LM: delta_ratio coef = −0.367, p = 1.42e-06 (source: mecp2_delta_ratio_lm_coefficients.tsv); GLM OR = 7.2e-07, p = 1.97e-05 (source: mecp2_delta_ratio_glm_odds_ratios.tsv)

**What this shows:** MeCP2 differential binding (CUT&RUN, NOT ChIP) is enriched at hypermethylated DMRs by peak co-location, but gene-level directional coupling is weak (stratified rho ≈ 0.02-0.08; binary-gain Fisher NS). Delta-ratio regression recovers a modest significant signal: greater TET impairment predicts higher MeCP2 occupancy.

**Figures:**
- `11a_mecp2_overlap` — overlap % by DMR direction + Fisher
- `11b_mecp2_fold_by_dmr_direction` — MeCP2 fold by hyper/hypo/NS
- `11c_mc_vs_mecp2_scatter` — gene-level mC% vs MeCP2 fold
- `11d_mecp2_coordinated_genes` — coordinated-gene box + top-20 bar
- `11e_mecp2_integration_heatmap` — 2×2 mC × MeCP2 O/E (figure-only)
- `11f_mecp2_delta_ratio_lm` / `11g_mecp2_delta_ratio_glm` — delta-ratio regression
- `11h_mecp2_binary_gain_fisher` — binary gain enrichment
- `11i_mc_vs_mecp2_stratified` — scatter by coordination status
