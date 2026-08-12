## Section 47: CTCF Anchor Methylation Overlay
**Key numbers:**
- mC hyper at dynamic CpG regions (shores+shelves), LOST vs GAINED CTCF anchors: 22.3% vs 8.1%, OR = 3.28 [2.57–4.20], p = 5.4e-24; hmC hypo OR = 2.08, p = 3.9e-6 (47a_fisher_results.tsv)
- CpG islands (null, constitutively unmethylated): mC OR = 0.84, p = 1 (47b_multiregion_fisher.tsv); shores OR = 3.01, shelves OR = 3.68, promoters OR = 3.32 (all significant)
- Coordinated mC-up/hmC-down at lost vs gained anchors: 33.1% vs 18.5%, OR = 2.18, p = 1.3e-18 (47c_coordinated_enrichment.tsv)
- Distance gradient: <200 kb OR = 11.21 (p = 5.5e-23) decaying to >1 Mb OR = 0.87 (NS); CMH common OR = 2.87, p = 1.05e-23 (47d_distance_stratified_results.tsv)
- Loop logFC vs anchor mC: Spearman rho = -0.244, p = 4.2e-9, n = 565; hmC rho = +0.216; distance-adjusted partial beta = -2.24, p = 7.4e-10 (47e_correlation_results.tsv)

**What this shows:** Lost CTCF loop anchors are ~3.3-fold enriched for mC hypermethylation specifically at their flanking dynamic CpG regions, while constitutively unmethylated CpG islands show no effect. The negative loop-strength correlation and the strong distance gradient (effect 11× at <200 kb, gone >1 Mb) implicate local CTCF-site hypermethylation in the loss of short-range loops — the methylation-side confirmation of the Flavahan/Bernstein IDH-glioma insulator-decay model, carrying the same coordinated mC-up/hmC-down signature.

**Figures:**
- 47a_cpgi_* (mc_direction_barchart, mc_violin, hmc_violin) — CpG-island null panels
- 47a_dynamic_* (mc_direction_barchart, mc_violin, hmc_violin) — dynamic-region core test
- 47b_multiregion_OR_forest, 47b_multiregion_pct_heatmap — Fisher OR across region types
- 47c_coordinated_barchart, 47c_scatter_mc_vs_hmc, 47c_combined_effect_violin — coordinated pattern
- 47d_distance_stratified_OR, 47d_mc_violin_by_distance — distance stratification + CMH
- 47e_logfc_vs_mc_scatter, 47e_logfc_vs_hmc_scatter, 47e_residual_partial_correlation — loop-strength vs methylation
