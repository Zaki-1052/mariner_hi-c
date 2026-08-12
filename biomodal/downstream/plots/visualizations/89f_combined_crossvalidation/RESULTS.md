## Section 89: Cross-Validation of CALDER2 Subcompartments vs Homer PC1 Eigenvalues

**Key numbers:**
- Spearman rho = 0.928 between CALDER2 continuous rank and Homer PC1 eigenvalue at 100kb (89_pc1_rank_correlation.tsv)
- 91.5% concordance on A/B classification, Cohen's kappa = 0.829 (89_concordance_matrix.tsv)
- A precision=89.0% / recall=92.5%; B precision=93.7% / recall=90.8% (89_concordance_matrix.tsv)
- All pairwise Wilcoxon and KS tests between adjacent subcompartments p < 2.2e-16 (89_pc1_by_subcompartment_summary.tsv)
- Switching concordance: Fisher's OR = 10.62, p < 2.2e-16 — CALDER2 label changes strongly predict Homer PC1 shifts (89_switching_concordance.tsv)
- 23,840 / 23,874 CALDER2 bins matched to Homer (99.9%); 98.3% had the expected 4 Homer sub-bins per 100kb window

**What this shows:** The CALDER2 hierarchical decomposition and Homer PC1 eigenvalue analysis — two completely independent methods — produce highly concordant compartment calls. The Spearman rho of 0.928 between CALDER2's continuous domain rank and Homer's PC1 eigenvalue confirms these are measuring the same underlying chromatin property. The 8.5% disagreement rate is concentrated at A/B boundaries (bins near PC1 = 0), not systematic misclassification. Switching concordance (OR = 10.62) validates that both tools identify the same compartment transitions between control and mutant. Uses adult/late (250402) timepoint.

**Figures:**
- 89a_pc1_violin_by_subcompartment/ — Homer PC1 eigenvalue distribution per CALDER2 subcompartment (Kruskal-Wallis)
- 89b_rank_vs_pc1_scatter/ — CALDER2 continuous rank vs Homer PC1 scatter (Spearman)
- 89c_concordance_heatmap/ — 2x2 confusion matrix with accuracy/kappa
- 89d_rank_density_by_subcompartment/ — CALDER2 continuous rank density overlays (KS tests)
- 89e_switching_concordance/ — compartment switching agreement (Fisher's test)
