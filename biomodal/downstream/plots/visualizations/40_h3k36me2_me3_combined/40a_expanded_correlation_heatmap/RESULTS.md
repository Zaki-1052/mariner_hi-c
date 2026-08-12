## Section 40: H3K36me2/me3 Combined Analysis and Conversion Dynamics
**Key numbers:**
- H3K36me2 vs H3K36me3 gene-level correlation: Spearman rho = -0.082 (marks essentially decoupled) (h3k36_expanded_correlations.tsv)
- H2AK119ub vs 5mC diff: rho = +0.410 — strongest positive correlate of hypermethylation in the 9-variable matrix (h3k36_expanded_correlations.tsv)
- H3K36me3 vs 5mC diff: rho = +0.018; H3K36me3 vs H2AK119ub: rho = -0.040 (both near zero) (h3k36_expanded_correlations.tsv)
- Combined profile: 20,915 genes; 10,301 with me3 fold, 2,178 with me2 fold, 1,517 with both (h3k36_combined_gene_profile.tsv)

**What this shows:** Analyzed jointly, H3K36me2 and me3 move nearly independently (rho=-0.082) — no single SETD2-driven me2→me3 conversion flux unites them. The clustered cross-mark correlation heatmap is the clearest single summary of the DNMT question: among all six histone/accessibility marks, only H2AK119ub correlates positively with mC gain (+0.41), while both H3K36 marks sit near zero or negative. The hypermethylation axis is ubiquitin/Polycomb-linked, not H3K36/transcription-linked.

**Figures:**
- 40a_expanded_correlation_heatmap — clustered 9-variable Spearman correlation heatmap
- 40b_me2_vs_me3_scatter — me2 vs me3 fold with conversion-shift/block quadrants (console-only rho/quadrant counts)
- 40c_me2_me3_ratio_delta — (me2 − me3) fold delta by DMR status (Kruskal-Wallis)
- 40d_three_way_venn — gene overlap of sig me2, sig me3, mC DMR sets
- 40e_direction_flow — direction patterns for triple-significant genes
- 40f_go_comparison — GO BP for me2-only vs me3-only vs shared genes
