## Section 28: Coordinated Q4 (mC↑/hmC↓) Gene Characterization
**Key numbers:**
- Quadrant breakdown: Q1 = 411, Q2 = 1,255, Q3 = 116, Q4 (coordinated) = 6,589, non-Q4 = 1,782; total co-significant = 8,371 (coordinated_gene_characteristics.tsv)
- Q4 is 78.7% of co-significant genes — the dominant mC↑/hmC↓ pattern (coordinated_gene_characteristics.tsv)
- Hi-C loop involvement: Q4 = 333/6,589 = 5.1% at loop anchors vs non-Q4 = 210/1,782 = 11.8% (coordinated_gene_characteristics.tsv)

**What this shows:** The coordinated Q4 class (mC up / hmC down) dominates genes significant in both marks, reaffirming a TET-mediated demethylation block as the prevailing direction of change. Across the 9 epigenomic dimensions profiled, Q4 genes are characterized in the 28a composite (effect size, expression, ATAC, K119ub, K27ac, chromatin state, MeCP2, loop). Notably Q4 genes are *less* likely to sit at loop anchors than non-Q4 genes (5.1% vs 11.8%), indicating the coordinated phenotype is primarily a gene-body chemical-conversion effect rather than a 3D-architecture effect. The 9-dimension Fisher/Wilcoxon statistics are console-only; only the master per-gene table is persisted.

**Figures:**
- 28a_coordinated_composite — 3×3 grid, Q4 vs non-Q4 across 9 dimensions
- 28a_panel_mc_diff / 28a_panel_hmc_diff — methylation effect-size boxplots
- 28a_panel_log2fc — RNA-seq log2FC
- 28a_panel_atac — net ATAC change
- 28a_panel_k119ub — K119ub gene-body signal
- 28a_panel_k27ac — H3K27ac gained/lost
- 28a_panel_chromatin — chromatin-state stacked bar
- 28a_panel_mecp2 — MeCP2 up/down
- 28a_panel_loop — Hi-C loop involvement
- 28b_mc_hmc_concordance_scatter — mc_diff vs hmc_diff colored by concordance
- 28c_mc_vs_expression_per_group — mC change vs RNA-seq log2FC
- 28d_go_enrichment_coordinated — GO BP dot plot
- 28e_all_quadrants_comprehensive — 4-quadrant epigenomic panel
