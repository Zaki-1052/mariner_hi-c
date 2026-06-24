## Section 60: MeCP2-Status Epigenetic Profiles
**Key numbers:**
- MeCP2-Up, K119ub BigWig mean = +0.630 (q = 1.77e-19, ***); MeCP2-Down = -0.039 (q = 0.182, ns) (source: 60_mecp2_status_stats.tsv)
- MeCP2-Up: 5mC = +0.049 (q = 9.2e-13), 5hmC = -0.018 (q = 3.0e-08), K27ac = -0.718 (q = 1.2e-09), K27me3 = +1.110 (q = 0.018) (source: 60_mecp2_status_stats.tsv)
- MeCP2-Down mirror: 5mC = -0.028, K27ac = +0.481, K27me3 = -0.337 (all q < 0.05) (source: 60_mecp2_status_stats.tsv)
- 35 MeCP2-Down genes exported (source: 60_mecp2_down_gene_table.tsv)

**What this shows:** MeCP2-Up and MeCP2-Down loci display mirror-image epigenetic shifts across 5mC, 5hmC, K27ac, and K27me3, but K119ub is significantly elevated only at MeCP2-Up loci and flat (ns) at MeCP2-Down loci. Up-vs-Down Mann-Whitney is significant for every mark, formalizing the cascade and its diagnostic K119ub asymmetry.

**Figures:**
- 60a_violin_marks_by_mecp2 — 5-facet violin+box with ggsignif brackets
- 60b_summary_heatmap — mean-change heatmap (5 marks x Up/Down) with stars
- 60c_effect_size_lollipop — mean +/- 95% CI for 6 marks
- 60d_down_genes_strip — per-gene strip for the 35 MeCP2-Down genes
- 60e_composite — patchwork composite
