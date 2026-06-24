## Section 39: H3K36me2 Polycomb Boundary Analysis (NSD/DNMT3A axis)
**Key numbers:**
- me2 peaks at H3K27me3 boundaries (±5kb): n=660, median Log2FC=+0.636; away: n=5,966, median=+0.524; Wilcoxon p=2.97e-47 (h3k36me2_k27me3_boundary_analysis.tsv)
- Significant H3K36me2 gene-level peaks: 2,909 genes (1,557 gained / 1,352 lost) (h3k36me2_gene_level_summary.tsv)
- me2 vs 5mC diff: Spearman rho = -0.185; me2 vs 5hmC diff: rho = +0.320 (h3k36_expanded_correlations.tsv)
- me2 vs H3K27me3: rho = -0.483 (strongest anticorrelation in the cross-mark matrix) (h3k36_expanded_correlations.tsv)

**What this shows:** H3K36me2 — the NSD-deposited, DNMT3A-recruiting, anti-PRC2 barrier mark — is broadly gained in BAP1-KO and gained slightly more at H3K27me3 boundary zones, consistent with expansion against Polycomb domains (strong me2 ↔ H3K27me3 anticorrelation). But me2 gain is negatively correlated with 5mC gain and positively with 5hmC, so the NSD/H3K36me2 → DNMT3A route is not driving gene-body hypermethylation either.

**Figures:**
- 39a_h3k36me2_volcano_annotation — DiffBind volcano + up-vs-down annotation bar
- 39b_me2_x_mc_oe_heatmap — me2 direction × 5mC DMR direction O/E
- 39c_me2_at_k27me3_boundary — me2 fold at vs away from H3K27me3 boundaries (Wilcoxon)
- 39d_me2_x_k27me3_oe_heatmap — me2 × H3K27me3 mutual-exclusivity O/E
- 39e_me2_vs_methylation_scatter — me2 fold vs 5mC and 5hmC change
- 39f_me2_by_chromatin_state — me2 gained/lost/none at mC DMRs by chromatin state
- 39g_me2_genic_vs_intergenic — me2 fold genic vs intergenic (Wilcoxon)
