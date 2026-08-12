## Section 20: Coordinated Gene RNA-seq Expression Integration
**Key numbers:**
- Q4 coordinated (mC up/hmC dn) = 6,589 of 8,371 co-significant genes = 78.7% (source: coordinated_rnaseq_expression.tsv)
- Tier 1 (padj<0.05) Q4 expression: Down = 787 (13.2%), Up = 157 (2.6%); Down/Up ratio = 5.01 (source: coordinated_rnaseq_expression.tsv, col expr_tier1)
- Tier 2 (padj<0.05 & |log2FC|>0.3) Q4: Down = 666, Up = 38; Down/Up ratio = 17.53 (source: coordinated_rnaseq_expression.tsv, col expr_tier2)
- Combined methylation effect vs log2FC across Q4: Spearman rho = -0.249 (n = 5,961 with RNA-seq) (recomputed from coordinated_rnaseq_expression.tsv)

**What this shows:** Among genes co-significant in both 5mC and 5hmC, the dominant coordinated quadrant (Q4: mC gain / hmC loss) is transcriptionally repressive. When expression changes significantly, Q4 genes are downregulated 5:1 (Tier 1) to 17.5:1 (Tier 2), and larger combined methylation change tracks with lower expression (rho = -0.249). Most Q4 genes (84%) are expression-unchanged at this timepoint, indicating the methylation lesion is largely subthreshold for transcription.

**Figures:**
- `20a_coordinated_expression_breakdown` — stacked bar of Up/Down/Unchanged for Coordinated vs Other mC DMR vs All genes (Tier 1 / Tier 2 facets, Fisher OR).
- `20b_methylation_vs_expression_scatter` — combined methylation effect vs log2FC scatter with lm fit and Spearman rho.
- `20c_log2fc_violin_comparison` — log2FC violin/box, Coordinated vs Other significant mC DMR genes (Wilcoxon).
- `20d_mc_expression_heatmap` — 2x2 5mC-direction x expression-direction enrichment heatmap.

_Note: prose docs (FIGURES.md Fig 32, concordant_vs_discordant_trends.md) report stale run-4 counts (Q4=5,708; 85%); see docs/results/section_20_21_rnaseq_quadrants.md for the reconciliation._
