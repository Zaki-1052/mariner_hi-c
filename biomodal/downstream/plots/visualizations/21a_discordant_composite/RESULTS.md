## Section 21: Discordant (Q2) vs Coordinated (Q4) mC/hmC Gene Characterization
**Key numbers:**
- Quadrant counts: Q1 = 411, Q2 (discordant) = 1,255, Q3 = 116, Q4 (coordinated) = 6,589 (source: discordant_gene_characteristics.tsv, col quadrant)
- Expression direction (padj<0.05): Q2 = 190 Up / 10 Down vs Q4 = 157 Up / 787 Down; Fisher sample OR ≈ 95, p ≈ 1.2e-104 (source: discordant_gene_characteristics.tsv)
- K119ub gene-body log2FC median: Q2 = -0.0385 (loss) vs Q4 = +0.0302 (gain) (source: discordant_gene_characteristics.tsv, col k119ub_gb_log2fc)
- Hi-C loop anchors: Q2 = 13.1% (165/1,255) vs Q4 = 5.1% (333/6,589); Fisher sample OR ≈ 2.84, p ≈ 1.1e-22 (source: discordant_gene_characteristics.tsv, col loop_involved)

**What this shows:** The minority discordant quadrant (Q2: mC loss / hmC gain) is the near-inverse of the dominant coordinated quadrant (Q4) across nine dimensions. Q2 genes are upregulated, lose K119ub, gain accessibility and H3K27ac, sit predominantly at Repressed_Promoter chromatin (67.4% vs Q4's 63.8% Active_Promoter), and are ~2.6x more likely to be at loop anchors. Q4 genes show the opposite: repression, K119ub gain, MeCP2-up bias (~5.9:1), and active-promoter chromatin. Q2 has smaller methylation effect sizes and no GO enrichment vs Q4 background.

**Figures:**
- `21a_discordant_composite` (+ `21a_panel_*` for mc_diff, hmc_diff, log2fc, atac, k119ub, k27ac, chromatin, mecp2, loop) — 3x3 Q2-vs-Q4 comparison panels.
- `21b_mc_hmc_concordance_scatter` — mC vs hmC scatter, Q4/Q2 highlighted over background.
- `21c_mc_vs_expression_per_group` — faceted mC-change vs log2FC with per-group Spearman.
- `21e_all_quadrants_comprehensive` — 4-panel log2FC / net ATAC / K119ub / chromatin across all quadrants.
- (`21d_go_enrichment_discordant` not produced — GO returned no significant terms.)

_Note: prose docs report stale run-4 counts (Q2=755, Q4=5,708); the numbers above are from the current run-5 TSV. See docs/results/section_20_21_rnaseq_quadrants.md._
