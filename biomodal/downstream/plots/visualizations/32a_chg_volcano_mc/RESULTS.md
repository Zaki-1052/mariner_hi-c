## Section 32: CHG Context Exploratory Analysis
**Key numbers:**
- 320 significant CHG 5mC gene-body DMRs (157 hyper / 163 hypo); 92 significant CHG 5hmC DMRs (source: chg_mc_significant_dmrs.tsv, chg_hmc_significant_dmrs.tsv)
- 206/320 (64.4%) overlap CG mC, but 126/206 (61.2%) are direction-DISCORDANT with CG mC (source: chg_mc_significant_dmrs.tsv)
- 244/320 (76.3%) of CHG mC DMRs are on chr7+chr8 (chr8=149, chr7=95) (source: chg_mc_significant_dmrs.tsv)
- 114 CHG-unique genes (not CG mC sig), 55 also in CG hmC; KEGG only 2 pathways at q<0.05, GO BP 0 at q<0.05 (source: chg_unique_genes.tsv, chg_enrichment_kegg.tsv, chg_enrichment_go_bp.tsv)

**What this shows:** CHG (non-CG) methylation produces a statistically detectable but biologically marginal signal under BAP1 loss. Effect sizes are ~100x smaller than CG (plotted in permille vs percent), the hits concentrate on two chromosomes, and the majority of genes significant in both contexts change in opposite directions — so CHG does not corroborate the CG-context findings and is best read as low-level, largely CG-independent background. (Note: section-script comments saying "70 genes" are stale; run-5 count is 320.)

**Figures:**
- 32a_chg_volcano_mc / 32b_chg_volcano_hmc — CHG 5mC and 5hmC volcano plots (permille effect axis)
- 32c_chg_direction_breakdown — hyper vs hypo counts for CHG mC and hmC
- 32d_chg_cg_venn_mc / 32e_chg_cg_venn_three_way — CHG-vs-CG significant-gene Venn diagrams
- 32f_chg_cg_concordance_scatter / 32g_chg_cg_concordance_bar — CHG vs CG direction concordance
- 32h_chg_top_dmr_genes / 32i_chg_unique_genes — ranked all-gene and CHG-unique gene charts
- 32j_chg_chromosome_distribution / 32k_chg_vs_cg_chr_comparison — chromosome distribution (chr7/chr8 enrichment)
- 32l_chg_enrichment_go / 32m_chg_enrichment_kegg — GO BP and KEGG enrichment (exploratory thresholds)
