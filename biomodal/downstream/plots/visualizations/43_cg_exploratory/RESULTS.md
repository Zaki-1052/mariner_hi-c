## Section 43: CG Context Exploratory Analysis
**Key numbers:**
- 12 figure panels: direction breakdown, mC/hmC Venn, per-chromosome enrichment (Fisher O/E), effect-size and direction by chromosome, chrX-removal asymmetry test
- Quantitative outputs (per-chromosome significant-DMR counts, Fisher O/E, Spearman gene-density correlation, Kruskal-Wallis effect-size) are computed at runtime from mc_dmr/hmc_dmr
- Exported TSVs: cg_mc_chromosome_distribution.tsv, cg_hmc_chromosome_distribution.tsv, cg_top50_mc_dmr_genes.tsv (not in this group's verified starter set)

**What this shows:** An exploratory profile of the primary CG significant DMRs across chromosomes. The chrX-removal panel (43m) demonstrates the global hyper/hypo direction asymmetry is autosomal and not an artifact of X-linked DMRs in the mixed-sex cohort; the gene-density panel (43k) controls for the trivial "gene-dense chromosomes have more hits" confound. (Specific per-chromosome counts are not cited here as they were not confirmed against the exported TSVs.)

**Figures:**
- 43a_cg_direction_breakdown — hyper/hypo bars for mC and hmC
- 43b_cg_mc_hmc_venn — significant-gene overlap (5mC vs 5hmC)
- 43c/43d_cg_*_chromosome_distribution — per-chr counts vs expected, Fisher labels
- 43e_cg_mc_vs_hmc_chr_comparison — normalized % DMRs per chr
- 43f_cg_chr_obs_exp_heatmap — observed/expected ratio heatmap (PNG/PDF only)
- 43g/43h_cg_*_chr_direction — stacked direction per chromosome
- 43i_cg_regional_chr_enrichment — chromosome enrichment across genomic regions
- 43j_cg_top_dmr_genes — top 50 mC DMR genes by q-value
- 43k_cg_mc_sig_rate_vs_gene_density — sig rate vs total genes (Spearman)
- 43l_cg_mc_effect_size_by_chr — effect-size boxplots per chromosome (Kruskal-Wallis)
- 43m_chrX_direction_comparison — direction asymmetry with vs without chrX
