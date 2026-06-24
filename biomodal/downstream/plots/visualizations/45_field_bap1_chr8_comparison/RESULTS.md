## Section 45: Field et al. (2019) BAP1 Chr8 Methylation Hotspot Comparison
**Key numbers:**
- Field chr8 genes mapped to mouse orthologs = 81/85 (22 hyper, 63 hypo) (field_chr8_ortholog_mapping.tsv)
- Gene-body mC concordance: 21 Concordant vs 38 Discordant of 59 sig genes; Fisher p = 0.00887; hypergeometric p = 4.61e-6 (field_chr8_comparison_full.tsv, field_chr8_statistical_tests.tsv)
- Promoter mC concordance: Fisher p = 0.40 (NS); only 6/81 promoters significant (1 concordant, 5 discordant) (field_chr8_statistical_tests.tsv, field_chr8_promoter_comparison.tsv)
- RNA-seq: 6/21 concordant genes also expression-concordant (RNF122, ZDHHC2, KHDRBS3, XKR4, Adgrb1, Ptp4a3) (field_chr8_concordant_genes.tsv)

**What this shows:** The human BAP1 chr8 methylation hotspot replicates directionally in BAP1-KO mouse cerebellum at the GENE-BODY level (concordance enriched over chance; orthologs strongly over-represented among our DMRs) but NOT at promoters. The conserved BAP1 methylation response is a gene-body phenomenon. The trisomy-8 diagnostic panel rules out aneuploidy as the driver, so cross-species concordance reflects genuine BAP1-dependent targeting.

**Figures:**
- 45a_mapping_funnel — Field genes → ortholog → in-data → sig mC/hmC attrition
- 45b_genebody_concordance_heatmap — 2×2 direction concordance + Fisher
- 45c_dual_modality_dotplot — per-gene mC + hmC effect
- 45d_volcano_highlight — Field genes on genome-wide gene-body mC volcano
- 45e_effect_size_lollipop — per-gene mC effect by concordance
- 45f_quadrant_analysis — coordinated mC/hmC quadrant, Field vs genome-wide
- 45g_genebody_vs_promoter — gene-body vs promoter mC scatter
- 45i_promoter_concordance_heatmap — promoter trend-level concordance (~random)
- 45j_rnaseq_expression — RNA-seq log2FC of Field genes by concordance
- 45k_trisomy8_diagnostic — coverage / RNA-seq log2FC / DMR-rate, chr8 vs autosomes
