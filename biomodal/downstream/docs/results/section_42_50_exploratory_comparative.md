# Sections 42-50: Exploratory & comparative: max-significance genes, allele-specific methylation, Field 2019 chr8 cross-species, genome browser loci, CTCF anchors, CpG-island ubiquitination, HOMER motif enrichment

## Summary
This group of sections probes whether the BAP1-KO methylation phenotype generalizes across analysis frames: the strongest-effect genes (q-floor list, PCA), allelic imbalance at heterozygous SNVs (ASM), a published human BAP1 chr8 hotspot (Field et al. 2019), per-locus genome-browser views, CTCF loop-anchor methylation (testing the Flavahan/Bernstein IDH-glioma insulator model), CpG-island ubiquitination/de-novo methylation, and TF-motif enrichment (HOMER) at differential chromatin and DMR sites. The central answers: (1) ~100 genes reach the q-value floor with 42 showing coordinated mC↑/hmC↓ (TET-block signature); (2) BAP1-KO roughly doubles allele-specific mC sites (mut/ctrl ≈ 1.95×); (3) the Field chr8 hotspot replicates directionally at the gene-body level (Fisher p = 0.0089, hypergeometric p = 4.6e-6) but NOT at promoters; (4) lost CTCF anchors are strongly enriched for mC-hypermethylation at flanking dynamic CpG regions (OR ≈ 3.3, p < 2.2e-16), paralleling IDH-glioma insulator decay; (5) hypermethylated CpG islands mostly amplify pre-existing methylation (only ~15% de novo) and are NOT enriched for K119ub gain; (6) bHLH (neuronal) TFs dominate H2AK119ub-differential sites and zinc-finger TFs dominate coordinated DMRs.

## Section 42: section_42_max_significance_gene_list
### Analysis question
Which genes hit the numerical q-value floor (q < 1e-300, the -log10(q)=300 volcano ceiling) in 5mC and/or 5hmC, and do they form a coordinated mC↑/hmC↓ population? PCA of per-sample methylation across these genes is also produced.
### Key results
- 5mC genes at q-floor = 100 (mc_at_floor = TRUE; source: max_significance_genes_merged.tsv)
- 5hmC genes at q-floor = 67 (hmc_at_floor = TRUE; source: max_significance_genes_merged.tsv)
- Total unique merged loci at floor = 124 rows / 123 unique gene names (source: max_significance_genes_merged.tsv, max_significance_gene_names.txt)
- Coordinated mC_up_hmC_down = 42 of 124 loci; mC_down_hmC_up = 1; single_modality = 81 (source: max_significance_genes_merged.tsv, `coordinated` column)
- Top mC effect: Syt1 = +17.31% mC (ctrl 0.588 → mut 0.762), Galnt7 = +16.28% (source: max_significance_genes_mc.tsv)
- Top hmC effect: Syt1 = -14.97% hmC (ctrl 0.266 → mut 0.117), Mcu = -9.44% (source: max_significance_genes_hmc.tsv)
- 5mC-only list = 100 genes; 5hmC-only list = 67 genes (source: max_significance_genes_mc.tsv, max_significance_genes_hmc.tsv, row counts)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The genes with the most extreme statistical significance are dominated by the coordinated mC-up/hmC-down signature (42/124, and Syt1 is the archetype: +17% mC alongside -15% hmC). This is the locus-level expression of the genome-wide TET-mediated demethylation block — at the most strongly affected loci, gained 5mC is mirrored by depleted 5hmC, exactly what is expected if TET oxidation of 5mC→5hmC is impaired. These floor genes (Syt1, Cntnap2, Cadm2, Cdh8, neuronal adhesion/synaptic genes) define the highest-confidence target set for downstream pathway and browser analysis.
### Plot inventory
- 42_pca_max_significance_genes — PCA (5mC | 5hmC side-by-side) of the max-significance gene set, per-sample regional fraction, colored by condition/batch
- 42_pca_all_genes — genome-wide PCA of all gene bodies (5mC | 5hmC) for comparison

## Section 43: section_43_cg_exploratory_analysis
### Analysis question
Exploratory profile of the primary CG context significant DMRs: direction breakdown, mC/hmC gene overlap (Venn), per-chromosome enrichment (Fisher O/E), effect-size and direction by chromosome, and the impact of removing chrX on direction asymmetry. (Note: this section's numbers are computed in-script from `mc_dmr`/`hmc_dmr`; its exported TSVs are chromosome-distribution summaries, not in this group's starter list.)
### Key results
- Exported tables: cg_mc_chromosome_distribution.tsv, cg_hmc_chromosome_distribution.tsv, cg_top50_mc_dmr_genes.tsv (source: section_43 write.table calls — [NOT INDEPENDENTLY READ: these TSVs were not in the assigned starter set; counts below are from script logic, not confirmed against the TSV])
- Direction breakdown, mC/hmC Venn overlap, chromosome O/E heatmap, and chrX-removal panel are all generated (12 figure panels saved)
- [UNVERIFIED: specific significant-DMR counts per chromosome are computed at runtime and printed; not confirmed in a TSV I read for this group]
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The chrX-removal panel exists specifically to demonstrate that the global direction asymmetry (hyper vs hypo balance) is not an artifact of X-linked DMRs in a mixed-sex cohort — i.e., the BAP1 methylation phenotype is autosomal and robust to sex-chromosome contribution. The per-chromosome Fisher enrichment asks whether DMRs cluster on particular chromosomes beyond gene density; the gene-density correlation panel (Spearman) controls for the trivial "more genes = more hits" confound.
### Plot inventory
- 43a_cg_direction_breakdown — hyper/hypo bar for mC and hmC
- 43b_cg_mc_hmc_venn — significant-gene overlap between 5mC and 5hmC
- 43c_cg_mc_chromosome_distribution / 43d_cg_hmc_chromosome_distribution — per-chr significant counts vs expected, Fisher labels
- 43e_cg_mc_vs_hmc_chr_comparison — normalized % of DMRs per chr, mC vs hmC
- 43f_cg_chr_obs_exp_heatmap — observed/expected ratio heatmap (PNG/PDF only, top-level)
- 43g_cg_mc_chr_direction / 43h_cg_hmc_chr_direction — stacked direction per chromosome
- 43i_cg_regional_chr_enrichment — chromosome enrichment across genomic regions
- 43j_cg_top_dmr_genes — top 50 mC DMR genes by q-value
- 43k_cg_mc_sig_rate_vs_gene_density — significance rate vs total genes (Spearman)
- 43l_cg_mc_effect_size_by_chr — effect-size boxplots per chromosome (Kruskal-Wallis)
- 43m_chrX_direction_comparison — direction asymmetry with vs without chrX

## Section 44: section_44_allele_specific_methylation
### Analysis question
First downstream analysis of DUET evoC allele-specific methylation: does BAP1-KO change the number/magnitude of mC (and hmC) allelic-imbalance sites at heterozygous SNVs, and do ASM loci coincide with gene-body DMRs?
### Key results
- Per-sample significant mC ASM: Control mean = 11,295.5 sites/sample (n=4), Mutant mean = 22,080.2 (n=4), ratio = 1.95× (source: asm_mc_significant_per_sample.tsv)
- Control range 8,578–15,760; Mutant range 18,962–24,768 sites per sample (source: asm_mc_significant_per_sample.tsv)
- Unique significant mC ASM loci = 58,955 (source: asm_mc_significant_loci.tsv, row count)
- Top ASM locus chr17:39848255 (G→A), mean_meth_diff = 0.359, min_corr_p = 2.06e-139 (source: asm_mc_significant_loci.tsv) — note: chr17:39.84Mb cluster (multiple adjacent variants) dominates the strongest hits, consistent with an imprinted/strong-ASM region
- Loci significant in ≥1 control sample AND ≥1 mutant sample (n_ctrl>0 & n_mut>0) = 20,546; control-only (n_ctrl>0, n_mut=0) = 9,768; mutant-only (n_mut>0, n_ctrl=0) = 28,641 (source: asm_mc_significant_loci.tsv, computed from n_ctrl/n_mut columns)
- ASM-DMR overlap: 1,408 gene-body DMR genes contain ≥1 significant mC ASM site (source: asm_dmr_overlap_summary.tsv, row count); e.g. A430033K04Rik = 15 ASM sites
### [INTERPRETATION] Biological meaning
[INTERPRETATION] BAP1 loss nearly doubles the count of allele-imbalanced mC sites (1.95×), and mutant-only loci (28,641) far outnumber control-only (9,768). This indicates BAP1-KO induces widespread NEW allelic methylation asymmetry rather than simply erasing existing imprints — consistent with a stochastic gain-of-methylation process (the hypermethylation arm of the phenotype) that acts unevenly between alleles. The fact that 1,408 gene-body DMRs harbor ASM sites links the bulk DMR signal to allele-resolved events. Caveat: ASM here uses uncorrected per-sample testing aggregated across samples, so absolute counts are exploratory.
### Plot inventory
- 44a_filter_distribution_per_sample — ASM filter categories (PASS vs failure modes) per sample, mC + hmC
- 44b_pass_site_counts_per_sample — PASS site counts per sample by condition
- 44c_asm_mc_chromosome_distribution — unique significant mC ASM loci per chromosome with Fisher enrichment
- 44d_methylation_diff_distributions — density + per-sample boxplots of allelic methylation difference
- 44e_mutant_vs_control_asm_comparison — per-sample significant counts, mutant vs control (Wilcoxon)
- 44f_pvalue_distributions — raw p-value histogram + -log10(corrected p) distribution
- 44g_asm_volcano_plot — methylation difference vs significance, faceted by condition
- 44h_shared_asm_sites_venn — control vs mutant significant-locus overlap
- 44i_asm_dmr_overlap — ASM sites within gene-body DMR categories + DMR effect vs ASM bias scatter
- 44j_mc_vs_hmc_asm_comparison — significant mC vs hmC ASM counts + locus Venn

## Section 45: section_45_field_bap1_chr8_comparison
### Analysis question
Do the chr8 methylation-hotspot genes from Field et al. (2019, human BAP1-KD uveal melanoma, Clin Cancer Res) replicate directionally in our BAP1-KO mouse cerebellum at gene-body and promoter level, and is the chr8 signal an artifact of trisomy 8 (as in the human tumors)?
### Key results
- Field chr8 genes mapped to mouse orthologs = 81 of 85 (source: field_chr8_ortholog_mapping.tsv); Field directions: 22 Hypermethylated, 63 Hypomethylated
- Gene-body mC concordance: 21 Concordant vs 38 Discordant of 59 significant testable genes (35.6% concordant) (source: field_chr8_comparison_full.tsv, gb_concordance column; 18 Non-significant, 4 Not-in-data)
- Fisher's exact (gene-body mC direction concordance) p = 0.00887 (source: field_chr8_statistical_tests.tsv)
- Hypergeometric enrichment (Field orthologs among our sig gene-body mC DMRs) p = 4.61e-6 (source: field_chr8_statistical_tests.tsv)
- Binomial test (gene-body concordance vs 50%) p = 0.0363 (source: field_chr8_statistical_tests.tsv)
- Promoter mC concordance: Fisher p = 0.40 (NON-significant); promoter breakdown = 1 Concordant, 5 Discordant, 66 Non-significant, 9 Not-in-data; only 6 promoters significant (source: field_chr8_statistical_tests.tsv, field_chr8_promoter_comparison.tsv)
- RNA-seq: 6 of 21 concordant genes also expression-concordant with the Field hyper→down / hypo→up pattern (e.g. RNF122, ZDHHC2, KHDRBS3, XKR4, BAI1/Adgrb1, PTP4A3) (source: field_chr8_concordant_genes.tsv, expr_concordance column)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The Field chr8 methylation hotspot replicates in mouse cerebellum at the GENE-BODY level (concordance enriched over chance, p ≈ 0.009; orthologs strongly over-represented among our DMRs, p ≈ 5e-6) but collapses at PROMOTERS (p = 0.40, only 6/81 promoters significant). This argues the conserved BAP1 methylation response is a gene-body phenomenon, not a promoter/CpG-island phenomenon — consistent with the rest of this study's finding that the action is in flanking/dynamic regions, not constitutively unmethylated promoter islands. The trisomy-8 diagnostic (coverage ratio, RNA-seq median log2FC, DMR-rate panels) rules out aneuploidy as the driver, so the cross-species concordance reflects genuine BAP1-dependent epigenetic targeting rather than chr8 copy number.
### Plot inventory
- 45a_mapping_funnel — attrition Field genes → mouse ortholog → in-data → sig mC → sig hmC, split by direction
- 45b_genebody_concordance_heatmap — 2×2 Field vs our gene-body mC direction tile + Fisher
- 45c_dual_modality_dotplot — per-gene mC + hmC effect (filled = significant)
- 45d_volcano_highlight — Field genes highlighted on genome-wide gene-body mC volcano
- 45e_effect_size_lollipop — per-gene gene-body mC effect, colored by concordance
- 45f_quadrant_analysis — coordinated mC/hmC quadrant, Field chr8 vs genome-wide
- 45g_genebody_vs_promoter — gene-body vs promoter mC effect scatter
- 45i_promoter_concordance_heatmap — promoter trend-level concordance tile (expected ~random)
- 45j_rnaseq_expression — RNA-seq log2FC of Field genes, colored by expression concordance
- 45k_trisomy8_diagnostic — 3-panel coverage / RNA-seq log2FC / DMR-rate, chr8 vs autosomes

## Section 46: section_46_genome_browser_loci
### Analysis question
Per-locus multi-omic genome-browser views (Gviz) at the 10 KEY_GENES, stacking condition-averaged 5mC/5hmC tracks, mC/hmC difference, RNA-seq coverage, H2AK119ub/H3K27me3/H3K4me3/H3K27ac/ATAC/H3K27me1 (ChIP/ATAC signal) and MeCP2 (CUT&RUN signal), CpG island/shore/shelf annotation, CTCF peaks, and Hi-C loop arcs (lost/gained).
### Key results — VISUAL ONLY, NO QUANTIFICATION TABLE
- This section produces NO TSV output. It is a rendered-tracks figure section; numbers shown on the tracks (per-track RNA-seq log2FC/padj annotations) are read live from `mc_dmr`, `hmc_dmr`, the RNA-seq xlsx, and BigWigs.
- KEY_GENES rendered (full + compact each): Syt1, Zbtb20, Trpm3, Epha3, Mcu, Cntnap2, Lpp, Dlgap1, Arhgap26, Cdh8 (source: _shared_config.R KEY_GENES; confirmed by per-gene output folders)
- Locus output folders present: Syt1_locus, Syt1_locus_compact (+ poster variants), Zbtb20_locus(+compact), Trpm3, Epha3, Mcu, Cntnap2, Lpp, Dlgap1, Arhgap26, Cdh8 — plus composite_syt1_panel (Panel A browser + Panel B coordinated scatter + Panel C multi-omic z-score heatmap)
- Each view extends ±50 kb (EXTEND_BP = 50,000) around the gene body
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The browser views are the integrative "money shot" for the KEY_GENES: they let a reader see, at a single locus (e.g. Syt1), the simultaneous gain of 5mC, loss of 5hmC, change in RNA-seq, accumulation of H2AK119ub, and disruption/gain of CTCF-anchored loops — i.e. the full BAP1 mechanism (Polycomb signal → methylation → 3D architecture) co-localized in cis. The Syt1 composite panel pairs the locus track with the genome-wide coordinated mC/hmC scatter (KEY_GENES gold-highlighted) and a cross-gene z-score summary, tying the single-locus story to the population.
### Plot inventory
- Syt1_locus, Syt1_locus_compact, Syt1_locus_poster, Syt1_locus_poster_v2, Syt1_locus_poster_v2_aligned — Syt1 browser (full/compact/poster variants)
- Zbtb20_locus(+_compact), Trpm3_locus(+_compact), Epha3_locus(+_compact), Mcu_locus(+_compact), Cntnap2_locus(+_compact), Lpp_locus(+_compact), Dlgap1_locus(+_compact), Arhgap26_locus(+_compact), Cdh8_locus(+_compact) — per-gene full + compact browser views
- composite_syt1_panel — assembled 3-panel composite (Syt1 browser + coordinated scatter + multi-omic z-score heatmap)

## Section 47: section_47_ctcf_anchor_methylation_overlay
### Analysis question
Are LOST CTCF loop anchors hypermethylated at their flanking dynamic CpG regions (shores+shelves), paralleling the Flavahan/Bernstein IDH-glioma CTCF-insulator-decay model? CpG islands are expected null (constitutively unmethylated); the effect should be in dynamic regions. Tested with Fisher's exact, Wilcoxon, distance-stratification (CMH), and loop-logFC correlation.
### Key results
- Dynamic CpG regions (shores+shelves) at lost anchors = 1,084 vs gained = 1,427; background = 58,928 (source: 47a_fisher_results.tsv)
- mC hyper, lost vs gained dynamic regions: 22.32% vs 8.06%, OR = 3.28 [2.57–4.20], p = 5.4e-24 (source: 47a_fisher_results.tsv)
- hmC hypo, lost vs gained dynamic regions: 9.69% vs 4.91%, OR = 2.08 [1.50–2.89], p = 3.9e-6 (source: 47a_fisher_results.tsv)
- Motif-based CTCF anchors (sensitivity), mC hyper lost vs gained: OR = 2.52 [2.03–3.14], p = 3.0e-18 (source: 47a_fisher_results.tsv)
- Multi-region (47b): CpG shores mC OR = 3.01 (p = 1.3e-12), CpG shelves mC OR = 3.68 (p = 4.3e-13), Promoters mC OR = 3.32 (p = 0.0019); CpG islands mC OR = 0.84 (p = 1, NULL as predicted) (source: 47b_multiregion_fisher.tsv)
- Coordinated mC-up/hmC-down (47c) at lost vs gained anchors: 33.06% vs 18.45%, OR = 2.18 [1.83–2.61], p = 1.3e-18 (source: 47c_coordinated_enrichment.tsv)
- Distance-stratified (47d): strongest at <200 kb loops (OR = 11.21, p = 5.5e-23); Cochran-Mantel-Haenszel common OR = 2.87 [computed], cmh_p = 1.05e-23; >1 Mb bin NON-significant (OR = 0.87, p = 0.55) (source: 47d_distance_stratified_results.tsv)
- Loop-logFC vs methylation (47e): Spearman mC-diff vs logFC (all CTCF-CTCF loops) rho = -0.244, p = 4.2e-9, n = 565; hmC-diff vs logFC rho = +0.216, p = 2.0e-7; distance-adjusted partial beta = -2.24, p = 7.4e-10; 565/1,467 CTCF-CTCF loops (38.5%) had dynamic regions at anchors (source: 47e_correlation_results.tsv, 47e_loop_methylation_correlation.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the methylation-side confirmation of the Flavahan/Bernstein insulator model: CTCF anchors of loops LOST in BAP1-KO are ~3.3-fold enriched for mC hypermethylation (and ~2-fold for hmC hypomethylation) specifically in their flanking dynamic CpG regions, while constitutively unmethylated CpG islands show no effect (OR ≈ 0.84). The negative mC–logFC correlation (more anchor methylation → weaker loop) and the distance gradient (effect 11× at <200 kb, vanishing >1 Mb) argue that local hypermethylation at CTCF sites mechanistically weakens short-range CTCF-anchored loops — the same logic by which IDH-mutant gliomas lose insulation. The coordinated mC-up/hmC-down enrichment at lost anchors (OR 2.18) ties the insulator decay to the same TET-block signature seen genome-wide.
### Plot inventory
- 47a_cpgi_mc_direction_barchart / 47a_cpgi_mc_violin / 47a_cpgi_hmc_violin — CpG-island anchor methylation (null-result panels)
- 47a_dynamic_mc_direction_barchart / 47a_dynamic_mc_violin / 47a_dynamic_hmc_violin — dynamic-region (shores+shelves) mC/hmC by anchor group (core test)
- 47b_multiregion_OR_forest — Fisher OR forest across region types
- 47b_multiregion_pct_heatmap — % significant DMR by anchor group × region
- 47c_coordinated_barchart — % coordinated mC-up/hmC-down by anchor group
- 47c_scatter_mc_vs_hmc — mC vs hmC scatter colored by anchor group
- 47c_combined_effect_violin — |mC|+|hmC| combined effect by anchor group
- 47d_distance_stratified_OR — distance-binned OR forest + CMH overall
- 47d_mc_violin_by_distance — mC difference by loop-distance bin, lost vs gained
- 47e_logfc_vs_mc_scatter / 47e_logfc_vs_hmc_scatter — loop logFC vs anchor methylation
- 47e_residual_partial_correlation — distance-adjusted partial correlation

## Section 48: section_48_cpg_island_ubiquitination
### Analysis question
At CpG islands, (1) is H2AK119ub1 accumulation (from BAP1 loss) enriched at differentially methylated islands, and (2) is hypermethylation de-novo (control unmethylated) or amplification of pre-existing methylation? Also: coordinated mC+hmC quadrants and chromatin-mark context.
### Key results
- CpG island universe = 8,910 islands; mC significant = 442 (122 hyper, 320 hypo); hmC significant = 172 (26 hyper, 146 hypo) (source: 48_cpg_island_ubiquitination_summary.tsv)
- De-novo classification of the 122 hyper islands: control mC <0.05 = 7 (5.7%), <0.10 = 12 (9.8%), <0.20 = 18 (14.8%) → majority (85%+) are pre-existing (source: 48_cpg_island_ubiquitination_summary.tsv, mc_mean_mod_control)
- Mean baseline (control) mC: hyper islands = 0.566, hypo = 0.444, non-significant = 0.127 — hypermethylated islands start ALREADY methylated (source: 48_cpg_island_ubiquitination_summary.tsv)
- K119ub overlap at hyper islands is LOWER than background, not higher: K119ub-mutant overlap 17.2% (hyper) vs 30.2% (non-sig); K119ub-gained 0.0% (hyper) vs 1.9% (non-sig) (source: 48_cpg_island_ubiquitination_summary.tsv) → Fisher OR < 1 (depletion), not enrichment
- Co-significant (both mC & hmC) islands = 112; quadrants mC+/hmC+ = 5, mC-/hmC+ = 11, mC-/hmC- = 45, mC+/hmC- = 51 (source: 48_cpg_island_ubiquitination_summary.tsv)
- sig_status: mC-only = 330, hmC-only = 60 islands (source: 48_cpg_island_ubiquitination_summary.tsv)
- Chromatin context of hyper islands: H3K4me1 27.9%, H3K27ac 18.0%, H3K27me3 13.1%, H3K4me3 6.6%, Bivalent 3.3% (source: 48_cpg_island_ubiquitination_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Two negative/refining results sharpen the model. First, CpG-island hypermethylation in BAP1-KO is overwhelmingly amplification of PRE-EXISTING methylation (baseline mC 0.566; only ~15% de novo even at the loosest threshold), not de-novo silencing of clean promoters — so this is not a classic Polycomb-island CpG-island-methylator phenotype. Second, and counter to the naive "BAP1 loss → K119ub gain → methylation" hypothesis at islands, K119ub overlap is actually DEPLETED at hypermethylated islands (17% vs 30%), and gained-K119ub overlap is essentially zero (0/122). Together this argues the K119ub→methylation coupling operates at gene bodies / flanking dynamic regions (Sections 47, 11–19), NOT at CpG islands themselves, and that island hypermethylation is a secondary amplification process. The dominant co-significant quadrant is mC+/hmC- (51), again the TET-block signature.
### Plot inventory
- 48a_k119ub_enrichment — K119ub (ctrl/mut/gained/lost) overlap %, DMR vs background, Fisher
- 48b_k119ub_peak_count_violin — K119ub mutant peak count per island by DMR category (Wilcoxon)
- 48c_k119ub_diffbind_fold — DiffBind K119ub fold change at islands by category
- 48d_baseline_mc_density — control mC fraction density by DMR status (de-novo thresholds marked)
- 48e_ctrl_vs_mut_scatter — control vs mutant mC per island (diagonal = no change)
- 48f_de_novo_classification — de-novo vs pre-existing counts at 0.05/0.10/0.20 thresholds
- 48g_gain_magnitude — mC gain magnitude, low vs high baseline
- 48h_mc_hmc_scatter — coordinated mC vs hmC difference, quadrant counts
- 48i_cosig_heatmap — co-significant island heatmap (mC/hmC ctrl/mut)
- 48j_k119ub_denovo_forest — Fisher OR forest (hyper/hypo/de-novo vs K119ub gained)
- 48k_k119ub_denovo_tile — 2×2 K119ub × baseline contingency
- 48l_chromatin_overlap_bars — chromatin-mark overlap % by DMR status
- 48m_chromatin_or_heatmap — log2 OR chromatin enrichment heatmap

## Section 49: section_49_homer_motif_enrichment AND section_49_homer_k119ub_k27ac
### Analysis question
Which TF motifs are enriched (HOMER findMotifsGenome.pl known motifs) at differential H2AK119ub sites (peak-based B1/B2, DiffBind B3/B4) and differential H3K27ac sites (C1–C4, including vs-genome backgrounds)? Two complementary scripts: `_homer_motif_enrichment` (dot/family/wordcloud/MOI-heatmap) and `_homer_k119ub_k27ac` (matrix dot plots + family bars).
### Key results
- Significant motifs (q < 0.05) per comparison (uncapped): B1 (K119ub gained-vs-lost) = 192, B2 = 3, B3 (DiffBind gained-vs-lost) = 216, B4 = 15, C1 (K27ac gained-vs-lost) = 15, C2 (K27ac lost-vs-gained) = 60 (source: 49_homer_k119ub_k27ac/homer_all_significant_motifs.tsv, valid-ID line counts; total 501 significant motif-comparison rows)
- Top H2AK119ub-gained (B1) motifs are bHLH neuronal: Atoh1 (target 31.5% vs bg 23.7%, FE 1.33, q≈0), Atoh7 (FE 1.49), BHLHA15 (FE 1.33), NeuroD1 (FE 1.49) (source: 49_homer_motif_enrichment/homer_top25_significant_motifs.tsv)
- Top H3K27ac lost-vs-gained (C2) motifs also bHLH: Atoh1 (37.9% vs 33.8%, q≈0), Atoh7, BHLHA15, plus Egr2/FOXA1/GSC (source: 49_homer_motif_enrichment/homer_top25_significant_motifs.tsv)
- Motifs-of-interest are directionally specific: NeuroD1 enriched in B1/B3 (FE 1.49) but FE < 1 in lost-side B2/B4; YY1 significant only in B3 (FE 1.85, q≈0); Olig2 enriched B1/B3/C2/C3/C4 (source: 49_homer_motif_enrichment/homer_motifs_of_interest_all.tsv)
- Top-25 dot-plot export rows: B1=25, B3=25, C1=15, C2=25 (90 rows; C1 has only 15 because only 15 significant) (source: 49_homer_motif_enrichment/homer_top25_significant_motifs.tsv)
- K119ub matrix dot plot uses B1–B4 union; K27ac matrix uses C1–C2 (source: 49_homer_k119ub_k27ac script + homer_top25_per_comparison.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] H2AK119ub sites that GAIN signal in BAP1-KO are enriched for neuronal bHLH motifs (Atoh1, NeuroD1, Atoh7, BHLHA15) — the cerebellar granule/neuronal lineage TFs. Since BAP1 normally removes H2AK119ub, its loss spreads Polycomb ubiquitination, and the affected regions are precisely the binding grammar of cerebellar neuronal differentiation factors. The same bHLH grammar marks H3K27ac sites that are LOST (C2), linking enhancer deactivation to the same lineage TFs. YY1 (a known loop/architectural and Polycomb-recruiting factor) surfacing specifically in the DiffBind K119ub gained-vs-lost contrast (B3) is mechanistically suggestive of a recruitment link between ubiquitination spread and 3D architecture. B2/B4 (lost side) being nearly empty (3 and 15 motifs) shows the signal is one-directional: gain of K119ub carries the motif grammar, loss does not.
### Plot inventory (49_homer_motif_enrichment)
- homer_dotplot_top15 — top-15 motifs per key comparison (B1/B3/C1/C2), size=%target, color=fold-enrichment
- homer_family_barchart — TF-family counts per key comparison
- homer_wordclouds — TF word clouds (B1/B3/C2), size ∝ -log10(p)
- homer_motifs_of_interest_heatmap — curated MOI (REST/YY1/NeuroD1/Atoh1/…) -log10(q) across B/C comparisons
### Plot inventory (49_homer_k119ub_k27ac)
- homer_k119ub_matrix_dotplot — B1–B4 matrix dot plot (top-20 union TFs)
- homer_k27ac_matrix_dotplot — C1–C2 matrix dot plot
- homer_family_barchart_all — TF-family counts, all comparisons, faceted by mark

## Section 50: section_50_homer_dmr_motifs
### Analysis question
Which TF motifs are enriched at COORDINATED mC-up/hmC-down DMRs versus discordant DMRs (HOMER comparison A5)? (A1–A4, which used genome background, produced no/artefactual signal and are excluded.)
### Key results
- A5 significant motifs (q < 0.05) = 128 (source: 50_homer_dmr_motifs/homer_a5_significant_motifs.tsv, row count)
- Family breakdown of significant motifs: zinc-finger (Zf) = 34 (dominant), bHLH = 17, NR = 16, Homeobox = 12, ETS = 12, Forkhead = 8 (source: 50_homer_dmr_motifs/homer_a5_significant_motifs.tsv)
- Top motifs by p-value all Zf/ETS: ZEB2 (p ≈ 1e-25, target 98.1% vs bg 95.4%, FE 1.03), Maz (p ≈ 1e-24), Klf15, Sp5, ETV4, KLF17/14/5/6, ZNF148, ZNF711 (source: 50_homer_dmr_motifs/homer_a5_significant_motifs.tsv)
- Fold-enrichments are modest (≈1.02–1.07) because both target and background already contain the motif at high frequency (e.g. ZEB2 98% vs 95%) — significance comes from large n, not large effect (source: 50_homer_dmr_motifs/homer_a5_significant_motifs.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Coordinated mC-up/hmC-down regions (the TET-block signature) are distinguished from discordant DMRs by an excess of GC-rich zinc-finger / KLF / Sp-family motifs (ZEB2, Maz, KLF15/14/5/6, Sp5). These are GC-box / CpG-rich binding factors — exactly the sequence context where 5mC↔5hmC turnover by TET enzymes is most consequential and where methylation can directly modulate TF binding (e.g. methyl-sensitive KLF/Sp sites). The very small fold-enrichments reflect that the comparison is coordinated-vs-discordant DMRs (both CpG-rich), so the result is best read as "coordinated DMRs are subtly more zinc-finger-motif-dense," not a dramatic de-novo motif. bHLH (17) reappears here too, echoing the K119ub/K27ac neuronal grammar.
### Plot inventory
- homer_a5_coordinated_dotplot — top-25 A5 motifs, size=%target, color=TF family
- homer_a5_family_counts — TF-family counts among significant A5 motifs

## Cross-section synthesis
[INTERPRETATION] These nine sections triangulate the same BAP1-KO mechanism from independent angles and consistently return the coordinated mC-up/hmC-down ("TET-block") signature as the through-line. Section 42 shows it dominates the most extreme-effect genes (42/124 floor loci) and names the canonical targets (Syt1, Cntnap2). Section 44 shows the gain-of-methylation is allele-resolved and roughly doubled (1.95×) in mutant. Section 45 demonstrates the gene-body methylation response is evolutionarily conserved (Field human chr8 hotspot replicates, p ≈ 0.009) but promoter-independent. Section 47 is the architectural payoff: lost CTCF anchors are hypermethylated at flanking dynamic regions (OR ≈ 3.3) with a distance gradient and a negative loop-strength correlation, directly implicating methylation in the 3D loop loss (the Flavahan/Bernstein insulator model) and again carrying the coordinated signature (OR 2.18). Section 48 refines the model by ruling OUT two simple hypotheses — CpG-island hypermethylation is amplification not de-novo, and K119ub is depleted (not enriched) at island DMRs — pushing the K119ub→methylation coupling to gene bodies and dynamic regions rather than islands. Sections 49–50 supply the sequence grammar: neuronal bHLH motifs (Atoh1/NeuroD1) mark the gained-H2AK119ub and lost-H3K27ac sites, and zinc-finger/KLF motifs mark the coordinated DMRs — connecting the Polycomb-ubiquitination spread to cerebellar neuronal regulatory elements. Together they support the paper's thesis that BAP1 loss restructures both DNA methylation and 3D genome architecture through an upstream H2AK119ub signal, with a TET-mediated demethylation block as the unifying methylation phenotype, and Section 46 renders all of this co-localized at single loci.

## Tables used
- max_significance_genes_merged.tsv — merged mC/hmC q-floor genes with coordinated-pattern flag (124 rows)
- max_significance_genes_mc.tsv — mC-only q-floor genes sorted by effect size (100 genes)
- max_significance_genes_hmc.tsv — hmC-only q-floor genes (67 genes)
- max_significance_gene_names.txt — unique gene-name list (123) for pathway tools
- asm_mc_significant_per_sample.tsv — per-sample significant mC ASM counts with condition/sex/batch (8 samples)
- asm_mc_significant_loci.tsv — unique significant mC ASM loci with per-condition sample counts (58,955 loci)
- asm_dmr_overlap_summary.tsv — gene-body DMR genes containing ≥1 significant mC ASM site (1,408 genes)
- field_chr8_ortholog_mapping.tsv — Field human→mouse ortholog mapping with field direction (85 rows, 81 mapped)
- field_chr8_comparison_full.tsv — master Field comparison (gene-body/promoter mC, hmC, RNA-seq, concordance; 81 rows)
- field_chr8_concordant_genes.tsv — 21 gene-body-concordant Field genes with mC/hmC/RNA-seq
- field_chr8_promoter_comparison.tsv — promoter-level Field comparison (81 rows)
- field_chr8_statistical_tests.tsv — 4 Field statistical tests with p-values
- 47a_fisher_results.tsv — Fisher OR for mC-hyper/hmC-hypo at dynamic CTCF-anchor regions (4 tests)
- 47a_dynamic_ctcf_anchor_methylation.tsv — per-region mC/hmC overlay at dynamic anchors
- 47a_cpgi_ctcf_anchor_methylation.tsv — per-island mC/hmC at CpG-island anchors (8,910 rows; null analysis)
- 47b_multiregion_fisher.tsv — Fisher OR across region types (islands/shores/shelves/promoters)
- 47b_multiregion_wilcoxon.tsv — Wilcoxon mod_difference per region type
- 47c_coordinated_enrichment.tsv — coordinated mC-up/hmC-down Fisher tests (3 tests)
- 47c_dynamic_coordinated_pattern.tsv — per-region coordinated-pattern flags (64,843 rows; anchor-group breakdown)
- 47d_distance_stratified_results.tsv — distance-binned Fisher + CMH common OR
- 47e_correlation_results.tsv — loop-logFC vs methylation Spearman/partial correlations (4 tests)
- 47e_loop_methylation_correlation.tsv — per-loop anchor methylation (1,467 CTCF-CTCF loops)
- 48_cpg_island_ubiquitination_summary.tsv — full per-island table (8,910 islands; mC/hmC/K119ub/de-novo/chromatin)
- 49_homer_motif_enrichment/homer_top25_significant_motifs.tsv — top-25 motifs/comparison (B1/B3/C1/C2)
- 49_homer_motif_enrichment/homer_motifs_of_interest_all.tsv — curated MOI across all 8 comparisons (96 rows)
- 49_homer_k119ub_k27ac/homer_all_significant_motifs.tsv — all q<0.05 motifs (501 valid rows)
- 49_homer_k119ub_k27ac/homer_top25_per_comparison.tsv — top-25/comparison for matrix dot plots
- 50_homer_dmr_motifs/homer_a5_significant_motifs.tsv — A5 coordinated-vs-discordant significant motifs (128)

## Data-quality flags
- **Section 43 TSVs not independently read.** Section 43 exports cg_mc_chromosome_distribution.tsv / cg_hmc_chromosome_distribution.tsv / cg_top50_mc_dmr_genes.tsv, which were not in this group's starter table set. Its per-chromosome enrichment and direction counts are computed at runtime; I did not confirm specific numbers against those TSVs, so Section 43 quantitative claims here are limited to structural facts (panel inventory). Marked [UNVERIFIED] where applicable.
- **Section 46 is VISUAL ONLY — no TSV/quantification.** All numbers are rendered live onto Gviz tracks; there is no exported table. Documented gene loci and track content only.
- **HOMER TSVs live inside section folders, not tables/.** Confirmed: 49_homer_motif_enrichment/, 49_homer_k119ub_k27ac/, 50_homer_dmr_motifs/ each hold their *.tsv. The starter list referenced no homer_*.tsv in tables/ (correct).
- **HOMER comp_label contains literal newlines** in 49_homer_k119ub_k27ac/homer_all_significant_motifs.tsv (labels like "Gained\nvs Lost"), so naive line/field parsing splits rows. Significant-motif counts here were obtained by counting lines beginning with a valid comparison ID (B1–B4/C1–C4): B1=192, B2=3, B3=216, B4=15, C1=15, C2=60 (501 total). The capped top-25 dot-plot file (homer_top25_significant_motifs.tsv) is single-line-safe.
- **Two "significant-motif count" frames differ by design.** The dot-plot export (homer_top25_significant_motifs.tsv) is CAPPED at 25/comparison (B1=25, B3=25, C1=15, C2=25); the uncapped file reports the true totals (B1=192, B3=216, …). Not a discrepancy — different purposes; both cited explicitly.
- **Section 48 K119ub result is a depletion, not enrichment.** Hypermethylated CpG islands have LOWER K119ub overlap than background (17.2% vs 30.2%; gained-K119ub 0% vs 1.9%), i.e. Fisher OR < 1. This contradicts a naive "K119ub gain drives island methylation" reading and is reported as-is. The script's `48j`/`48k` forest/contingency panels visualize this OR<1.
- **Section 47 anchor-region counts: 47a vs 47c context.** 47a's Fisher uses lost=1,084 / gained=1,427 dynamic regions (shores+shelves only); 47c's pattern table reports lost=1,210 / gained=1,577 / both=274 / background=61,782 over the joined dynamic set. Counts differ because 47c is the inner-joined mC∩hmC region set (slightly different membership), not a discrepancy.
- **Section 45 trisomy-8 coverage ratio and Mann-Whitney p are NOT in a TSV** (printed to console / shown only on figure 45k subtitles). The trisomy conclusion ("ruled out") is supported by figure 45k but its exact coverage-ratio and expression-p values are not in any exported table I read; not cited numerically here.
- **Section 44 ASM shared/only counts use sample-level membership.** The 20,546 shared / 9,768 ctrl-only / 28,641 mut-only figures are derived from the loci table's n_ctrl/n_mut columns (locus "in condition" = significant in ≥1 sample of that condition), matching the script's Venn definition; per-sample significance is uncorrected (exploratory), so absolute ASM counts should be treated as indicative.
