# Sections 59-63: MeCP2 mechanistic model: log2 quadrant scatters, MeCP2-status epigenetic profiles, stoichiometry test, multifeature chromatin regression, master heatmap

## Summary
This group builds and tests the BAP1-KO MeCP2 mechanistic model by integrating H2AK119ub, MeCP2 CUT&RUN, CG 5mC/5hmC, and seven histone/ATAC marks at gene-body and peak resolution. The central questions are: (1) do MeCP2-Up loci show a coordinated cascade (K119ub up -> 5mC up, 5hmC down, K27ac down, K27me3 up) with a mirror-image at MeCP2-Down loci; (2) is the 5mC up / 5hmC down shift driven by direct DNMT3A dehydroxymethylase (stoichiometric, slope ~ -1, total methylation conserved) or by independent TET inhibition plus de novo DNMT3A; and (3) can chromatin features explain MeCP2 binding *beyond* CG methylation. The answers: MeCP2 status defines mirror-image epigenetic profiles with K119ub flat at MeCP2-Down loci; the genome-wide Delta5mC~Delta5hmC slope is -0.959 but total methylation is *not* conserved (it rises at MeCP2-Up loci), favoring TET inhibition over pure dehydroxymethylase; and chromatin marks explain MeCP2 binding far better than CG methylation alone (binding-level R2 0.246 chromatin-only vs 0.017 CG-only), with K119ub the dominant predictor.

## Section 59: quadrant_log2_scatter
### Analysis question
Build the master gene-level integration table and test pairwise log2-fold relationships between K119ub, MeCP2 binding, 5mC, and 5hmC via quadrant scatter plots with Spearman correlations.
### Key results
- Master table size = 23,150 genes (source: 59_quadrant_master.tsv)
- K119ub-quantifiable genes = 21,604 (source: 59_quadrant_master.tsv, gb_signal_class == "quantifiable")
- MeCP2-Up genes (mecp2_nearest_fdr < 0.05 & mecp2_mean_fold > 0) = 79 (source: 59_quadrant_master.tsv)
- MeCP2-Down genes (mecp2_nearest_fdr < 0.05 & mecp2_mean_fold < 0) = 34 (source: 59_quadrant_master.tsv)
- Genes carrying methylation data (mc_diff present) = full master join; MeCP2 fold present for the 113 significant genes plus continuous folds (source: 59_quadrant_master.tsv)
- 5mC vs 5hmC reciprocal structure is captured downstream by the stoichiometry slope (see section 61); the quadrant labels in 59a/59d/59e are derived directly from these columns.
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The master table operationalizes the mechanistic model: a positive K119ub vs MeCP2 association (59a top-right quadrant) supports H2AK119ub accumulation recruiting/stabilizing MeCP2, while the K119ub vs 5mC (59d, positive) and K119ub vs 5hmC (59e, negative) panels encode the predicted "ubiquitin -> methylation gain / hydroxymethylation loss" cascade. The euchromatin (K27ac, 59b) vs heterochromatin (K27me3, 59c) split tests whether MeCP2 redistribution tracks Polycomb gain rather than active-mark loss.
### Plot inventory
- 59a_k119ub_vs_mecp2 — gene-body K119ub log2FC vs concentration-weighted MeCP2 DiffBind fold, colored by MeCP2 status
- 59a2_k119ub_vs_mecp2_peaklevel — same comparison at MeCP2 peak resolution (every peak a point)
- 59b_mecp2_vs_k27ac / 59b_mecp2_euchromatin_k27ac — MeCP2 fold vs H3K27ac fold (euchromatin)
- 59c_mecp2_vs_k27me3 / 59c_mecp2_heterochromatin_k27me3 — MeCP2 fold vs H3K27me3 fold (heterochromatin)
- 59bc_chromatin_composite — euchromatin + heterochromatin side-by-side
- 59d_k119ub_vs_5mc — K119ub log2FC vs Delta5mC, colored by 5mC DMR status
- 59e_k119ub_vs_5hmc — K119ub log2FC vs Delta5hmC, colored by 5hmC DMR status
- 59f_mecp2_vs_5mc / 59f_composite — MeCP2 fold vs Delta5mC (plus composite spacer)
- 59g_mecp2_vs_5hmc — MeCP2 fold vs Delta5hmC
- 59h_composite — all 7 quadrant panels stacked

## Section 60: mecp2_status_epigenetic_profiles
### Analysis question
Formalize the mirror-image test: do MeCP2-Up and MeCP2-Down loci show opposite-direction epigenetic shifts across K119ub (BigWig + DiffBind), 5mC, 5hmC, K27ac, K27me3, with K119ub specifically *flat* at MeCP2-Down loci?
### Key results
- MeCP2-Up, K119ub BigWig: mean = +0.6298, wilcox_q = 1.77e-19 (***), n = 78 (source: 60_mecp2_status_stats.tsv)
- MeCP2-Down, K119ub BigWig: mean = -0.0391, wilcox_q = 0.182 (ns), n = 35 — K119ub is flat at MeCP2-Down loci (source: 60_mecp2_status_stats.tsv)
- MeCP2-Up: 5mC mean = +0.0491 (q = 9.2e-13, ***) and 5hmC mean = -0.0175 (q = 2.98e-08, ***), n = 78 (source: 60_mecp2_status_stats.tsv)
- MeCP2-Up: K27ac mean = -0.7181 (q = 1.23e-09, ***, n = 61) and K27me3 mean = +1.1101 (q = 0.018, *, n = 17) (source: 60_mecp2_status_stats.tsv)
- MeCP2-Down mirror: 5mC = -0.0278 (q = 7.2e-04), K27ac = +0.4814 (q = 1.58e-03), K27me3 = -0.3366 (q = 0.018) (source: 60_mecp2_status_stats.tsv)
- Up-vs-Down Mann-Whitney is significant for every mark, e.g. K119ub BigWig mw_updown_q = 5.42e-12, 5mC mw_updown_q = 5.42e-12, K27ac mw_updown_q = 1.49e-08 (source: 60_mecp2_status_stats.tsv)
- MeCP2-Down gene table = 35 genes exported (source: 60_mecp2_down_gene_table.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] MeCP2-Up loci show the full predicted cascade (K119ub gain drives 5mC up, 5hmC down, K27ac loss, K27me3 gain), consistent with Polycomb-domain compaction recruiting MeCP2. The MeCP2-Down loci show the mirror-image for every mark *except* K119ub, which stays flat — this asymmetry argues that MeCP2 loss is not simply the reverse of K119ub gain but reflects a separate redistribution (e.g., MeCP2 vacating active enhancers that gain accessibility). The flat-K119ub-at-Down signature is the cleanest single discriminator between the two mechanisms.
### Plot inventory
- 60a_violin_marks_by_mecp2 — 5-facet violin+box (K119ub, 5mC, 5hmC, K27ac, K27me3) by MeCP2 status with ggsignif brackets
- 60b_summary_heatmap — pheatmap of mean change (5 marks x Up/Down) with BH-Wilcoxon stars
- 60c_effect_size_lollipop — mean +/- 95% CI dot-and-whisker for all 6 marks
- 60d_down_genes_strip — per-gene strip chart for the 35 MeCP2-Down genes
- 60e_composite — patchwork of 60a + 60c + 60d

## Section 61: stoichiometry_mechanism
### Analysis question
Discriminate DNMT3A dehydroxymethylase (direct 5hmC->5mC, total methylation conserved, Delta5mC~Delta5hmC slope = -1) from TET inhibition + de novo DNMT3A (total methylation rises, slope != -1); test whether K27me3->ATAC accessibility mediates TET-efficiency change; define a neuronal gene set; and compare BAP1-KO vs TET-KO.
### Key results
- Global Delta5mC~Delta5hmC slope = -0.9589 [95% CI -0.9778, -0.9401], R2 = 0.3218, differs from -1 = TRUE (source: 61_stoichiometry_slopes.tsv)
- Total methylation (5mC+5hmC) change: All genes mean Delta = -0.00139 (q = 2.48e-112, ***), but MeCP2-Up mean Delta = +0.0316 (q = 3.19e-11, ***) and Coordinated mean Delta = +0.00697 (q = 2.12e-178, ***) — total is NOT conserved at the relevant loci (source: 61_stoichiometry_summary.tsv)
- Per-group slopes: MeCP2-Up = -0.894 [-1.291, -0.496] (CI includes -1), Coordinated = -1.288 [-1.323, -1.252], Neuronal = -0.721 [-0.805, -0.637] (source: 61_stoichiometry_slopes.tsv)
- Chromatin-accessibility model R2 for delta_ratio: ATAC only = 0.0949, ATAC+K27me3 = 0.1887, ATAC+K27me3+K119ub = 0.2290 (n = 3,434), +K27ac = 0.3033 (n = 1,379) (source: 61_delta_ratio_chromatin_models.tsv)
- Mediation K27me3->ATAC->delta_ratio: indirect effect (a*b) = -0.00279 [95% bootstrap CI -0.00337, -0.00226], proportion mediated = 18.0% (source: 61_mediation_results.tsv)
- BAP1-KO vs TET-KO delta_ratio Spearman: All genes rho = 0.220 (p = 1.18e-226, n = 20,851); Neuronal rho = 0.098 (p = 8.6e-04, n = 1,149); Coordinated rho = 0.176 (p = 1.9e-43, n = 6,069) (source: 61_tetko_comparison.tsv)
- Neuronal gene set = 1,149 genes; overlap with coordinated = 999, overlap with MeCP2-Up = 16 (source: 61_neuronal_gene_set.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The genome-wide slope of -0.959 superficially looks stoichiometric, but it is statistically distinguishable from -1 and, more decisively, total methylation *rises* at MeCP2-Up (+0.032) and coordinated (+0.007) loci rather than being conserved. That excess mass cannot come from pure 5hmC->5mC conversion, so the data favor TET inhibition plus continued de novo DNMT3A deposition. The Coordinated subset slope of -1.29 (steeper than -1) is consistent with added de novo 5mC on top of conversion. Mediation shows K27me3 gain partly closes chromatin (ATAC down), which accounts for ~18% of the TET-efficiency drop — chromatin compaction is a contributing, not dominant, route. The modest BAP1-KO/TET-KO correlation (rho 0.22 genome-wide) supports a shared TET-pathway component without full phenocopy.
### Plot inventory
- 61a_total_methylation_conservation — violin of Delta(5mC+5hmC) by gene group
- 61b_stoichiometry_scatter — Delta5mC vs Delta5hmC with OLS slope and reference -1 line
- 61c_chromatin_model_r2 — R2 bar for nested chromatin-accessibility models of delta_ratio
- 61d_mediation_panels — 3-panel mediation (K27me3->ATAC, ATAC->delta_ratio, bootstrap histogram)
- 61e_neuronal_venn — 3-way Venn (Neuronal x Coordinated x MeCP2-Up)
- 61f_bap1_vs_tetko_scatter — BAP1-KO vs TET-KO delta_ratio scatter
- 61g_composite — composite of 61a/61b/61c/61e

## Section 61h: mecp2_up_quadrant_enrichment
### Analysis question
Extract the top-right quadrant genes (MeCP2 significantly up AND K119ub log2FC > 0) from the 59a scatter and run GO BP / KEGG enrichment to test whether the consistent MeCP2-Up + K119ub-Up subset is neuronally enriched; repeat under relaxed (no FDR) and genome-wide-background settings.
### Key results
- Strict top-right quadrant (FDR < 0.05 + K119ub > 0) = 72 genes (source: 61h_mecp2_up_k119ub_up_genes.tsv)
- Strict + custom background: 1 significant GO BP term (q < 0.05) — "basement membrane assembly", q = 0.0209, ratio 3/72 (source: 61h_quadrant_go_bp.tsv)
- Strict + genome-wide background: 1 significant GO BP term (source: 61h_strict_genomewide_go_bp.tsv)
- Relaxed (fold > 0 + K119ub > 0) + custom background: 50 significant GO BP terms (source: 61h_relaxed_go_bp.tsv)
- Relaxed + genome-wide background: 1,117 significant GO BP terms; top terms are "regulation of synapse structure or activity" (q = 6.95e-13), "synapse assembly" (q = 1.88e-12), "axon guidance" (q = 3.06e-11); 13/25 top terms neuronal (source: 61h_relaxed_genomewide_go_bp.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The strict 72-gene set is too small and too tightly background-matched to yield neuronal enrichment (only a single, non-neuronal term survives) — this is a power/background artifact, not evidence against neuronal involvement. Relaxing the MeCP2 FDR cut and using a genome-wide background recovers a strong synaptic/axon-guidance signature, indicating that the MeCP2-Up + K119ub-Up program is broadly neuronal. The dependence of the neuronal signal on the background choice is the key caveat carried into 61i-61k.
### Plot inventory
- 61h_relaxed_go_dotplot — GO BP dotplot, relaxed quadrant, custom background
- 61h_relaxed_genomewide_dotplot — GO BP dotplot, relaxed quadrant, genome-wide background
- (no strict dotplots: strict analyses had <3 significant terms, below the plotting gate)

## Section 61i: quadrant_characterization
### Analysis question
Characterize the 72 strict quadrant genes by chromatin state and compare neuronal vs non-neuronal members across all marks.
### Key results
- Quadrant genes = 72; neuronal = 15 (21%), non-neuronal = 57 (79%) (source: 61i_quadrant_genes_full.tsv; 61i_neuronal_vs_nonneuronal_comparison.tsv)
- Chromatin-state composition: Active_Promoter = 35, Other = 30, Bivalent_Promoter = 4, Repressed_Promoter = 3 (source: 61i_quadrant_genes_full.tsv)
- Only K119ub log2FC differs significantly between neuronal and non-neuronal members: neuronal mean = 0.411 vs non-neuronal mean = 0.772, wilcox_q = 0.0111 (*) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)
- MeCP2 fold does not differ: neuronal mean = 0.862 vs non-neuronal mean = 0.684, wilcox_q = 0.497 (ns) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)
- Delta ratio trend: neuronal mean = -0.0525 vs non-neuronal -0.0315, wilcox_q = 0.083 (ns) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)
- 5hmC diff: neuronal mean = -0.0306 vs non-neuronal -0.0143, wilcox_q = 0.116 (ns) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Within the consistent quadrant, neuronal and non-neuronal genes are largely indistinguishable across MeCP2, 5mC, 5hmC, K27me3, K27ac, and ATAC — the only separation is that non-neuronal members carry somewhat higher K119ub gain. The dominance of Active_Promoter (35/72) and "Other" (30/72) states indicates the quadrant program operates mostly at active/unmarked promoters rather than canonical Polycomb domains, tempering a purely Polycomb-centric reading of the MeCP2-Up cascade at the gene-body level.
### Plot inventory
- 61i_chromatin_state_bar — chromatin-state bar chart of the 72 quadrant genes
- 61i_neuronal_vs_nonneuronal — 4-panel boxplots (MeCP2 fold, K119ub log2FC, 5mC, 5hmC) by gene class
- 61i_pergene_strip — per-gene multi-mark strip (6 marks)
- 61i_composite — state bar + neuronal-vs-non composite

## Section 61j: peak_level_enrichment
### Analysis question
Move from gene-level to peak-level: annotate the 7,686 MeCP2-Up DiffBind peaks (and the 6,849 MeCP2-Up + K119ub-Up peaks) to nearest genes and run GO BP enrichment under custom vs genome-wide backgrounds, with MeCP2-Down as a contrast.
### Key results
- MeCP2-Up peaks = 7,686; MeCP2-Down peaks = 1,200; MeCP2-Up + K119ub-Up peaks = 6,849 (source: 62_mecp2_peak_chromatin_signal.tsv, FDR/Fold/k119ub_log2fc columns)
- MeCP2-Up peaks -> 2,107 unique genes; MeCP2-Up + K119ub-Up -> 1,897 unique genes (source: 61j_mecp2up_peak_genes.tsv; 61j_mecp2up_k119up_peak_genes.tsv)
- MeCP2-Up, custom background: 129 significant GO BP terms (source: 61j_mecp2up_custom_go_bp.tsv)
- MeCP2-Up, genome-wide background: 818 significant GO BP terms; top = "synapse assembly" (q = 5.34e-22), "regulation of synapse organization" (q = 1.16e-21); 16/25 top terms neuronal (source: 61j_mecp2up_genomewide_go_bp.tsv)
- MeCP2-Up + K119ub-Up: 131 (custom bg) and 742 (genome-wide bg) significant GO BP terms (source: 61j_mecp2up_k119up_custom_go_bp.tsv; 61j_mecp2up_k119up_genomewide_go_bp.tsv)
- MeCP2-Down, genome-wide background: 517 significant GO BP terms; top = "regulation of synapse structure or activity" (q = 9.73e-11) (source: 61j_mecp2down_genomewide_go_bp.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] At peak resolution the neuronal/synaptic enrichment is robust and no longer dependent on relaxing thresholds (818 terms for MeCP2-Up genome-wide, top terms all synaptic). Both MeCP2-Up and MeCP2-Down peaks land near synaptic genes, so MeCP2 *redistribution* (not just gain) is concentrated at the neuronal regulatory landscape. Adding the K119ub-Up requirement barely changes the gene set (1,897 vs 2,107 genes) or its enrichment, indicating most MeCP2-Up peaks already coincide with K119ub gain — consistent with the ubiquitin-driven model.
### Plot inventory
- 61j_mecp2up_custom_dotplot / 61j_mecp2up_genomewide_dotplot — MeCP2-Up GO dotplots (custom / genome-wide bg)
- 61j_up_k119up_custom_dotplot / 61j_up_k119up_genomewide_dotplot — MeCP2-Up + K119ub-Up GO dotplots
- 61j_mecp2down_genomewide_dotplot — MeCP2-Down GO dotplot
- 61j_composite — combined bar-panel ORA figure (from section_61jk_composites.R)

## Section 61jk: composites (covered with 61j/61k)
### Analysis question
Render the multi-panel composite figures for the 61j peak-level ORA and the 61k GSEA results (no new statistics; reads the 61j/61k TSVs).
### Key results
- Composite reads the five 61j ORA tables and both 61k GSEA tables already documented above; no new numeric outputs are written by this script.
- 61j composite annotation = "7,686 MeCP2-Up peaks -> 2,107 genes" (matches 61j_mecp2up_peak_genes.tsv = 2,107 and 62 cache MeCP2-Up = 7,686).
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Purely presentational; consolidates the peak-level ORA and GSEA evidence into two figures for the paper.
### Plot inventory
- 61j_peak_level_enrichment/61j_composite — peak-level ORA bar composite (red = neuronal terms)
- 61k_gsea_mecp2_k119ub/61k_composite — GSEA NES bar composite (MeCP2 ranking | K119ub ranking)

## Section 61k: gsea_mecp2_k119ub
### Analysis question
Threshold-free test: rank ALL genes by continuous MeCP2 fold and by K119ub log2FC, then run GSEA (gseGO BP) to ask whether neuronal/synaptic gene sets concentrate among the largest MeCP2 and K119ub increases, with proper weighting.
### Key results
- MeCP2 ranking GSEA: 1 significant GO BP term (q < 0.05) — "synapse assembly", NES = +1.68, setSize = 299, p.adjust = 0.0245 (source: 61k_gsea_mecp2_go_bp.tsv)
- K119ub ranking GSEA: 115 significant GO BP terms (61 positive NES, 54 negative NES) (source: 61k_gsea_k119ub_go_bp.tsv)
- K119ub top positive-NES terms are RNA/tRNA metabolism: "tRNA processing" NES = +2.17 (q = 9.2e-06), "tRNA modification" NES = +2.15 (q = 8.4e-05), "RNA modification" NES = +2.15 (q = 6.1e-06) (source: 61k_gsea_k119ub_go_bp.tsv)
- K119ub neuronal/synaptic significant terms = 3: "translation at presynapse" NES = +1.94 (q = 0.0163), "translation at synapse" NES = +1.89 (q = 0.0257), "translation at postsynapse" NES = +1.89 (q = 0.0257) (source: 61k_gsea_k119ub_go_bp.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Threshold-free GSEA gives a *weaker* neuronal signal than the peak-level ORA: only synapse assembly survives for the MeCP2 ranking, and the K119ub ranking's strongest enrichments are RNA/tRNA-processing rather than synaptic. This means the strong synaptic ORA in 61j is partly an artifact of where MeCP2/K119ub peaks fall (gene-dense neuronal loci) rather than the genes with the very largest fold changes. The honest reading is that neuronal involvement is real at the binding-site level but should not be overstated as a fold-change-driven program — directly anticipating the section 78 correction for neuronal-gene-set bias.
### Plot inventory
- 61k_gsea_mecp2_dotplot — MeCP2-ranked GSEA dotplot (split by sign)
- 61k_gsea_k119ub_dotplot — K119ub-ranked GSEA dotplot
- 61k_mecp2_running_synapse_assembly — running-score plot for synapse assembly (MeCP2 ranking)
- 61k_k119ub_running_translation_at_presynapse — running-score plot for translation at presynapse (K119ub ranking)
- 61k_composite — side-by-side NES bar composite (from section_61jk_composites.R)

## Section 62: multifeature_chromatin_regression
### Analysis question
Core test of whether chromatin features explain MeCP2 binding better than CG methylation alone. Extract continuous BigWig signal for 7 histone/ATAC marks + CG 5mC/5hmC at all MeCP2 DiffBind peaks, then fit nested CG-only / Chromatin-only / Full regressions for both binding level (where MeCP2 binds) and differential (where MeCP2 changes), plus variance partitioning, LASSO, standardized coefficients, and chromatin-state stratification.
### Key results
- Binding-level R2: CG only = 0.0168, Chromatin only = 0.2456, Full = 0.2603 (n = 202,574 peaks) — chromatin >> CG (source: 62_model_comparison_summary.tsv)
- Differential R2: CG only = 0.0128, Chromatin only = 0.1704, Full = 0.1732 (n = 202,574) — MeCP2 redistribution tracks chromatin remodeling, not methylation (source: 62_model_comparison_summary.tsv)
- Variance partition (binding level): Chromatin-unique = 24.3%, CG-unique = 1.46%, Shared = 0.22%, Unexplained = 74.0% (source: 62_variance_partition.tsv)
- Variance partition (differential): Chromatin-unique = 16.0%, CG-unique = 0.28%, Shared = 1.00%, Unexplained = 82.7% (source: 62_variance_partition.tsv)
- Largest standardized betas (binding-level full model): K119ub = +0.199 (p < 2.2e-16), ATAC = +0.114, CG 5mC = +0.089, K27me3 = +0.061; CG 5hmC = -0.0177 (source: 62_standardized_coefficients.tsv)
- LASSO (lambda.1se, binding) retains K119ub (+0.187), K36me3 (+0.071), ATAC (+0.039), CG 5mC (+0.037), K27me3 (+0.017), K27me1 (+0.0003); drops CG 5hmC, K27ac, K4me3 (source: 62_lasso_selected_features.tsv)
- Chromatin-state-stratified R2 gain (Full - CG only) is largest at Active_Enhancer (+0.376, n = 7,795) and Poised_Enhancer (+0.210, n = 24,954); smallest at Repressed_Promoter (+0.008) and Polycomb (+0.027) (source: 62_chromatin_state_r2.tsv)
- CG-unexplained MeCP2 residual correlates most with K119ub change (rho = 0.209) and K27me3 change (rho = 0.127), both p ~ 0 (source: 62_residual_chromatin_regression.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the decisive quantitative answer to Jai's hypothesis: CG methylation alone explains <2% of MeCP2 binding variance, while chromatin marks explain ~25% (binding) and ~17% (differential). K119ub is by far the dominant predictor both as a standardized coefficient and as the strongest correlate of the CG-unexplained residual, directly supporting H2AK119ub (not CG methylation) as the upstream signal positioning MeCP2 in BAP1-KO. The enhancer-biased R2 gain shows the chromatin advantage is concentrated where MeCP2 is being redistributed, not at static Polycomb domains.
### Plot inventory
- 62a_model_r2_comparison — R2 bar for binding-level and differential CG/Chromatin/Full models
- 62b_variance_partition — stacked variance partition (CG-unique | Shared | Chromatin-unique | Unexplained)
- 62c_lasso_path — LASSO coefficient path (binding-level)
- 62d_standardized_coefficients — standardized-beta forest plot (full binding model)
- 62e_chromatin_state_r2 — CG-only vs Full R2 by chromatin state
- 62f_residual_vs_k27me3 — CG-unexplained MeCP2 residual vs H3K27me3 change scatter
- 62g_composite — composite of 62a/62b/62d/62f

## Section 63: mecp2_master_heatmap
### Analysis question
Cluster all significant MeCP2 peaks by multi-mark ctrl signal (Z-scored, ward.D2, k=4) to reveal whether MeCP2-Up peaks partition into distinct, e.g. Polycomb-enriched, chromatin profiles, annotated by MeCP2 direction and chromatin state.
### Key results
- Significant MeCP2 peaks total = 8,886 (7,686 Up + 1,200 Down, from the 62 cache); heatmap subsampled to 5,000 peaks (source: 62_mecp2_peak_chromatin_signal.tsv; 63_cluster_assignments.tsv)
- Subsampled set: 4,306 MeCP2-Up and 694 MeCP2-Down peaks (source: 63_cluster_assignments.tsv)
- 4 clusters: Cluster 1 = 1,874 peaks, Cluster 2 = 1,692, Cluster 3 = 1,064, Cluster 4 = 370 (source: 63_cluster_assignments.tsv)
- Every cluster's most common chromatin state is "Unmarked" (Cluster 1: 843 Unmarked / 570 Poised_Enhancer / 426 Active_Enhancer; Cluster 3: 711 Unmarked / 168 Polycomb / 118 Poised_Enhancer) (source: 63_cluster_assignments.tsv)
- Peak signal matrix exported = 5,000 peaks x 10 marks (MeCP2, H3K27me3, H2AK119ub, H3K27me1, H3K27ac, ATAC, H3K4me3, H3K36me3, CG 5mC, CG 5hmC) (source: 63_peak_signal_matrix.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] MeCP2 binding sites do not collapse into one homogeneous Polycomb block; they split into four chromatin-defined clusters, with enhancer-flavored clusters (1, with substantial Poised/Active_Enhancer) distinct from a more Polycomb-leaning cluster (3, the only one with appreciable Polycomb membership). The persistence of "Unmarked" as the modal state across clusters echoes 61i and indicates much of MeCP2 occupancy is at regions without strong canonical histone marks, where CG mC and K119ub carry the signal.

## Cross-section synthesis
Sections 59-63 assemble the paper's MeCP2 mechanistic spine. Section 59 builds the integrated gene-body table; section 60 shows MeCP2 status defines mirror-image epigenetic profiles with the diagnostic flat-K119ub at MeCP2-Down loci; section 61 uses the stoichiometry slope (-0.959, distinguishable from -1) together with non-conserved total methylation (rising at MeCP2-Up/coordinated loci) to favor TET inhibition + de novo DNMT3A over pure dehydroxymethylase, while showing K27me3->ATAC compaction mediates only ~18% of the TET-efficiency drop; sections 61h-61k probe whether the consistent MeCP2-Up + K119ub-Up program is neuronal, converging on a robust peak-level synaptic enrichment that GSEA shows is weaker when weighted by fold change (flagging neuronal-set bias); section 62 delivers the quantitative verdict that chromatin marks — K119ub foremost — explain MeCP2 binding an order of magnitude better than CG methylation; and section 63 shows MeCP2 sites are chromatin-heterogeneous rather than a single Polycomb block. Together they support the thesis that BAP1 loss drives H2AK119ub accumulation, which (not CG methylation per se) repositions MeCP2 and couples to a TET-block methylation shift across the neuronal regulatory genome.

## Tables used
- 59_quadrant_master.tsv — master gene-level integration (23,150 genes x 28 cols): K119ub gb_log2fc, MeCP2 fold/FDR, 5mC/5hmC, K27ac/K27me3/K119ub DiffBind folds, chromatin state
- 60_mecp2_status_stats.tsv — per-mark one-sample and Up/Down/NS two-sample test stats by MeCP2 status (18 rows)
- 60_mecp2_down_gene_table.tsv — the 35 MeCP2-Down genes with all-mark values
- 61_stoichiometry_summary.tsv — total-methylation conservation Wilcoxon by gene group
- 61_stoichiometry_slopes.tsv — Delta5mC~Delta5hmC OLS slopes (global + per group) vs reference -1
- 61_delta_ratio_chromatin_models.tsv — nested chromatin-accessibility model R2 for delta_ratio
- 61_mediation_results.tsv — K27me3->ATAC->delta_ratio mediation paths + bootstrap indirect effect
- 61_tetko_comparison.tsv — BAP1-KO vs TET-KO delta_ratio Spearman by group
- 61_neuronal_gene_set.tsv — 1,149 neuronal genes with coordinated/MeCP2-Up flags (NARROW DMR-enriched set; see flags)
- 61h_mecp2_up_k119ub_up_genes.tsv — 72 strict top-right quadrant genes
- 61h_quadrant_go_bp.tsv / 61h_strict_genomewide_go_bp.tsv — strict GO (1 sig term each)
- 61h_relaxed_go_bp.tsv / 61h_relaxed_genomewide_go_bp.tsv — relaxed GO (50 / 1,117 sig terms)
- 61i_quadrant_genes_full.tsv — 72 quadrant genes with all marks + chromatin state + neuronal flag
- 61i_neuronal_vs_nonneuronal_comparison.tsv — neuronal vs non-neuronal Wilcoxon per metric
- 61j_mecp2up_peak_genes.tsv (2,107) / 61j_mecp2up_k119up_peak_genes.tsv (1,897) — peak-annotated gene lists
- 61j_mecp2up_*_go_bp.tsv / 61j_mecp2up_k119up_*_go_bp.tsv / 61j_mecp2down_genomewide_go_bp.tsv — peak-level ORA tables
- 61k_gsea_mecp2_go_bp.tsv (1 sig) / 61k_gsea_k119ub_go_bp.tsv (115 sig) — threshold-free GSEA
- 62_model_comparison_summary.tsv — binding/differential CG/Chromatin/Full R2 (n = 202,574)
- 62_variance_partition.tsv — CG-unique / Shared / Chromatin-unique / Unexplained
- 62_standardized_coefficients.tsv — standardized betas, full binding model
- 62_lasso_selected_features.tsv — LASSO coefficients (binding 1se/min, differential 1se)
- 62_chromatin_state_r2.tsv — per-state CG-only vs Full R2 and gain
- 62_residual_chromatin_regression.tsv — per-mark Spearman of CG-unexplained residual
- 62_mecp2_peak_chromatin_signal.tsv — 202,650-peak signal cache (peak counts, log2FCs)
- 63_cluster_assignments.tsv — 5,000 subsampled peaks with cluster, MeCP2 class, chromatin state
- 63_peak_signal_matrix.tsv — Z-scored 5,000 x 10 mark matrix for the heatmap

## Data-quality flags
- NARROW neuronal gene set bias (group note): section 61's neuronal set (61_neuronal_gene_set.tsv, 1,149 genes) is extracted from GO terms matching "synap|neuron|axon" in an enrichment_go_bp.tsv that is itself DMR-enriched; section 78 later shows this set's directionality is biased. The neuronal stoichiometry slope (-0.721) and neuronal/coordinated overlaps (999 of 1,149) inherit this bias and should not be read as unbiased estimates.
- Background-dependent neuronal enrichment: 61h strict analyses yield only 1 GO term (and that one non-neuronal, "basement membrane assembly"), whereas relaxed + genome-wide background yields 1,117 terms dominated by synaptic categories. The neuronal claim is real at peak level (61j) and weaker under threshold-free GSEA (61k, MeCP2 ranking = 1 term). State the background and weighting when citing neuronal enrichment.
- GSEA vs ORA discordance: 61k K119ub ranking's strongest enrichments are RNA/tRNA-processing terms (NES +2.17 etc.), not synaptic, and the MeCP2 ranking yields a single significant term. This is a genuine source-vs-source nuance versus 61j's strong ORA neuronal signal, not an error — it reflects ORA peak-location bias vs GSEA fold-weighting.
- MeCP2-Up count: section 59/60 classification gives 79 MeCP2-Up / 34 MeCP2-Down genes at gene level (from 59_quadrant_master.tsv FDR<0.05), but per-mark test n in 60_mecp2_status_stats.tsv is 78 Up / 35 Down (K119ub BigWig row) because the gb_signal_class=="quantifiable" and per-mark non-NA filters drop/shift one gene each. Not a discrepancy in the data, but the headline "n" differs by mark.
- Peak-level K119ub-Up count: 61j composite cites 6,849 MeCP2-Up+K119ub-Up peaks consistent with the 62 cache; the gene list collapses to 1,897 genes. Peak vs gene counts must not be conflated.
- FIGURES.md does not yet contain entries for sections 59-63 (its "quadrant" references are sections 20-21); no prose-vs-table conflict exists to reconcile, but FIGURES.md is incomplete for this group.
- Section 63 heatmap is subsampled (5,000 of 8,886 significant peaks, seed=42); cluster sizes and per-cluster state counts are from the subsample, not the full set. The 63a/63d pheatmaps are visual-only (no per-cluster numeric table beyond cluster assignments).
- All numbers reflect run-5 (deep-seq, 8 samples, sex covariate) per project primary-run status; the 62 BigWig signal cache and 61/63 derivatives depend on the external BigWigs in /Users/zakiralibhai/sdsc/bigwigs/.
