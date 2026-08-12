# Sections 32-37: CHG context exploratory, multi-mark DiffBind quantitative integration, genomic-interval & gene-level permutation validation

## Summary
This group answers two questions. (1) **Is there any signal in the non-CG (CHG) context, and does it agree with the primary CG findings?** Section 32 calls 320 significant gene-body CHG 5mC DMRs but they sit at effect sizes ~100x smaller than CG (permille vs percent), are concentrated on chr7+chr8 (76%), and are direction-**discordant** with CG mC in 61% of overlapping genes — i.e. CHG is largely independent low-level noise, not a clean corroboration of CG. (2) **Does the BAP1 mechanistic cascade (K119ub gain → K27me3 gain → accessibility/K27ac loss → methylation change) hold quantitatively, and does it survive a null-model that accounts for genomic non-independence?** Section 33 upgrades binary peak overlap to quantitative DiffBind log2FC integration across 4 marks and shows H2AK119ub is by far the strongest predictor of mC hypermethylation (OR≈4.7 per unit log2FC; 4-mark logistic AUC=0.818). Sections 34-37 are confirmatory permutation tests: the gene-level label-shuffle test (section 37) cleanly validates **all 15** prior Fisher/O-E enrichments, while the interval-based regioneReloaded tests (34-36) completed the permutation computation and heatmaps but crashed before writing their Fisher-vs-permutation comparison tables.

## Section 32: section_32_chg_exploratory_analysis
### Analysis question
Characterize significant CHG (non-CG) 5mC and 5hmC DMRs from run-5, and ask whether CHG changes recapitulate, extend, or contradict the primary CG-context signal.
### Key results
- CHG mC significant DMRs = 320 genes (source: chg_mc_significant_dmrs.tsv)
- CHG mC direction split = 157 hypermethylated / 163 hypomethylated (source: chg_mc_significant_dmrs.tsv)
- CHG mC genes also significant in CG mC = 206 / 320 = 64.4% overlap (source: chg_mc_significant_dmrs.tsv, `cg_mc_significant` column)
- Direction discordance with CG mC = 126 / 206 overlapping genes discordant = 61.2% (source: chg_mc_significant_dmrs.tsv, `concordant_with_cg_mc` = FALSE among `cg_mc_significant`=TRUE)
- CHG-unique genes (significant CHG mC, not CG mC) = 114, of which 55 also significant in CG hmC (source: chg_unique_genes.tsv)
- chr7+chr8 concentration = 244 / 320 = 76.3% (chr8=149, chr7=95) (source: chg_mc_significant_dmrs.tsv, `chr` column)
- CHG hmC significant DMRs = 92 genes (56 hyper, 36 hypo) (source: chg_hmc_significant_dmrs.tsv)
- Top CHG mC genes by q-value: Rbpms (chr8, q=1.05e-305), Dlc1 / Nrg1 / Sgcz (chr8, q=1.21e-304), Mtus1 (chr8, q=8.42e-268) (source: chg_mc_significant_dmrs.tsv)
- KEGG enrichment: only 2 pathways at q<0.05 — "Folate transport and metabolism" (q=0.0483, count=4, FoldEnrichment=12.7) and "Mucin type O-glycan biosynthesis" (q=0.0483, count=4, FE=10.7) (source: chg_enrichment_kegg.tsv)
- GO BP enrichment: 0 terms at q<0.05, 27 terms at q<0.1; top hits all glycosylation/glycoprotein processes (q≈0.053, e.g. glycoprotein metabolic process count=14) (source: chg_enrichment_go_bp.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The CHG "signal" is best read as low-amplitude, mostly CG-independent background rather than a second corroborating epigenetic layer. Three observations support this: effect sizes are ~100x smaller than CG (the script plots CHG in permille, CG in percent); 61% of the genes that are significant in both contexts move in *opposite* directions; and the hits pile up on two chromosomes (chr7+chr8), the hallmark of a regional/structural artifact or a small number of large gene-dense domains dominating the multiple-testing landscape. The recurrence of glycosylation/O-glycan and folate-metabolism terms (folate is the methyl-donor pathway) is biologically suggestive but rests on tiny gene counts (4 per KEGG pathway) and does not survive a strict q<0.05 GO threshold, so it should be treated as hypothesis-generating only. The honest conclusion is that BAP1 loss does not produce a coherent, CG-concordant CHG remethylation program.
### Plot inventory
- 32a_chg_volcano_mc — volcano of gene-body CHG 5mC DMRs (permille effect axis), top 15 genes labelled
- 32b_chg_volcano_hmc — volcano of gene-body CHG 5hmC DMRs
- 32c_chg_direction_breakdown — bar chart of hyper vs hypo counts for CHG mC and hmC
- 32d_chg_cg_venn_mc — 2-set Venn: CG mC vs CHG mC significant gene overlap
- 32e_chg_cg_venn_three_way — 3-set Venn: CG mC, CG hmC, CHG mC
- 32f_chg_cg_concordance_scatter — CG mC (%) vs CHG mC (permille) per overlapping gene, coloured by concordance
- 32g_chg_cg_concordance_bar — concordant/discordant summary bars for mC and hmC
- 32h_chg_top_dmr_genes — all 320 CHG mC genes ranked by q-value with CG-overlap shape annotation
- 32i_chg_unique_genes — lollipop of the 114 CHG-unique genes
- 32j_chg_chromosome_distribution — per-chromosome significant-DMR counts vs genome-wide expectation, with per-chr Fisher p
- 32k_chg_vs_cg_chr_comparison — normalized chromosome distribution, CHG vs CG
- 32l_chg_enrichment_go — GO BP dotplot (exploratory q<0.1)
- 32m_chg_enrichment_kegg — KEGG pathway dotplot (exploratory)

## Section 33: section_33_multi_mark_diffbind_integration
### Analysis question
Replace binary peak-overlap integration with quantitative DiffBind fold-change integration across 4 chromatin marks (ATAC, H3K27ac, H3K27me3, H2AK119ub), and test whether marks predict and converge on BAP1-driven methylation change in the cascade-predicted directions.
### Key results
- Per-mark DiffBind peaks (sig up / sig down at FDR<0.05): ATAC 75,867 peaks (9,263 up / 4,159 down); H3K27ac 25,669 (5,077 / 6,706); H3K27me3 18,324 (2,293 / 4,811); H2AK119ub 41,392 (16,715 / 5,097) (source: run.txt log §Section 33; per-peak data underlies diffbind_gene_level_all_marks.tsv)
- Multi-mark gene profile = 20,915 genes; complete cases across all 4 marks = 1,379 (6.6%) (source: diffbind_gene_level_all_marks.tsv; coverage % from run.txt log)
- 4-mark logistic regression for mC hypermethylation: H2AK119ub OR=4.707 [3.326, 6.663] p=2.4e-18; H3K27me3 OR=1.443 [1.188, 1.754] p=2.2e-04; H3K27ac OR=0.480 [0.363, 0.636] p=2.8e-07; ATAC-seq OR=0.221 [0.125, 0.391] p=2.2e-07 (source: diffbind_logistic_model_coefficients.tsv)
- Model performance: 4-mark AUC=0.818 [0.794, 0.841], McFadden R2=0.235 (N=1,217); extended 4-mark+baseline mC/hmC AUC=0.869 [0.850, 0.889], McFadden R2=0.348 (source: run.txt log §Section 33, lines 5202-5204 — not in a TSV)
- Cross-mark Spearman: H2AK119ub fold vs mC difference rho=+0.410; H2AK119ub fold vs hmC difference rho=-0.300; H3K27me3 vs H3K27ac rho=-0.497; mC vs hmC difference rho=-0.685 (source: diffbind_cross_mark_correlations.tsv)
- Quantitative O/E, cascade-confirming quadrants: mC Up + K119ub Gained O/E=1.21 (obs 1779, Fisher OR=14.75, p=7.5e-175); mC Down + K119ub Lost O/E=2.38; mC Up + K27me3 Gained O/E=1.64 (Fisher OR=14.12, p=8.6e-69); mC Up + K27ac Lost O/E=1.26 and mC Down + K27ac Gained O/E=2.17 (Fisher OR=0.062, p=4.1e-112) (source: diffbind_quantitative_oe_comparison.tsv)
- Convergence: 178 genes have 3+ marks changing in the predicted cascade direction (158 with exactly 3, 20 with all 4); of these 128 are mC-hyper DMRs, 14 mC-hypo, 36 non-DMR (source: diffbind_convergent_genes_3plus.tsv)
- Convergence distribution across all 20,915 genes: 16,828 with 0 concordant marks, 3,094 with 1, 815 with 2, 158 with 3, 20 with 4 (source: diffbind_convergence_per_gene.tsv, `n_concordant`)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the quantitative backbone of the paper's mechanism. The logistic model ranks the four marks exactly as the BAP1 model predicts: H2AK119ub gain is the dominant positive predictor of mC hypermethylation (OR≈4.7, an order of magnitude above K27me3), while open-chromatin/active marks (ATAC, K27ac) are protective (OR<1). Because BAP1 is the H2AK119ub deubiquitinase, its loss raising K119ub being the top predictor of remethylation places K119ub upstream in the causal chain rather than as a passenger. The cross-mark correlations independently reproduce the cascade wiring (K119ub↑ tracks mC↑ and hmC↓; K27me3 opposes K27ac). The O/E quadrants show the cascade is strongest at K119ub and K27me3 (Fisher OR 14-15) and the convergence count (178 genes with ≥3 concordant marks, 72% of them hyper-DMRs) identifies the highest-confidence BAP1-responsive loci where multiple epigenetic layers move together. The AUC jump from 0.818 (marks only) to 0.869 (adding baseline mC/hmC) indicates baseline methylation state adds independent predictive information beyond the dynamic mark changes.
### Plot inventory
- 33a_diffbind_volcano_plots — 2x2 grid of per-mark DiffBind volcanoes (ATAC, K27ac, K27me3, K119ub)
- 33b_cross_mark_correlation_heatmap — 7x7 Spearman correlation heatmap (4 marks + mC + hmC + delta-ratio), clustered
- 33c_quantitative_oe_dotplot — O/E enrichment per mark for mC and hmC perspectives, faceted by mark
- 33d_methylation_vs_mark_scatters — 4x2 grid of methylation difference vs mark log2FC with loess fits
- 33e_logistic_regression_forest — forest plot of 4-mark logistic ORs predicting mC hypermethylation
- 33f_convergence_analysis — cascade-concordance-by-DMR-status stacked bar + mark-combination intersection bars

## Section 34: section_34_permutation_dmr_chromatin_marks
### Analysis question
Confirmatory. Use regioneReloaded `crosswisePermTest` (interval-shuffling, per-chromosome) to test whether the spatial co-occurrence of mC Hyper/Hypo DMR intervals with 8 chromatin-mark peak sets survives a genomic null, validating the Fisher tests from sections 12a, 14a, 14b, 19f.
### Key results
- Design: 2 DMR direction sets (mC Hyper = 7,522 regions; mC Hypo = 3,267) x 8 mark sets (ATAC Up/Down, K119ub Ctrl/Mut/Gained/Lost, H3K27ac Ctrl/Mut) = 16 pairwise permutation tests at ntimes=100,000 (source: biomodal_perm_47851857.out)
- Permutation completed: crosswise z-score range across the 16 cells = [-8.82, 65.14] (source: biomodal_perm_47851857.out, line 112)
- Recomputed parent Fisher tests (inline): 12a ATAC OR=0.07 p=6.92e-179; 14a K119ub OR=0.70 p=3.93e-16; 14b K119ub-differential OR=4.45 p=4.70e-97; 19f H3K27ac OR=1.36 p=1.32e-05 (source: biomodal_perm_47851857.out, lines 128-131)
- Peak set sizes after chr filter: K119ub Gained (mut-only) 6,172; K119ub Lost (ctrl-only) 6,096; ATAC Up 7,618; ATAC Down 3,744; H3K27ac Mut 15,657; H3K27ac Ctrl 13,629 (source: biomodal_perm_47851857.out, lines 90-97)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The permutation computation itself succeeded and the extreme positive z-scores (up to 65) confirm that DMR intervals co-localize with chromatin-mark peaks far more than chromosome-matched random intervals would — the spatial associations underlying the K119ub/ATAC/K27ac Fisher tests are not artifacts of genomic clustering. However, this section never produced its final concordance table (see data-quality flags), so the formal "Confirmed/Weakened" labels for the four parent tests live only in the gene-level analogue (section 37, tests 37-01/37-04/37-05), where they all pass.
### Plot inventory
- 34a_crosswise_dmr_x_marks — 2x8 crosswise z-score association heatmap (DMR direction x marks). Only this panel was produced; 34b (Fisher-vs-permutation forest) and 34c (local z-score curve) were not generated because the script crashed at the comparison-table step.

## Section 35: section_35_permutation_atac_loops
### Analysis question
Confirmatory. Three regioneReloaded sub-analyses validating sections 13b/13c/13d/27/31: (35A) ATAC peaks x 6 ChIP marks, (35B) ATAC peaks x 7 chromatin states, (35C) Hi-C loop anchors x 6 chromatin features.
### Key results
- Total design = 38 pairwise permutation tests at ntimes=100,000 (35A: 2x6=12; 35B: 2x7=14; 35C: 2x6=12) (source: biomodal_perm_47851858.out, lines 125-129)
- Loop input: 2,910 loops → 2,327 unique gained-loop anchors (from 1,723 up_in_mutant loops) and 1,520 lost-loop anchors (from 1,187 down_in_mutant loops) (source: biomodal_perm_47851858.out, lines 95-97)
- ChIP peak sets: CTCF 32,487; H3K27ac 15,105; H3K27me3 15,809; H3K4me1 113,781; H3K4me3 6,581; Bivalent 318 (source: biomodal_perm_47851858.out, lines 79-84)
- MeCP2 differential peaks (35C feature set): MeCP2 Up 7,686; MeCP2 Down 1,200 (source: biomodal_perm_47851858.out, lines 87-88)
- Chromatin-state region counts feeding 35B: Active_Promoter 6,399; Repressed_Promoter 3,246; Other 10,994; Bivalent_Promoter 152; Polycomb 57; Active_Enhancer 50; Poised_Enhancer 93 (source: biomodal_perm_47851858.out, lines 115-121)
- All three crosswisePermTest objects completed and were cached to permutation_35_atac_loops.rds (source: biomodal_perm_47851858.out, line 162)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The asymmetric anchor counts (2,327 gained vs 1,520 lost anchors; 1,723 vs 1,187 loops) restate the central architectural finding that BAP1 loss gains more loops than it loses, consistent with Polycomb-domain compaction. The permutation framework was set up to ask whether gained vs lost anchors differentially co-localize with accessibility loss, MeCP2 redistribution, and DMRs — the spatial computations ran to completion, but as with section 34 the concordance table was never emitted, so the quantitative anchor-vs-feature confirmation must be read from the gene-level equivalents in section 37 (tests 37-09, 37-11) which both pass.
### Plot inventory
- 35a_crosswise_atac_x_chip — 2x6 heatmap: ATAC direction x 6 ChIP marks
- 35b_crosswise_atac_x_chromatin_state — 2x7 heatmap: ATAC direction x chromatin states
- 35c_crosswise_loops_x_features — 2x6 heatmap: loop-anchor direction x chromatin features. Panels 35d (Fisher-vs-permutation table) and 35e (local z-score) were not generated — the script crashed at `extract_perm_results` before figure 35d.

## Section 36: section_36_permutation_domains
### Analysis question
Confirmatory. Domain-scale regioneReloaded tests validating sections 29 and 30: (36A) 4 DMR direction sets x A/B compartments + compartment shifts; (36B) 4 DMR direction sets x Polycomb domains (H3K27me3, Bivalent).
### Key results
- Design = 24 pairwise tests at ntimes=100,000 (36A: 4x4=16; 36B: 4x2=8) (source: biomodal_perm_47851859.out, lines 97-100)
- DMR direction sets: mC Hyper 7,522; mC Hypo 3,267; hmC Hypo 9,532; hmC Hyper 1,967 regions (source: biomodal_perm_47851859.out, lines 75-78)
- Compartment bins (mm10 standard chr): A compartment 47,241; B compartment 54,480; B→A shift 5,430; A→B shift 2,641 (from 101,721 filtered bins) (source: biomodal_perm_47851859.out, lines 85-89)
- Polycomb domains: H3K27me3 15,809 peaks; Bivalent 318 peaks (source: biomodal_perm_47851859.out, lines 92-93)
- Both crosswise objects (36A, 36B) plus the 36c local z-score computed and cached; figures 36a, 36b, 36c all saved (source: biomodal_perm_47851859.out, lines 116-163)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This section reached further than 34/35 — it produced the compartment heatmap (36a), Polycomb heatmap (36b), AND the local z-score curve (36c) for mC-hyper DMRs at the A compartment, only failing at the final Fisher-comparison dot-plot. The large, well-powered B→A (5,430) and A→B (2,641) shift sets give the permutation test real power, and the companion gene-level tests (section 37, tests 37-10a/37-10b) show the compartment-shift × mC-direction relationship is Confirmed: B→A shifts are strongly anti-associated with hypermethylation (37-10a perm z=-25.8) while A→B shifts associate with it (37-10b perm z=+5.1). [INTERPRETATION] Mechanistically this couples 3D-compartment remodeling to the methylation change — regions that open into the active A compartment lose the repressive remethylation signal, and vice versa.
### Plot inventory
- 36a_crosswise_dmr_x_compartment — 4x4 heatmap: 4 DMR direction sets x A/B compartment + shifts
- 36b_crosswise_dmr_x_polycomb — 4x2 heatmap: DMR direction x H3K27me3 / Bivalent
- 36c_local_zscore_compartment — local z-score curve, mC Hyper DMRs at A compartment, ±50kb. Panel 36d (Fisher-vs-permutation) was not generated — crashed at `extract_perm_results`.

## Section 37: section_37_permutation_gene_level
### Analysis question
Confirmatory and the load-bearing permutation section. Chromosome-stratified label-shuffle permutation of 15 gene-level 2x2 contingency tables (the unit is a gene with two categorical labels, not an interval), generating an empirical null OR distribution to give honest p-values for the Fisher/O-E enrichments in sections 12e, 15a-c, 19g, 20d, 27c, 29, 31b, and 33c.
### Key results
- 15 / 15 gene-level tests Confirmed; 0 Weakened, 0 Strengthened, 0 Concordant-NS (source: permutation_37_summary.tsv, `concordance` column)
- Strongest effects (|perm z|): mC direction x expression direction perm z=-31.93 (test 37-08, Fisher OR=0.012, p=2.8e-162, n=1,177); mC direction x ATAC direction perm z=-30.12 (37-01, OR=0.095, n=2,958); hmC direction x ATAC direction perm z=+28.89 (37-03, OR=10.48, n=2,973); mC direction x K119ub DiffBind perm z=+28.75 (37-07d, OR=14.75, n=2,857) (source: permutation_37_summary.tsv)
- DiffBind cascade tests (validating section 33c) all Confirmed: 37-07a mC x ATAC OR=0.059 z=-14.67; 37-07b mC x K27ac OR=0.062 z=-23.96; 37-07c mC x K27me3 OR=14.12 z=+19.33; 37-07d mC x K119ub OR=14.75 z=+28.75 (source: permutation_37_summary.tsv)
- Compartment-shift tests (validating section 29): 37-10a B→A shift x mC direction OR=0.030 z=-25.77; 37-10b A→B shift x mC direction OR=2.73 z=+5.09 (source: permutation_37_summary.tsv)
- Loop/MeCP2 anchor tests: 37-09 K119ub direction x Hyper-at-loop-anchors OR=5.93 z=+12.45 (validates 27c); 37-11 Loop direction x MeCP2 Sig-Up at anchors OR=0.336 z=-4.30 (validates 31b) (source: permutation_37_summary.tsv)
- Weakest-but-still-significant effect: hmC direction x MeCP2 direction (37-02, validates 15a) OR=1.27, perm z=+4.19, empirical p=2.0e-04 (source: permutation_37_summary.tsv)
- Empirical p-values: 13/15 tests at the floor p≈9.999e-05 (= 1/(ntimes+1) at 100,000 permutations); the two weakest (37-02, 37-11) at p=2.0e-04 and 8.0e-04 (source: permutation_37_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the validation that actually closes the loop for the paper. Every key cross-tabulation that the discordance/cascade/compartment story rests on — mC↔expression anticorrelation, mC↔ATAC anticorrelation, hmC↔ATAC positive correlation, the four DiffBind cascade directions, the compartment-shift coupling, and the loop-anchor/MeCP2 relationships — survives a null model that preserves per-chromosome marginal counts and genomic neighborhood structure. Because label-shuffling within chromosomes removes any inflation from gene-density clustering, the uniformly "Confirmed" verdict means the associations are driven by genuine co-assignment of labels at the same genes, not by chromosomal co-residence. The empirical p-values hitting the 1/100,001 resolution floor for 13 of 15 tests indicates the observed odds ratios fall entirely outside 100,000 permuted nulls. This gene-level test is also the only one of the four permutation sections that ran to completion, so it carries the formal confirmatory weight that sections 34-36 were designed to provide at the interval level.
### Plot inventory
- 37a_zscore_forest_plot — forest plot of all 15 permutation z-scores, coloured by source section
- 37b_null_distribution_top4 — 2x2 null-distribution histograms for the 4 strongest effects with observed OR marked
- 37c_fisher_vs_permutation_table — formatted Fisher-OR vs permutation-z comparison table
- 37d_log2or_vs_zscore_scatter — scatter of observed log2(OR) vs permutation z-score with trend line

## Cross-section synthesis
The group moves from "is there a competing signal?" to "is the main signal robust?" Section 32 rules out CHG as a meaningful parallel mechanism — its 320 DMRs are tiny, chr7/chr8-skewed, and 61% direction-discordant with CG, so the paper's CG-centric narrative is not undercut by non-CG methylation. Section 33 then builds the quantitative spine of that CG narrative: across 20,915 genes, raised H2AK119ub is the single strongest quantitative predictor of mC hypermethylation (OR≈4.7, logistic AUC=0.818), with K27me3, K27ac and ATAC arranged exactly as the BAP1 cascade predicts, and 178 genes showing ≥3 marks converging in-direction. Sections 34-37 stress-test those enrichments against genomic non-independence: the gene-level label-shuffle (section 37) confirms all 15 underlying Fisher/O-E tests — including the four DiffBind cascade directions from section 33c and the compartment-shift coupling — while the interval-level regioneReloaded tests (34-36) completed their permutation computations and heatmaps (with z-scores up to 65 in section 34) but did not emit their comparison tables. Together the group delivers the paper's core claim with quantitative effect sizes and a passing null model: [INTERPRETATION] BAP1 loss raises H2AK119ub, which is the upstream driver of coordinated chromatin and DNA-methylation remodeling.

## Tables used
- chg_mc_significant_dmrs.tsv — 320 significant CHG 5mC gene-body DMRs with CG mC/hmC cross-reference and concordance flag
- chg_hmc_significant_dmrs.tsv — 92 significant CHG 5hmC gene-body DMRs
- chg_unique_genes.tsv — 114 genes significant in CHG mC but not CG mC, with CG hmC overlap flag
- chg_enrichment_go_bp.tsv — GO Biological Process enrichment of the 320 CHG mC genes (exploratory, 0 at q<0.05)
- chg_enrichment_kegg.tsv — KEGG pathway enrichment of the 320 CHG mC genes (2 at q<0.05)
- diffbind_gene_level_all_marks.tsv — 20,915-gene profile: methylation + delta-ratio + chromatin state + 4-mark DiffBind fold/FDR/n_peaks
- diffbind_cross_mark_correlations.tsv — 7x7 Spearman correlation matrix (4 marks + mC + hmC + delta-ratio)
- diffbind_quantitative_oe_comparison.tsv — 16 enriched O/E quadrants (4 marks x mC/hmC perspectives) with Fisher OR/p
- diffbind_logistic_model_coefficients.tsv — 4-mark logistic regression ORs (95% CI, p) for predicting mC hypermethylation
- diffbind_convergence_per_gene.tsv — per-gene cascade-concordance score (0-4 marks) for all 20,915 genes
- diffbind_convergent_genes_3plus.tsv — 178 genes with ≥3 concordant marks, by DMR status
- permutation_37_summary.tsv — 15 gene-level label-shuffle permutation tests with Fisher OR, perm z, empirical p, concordance (all Confirmed)

## Data-quality flags
- **Sections 34, 35, 36 never wrote their comparison TSVs.** `permutation_34_comparison.tsv`, `permutation_35_comparison.tsv`, and `permutation_36_comparison.tsv` are absent from the tables directory. The SLURM run-5 jobs (biomodal_perm_47851857/858/859) each crashed at the Fisher-vs-permutation table step: section 34 with `replacement has length zero` at `perm_results$n_hits[match_idx]` (column-name mismatch in the cached object), and sections 35/36 with an `extract_perm_results -> do.call ... differing number of rows` error. Consequently only the `.rds` caches and the heatmap panels (34a; 35a/b/c; 36a/b/c) exist. Their quantitative confirmation must be read from the inline Fisher recomputations in the logs (section 34 only) and from the gene-level analogues in section 37. These three sections are effectively **visual-only / partially-complete** in run-5.
- **`.rds` caches are not directly readable here** (permutation_34_dmr_marks.rds, permutation_35_atac_loops.rds, permutation_36_domains.rds); all section-34/35/36 numbers above come from the section scripts and the run-5 SLURM logs (downstream/logs/5/biomodal_perm_4785185{7,8,9}.out), not from TSVs.
- **AUC / McFadden R2 for section 33 are not in any TSV** — only in the captured run log (downstream/logs/down/run.txt, lines 5202-5204). The logistic *coefficients* (ORs) are in diffbind_logistic_model_coefficients.tsv, but the AUC=0.818 / extended AUC=0.869 are log-only.
- **Stale "70" in section 32 script comments.** The script header/figure comments ("All 70 CHG DMR genes", "all 70 CHG mC significant genes") predate run-5; the actual run-5 count is 320 (confirmed in chg_mc_significant_dmrs.tsv and run.txt line 4934). Figure titles are generated dynamically from `nrow()`, so the plots say 320 — only the source-code comments are stale.
- **FIGURES.md says CHG has "No Signal" / "No significant changes in non-CpG methylation"** (lines 35-37) and the project CLAUDE.md calls non-CG "noise." That is a global-summary statement; section 32 nonetheless reports 320 significant CHG mC DMRs under the modality GLM. This is not a contradiction in the numbers (the CHG effect sizes are ~100x smaller and the global methylation level is ~0.6%), but the prose "No Signal" framing and the section-32 "320 significant DMRs" framing are in tension and should be reconciled in the manuscript (recommend describing CHG hits as statistically significant but biologically negligible and CG-discordant).
- **`extended_logistic_regression_6mark.tsv` is out of this group's scope.** It was in the A6 starting table list but is written by `section_41_dnmt3a_vs_dnmt3b_discrimination.R` (adds H3K36me2/me3 to the model), not by any section 32-37 script. It is not cited as a 32-37 result here. For reference its 6-mark model still ranks H2AK119ub top (OR=10.29, p=1.0e-10).
- **Section 37 log vs TSV minor numeric drift.** A few perm-z values differ slightly between the SLURM log and the final TSV (e.g. test 37-08 z=-31.78 in log vs -31.93 in permutation_37_summary.tsv; 37-11 z=-4.29 vs -4.30) due to permutation stochasticity between runs. The TSV is treated as canonical; all section-37 numbers above are from the TSV.
