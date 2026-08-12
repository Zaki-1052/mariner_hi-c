# Sections 27-31: Hi-C 3D genome integration: loop anchors, A/B compartments, Polycomb targets, MeCP2-loop association

## Summary
This group connects BAP1-KO DNA-methylation changes to 3D genome architecture, asking whether differential Hi-C loop anchors, A/B compartments, Polycomb domains, and the methyl-reader MeCP2 spatially organize the mC↑/hmC↓ phenotype. It integrates modality DMR calls (`mc_dmr`, `hmc_dmr`), Mariner differential loops (late timepoint), HOMER A/B compartment PC1, the 7-category chromatin-state system, H2AK119ub peaks, and MeCP2 CUT&RUN differential binding. The central answers: (1) methylation associates with loop direction, not loop presence — lost-loop anchors are ~2.6x more likely hypermethylated than gained (OR=2.59) and loop direction survives multivariate adjustment; (2) hypermethylation is strongly enriched in the A (euchromatin) compartment (OR=13.64), matching the TET-KO DNMT3A-redistribution signature; (3) Polycomb/H3K27me3 targets are depleted from hypermethylation (OR=0.063) and enriched in hypomethylation (OR=9.80), confirming that active genes (not repressed heterochromatin) are the de-novo-methylation targets; (4) MeCP2 significant gain is modestly enriched at lost vs gained loops (enrichment 1.51 vs 0.65), consistent with more 5mC drawing more MeCP2 at disrupted anchors.

## Section 27: methylation_hic_loop_anchor_integration
### Analysis question
Are genes at differential loop anchors enriched for coordinated mC↑/hmC↓, and does methylation direction track loop direction (lost vs gained), K119ub gain, and shared-anchor status?

### Key results
- Pooled loop-anchor genes are LESS coordinated than background, GREAT-style: OR = 0.672, p = 5.33e-23; anchor 43.37% vs background 53.26% coordinated (source: methylation_loop_anchor_coordinated_enrichment.tsv).
- Nearest-gene method agrees: OR = 0.683, p = 2.55e-17; anchor 43.41% vs background 52.90% (source: methylation_loop_anchor_coordinated_enrichment.tsv).
- Directional split — Lost-vs-Gained hypermethylation: lost 46.86% hyper vs gained 25.40% hyper, Fisher OR = 2.589, p = 1.10e-35 (source: methylation_direction_by_loop_direction.tsv).
- Lost-vs-Background hyper OR = 1.569 (p = 6.13e-15); Gained-vs-Background hyper OR = 0.606 (p = 2.41e-20); n_lost = 1,321, n_gained = 1,827 anchor genes in DMR universe (source: methylation_direction_by_loop_direction.tsv).
- Delta-ratio (Lost vs Gained) Wilcoxon p = 1.07e-61; median delta-ratio Lost = -0.01413 vs Gained = +0.00105 (more negative = greater TET impairment) (source: methylation_direction_by_loop_direction.tsv).
- K119ub-gained-at-anchor × hypermethylated convergence: OR = 1.791, p = 2.05e-06; contingency a=147, b=862, c=164, d=1723; triple convergence (K119ub gained + hyper + lost loop) = 119 genes (source: k119ub_loop_methylation_convergence.tsv).
- Logistic regression — loop direction is the dominant predictor: loop_dir_binaryLost OR = 2.290, p = 5.53e-32; anchor_typeActive_Enhancer OR = 1.660 (p = 1.91e-03, vulnerable); anchor_typeRepressed_Promoter OR = 0.180, anchor_typePolycomb OR = 0.232 (protected); distance terms NS (log_loop_distance p = 0.79, log_dist_tss p = 0.66) (source: logistic_regression_methylation_loop.tsv).
- Linear model of delta-ratio: loop_dir_binaryLost coefficient = -0.01762, p = 3.00e-63 (lost loops have more negative delta-ratio); Polycomb/Repressed_Promoter anchors have positive coefficients (less TET impairment) (source: linear_model_delta_ratio_loop.tsv).
- Shared-anchor hubs (n=236 genes) show LOWER coordinated rate (34.32%) than non-shared anchors (44.17%) and background (53.26%); shared-vs-BG OR = 0.459 (p = 7.74e-09), shared-vs-non-shared OR = 0.661 (p = 3.92e-03) (source: shared_anchor_methylation_profile.tsv).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The pooled depletion of coordination at anchors (OR<1) is not a null result but a cancellation artifact: lost and gained loops carry opposite methylation signatures that average out when pooled. Stratifying by direction reveals the real biology — lost loops (BAP1-KO down_in_mutant) mark genes undergoing TET-impaired hypermethylation, whereas gained loops mark genes trending hypomethylated. Because loop direction stays significant after controlling for distance-to-TSS, loop distance, and anchor chromatin state, loop disruption and methylation are coupled rather than both being passive readouts of promoter proximity. The K119ub→hyper convergence and the 119 triple-convergence genes support a causal chain in which ectopic H2AK119ub accrual at lost-loop anchors blocks TET and drives 5mC gain. Active_Enhancer anchors being the most vulnerable (OR=1.66) and Polycomb/Repressed_Promoter anchors being protected (OR≈0.2) prefigures the Polycomb-depletion result in section 30. The shared-anchor exception suggests structural hub anchors are buffered against methylation perturbation.

### Plot inventory
- 27a_coordinated_enrichment_at_loop_anchors — grouped bar, % coordinated at anchors vs background (GREAT + nearest-gene).
- 27b_mc_direction_by_loop_direction — bar, hypermethylation rate Lost vs Gained vs Background.
- 27b_mc_diff_violin_by_loop_direction — violin/box of mC difference by loop direction.
- 27b_delta_ratio_violin_by_loop_direction — violin/box of delta demethylation ratio by loop direction.
- 27c_k119ub_methylation_loop_convergence_heatmap — O/E heatmap, K119ub status × methylation status.
- 27c_triple_convergence_summary — bar, multi-layer convergence gene counts.
- 27d_logistic_regression_forest_plot — forest plot of odds ratios for predictors of hypermethylation.
- 27d_linear_model_coefficients — coefficient plot for delta-ratio linear model.
- 27e_shared_anchor_coordinated_rate — bar, coordinated rate at shared vs non-shared anchors vs background.
- 27e_shared_anchor_delta_ratio_violin — violin/box of delta-ratio across shared/non-shared/background.

## Section 28: coordinated_mc_hmc_analysis
### Analysis question
Characterize the majority coordinated Q4 (mC↑/hmC↓) genes against non-coordinated (Q1+Q2+Q3) genes across 9 epigenomic dimensions, including Hi-C loop involvement.

### Key results
- Quadrant breakdown (genes significant in both mC and hmC): Q1 (mC dn/hmC dn) = 411; Q2 (mC dn/hmC up) = 1,255; Q3 (mC up/hmC up) = 116; Q4 (mC up/hmC dn, COORDINATED) = 6,589; non-Q4 reference = 1,782; total = 8,371 (source: coordinated_gene_characteristics.tsv).
- Loop involvement: Q4 coordinated genes are LESS likely at loop anchors — 333/6,589 = 5.1% vs non-Q4 210/1,782 = 11.8% (source: coordinated_gene_characteristics.tsv).
- Q4 fraction of co-significant genes = 6,589/8,371 = 78.7% (the dominant pattern) (source: coordinated_gene_characteristics.tsv).
- Table provides per-gene mc_diff, hmc_diff, combined_effect, RNA-seq log2FC/padj, net_atac, k119ub_gb_log2fc, k27ac_status, chromatin_state, mecp2_status, loop_involved for all 8,371 genes (source: coordinated_gene_characteristics.tsv).
- Note: the 9-dimension Fisher/Wilcoxon p-values (|mC diff|, log2FC, ATAC, K119ub, K27ac, chromatin state, MeCP2, loop) are emitted to console only and are not written to a TSV (visual-only in 28a composite); only the master gene table is persisted.

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The coordinated Q4 class dominates the co-significant set, reaffirming the project's headline that 5hmC→5mC conversion (TET block) is the prevailing direction of change. Crucially, Q4 genes are depleted at loop anchors (5.1% vs 11.8%), which dovetails with section 27: the strongest coordinated methylation occurs at non-structural genes, while loop-anchor genes split by loop direction and partly cancel. This positions the coordinated phenotype as primarily a gene-body chemical-conversion effect rather than a direct 3D-architecture effect, with the 3D coupling concentrated specifically at the lost-loop subset.

### Plot inventory
- 28a_coordinated_composite — 3×3 grid: Q4 vs non-Q4 across all 9 dimensions.
- 28a_panel_mc_diff, 28a_panel_hmc_diff — |mC| / |hmC| effect-size boxplots.
- 28a_panel_log2fc — RNA-seq log2FC boxplot.
- 28a_panel_atac — net ATAC change boxplot.
- 28a_panel_k119ub — K119ub gene-body signal boxplot.
- 28a_panel_k27ac — H3K27ac gained/lost bar.
- 28a_panel_chromatin — chromatin-state stacked bar.
- 28a_panel_mecp2 — MeCP2 up/down bar.
- 28a_panel_loop — Hi-C loop involvement stacked bar.
- 28b_mc_hmc_concordance_scatter — mc_diff vs hmc_diff scatter colored by concordance.
- 28c_mc_vs_expression_per_group — mC change vs RNA-seq log2FC per group.
- 28d_go_enrichment_coordinated — GO BP dot plot (Q4 vs non-Q4 background).
- 28e_all_quadrants_comprehensive — 4-quadrant epigenomic panel.

## Section 29: ab_compartment_methylation_mapping
### Analysis question
Do BAP1-KO methylation changes follow the Lopez-Moyado TET-KO signature — hypermethylation in A (euchromatin) and hypomethylation in B (heterochromatin) — indicating DNMT3A redistribution?

### Key results
- mC hyper → A enriched: OR = 13.642 (95% CI 11.92-15.67), p = 0; 48.76% of A-compartment genes hyper vs 6.52% of B (a=6,022, b=244, c=6,329, d=3,499) (source: compartment_fisher_tests.tsv).
- mC hypo → B enriched: OR = 1.726 (1.58-1.89), p = 8.33e-32; 24.21% hypo in B vs 15.61% in A (source: compartment_fisher_tests.tsv).
- hmC loss → A enriched: OR = 9.774 (8.83-10.83), p = 0; 60.06% of A genes lose hmC vs 13.33% of B (source: compartment_fisher_tests.tsv).
- hmC gain → B enriched: OR = 3.030 (2.73-3.36), p = 4.30e-95 (reciprocal test) (source: compartment_fisher_tests.tsv).
- Compartment composition: A = 12,351 genes (75.7% hyper rate among mC-sig), B = 3,743 genes (21.2% hyper rate); A median mc_diff = +0.01073, B = -0.00417 (source: compartment_methylation_summary.tsv).
- B→A shifted bins are DEPLETED for mC hyper / enriched for hypo: B→A-shift × mC-hyper OR = 0.0625, p = 1.37e-66 (B→A shift n=427 genes, only 17 hyper vs 239 hypo); A→B-shift × mC-hypo OR = 0.454, p = 5.53e-04 (A→B n=202, 117 hyper vs 18 hypo) (source: compartment_fisher_tests.tsv, compartment_methylation_summary.tsv).
- Gene matching: 16,094 genes assigned to compartments and matched to DMRs (source: compartment_gene_assignment.tsv, 16,095 lines incl. header).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The strong A-compartment enrichment of both hypermethylation (OR=13.6) and hmC loss (OR=9.8) is the cleanest cross-modality confirmation that BAP1 loss phenocopies the TET-KO DNMT3A-redistribution signature: when TET is impaired in euchromatin, DNMT3A re-deposits 5mC in the open A compartment where it normally counteracts methylation. The static A/B picture supports convergent mechanisms with TET-KO. The dynamic compartment-shift result is opposite in sign and equally informative — bins actively switching B→A in the mutant are losing heterochromatic character and are enriched for hypomethylation (not hyper), i.e. they are gaining activity/accessibility rather than acquiring de novo methylation. So the de novo methylation targets are stable-A euchromatin, not freshly opened regions.

### Plot inventory
- 29a_mc_violin_by_compartment — mC difference violin, A vs B.
- 29b_hmc_violin_by_compartment — hmC difference violin, A vs B.
- 29c_mc_violin_by_shift — mC difference violin by compartment shift (B→A, A→B, Stable).
- 29d_hmc_violin_by_shift — hmC difference violin by shift.
- 29e_dmr_direction_stacked_bar — mC DMR direction proportions by compartment and shift.
- 29f_pc1_vs_mc_scatter — control PC1 vs mC difference scatter with loess + Spearman.
- 29g_composite_compartment_summary — composite (29a+29b+29e).

## Section 30: polycomb_target_enrichment
### Analysis question
Are classic Polycomb/H3K27me3 targets the primary hypermethylation targets, or — as the dual-mechanism model predicts — are they protected because heterochromatin is inaccessible to DNMT3A?

### Key results
- Chromatin-State Polycomb (Repressed_Promoter + Polycomb + Bivalent, n=3,447) DEPLETED from hypermethylation: OR = 0.0633 (95% CI 0.0533-0.0748), p = 0, q = 0; 4.41% Polycomb hyper rate vs 42.14% non-Polycomb (source: polycomb_fisher_tests.tsv).
- Chromatin-State Polycomb ENRICHED in hypomethylation: OR = 9.797 (8.996-10.673), p = 0, q = 0; 49.14% Polycomb hypo vs 8.98% non-Polycomb (source: polycomb_fisher_tests.tsv).
- Result is robust across all three Polycomb definitions: Strict (no Bivalent) × mC-hyper OR = 0.0422; H3K27me3-overlap × mC-hyper OR = 0.1027; both q = 0 (source: polycomb_fisher_tests.tsv).
- Per-chromatin-state hypermethylation rate: Active_Promoter 71.77% hyper (OR=10.00, q=0, n=6,356) and Active_Enhancer 65.31% (OR=3.37, q=3.93e-05) are the hyper targets; Polycomb 1.79% (OR=0.0323, q=1.15e-09) and Repressed_Promoter 2.99% (OR=0.0427, q=0) are protected; Bivalent_Promoter 36.00% is at background (OR=1.003, ns) (source: polycomb_per_chromatin_state_enrichment.tsv).
- Magnitude is also smaller at Polycomb among hyper genes: MC hyper Polycomb median |mc_diff| = 0.01886 (n=152) vs non-Polycomb = 0.02904 (n=7,361), Wilcoxon q = 1.73e-10 (source: polycomb_methylation_magnitude.tsv).
- hmC also depleted at Polycomb: Chromatin-State × hmC-hypo OR = 0.131 (q=0); × hmC-hyper OR = 11.39 (q=0) (Polycomb regions gain rather than lose hmC) (source: polycomb_fisher_tests.tsv).
- Universe = 20,915 genes with chromatin state + DMR data (source: polycomb_fisher_tests.tsv n_universe; polycomb_gene_classification.tsv = 20,915 genes).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the decisive falsification of the naive "H2AK119ub marks Polycomb, so Polycomb genes get hypermethylated" expectation. Polycomb-marked, H3K27me3-repressed chromatin is essentially excluded from de novo methylation (hyper OR≈0.06) and instead drifts hypomethylated (OR≈9.8) — exactly what the dual-mechanism model predicts, because compact heterochromatin is inaccessible to DNMT3A. The hypermethylation burden falls on normally active genes (Active_Promoter ~72%, Active_Enhancer ~65%) that acquire ectopic H2AK119ub and lose TET activity. This explains the section-27 logistic result (Active_Enhancer anchors vulnerable, Polycomb anchors protected) and the section-29 A-compartment enrichment (active genes live in A) as one coherent mechanism.

### Plot inventory
- 30a_polycomb_vs_non_polycomb_stacked_bar — mC DMR status, Polycomb vs non-Polycomb.
- 30b_fisher_forest_plot — forest plot of ORs across Polycomb definitions and directions.
- 30c_mc_magnitude_violin — |mc_diff| violins, Polycomb vs non-Polycomb by direction.
- 30d_per_state_hypermethylation_rate — % hyper per chromatin state with genome-wide reference line.
- 30e_hmc_magnitude_violin — |hmc_diff| violins, Polycomb vs non-Polycomb by direction.
- 30f_composite_polycomb_summary — composite (30a+30b+30d).

## Section 31: mecp2_loop_anchor_integration
### Analysis question
Does MeCP2 (CUT&RUN) differential binding associate with Hi-C loop anchor positions and loop direction, connecting the 5mC-preferring reader to 3D architecture?

### Key results
- MeCP2 significant gain × loop direction (gene-level Fisher): OR = 0.424, p = 1.32e-03; MeCP2-Sig-Up is enriched at Lost loops (Observed 37 vs Expected 24.5, enrichment 1.51) and depleted at Gained loops (Observed 23 vs Expected 35.5, enrichment 0.65) (source: mecp2_loop_direction_fisher.tsv).
- Per-loop direct overlap: 2,910 differential loops scored for any-anchor MeCP2 overlap; export table has 2,910 loops (source: mecp2_loop_anchor_overlap.tsv, 2,911 lines incl. header).
- Gene-level MeCP2 fold by loop group: Gained n=2,246, Lost n=1,553, Background (no loop) n=17,701 anchor genes with MeCP2 data (source: mecp2_fold_at_loop_anchor_genes.tsv group counts).
- Loop-level scatter: 2,645 loops have MeCP2 gene-level data for the loop-logFC vs mean-MeCP2-fold correlation (source: mecp2_loop_gene_level_scatter_data.tsv, 2,646 lines incl. header).
- The 31a Fisher ORs (gained/MeCP2-Up, lost/MeCP2-Down), 31c Wilcoxon p-values (Gained/Lost/Background), and 31d Spearman rho are emitted to console only and are NOT written to a persisted TSV — the saved tables are per-loop/per-gene data, not the test statistics.

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The directional MeCP2 result is consistent with the Mellen-et-al reader logic: MeCP2 binds 5mC with high affinity and 5hmC with low affinity, so the lost-loop anchors that hypermethylate in BAP1-KO (section 27) draw more MeCP2 (gain enrichment 1.51 at lost loops), while gained loops — which trend hypomethylated — are depleted of MeCP2 gain (0.65). The gene-level Fisher OR<1 reflects that the contingency is built as Gained-vs-Lost, so the <1 value encodes "MeCP2 gain favors Lost over Gained," i.e. the biologically expected direction. This closes the loop from chemistry (5hmC→5mC) to reader occupancy (MeCP2) to 3D structure (lost loops), the core mechanistic narrative of the paper's title.

### Plot inventory
- 31a_mecp2_peak_overlap_at_loop_anchors — grouped bar, % loops with anchor MeCP2 overlap by loop direction.
- 31b_mecp2_loop_direction_fisher_heatmap — Obs/Exp/Enrichment heatmap, MeCP2-Sig-Up × loop direction.
- 31c_mecp2_fold_by_loop_direction — violin/box of MeCP2 fold change, Gained/Lost/Background.
- 31d_loop_logfc_vs_mecp2_fold_scatter — loop logFC vs mean MeCP2 fold scatter with gene labels.

## Cross-section synthesis
These five sections build a single coherent 3D-genome argument that supports the paper's thesis that BAP1 loss restructures both methylation and the MeCP2-readable chromatin landscape. Section 28 establishes that the dominant coordinated mC↑/hmC↓ phenotype is largely a gene-body chemical effect that is depleted, not enriched, at loop anchors (5.1% vs 11.8%). Section 27 resolves why: at anchors the effect splits by loop direction — lost loops hypermethylate (OR=2.59 vs gained), gained loops trend hypomethylated — and loop direction predicts methylation after multivariate adjustment, with Active_Enhancer anchors vulnerable and Polycomb anchors protected. Section 30 generalizes that protection genome-wide: Polycomb/H3K27me3 targets are excluded from hypermethylation (OR=0.063) and the de novo methylation lands on active genes, while section 29 places those active genes spatially in the A compartment (hyper OR=13.6), matching the TET-KO DNMT3A-redistribution signature. Section 31 completes the causal chain by showing the 5mC reader MeCP2 gains preferentially at the same hypermethylating lost-loop anchors (enrichment 1.51 vs 0.65). Together: ectopic H2AK119ub → TET block → 5mC gain at active, A-compartment, lost-loop genes → increased MeCP2 → disrupted 3D architecture.

## Tables used
- methylation_loop_anchor_coordinated_enrichment.tsv — 27a Fisher OR for coordinated-gene enrichment at anchors (GREAT + nearest-gene).
- methylation_direction_by_loop_direction.tsv — 27b hypermethylation rates and delta-ratio by loop direction (Lost/Gained/Background).
- k119ub_loop_methylation_convergence.tsv — 27c K119ub-gain × hyper Fisher and triple-convergence counts.
- logistic_regression_methylation_loop.tsv — 27d logistic ORs for predictors of hypermethylation.
- linear_model_delta_ratio_loop.tsv — 27d linear-model coefficients for delta demethylation ratio.
- shared_anchor_methylation_profile.tsv — 27e coordinated rate and delta-ratio at shared vs non-shared anchors.
- anchor_gene_associations_great.tsv — 27 GREAT-style anchor→gene map (6,795 associations); input to section 31.
- coordinated_gene_characteristics.tsv — 28 master per-gene table (8,371 genes, 22 columns) with quadrant and loop_involved.
- compartment_fisher_tests.tsv — 29 four directional + two shift Fisher tests for A/B enrichment.
- compartment_methylation_summary.tsv — 29 per-compartment / per-shift methylation summary statistics.
- compartment_gene_assignment.tsv — 29 per-gene compartment, shift, PC1, mC/hmC (16,094 genes).
- polycomb_fisher_tests.tsv — 30 Polycomb × DMR Fisher tests across 3 definitions + per-state (19 tests).
- polycomb_methylation_magnitude.tsv — 30 Wilcoxon of |diff| magnitude, Polycomb vs non-Polycomb.
- polycomb_per_chromatin_state_enrichment.tsv — 30 per-state hyper rates and Fisher ORs.
- polycomb_gene_classification.tsv — 30 per-gene Polycomb flags + DMR status (20,915 genes).
- mecp2_loop_direction_fisher.tsv — 31b MeCP2-gain × loop-direction Obs/Exp/Enrichment + Fisher OR.
- mecp2_loop_anchor_overlap.tsv — 31a per-loop MeCP2 anchor-overlap booleans (2,910 loops).
- mecp2_fold_at_loop_anchor_genes.tsv — 31c per-gene MeCP2 nearest-fold by loop group.
- mecp2_loop_gene_level_scatter_data.tsv — 31d per-loop logFC vs mean MeCP2 fold (2,645 loops).

## Data-quality flags
- STALE PROSE vs TSV (section 27): TODO.md row 27 and FIGURES.md state Lost-vs-Gained "OR=2.54" and logistic "OR=2.42". Current TSVs show Lost-vs-Gained Fisher OR = 2.589 (methylation_direction_by_loop_direction.tsv) and logistic loop_dir_binaryLost OR = 2.290 (logistic_regression_methylation_loop.tsv). The "OR=2.54/2.42" values are not confirmed in the current tables — treat as prior-run numbers. Hyper rates also differ: prose says lost 43.0%/gained 22.9%/bg 31.6%, TSV says 46.86%/25.40%/35.98%.
- STALE PROSE vs TSV (section 27 shared anchors): TODO.md says shared 35.2% / non-shared 45.7% / background 54.6% coordinated with shared OR=0.45; current TSV says 34.32% / 44.17% / 53.26% with shared-vs-BG OR=0.459. Same conclusion, slightly different numbers (prior run).
- STALE PROSE vs TSV (section 28): TODO.md/FIGURES.md state "Q4 (5,708) vs non-Q4 (1,042)" and loop involvement "5.1% Q4". Current coordinated_gene_characteristics.tsv shows Q4 = 6,589 and non-Q4 = 1,782 (Q1=411, Q2=1,255, Q3=116); loop involvement Q4 = 5.1% (matches) but non-Q4 = 11.8% (prose implies the discordant 15.5% from concordant_vs_discordant_trends.md, which uses a different Q2-only comparison). The "5,708/1,042" counts are NOT confirmed in the current table — prior-run numbers.
- STALE PROSE vs TSV (section 29): TODO.md/FIGURES.md state mC-hyper→A "OR=14.71" and hmC-loss→A "OR=9.35". Current compartment_fisher_tests.tsv shows OR = 13.642 and 9.774 respectively. Compartment composition prose (A=12,347 / B=3,739; 16,086 matched) vs TSV (A=12,351 / B=3,743; 16,094 matched) — minor prior-run drift. Group-brief note "Hyper enriched in A (OR~14.71)" is therefore [UNVERIFIED: 14.71 per TODO.md, not confirmed in tables — actual is 13.64].
- STALE PROSE vs TSV (section 30): TODO.md/FIGURES.md state Polycomb hyper "OR=0.064" and hypo "OR=8.71", and "Active_Promoter 65.2% hyper vs Repressed_Promoter 2.3%". Current TSVs show hyper OR = 0.0633, hypo OR = 9.797 (polycomb_fisher_tests.tsv); Active_Promoter = 71.77% and Repressed_Promoter = 2.99% (the 65.3% figure is Active_Enhancer, not Active_Promoter) (polycomb_per_chromatin_state_enrichment.tsv). The group-brief "OR~0.064 / OR~8.71" matches the prose but NOT the current tables — actuals are 0.0633 and 9.797.
- VISUAL/CONSOLE-ONLY tests: Section 28's nine cross-dimension Fisher/Wilcoxon statistics, and Section 31's 31a Fisher ORs, 31c Wilcoxon p-values, and 31d Spearman rho, are printed to the R console but not written to any TSV. Only per-gene/per-loop data tables are persisted, so those specific test statistics cannot be re-verified from tables and were intentionally not quoted as confirmed numbers above.
- Loop universe consistency: anchor_gene_associations_great.tsv contains 6,795 anchor-gene associations from the 2,910-loop late set; mecp2_loop_anchor_overlap.tsv covers all 2,910 loops; the slight reductions (2,645 loops in 31d, 16,094 genes in 29) are expected from inner-joins to MeCP2 / compartment data and are not errors.
