# Sections 38-41: H3K36me2/me3 analysis (NSD/SETD2 axis) and DNMT3A vs DNMT3B mechanistic discrimination

## Summary
This group tests whether BAP1-KO gene-body hypermethylation can be attributed to a specific de novo DNA methyltransferase pathway, using the established histone-to-DNMT recruitment logic: H3K36me2 (deposited by NSD1/NSD2) recruits **DNMT3A**, whereas H3K36me3 (deposited by SETD2 in transcribed gene bodies) recruits **DNMT3B**. Sections 38-39 profile H3K36me3 and H3K36me2 DiffBind changes against the mC/hmC DMR landscape; section 40 analyzes the two marks jointly (conversion dynamics, cross-mark correlations, GO); section 41 integrates H2AK119ub + H3K36me2/me3 + methylation into a logistic-regression and pathway-attribution framework. The main answer: of the 7,513 hypermethylated genes, only **26.8% receive any DNMT-pathway attribution**, and among those, the **DNMT3A axis (K119ub + H3K36me2) outnumbers the DNMT3B axis (H3K36me3) by ~19:1** (1,630 vs 84 genes). In the 6-mark logistic model, H2AK119ub is by far the dominant positive predictor of hypermethylation (OR = 10.29), H3K36me3 is non-significant (OR = 1.00), and H3K36me2 is weakly *negative* (OR = 0.62) — so H3K36me3/DNMT3B is essentially ruled out as a driver, consistent with the broader thesis that hypermethylation is a DNMT3A-redistribution / TET-block phenomenon rather than a SETD2/DNMT3B gene-body program.

## Section 38: section_38_h3k36me3_gene_body_analysis
### Analysis question
Does the SETD2-deposited gene-body mark H3K36me3 (the DNMT3B-recruiting signal) change coordinately with mC/hmC at the same genes in BAP1-KO? If DNMT3B drives gene-body hypermethylation, me3-gained genes should be enriched for mC-up / hmC-down.
### Key results
- me3 × mC direction concordance: Fisher OR = 4.30, p = 1.24e-04, n = 207 genes with both a significant mC DMR and significant me3 peak (source: h3k36me3_direction_concordance.tsv)
- me3 × hmC direction concordance: Fisher OR = 0.163, p = 1.60e-05, n = 184 genes (source: h3k36me3_direction_concordance.tsv) — me3 gain associates with hmC retention / me3 loss with hmC loss
- me3 fold vs mC difference (gene-level Spearman): rho = +0.018, p = 0.0670 (non-significant) (source: h3k36me3_direction_concordance.tsv)
- me3 fold vs hmC difference (gene-level Spearman): rho = -0.034, p = 5.65e-04 (source: h3k36me3_direction_concordance.tsv)
- Genes with a significant H3K36me3 gene-level peak (me3_fdr < 0.05): 333 total — 285 gained, 48 lost (source: h3k36me3_gene_level_summary.tsv)
- Genes carrying any H3K36me3 peak (gene-level table size): 11,952 (source: h3k36me3_gene_level_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The significant me3 × mC Fisher enrichment (OR = 4.30) looks superficially supportive of DNMT3B-coupled hypermethylation, but it rests on only 207 genes and the *continuous* gene-level correlation between me3 fold and mC difference is essentially zero (rho = +0.018, ns). The much stronger and opposite hmC signal (me3 × hmC OR = 0.163; rho = -0.034, p = 5.6e-04) indicates that me3 tracks hydroxymethylation/active transcription rather than the hypermethylation program. Taken together, H3K36me3 is a weak and inconsistent correlate of the mC changes — DNMT3B is not behaving like the primary effector. The marked imbalance of gained (285) over lost (48) me3 genes reflects a general transcription-coupled shift, not a methylation-specific one.
### Plot inventory
- 38a_h3k36me3_volcano_annotation — H3K36me3 DiffBind volcano (up/down/NS) + significant-peak genomic-annotation pie
- 38b_me3_x_mc_oe_heatmap — 2×2 observed/expected enrichment of me3 direction vs 5mC DMR direction (Fisher)
- 38c_me3_x_hmc_oe_heatmap — 2×2 O/E enrichment of me3 direction vs 5hmC DMR direction
- 38d_me3_vs_methylation_scatter — gene-level me3 fold vs 5mC and vs 5hmC change scatters (Spearman)
- 38e_me3_fold_coordinated_violin — me3 fold at coordinated (mC↑/hmC↓) vs other DMR vs non-DMR genes (Wilcoxon)
- 38f_dmr_me3_bed_overlap — coordinate-level overlap of mC hyper/hypo DMRs with me3 gained/lost peak BEDs (Fisher)
- 38g_chromosome_distribution — chromosome distribution of significant me3 peaks, chr13 histone cluster highlighted
- 38h_me3_fold_by_chromatin_state — me3 fold by 7-category chromatin state (Kruskal-Wallis)

## Section 39: section_39_h3k36me2_boundary_analysis
### Analysis question
Is the NSD-deposited H3K36me2 mark — which recruits DNMT3A and acts as an anti-PRC2 barrier — altered in BAP1-KO, and does it correlate with methylation change and with H3K27me3 boundary breach?
### Key results
- H3K36me2 peaks at H3K27me3 boundary zones (±5 kb of peak edges): n = 660, median me2 Log2FC = +0.636 (source: h3k36me2_k27me3_boundary_analysis.tsv)
- H3K36me2 peaks away from H3K27me3: n = 5,966, median me2 Log2FC = +0.524 (source: h3k36me2_k27me3_boundary_analysis.tsv)
- Boundary vs away Wilcoxon p = 2.97e-47 (source: h3k36me2_k27me3_boundary_analysis.tsv)
- Genes with a significant H3K36me2 gene-level peak: 2,909 — 1,557 gained, 1,352 lost (source: h3k36me2_gene_level_summary.tsv)
- me2 vs mC correlation (cross-mark matrix): Spearman rho = -0.185 (5mC diff vs H3K36me2) (source: h3k36_expanded_correlations.tsv)
- me2 vs hmC correlation: Spearman rho = +0.320 (5hmC diff vs H3K36me2) (source: h3k36_expanded_correlations.tsv)
- me2 vs H3K27me3 correlation: Spearman rho = -0.483 (strongest anticorrelation in the matrix) (source: h3k36_expanded_correlations.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] H3K36me2 is broadly *gained* in BAP1-KO (positive median folds everywhere; 1,557 gained vs 1,352 lost genes), and it is gained slightly more at H3K27me3 boundary zones (+0.636 vs +0.524, p = 3e-47), consistent with H3K36me2 expanding against Polycomb domains. The strong me2 ↔ H3K27me3 anticorrelation (rho = -0.483) supports the mutual-exclusivity / barrier model. Crucially for the DNMT story, me2 change is *negatively* correlated with mC change (rho = -0.185) and *positively* with hmC change (rho = +0.320) — i.e. me2-gained genes tend to lose, not gain, 5mC. This argues that the NSD/H3K36me2 → DNMT3A route is not driving the gene-body hypermethylation either; me2 gain co-travels with the active/hydroxymethylated compartment.
### Plot inventory
- 39a_h3k36me2_volcano_annotation — H3K36me2 DiffBind volcano + up-vs-down annotation bar
- 39b_me2_x_mc_oe_heatmap — 2×2 O/E of me2 direction vs 5mC DMR direction (Fisher)
- 39c_me2_at_k27me3_boundary — me2 fold at H3K27me3 boundary vs non-boundary (violin, Wilcoxon)
- 39d_me2_x_k27me3_oe_heatmap — me2 direction vs H3K27me3 direction mutual-exclusivity O/E
- 39e_me2_vs_methylation_scatter — me2 fold vs 5mC and vs 5hmC change scatters
- 39f_me2_by_chromatin_state — me2 gained/lost/none status at mC DMRs by chromatin state
- 39g_me2_genic_vs_intergenic — me2 fold in genic vs intergenic compartments (Wilcoxon)

## Section 40: section_40_h3k36me2_me3_combined
### Analysis question
Analyzed jointly, do H3K36me2 and me3 show evidence of an altered me2→me3 conversion balance (SETD2 activity), and how do both marks correlate with each other and with the wider chromatin/methylation landscape?
### Key results
- Cross-mark correlation matrix is a 9×9 Spearman table over {5mC diff, 5hmC diff, Delta ratio, H3K36me3, H3K36me2, H2AK119ub, H3K27me3, H3K27ac, ATAC-seq} (source: h3k36_expanded_correlations.tsv)
- H3K36me2 vs H3K36me3 correlation: Spearman rho = -0.082 (essentially decoupled) (source: h3k36_expanded_correlations.tsv)
- H2AK119ub vs 5mC diff: rho = +0.410 — the strongest positive correlate of hypermethylation in the matrix (source: h3k36_expanded_correlations.tsv)
- H3K36me3 vs H2AK119ub: rho = -0.040; H3K36me3 vs 5mC diff: rho = +0.018 (both near zero) (source: h3k36_expanded_correlations.tsv)
- Combined gene profile size: 20,915 genes; 10,301 with an me3 fold, 2,178 with an me2 fold, 1,517 with both me2 and me3 (source: h3k36_combined_gene_profile.tsv)
- Delta ratio vs 5hmC diff: rho = +0.917; vs 5mC diff: rho = -0.810 (confirms delta-ratio captures the mC↓/hmC↑ axis) (source: h3k36_expanded_correlations.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] H3K36me2 and me3 are nearly uncorrelated at the gene level (rho = -0.082), meaning the two marks move independently rather than through a single SETD2-driven me2→me3 conversion flux; there is no clean conversion-block or conversion-shift signature uniting them. The cross-mark matrix is the cleanest single summary of the whole DNMT question: H2AK119ub is the only histone mark with a strong positive correlation to mC gain (rho = +0.41), while both H3K36 marks sit near zero (me3 +0.018) or negative (me2 -0.185) against 5mC. This independently corroborates that the hypermethylation axis is ubiquitin-/Polycomb-linked, not H3K36/transcription-linked.
### Plot inventory
- 40a_expanded_correlation_heatmap — clustered 9-variable Spearman correlation heatmap (methylation + 6 marks)
- 40b_me2_vs_me3_scatter — gene-level me2 vs me3 fold scatter with conversion-shift/block quadrant counts
- 40c_me2_me3_ratio_delta — (me2 fold − me3 fold) delta by DMR status (Kruskal-Wallis)
- 40d_three_way_venn — gene overlap among significant me2, significant me3, and mC DMR sets
- 40e_direction_flow — direction-pattern bar for triple-significant (me2+me3+mC) genes
- 40f_go_comparison — GO BP enrichment for me2-only vs me3-only vs shared differential genes

## Section 41: section_41_dnmt3a_vs_dnmt3b_discrimination
### Analysis question
Integrating H2AK119ub, H3K36me2 and H3K36me3, can gene-body hypermethylation be attributed to DNMT3A (via H2AK119ub or via H3K36me2) versus DNMT3B (via H3K36me3), and which signal best predicts hypermethylation?
### Key results
- Total hypermethylated genes classified (mc_sig & mc_diff > 0): 7,513 (source: dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- Pathway attribution (41g pie): Unknown Mechanism = 5,498 (73.2%); DNMT3A via H2AK119ub = 1,478 (19.7%); Convergent (K119ub + me2) = 268 (3.6%); DNMT3A via H3K36me2 = 152 (2.0%); DNMT3B via H3K36me3 = 84 (1.1%); Convergent (K119ub + me3) = 33 (0.4%) (source: dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- Collapsed axis totals: DNMT3A-exclusive (K119ub + me2) = 1,630 (21.7%) vs DNMT3B-exclusive (me3) = 84 (1.1%) → DNMT3A:DNMT3B ≈ 19.4:1 among single-enzyme calls (source: dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- Among the 7,513 hyper genes: 1,779 K119ub-gained, 428 me2-gained, 563 me2-lost, 117 me3-gained, 12 me3-lost, 6,589 coordinated (mC↑/hmC↓) (source: dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- 6-mark logistic regression odds ratios for predicting hypermethylation: H2AK119ub OR = 10.29 (p = 1.04e-10, ***); H3K27me3 OR = 2.11 (p = 8.46e-05, ***); H3K27ac OR = 0.44 (p = 2.27e-04, ***); H3K36me2 OR = 0.62 (p = 0.0221, *); ATAC-seq OR = 0.50 (p = 0.125, ns); H3K36me3 OR = 1.00 (p = 0.996, ns) (source: extended_logistic_regression_6mark.tsv)
- Logistic coefficients (log-odds / estimate): k119ub_fold = +2.331, k27me3_fold = +0.745, k27ac_fold = -0.818, atac_fold = -0.696, me2_fold = -0.477, me3_fold = -0.0019 (source: extended_logistic_regression_6mark.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The attribution is decisive against a DNMT3B/SETD2 mechanism: H3K36me3-attributable genes are the smallest category (84 genes, 1.1%) and H3K36me3 has zero predictive value in the multivariate model (OR = 1.00, p = 0.996). The DNMT3A axis dominates the *attributable* fraction — almost entirely through H2AK119ub (1,478 + 301 convergent), with H3K36me2 contributing little (152) and even appearing as a weak negative predictor (OR = 0.62). H2AK119ub is overwhelmingly the strongest positive predictor (OR = 10.29), reinforcing the Gretarsson/Chen DNMT3A-UDR recruitment model where DNMT3A is targeted to H2AK119ub-marked nucleosomes. However, the largest single bucket is "Unknown Mechanism" (73.2%): most hypermethylated genes carry no significant gain in any of these three recruiting marks, which is more consistent with passive hypermethylation from a TET-mediated demethylation block (the 6,589/7,513 coordinated mC↑/hmC↓ genes) than with active de novo DNMT recruitment. The pie therefore should be read as "where recruitment evidence exists, it is DNMT3A/K119ub, not DNMT3B/me3 — but for most genes the driver is upstream of de novo methylation altogether."
### Plot inventory
- 41a_me3_at_hyper_k119ub_stratified — me3 fold at hyper genes stratified by K119ub-gain status (pairwise Wilcoxon)
- 41b_me2_me3_independent_pathway — me2/me3 fold at hyper genes lacking K119ub gain (ubiquitin-independent test)
- 41c_logistic_regression_forest — forest plot of 6-mark logistic ORs predicting hypermethylation, with 4-mark vs 6-mark AUC
- 41d_pathway_hypermethylation_rate — % hypermethylated by pathway group (K119ub-only / H3K36-only / Convergent / Neither; chi-square)
- 41e_decision_matrix_heatmap — me3 status × K119ub status O/E decision matrix with DNMT3A/DNMT3B/Both annotations
- 41f_k119ub_vs_me3_scatter — H2AK119ub vs H3K36me3 fold scatter colored by mC status, with mechanism quadrants
- 41g_pathway_attribution_pie — pie of hyper-gene DNMT pathway attribution (the headline 6-category breakdown)
- 41h_top50_summary_heatmap — top-50 hypermethylated genes × multi-mark z-scored profile, annotated by pathway

## Cross-section synthesis
Sections 38-41 form a single argument that runs from descriptive to discriminative. Section 38 shows H3K36me3 (the DNMT3B signal) is at best a weak, inconsistent correlate of methylation change; section 39 shows H3K36me2 (the DNMT3A/NSD signal) is broadly gained and acts against Polycomb boundaries but is *anti*-correlated with 5mC gain; section 40's cross-mark matrix crystallizes why — among all six histone/accessibility marks, only H2AK119ub correlates positively with hypermethylation (rho = +0.41), while both H3K36 marks sit near zero or negative. Section 41 then formalizes this with logistic regression (H2AK119ub OR = 10.29 vs H3K36me3 OR = 1.00) and explicit pathway attribution (DNMT3A axis 1,630 genes vs DNMT3B axis 84 genes, ~19:1), while flagging that 73% of hypermethylated genes carry no recruiting-mark evidence at all. For the paper's thesis, this group cleanly excludes a SETD2/DNMT3B gene-body methylation program, points the residual de novo activity to the DNMT3A/H2AK119ub (BAP1-loss → Polycomb-ubiquitin) axis, and leaves the bulk of hypermethylation to the upstream TET-block mechanism documented elsewhere (the 6,589 coordinated mC↑/hmC↓ genes).

## Tables used
- h3k36me3_direction_concordance.tsv — Fisher OR/p and Spearman rho/p for me3 direction vs mC and vs hmC (2 rows)
- h3k36me3_gene_level_summary.tsv — per-gene nearest-TSS H3K36me3 fold/FDR joined to mC/hmC DMR status and chromatin state (11,952 genes)
- h3k36me2_k27me3_boundary_analysis.tsv — me2 peak counts and median fold at vs away from H3K27me3 boundaries, with Wilcoxon p
- h3k36me2_gene_level_summary.tsv — per-gene significant H3K36me2 fold/FDR joined to mC/hmC DMR status (2,909 genes)
- h3k36_expanded_correlations.tsv — 9×9 Spearman cross-mark correlation matrix (methylation + 6 marks)
- h3k36_combined_gene_profile.tsv — merged per-gene me2 + me3 + multi-mark + methylation profile (20,915 genes)
- dnmt3a_vs_dnmt3b_pathway_attribution.tsv — per-hyper-gene boolean mark flags + DNMT pathway label (7,513 genes)
- extended_logistic_regression_6mark.tsv — 6-mark logistic regression term/estimate/OR/CI/p for predicting hypermethylation

## Data-quality flags
- **Embedded-newline parsing hazard.** `dnmt3a_vs_dnmt3b_pathway_attribution.tsv` has literal `\n` inside the `pathway` field (the R `case_when` labels contain newlines), so every data record spans two physical lines (15,026 physical lines = 7,513 records). Naive `wc -l` or column counting overcounts; numbers here were reconstructed by joining lines to 13 fields. Anyone re-reading this table must account for the embedded newlines.
- **Do NOT confuse with section 24/25 tables.** The `dnmt3a_model_coefficients.tsv`, `dnmt3a_feature_importance.tsv`, `dnmt3a_model_comparison.tsv`, `dnmt3a_cv_results.tsv`, `dnmt3a_stratified_comparison.tsv`, `dnmt3a_interaction_results.tsv`, and `dnmt3a_exclusive_model_comparison.tsv` files (with the TET AUC=0.793 vs DNMT3A AUC=0.693 result) are written by `section_24_dnmt3a_prediction.R` / `section_25_delta_ratio_models.R`, NOT by sections 38-41. Only `extended_logistic_regression_6mark.tsv` and `dnmt3a_vs_dnmt3b_pathway_attribution.tsv` belong to section 41.
- **AUC / quadrant / Venn numbers are console-only.** Section 41's 4-mark and 6-mark AUCs (printed to stdout via `cat`), section 41a/41b Wilcoxon p-values, section 40b's me2/me3 quadrant counts and me2-vs-me3 Spearman rho, and section 40d's triple-overlap Venn counts are NOT persisted to any TSV. They cannot be verified from tables and are omitted here rather than quoted. If needed, they require re-running the scripts and capturing stdout.
- **Plot format note.** All 29 panels (38a-h, 39a-g, 40a-f, 41a-h) rendered as JPG + PDF + SVG; no PNG was produced in this batch. All panel folders are populated (no empty/missing figures).
- **me3 master volcano counts not in a table.** Section 38a's peak-level up/down/NS counts and 39a's me2 peak-level counts live only in the DiffBind master files and console output. The gene-level summary tables give *gene*-level significant counts instead (me3: 333 genes; me2: 2,909 genes), which are what is reported here.
- **No source-vs-source numeric discrepancies** were found among the eight section-38-41 tables; prose docs (FIGURES.md/TODO.md) contain no section-38-41 numbers to reconcile against (TODO §22 is a conceptual placeholder predating the data).
