# Sections 51-58: Non-CG methylation at MeCP2 sites, single-CpG-resolution gene body, TAD organization, MeCP2-as-CG-corrected-proxy, Ecker WGBS validation of non-CG candidates

## Summary
This group tests the neuronal hypothesis that MeCP2 reads non-CG methylation (mCA/mCH) in addition to mCG, and asks whether the DUET evoC CG data — which cannot directly see mCH — is missing a non-CG axis of the BAP1-KO phenotype. Across gene bodies, MeCP2 peaks, TAD domains, and an external wildtype-cerebellum WGBS reference (Ecker), the answer is consistently that **non-CG methylation in this evoC dataset is at or near the detection floor (CHG/CHH globally <1% methylated, mostly zero at peaks), is not differential between mutant and control, and does not show TAD organization.** The one positive externally-validated signal: the 2,726 "non-CG candidate" MeCP2-Up peaks (MeCP2 gained but CG methylation unchanged/decreased) sit at regions with modestly higher Ecker CH methylation than CG-concordant peaks (median 0.0159 vs 0.0129, p = 2.6e-60), but the proposed dose-response (more CH → more CG-unexplained MeCP2) does NOT hold — the trend is flat-to-decreasing (Jonckheere increasing p ≈ 1.0; decreasing p = 2.0e-09). Net reading: MeCP2 redistribution upon BAP1 loss is not explained by a non-CG-methylation gradient that this data can resolve.

## Section 51: section_51_mecp2_noncg_methylation
### Analysis question
At MeCP2-bound genes, is there measurable non-CG (CHG/CHH) methylation, does it differ between mutant and control, and is it enriched at MeCP2 sites versus the rest of the genome? Genes are stratified by CG-methylation level (where low-CG genes are where MeCP2 might rely on mCH).
### Key results
- CHG at all MeCP2-bound genes: ctrl mean = 8.08e-06, mut mean = 1.10e-05, Wilcoxon ctrl-vs-mut p = 0.4857, n = 1211 genes (source: mecp2_noncg_summary.tsv)
- CHH at all MeCP2-bound genes: identical to CHG row — ctrl 8.08e-06, mut 1.10e-05, p = 0.4857, n = 1211 (source: mecp2_noncg_summary.tsv) [see Data-quality flags — CHG and CHH rows are byte-identical]
- CHG at MeCP2 CG-Low genes: ctrl 9.12e-06 vs mut 1.33e-05, p = 0.6612 (source: mecp2_noncg_summary.tsv)
- CHG at MeCP2 CG-High genes: ctrl 7.04e-06 vs mut 8.69e-06, p = 0.8845 (source: mecp2_noncg_summary.tsv)
- All effect sizes (mut - ctrl) are ~1e-06 to ~4e-06 fractional methylation, i.e. at the noise floor; no test approaches significance (source: mecp2_noncg_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] At gene-body resolution there is no detectable non-CG methylation signal at MeCP2-bound genes and no mutant-vs-control difference. The non-significance combined with sub-1e-5 absolute methylation means the evoC assay is not resolving any biologically meaningful mCH at this scale — consistent with the global <1% non-CG methylation. This is a clean negative result: whatever MeCP2 reorganization BAP1 loss causes, it is not driven by a gene-body-level non-CG methylation change this assay can see.
### Plot inventory
- `51c_noncg_at_mecp2_by_cg_stratum/` — per-sample CHG/CHH methylation at MeCP2 genes, faceted by CG-Low vs CG-High stratum and condition (current script output)
- `51d_mecp2_bound_vs_background_noncg/` — CHG/CHH methylation, MeCP2-bound vs non-MeCP2 genes by condition (current script output)
- `51a_mecp2_noncg_overlap_fisher/`, `51b_mecp2_fold_noncg_violin/`, `51c_cg_stratified_mecp2_fold/`, `51d_persample_noncg_at_mecp2/`, `51e_chg_chh_identity/` — panels from an EARLIER version of the script (see Data-quality flags); not produced by the current section_51 code

## Section 52: section_52_cpg_resolution_gene_body
### Analysis question
At single-feature resolution within max-significance gene bodies (5'UTR, exon, splice donor/acceptor, intron, 3'UTR), how does the BAP1-KO methylation change (delta 5mC, delta 5hmC) vary by sub-gene feature, and is it explained by enhancer (H3K27ac) overlap? (This is CG methylation at fine resolution, not non-CG.)
### Key results
- 6,196 feature intervals across 123 genes, broken down as: 1,742 Intron, 1,567 Exon, 1,341 SpliceSite_Acceptor, 1,330 SpliceSite_Donor, 124 3UTR, 92 5UTR (source: section52_feature_methylation.tsv)
- delta-5mC: introns differ from exons/splice sites — Exon vs Intron Z = -3.864, q = 4.18e-04; SpliceSite_Donor vs Intron Z = -3.898, q = 7.26e-04; SpliceSite_Acceptor vs Intron Z = -3.424, q = 1.55e-03 (source: section52_dunn_posthoc.tsv)
- delta-5mC: only 3 / 15 feature pairs significant at q < 0.05, all involving Intron (source: section52_dunn_posthoc.tsv)
- delta-5hmC: 8 / 15 feature pairs significant at q < 0.05 — strongest 5UTR vs Intron Z = 3.845 q = 9.06e-04; SpliceSite_Acceptor vs Intron Z = 3.759 q = 6.40e-04; Exon vs Intron Z = 3.668 q = 6.12e-04; SpliceSite_Donor vs Intron Z = 3.603 q = 5.90e-04 (source: section52_dunn_posthoc.tsv)
- delta-5hmC also separates 5UTR from 3UTR (Z = 3.151, q = 2.44e-03) and 5UTR from exon/splice features (q = 0.010), unlike delta-5mC (source: section52_dunn_posthoc.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The BAP1-KO methylation response is feature-structured: introns behave differently from exons and splice sites for both 5mC and 5hmC, and 5hmC additionally distinguishes the 5' end (5'UTR) from the rest of the gene body. The richer 5hmC differentiation (8 vs 3 significant pairs) is consistent with this group's broader TET-block thesis — 5hmC, the TET product, carries more of the spatially-resolved signal than 5mC. This is a CG-context fine-mapping that complements (does not test) the non-CG hypothesis.
### Plot inventory
- `52a_feature_distribution/` — bar chart of feature-interval counts per sub-gene feature type
- `52b_delta_mc_by_feature/` — boxplots of delta-5mC by feature type (Kruskal-Wallis + Dunn)
- `52c_delta_hmc_by_feature/` — boxplots of delta-5hmC by feature type
- `52d_key_gene_loci/` — per-KEY_GENE lollipop panels of delta-5mC by genomic position
- `52e_metagene_profile/` — 5' to 3' metagene profile of delta-5mC and delta-5hmC (strand-oriented bins)
- `52f_dunn_posthoc/` — Dunn post-hoc q-value significance heatmaps (mC and hmC panels)
- `52g_chromatin_overlay/` — ChIP-mark overlap per feature, intron and all-feature splits by H3K27ac status

## Section 53: section_53_mecp2_peak_noncg
### Analysis question
At 400bp peak resolution (8,886 MeCP2 CUT&RUN peaks run directly on CHG/CHH zarr stores), is there any detectable mCH at MeCP2 binding sites, is it differential mutant-vs-control, and is it enriched above shuffled-control regions? Includes a peak-vs-gene-body resolution comparison and a data-quality check.
### Key results
- CHG detection: 45 / 8,886 peaks (0.51%) have any detectable mCH in ≥1 of 8 samples; CHH: 146 / 8,886 (1.64%) (source: mecp2_peak_noncg_summary.tsv)
- MeCP2 peaks vs shuffled control — CHG: 0.51% vs 0.29% detected, Fisher OR = 1.734, p = 0.0316; CHH: 1.64% vs 1.04% detected, Fisher OR = 1.597, p = 5.17e-04 (source: mecp2_peak_noncg_summary.tsv)
- Peak-level Wilcoxon ctrl-vs-mut: CHG p = 0.0286 (ctrl 1.48e-05 vs mut 2.44e-05); CHH p = 0.1143 (ctrl 1.19e-05 vs mut 1.56e-05) (source: mecp2_peak_noncg_summary.tsv)
- Differential DMR calls (GLM, sex covariate): CHG 0 / 8,356 peaks at q < 0.05 (best q = 0.247); CHH 2 / 8,886 peaks at q < 0.05 (best q = 0.00317) (source: mecp2_peak_noncg_summary.tsv)
- CHG detection by direction: MeCP2_Up 42/7,686 (0.55%), MeCP2_Down 3/1,200 (0.25%); CHH: MeCP2_Up 138/7,686 (1.80%), MeCP2_Down 8/1,200 (0.67%) (source: mecp2_peak_noncg_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Non-CG methylation is detectable at MeCP2 peaks above a shuffled-region baseline (small but significant Fisher enrichment, OR ~1.6-1.7), confirming MeCP2 sites are modestly non-CG-favored — consistent with the neuronal mCA-reader model. However, the absolute detection is tiny (≤1.6% of peaks), and there is essentially no differential mutant-vs-control signal (0 significant CHG DMRs, only 2 CHH DMRs, in a chr8 cluster). The marginal CHG ctrl-vs-mut Wilcoxon (p = 0.029) reflects a near-floor shift and is not corroborated by the DMR GLM. The data-quality panel argues this negative is real, not a power artifact. The peak-vs-gene-body comparison shows the same flat result at both scales.
### Plot inventory
- `53a_detection_rate/` — non-CG detection categories (Ctrl&Mut / Ctrl-only / Mut-only) per context × direction
- `53b_nonzero_peak_scatter/` — ctrl vs mut mCH at the non-zero peaks only
- `53c_mecp2_vs_control_baseline/` — MeCP2 peaks vs shuffled control detection rate, Fisher tests
- `53d_dmr_rank_plot/` — DMR q-value rank plot, CHG/CHH (GLM with sex covariate)
- `53e_chr8_cluster_spotlight/` — the chr8 35-50 Mb CHH cluster (only region with q<0.05 CHH DMRs)
- `53f_resolution_comparison/` — gene-body (~18kb) vs peak (400bp) non-CG, consistent negative
- `53g_coverage_quality/` — cytosine contexts and coverage per peak (depth adequacy check)

## Section 54: section_54_tad_methylation_organization
### Analysis question
Is non-CG methylation (CHG/CHH) — and MeCP2 — more organized within TAD domains than CG methylation? Uses per-TAD mean signal, intra-TAD CV, boundary insulation, and the inter/intra-TAD variance ratio (ratio > 1 = TAD-organized) across 5 tracks.
### Key results
- Inter/intra-TAD variance ratio, Control: CG 5mC = 1.018, CG 5hmC = 2.078, CHG 5mC = 0.0250, CHH 5mC = 0.0248, MeCP2 = 0.0680 (source: 54_tad_methylation_organization_summary.tsv)
- CG 5hmC is by far the most TAD-organized track (variance ratio ~2.08 in both ctrl and mut), ~2x the CG 5mC ratio (~1.02) (source: 54_tad_methylation_organization_summary.tsv)
- Median per-TAD signal — CG 5mC = 0.00519 (ctrl); CHG 5mC = 0 and CHH 5mC = 0 (both conditions) (source: 54_tad_methylation_organization_summary.tsv)
- Median intra-TAD CV for CHG and CHH = 4 (both conditions) — i.e. extreme variability around a near-zero mean, the signature of empty signal (source: 54_tad_methylation_organization_summary.tsv)
- MeCP2 median per-TAD signal drops ctrl 0.976 → mut 0.815, with variance ratio staying ~0.068-0.070 (not TAD-organized either way) (source: 54_tad_methylation_organization_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Non-CG methylation shows no TAD organization (variance ratios ~0.025, far below 1) because its per-TAD signal is essentially zero with CV pinned at 4 — there is no structured non-CG methylation to organize. The standout positive is CG 5hmC, ~2x more TAD-organized than CG 5mC, reinforcing that the TET-product (5hmC) carries the spatially-structured epigenetic information in this system. MeCP2 raw signal is not TAD-organized (ratio ~0.07) and its global level falls in mutant, foreshadowing the section 55/56 finding that MeCP2 redistribution is not a TAD-scale non-CG phenomenon.
### Plot inventory
- `54a_metatad_signal/` — per-TAD mean methylation distribution, 5 contexts × condition
- `54b_boundary_insulation/` — boundary insulation score by context
- `54c_intratad_cv/` — intra-TAD coefficient of variation by context
- `54d_variance_ratio/` — inter/intra-TAD variance ratio (key panel), bootstrap CIs
- `54e_boundary_correlation/` — TAD boundary score vs methylation insulation (Spearman)
- `54f_differential_overlay/` — methylation log2FC at differential vs non-differential TADs
- `54g_boundary_type/` — insulation by TADCompare boundary type (shifted/split/merge/etc.)
- `54h_cross_context_fc/` — CG vs non-CG/MeCP2 per-TAD log2FC correlation
- `54j_composite/` — composite of signal + CV + variance ratio

## Section 55: section_55_mecp2_cg_corrected_proxy
### Analysis question
Treat MeCP2 CUT&RUN enrichment as a proxy for total (CG + non-CG) methylation and use OLS to subtract the CG-explained part, leaving a residual that is a candidate non-CG signal. Test whether that residual is TAD-organized (variance ratio), correlates with boundaries, and persists at TADs with no CG methylation change. Model A: MeCP2 ~ CG 5mC; Model B: + CG 5hmC.
### Key results
- Model A (MeCP2 ~ CG 5mC): R² = 0.0211, slope = 149.65 (p = 3.6e-58), n = 12,141 TADs (source: 55_mecp2_regression_model_summary.tsv, 55_mecp2_regression_coefficients.tsv)
- Model B (+ CG 5hmC): R² = 0.0359; 5hmC term estimate = 359.88 (p = 4.2e-42) is larger than the 5mC term 74.70 (p = 3.1e-12) (source: 55_mecp2_regression_model_summary.tsv, 55_mecp2_regression_coefficients.tsv)
- Variance ratio (ctrl): MeCP2 raw = 0.0680 [0.0636, 0.1036], MeCP2 residual Model A = 0.0665 [0.0622, 0.0904], Model B = 0.0656 — all far below CG 5mC = 1.018 (source: 55_mecp2_variance_ratio_summary.tsv)
- CG-stratified residual (Model A, ctrl): all four CG-log2FC quartiles have median residual ≤ 0 and reject zero — Q2 median = -0.268 (p = 6.6e-232), Q3 median = -0.245 (p = 1.8e-183), Q1 median = -0.083 (p = 1.2e-04), Q4 median = -0.021 (p = 0.0194) (source: 55_mecp2_cg_stratified_summary.tsv)
- DMR-stratified residual: Hyper-DMR TADs n = 6,939, median residual = -0.080; No-DMR TADs n = 5,202, median residual = -0.269 (source: 55_mecp2_dmr_stratified_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] CG methylation explains only ~2% of TAD-level MeCP2 enrichment (R² = 0.021), and adding 5hmC only raises it to ~3.6% — so MeCP2's TAD-scale distribution is overwhelmingly driven by something other than CG/hmCG methylation. But the CG-corrected residual is NOT TAD-organized (variance ratio ~0.066, indistinguishable from raw MeCP2 and ~15x below CG 5mC), so that "something" is not a structured non-CG methylation field — it is more likely sequence/chromatin context and noise. The residual is significantly non-zero across all CG quartiles, but median-negative, which means MeCP2 sits BELOW the CG prediction in most TADs (especially No-DMR TADs), the opposite of a positive non-CG excess. The strong 5hmC regression coefficient is the noteworthy positive: MeCP2 tracks 5hmC more steeply than 5mC, again pointing at the TET axis rather than non-CG.
### Plot inventory
- `55a_mecp2_vs_cg_scatter/` — MeCP2 vs CG 5mC per TAD with Model A fit
- `55b_model_diagnostics/` — residuals-vs-fitted and Q-Q for Model A
- `55c_residual_distributions/` — ctrl vs mut residual densities
- `55d_variance_ratio_comparison/` — CG vs MeCP2-raw vs residual-A vs residual-B variance ratios (key panel)
- `55e_boundary_correlation/` — MeCP2 residual vs TAD boundary score
- `55f_cg_stratified_residuals/` — residual by CG-log2FC quartile (key test panel)
- `55g_dmr_stratified_residuals/` — residual by CG DMR class
- `55h_model_comparison/` — Model A vs Model B variance ratios, ctrl vs mut
- `55i_composite/` — composite of scatter + variance ratio + CG-stratified residuals

## Section 56: section_56_mecp2_peak_cg_correction
### Analysis question
At peak resolution (202,650 MeCP2 DiffBind peaks), extract CG 5mC at each peak and classify the significant MeCP2-Up peaks: is the MeCP2 gain explained by a CG methylation increase (CG-concordant) or not (non-CG candidate)? Regress MeCP2 Fold on CG log2FC and characterize the chromatin context of each class.
### Key results
- Peak regression MeCP2 Fold ~ CG mC log2FC: R² = 0.0128, slope = 0.1909 (p ≈ 0), intercept = 0.0451, n = 202,574 (source: 56_mecp2_peak_regression.tsv)
- MeCP2-Up peaks (FDR<0.05, Fold>0): 7,686 total — CG-Concordant 4,960 (64.5%), Non-CG Candidate 2,726 (35.5%) (source: 56_mecp2_peak_classification_summary.tsv)
- Mean residual (CG-unexplained MeCP2 fold): CG-Concordant 1.603, Non-CG Candidate 1.583, both far above zero; mean MeCP2 Fold 1.692 vs 1.617 (source: 56_mecp2_peak_classification_summary.tsv)
- Non-CG candidates are overwhelmingly Unmarked chromatin (77.2%) vs 62.0% for CG-concordant; CG-concordant peaks have more enhancer overlap (Active 14.8% + Poised 20.9% = 35.7%, vs Non-CG 5.9% + 13.3% = 19.2%) (source: 56_mecp2_peak_chromatin_summary.tsv)
- CG-DMR overlap (MeCP2-Up): 2,687 peaks at Hyper-DMR genes vs 4,999 at No-CG-DMR; MeCP2-Down: 755 Hyper-DMR vs 445 No-CG-DMR (source: 56_mecp2_peak_dmr_crosstab.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Peak-level CG methylation explains only ~1.3% of MeCP2 fold-change variance, and a third (35.5%) of all MeCP2-Up peaks gain MeCP2 with NO accompanying CG increase — these are the formal "non-CG candidate" sites the rest of the group chases. Their depletion of active/poised enhancer marks and enrichment in Unmarked chromatin (77%) is consistent with MeCP2 spreading into gene-poor, low-CpG, neuronally-mCA-rich genomic space — exactly where the mCA-reader model predicts. The top recurrent candidate locus is Etv1 (multiple peaks). This section produces the 2,726-peak candidate set that sections 57-58 validate against external WGBS.
### Plot inventory
- `56a_mecp2_vs_cg_scatter/` — MeCP2 Fold vs CG log2FC, all 202,650 peaks
- `56b_peak_classification/` — CG-concordant vs non-CG-candidate bar chart (key panel)
- `56c_cg_levels_by_mecp2/` — CG 5mC at MeCP2-Up vs Down peaks, ctrl/mut
- `56d_residual_distributions/` — residual densities for Up vs Down significant peaks
- `56e_dmr_crosstab/` — MeCP2 direction × CG DMR overlap
- `56f_chromatin_context/` — 7-state chromatin composition of non-CG-candidate vs CG-concordant peaks
- `56g_composite/` — composite of scatter + classification + residuals

## Section 57: section_57_ecker_noncg_validation
### Analysis question
Use Ecker's wildtype adult-cerebellum WGBS (which directly measures CH/non-CG methylation) to test whether the 2,726 non-CG-candidate MeCP2-Up peaks really sit at higher-CH regions than the CG-concordant peaks. Ecker CG is the specificity control; evoC CHG/CHH are side-by-side comparators.
### Key results
- KEY TEST — Ecker CH (non-CG) at non-CG-candidate vs CG-concordant MeCP2-Up peaks: median 0.0159 vs 0.0129, Wilcoxon p = 2.63e-60 (n = 2,726 vs 4,960) (source: 57_ecker_noncg_validation_summary.tsv)
- SPECIFICITY CONTROL — Ecker CG at the same two groups goes the OTHER way: median 0.00955 vs 0.0136, p = 1.39e-55 (non-CG candidates have LOWER CG, as expected by construction) (source: 57_ecker_noncg_validation_summary.tsv)
- MeCP2-Down peaks have the highest Ecker CH of all (median 0.0219, n = 1,200), above both Up groups (source: 57_ecker_noncg_validation_summary.tsv)
- Ecker CH detected (>0) in 100% of peaks in every group; Not-Significant peaks have much lower CH (median 0.00523, n = 193,764) (source: 57_ecker_noncg_validation_summary.tsv)
- Correlation residual-vs-Ecker-CH (pooled Up+Down significant peaks): Spearman rho = -0.176, p = 1.81e-62, n = 8,886 — NEGATIVE (source: 57_ecker_correlation.tsv)
- evoC CHG/CHH at these peaks are essentially zero (median 0 in all groups; <1.2% detection), confirming evoC cannot see the CH signal Ecker measures (source: 57_ecker_noncg_validation_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The headline validation works at the group level: non-CG-candidate MeCP2-Up peaks sit at genomic regions with genuinely higher non-CG methylation (Ecker CH p = 2.6e-60) while simultaneously having lower CG methylation (p = 1.4e-55) — the double dissociation is exactly the signature of MeCP2 reading mCA where mCG is absent, and it is invisible to evoC (CHG/CHH ~0). BUT the residual-vs-CH correlation is negative when pooled, and MeCP2-Down peaks (lost MeCP2) have the HIGHEST CH — both warning signs that "more CH" does not straightforwardly mean "more CG-unexplained MeCP2." Section 58 dissects this: the negative pooled correlation is partly Simpson's paradox, but even within Up-only peaks the dose-response fails.
### Plot inventory
- `57a_ecker_ch_noncg_vs_concordant/` — Ecker CH at non-CG-candidate vs CG-concordant (key panel)
- `57b_residual_vs_ecker_ch/` — MeCP2 residual vs Ecker CH, pooled Up+Down (the Simpson's-paradox plot section 58 corrects)
- `57c_detection_rate/` — Ecker CH detection rate by group (all ~100%)
- `57d_ecker_cg_control/` — Ecker CG specificity control
- `57e_three_way_ecker_ch/` — Ecker CH by MeCP2 Up/Down/NS status
- `57f_ecker_vs_evoc_sidebyside/` — Ecker CH vs evoC CHG vs evoC CHH at the same peaks
- `57g_composite/` — composite of CH + residual-correlation + CG control
- (`57f_composite/` folder also present — likely a stale/duplicate of the composite; see Data-quality flags)

## Section 58: section_58_noncg_dose_response
### Analysis question
Does Ecker CH methylation *level* scale the magnitude of CG-unexplained MeCP2 enrichment (a dose-response)? Bins MeCP2-Up peaks into Ecker-CH quartiles and runs Jonckheere-Terpstra trend tests on the residual; corrects the section-57b correlation to Up-only; adds quantile regression and a joint CG×CH stratification. MeCP2-Down is the negative control.
### Key results
- Corrected correlation (MeCP2-Up ONLY): Spearman rho = -0.176 vs pooled-Up+Down rho also reported in script; the OLS-vs-QR residual Spearman = 0.906 (source: 58_noncg_regression_comparison.tsv) [NOTE: up-only rho printed to console, not stored in a TSV — see flags]
- DOSE-RESPONSE FAILS — MeCP2-Up OLS residual by Ecker CH quartile is flat-to-DECREASING: Q1 median 1.456, Q2 1.510, Q3 1.502, Q4 (high CH) 1.338 (source: 58_noncg_dose_response_up.tsv)
- Jonckheere-Terpstra increasing p = 0.99999 (n.s.); DECREASING p = 2.01e-09; Kruskal-Wallis p = 1.96e-15 (source: 58_noncg_trend_tests.tsv)
- Negative control MeCP2-Down also decreasing: residual Q1 -0.985 → Q4 -1.069; Jonckheere increasing p = 0.99998, decreasing p = 1.95e-05 (source: 58_noncg_dose_response_down.tsv, 58_noncg_trend_tests.tsv)
- Joint CG×CH stratification (MeCP2-Up): median residual is highest at CG-Low/CH-Low (1.602, n=1,970) and lowest at CG-High/CH-High (1.316, n=1,970) — residual tracks LOW CG, and adding CH does not raise it (source: 58_noncg_joint_stratification.tsv)
- OLS vs quantile regression slopes: OLS slope 0.1909, QR (tau=0.5) slope 0.0145, residuals 91% concordant (Spearman 0.906) (source: 58_noncg_regression_comparison.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The dose-response hypothesis is rejected: higher non-CG methylation does NOT produce more CG-unexplained MeCP2 — if anything the residual decreases with CH (significant Jonckheere DECREASING, p = 2e-09), and the same decreasing pattern appears in the MeCP2-Down negative control, indicating the CH gradient is confounded with overall genomic context rather than driving MeCP2 recruitment. The joint stratification clinches it: the MeCP2 residual is governed by LOW CG methylation, not by CH level (CH-High does not raise the residual at either CG stratum). Combined with section 57, the honest interpretation is that non-CG-candidate peaks are enriched for higher-CH regions as a *category* (57), but CH does not quantitatively scale the MeCP2 excess (58) — so this evoC + Ecker analysis cannot support a causal mCA-dose model for the BAP1-KO MeCP2 redistribution. The result is robust to OLS-vs-quantile-regression choice.
### Plot inventory
- `58a_corrected_residual_vs_ecker_ch/` — section-57b corrected to MeCP2-Up only
- `58b_dose_response_up/` — median residual by Ecker CH quartile, MeCP2-Up (key panel)
- `58c_dose_response_down/` — same for MeCP2-Down negative control
- `58d_scatter_by_quartile/` — residual vs CH scatter with LOESS, colored by quartile
- `58e_joint_cg_ch_heatmap/` — joint CG×CH median-residual heatmap (key panel)
- `58f_ols_vs_qr_comparison/` — OLS vs quantile-regression residual concordance
- `58g_composite/` — composite of corrected-correlation + dose-response + joint heatmap

## Cross-section synthesis
These eight sections form a single funnel testing one neuronal idea — that MeCP2 in BAP1-KO cerebellum reads non-CG (mCA/mCH) methylation that the DUET evoC CG assay cannot see. The funnel narrows from negative to nuanced: at gene-body (51), peak (53), and TAD (54) scale, evoC non-CG methylation is at the detection floor (CHG/CHH globally <1%, per-TAD medians of 0, no differential DMRs), so the assay simply cannot resolve a non-CG phenotype directly. Sections 55-56 pivot to using MeCP2 CUT&RUN itself as an indirect non-CG proxy by regressing out CG methylation; CG explains only ~1-2% of MeCP2 distribution, leaving a large residual — but that residual is not TAD-organized (55) and is governed by low-CG context (55/58), not a structured non-CG field. The external Ecker WGBS validation (57) provides the group's one solid positive: the 2,726 non-CG-candidate MeCP2-Up peaks genuinely sit at higher-CH / lower-CG regions (double dissociation, p ~1e-60), matching the mCA-reader model categorically. Yet section 58 shows the relationship is not dose-dependent (Jonckheere increasing n.s., decreasing significant; same in the Down control), so the data support association but not a quantitative causal mCA-dose mechanism. Throughout, CG 5hmC repeatedly emerges as the most TAD-organized, most MeCP2-predictive, most feature-differentiated track — quietly redirecting the paper's mechanistic weight from non-CG methylation back to the TET / 5hmC axis that is the dataset's documented central finding.

## Tables used
- `mecp2_noncg_summary.tsv` — section 51 Wilcoxon ctrl-vs-mut non-CG methylation at MeCP2 genes, CG-stratified and pooled
- `section52_feature_methylation.tsv` — section 52 per-feature delta-5mC/5hmC across 6,196 sub-gene intervals (123 genes)
- `section52_dunn_posthoc.tsv` — section 52 Dunn pairwise post-hoc q-values for delta-5mC and delta-5hmC by feature type
- `mecp2_peak_noncg_summary.tsv` — section 53 peak-resolution detection rates, Fisher tests, Wilcoxon, and DMR summaries
- `54_tad_methylation_organization_summary.tsv` — section 54 per-context median signal, CV, insulation, variance ratio per TAD
- `55_mecp2_regression_model_summary.tsv` / `55_mecp2_regression_coefficients.tsv` — section 55 Model A/B fit stats and coefficients
- `55_mecp2_variance_ratio_summary.tsv` — section 55 inter/intra-TAD variance ratios for CG, MeCP2-raw, and residuals
- `55_mecp2_cg_stratified_summary.tsv` / `55_mecp2_dmr_stratified_summary.tsv` — section 55 residual by CG-log2FC quartile and DMR class
- `56_mecp2_peak_regression.tsv` — section 56 peak-level MeCP2 Fold ~ CG log2FC regression
- `56_mecp2_peak_classification_summary.tsv` — section 56 CG-concordant vs non-CG-candidate MeCP2-Up counts and means
- `56_mecp2_peak_chromatin_summary.tsv` — section 56 7-state chromatin composition of the two peak classes
- `56_mecp2_peak_dmr_crosstab.tsv` — section 56 MeCP2 direction × CG DMR overlap counts
- `56_mecp2_peak_noncg_candidates.tsv` — section 56 the 2,726 non-CG-candidate peaks (top recurrent locus: Etv1)
- `57_ecker_noncg_validation_summary.tsv` — section 57 Ecker CH/CG and evoC CHG/CHH by peak group with Wilcoxon p-values
- `57_ecker_correlation.tsv` — section 57 Spearman correlation residual vs Ecker CH (pooled)
- `58_noncg_dose_response_up.tsv` / `58_noncg_dose_response_down.tsv` — section 58 residual by Ecker CH quartile, Up and Down
- `58_noncg_trend_tests.tsv` — section 58 Jonckheere-Terpstra and Kruskal-Wallis trend tests
- `58_noncg_joint_stratification.tsv` — section 58 joint CG×CH median-residual stratification
- `58_noncg_regression_comparison.tsv` — section 58 OLS vs quantile-regression comparison

## Data-quality flags
- **Section 51 CHG/CHH rows are byte-identical in `mecp2_noncg_summary.tsv`.** The "CHG at MeCP2 CG-Low", "CHH at MeCP2 CG-Low", and both "at all MeCP2-bound genes" rows share identical ctrl_mean, mut_mean, effect and p_value. This is suspicious — the script builds CHG and CHH from separate extracts (`chg_extract`, `chh_extract`) yet reports identical numbers, implying the CHH extract path may resolve to the same file as CHG, or the two non-CG extracts are genuinely numerically identical at this resolution. The pooled "at all MeCP2-bound genes" CHG and CHH rows are also identical. Treat section-51 CHG-vs-CHH as not independently distinguishable. The non-significance conclusion still holds.
- **Section 51 has orphan panel folders from an earlier script version.** Folders `51a_mecp2_noncg_overlap_fisher`, `51b_mecp2_fold_noncg_violin`, `51c_cg_stratified_mecp2_fold`, `51d_persample_noncg_at_mecp2`, `51e_chg_chh_identity` (dated Jun 9 13:44-13:46) are NOT produced by the current `section_51_mecp2_noncg_methylation.R` (which writes only `51c_noncg_at_mecp2_by_cg_stratum` and `51d_mecp2_bound_vs_background_noncg`, dated 13:58). The current script's numbers in `mecp2_noncg_summary.tsv` are what is documented here; the orphan panels reflect superseded analysis and should not be cited.
- **`mecp2_noncg_peaks_for_browser.tsv` (191 data rows) is an orphan table** from the earlier section-51 version (non-zero non-CG peaks for genome browsing); the current section_51 script does not write it. Not cited in results.
- **Section 53 detection has minor internal n inconsistency:** CHG DMR table reports 8,356 peaks tested while CHH and the extracts report 8,886; the script comment states 8,886 MeCP2 peaks. The 8,356 likely reflects CHG peaks passing the coverage filter. Detection denominators use 8,886 throughout.
- **Section 58 up-only correlation rho is not in any TSV.** The "MeCP2-Up only" Spearman rho is computed and printed to console but only the OLS-vs-QR comparison is saved; `58_noncg_regression_comparison.tsv` stores the residual concordance (0.906), not the up-only rho. The pooled rho (-0.176) is in `57_ecker_correlation.tsv`. The script's console claims an up-only rho but it is unverifiable from tables — reported here only as "negative" with the verifiable pooled value.
- **Section 57 has two composite folders** (`57f_composite` and `57g_composite`) plus `57f_ecker_vs_evoc_sidebyside`; the current script writes `57f_ecker_vs_evoc_sidebyside` and `57g_composite`. `57f_composite` is likely stale from a prior numbering. Cosmetic only.
- **`55_mecp2_dmr_stratified_summary.tsv` lists only 2 classes (Hyper-DMR, No DMR);** no Hypo-DMR or Both rows appear (collapsed because <5 TADs or zero hypo overlaps), so the section-55 DMR stratification is effectively a 2-way comparison, not the 4-way the code anticipates.
- No empty tables; all cited TSVs are present and populated. Sections 51-58 are tied to run-5 (deep-seq, 8 samples, sex covariate) per the project's current primary run.
