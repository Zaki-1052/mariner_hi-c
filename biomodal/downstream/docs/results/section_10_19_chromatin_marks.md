# Sections 10-19: 7-category chromatin state classification + CUT&RUN/ChIP mark integration (MeCP2, ATAC, H2AK119ub, H3K27ac) at DMRs

## Summary
This group asks where BAP1-KO DMRs fall in the chromatin landscape and whether the directional methylation changes (mC↑/hmC↓) are coupled to changes in four chromatin readouts: MeCP2 (CUT&RUN), ATAC-seq accessibility, H2AK119ub (the BAP1 substrate), and H3K27ac. It uses run-5 significant CpG DMRs (10,775 significant mC DMRs: 7,513 hypermethylated, 3,262 hypomethylated; and 11,484 significant hmC DMRs dominated by hmC-down: 9,521 decreased, 1,963 increased), classified into a 7-category chromatin-state system via histone-mark peak overlap, then intersected with condition-specific ChIP/CUT&RUN peak sets and continuous BigWig signal. The headline answers: (1) hypermethylation overwhelmingly lands at active promoters/enhancers (Active_Promoter is 93% hypermethylated) while hypomethylation lands at Repressed_Promoter/Polycomb; (2) ATAC and K119ub show strong, directionally-correct coupling at the DMR level (ATAC OR=0.068 i.e. hyper→ATAC-down; K119ub differential OR=4.46 hyper→gained), but the gene-level signal is much weaker once you account for the background gain rate (Section 17 "honest" re-assessment) and continuous signal (Section 18, Cliff's delta ≈ 0.10); (3) MeCP2 shows significant DMR-level overlap (OR=5.13) but essentially null gene-level directional coupling (OR≈0.92, p≈0.08).

## Section 10: section_10_chromatin_state
### Analysis question
Where do significant mC DMRs fall across a 7-category chromatin-state system (Active/Repressed/Bivalent promoter, Polycomb, Active/Poised enhancer, Unmarked) built from histone-mark peak overlaps + TSS distance, and does the distribution differ by methylation direction?
### Key results
- Total significant mC DMRs classified = 10,775 (7,513 hyper + 3,262 hypo) (source: chromatin_state_summary.tsv, summed)
- Active_Promoter = 4,906 DMRs total = 45.5% of all significant DMRs; among hyper it is 4,562/7,513 = 60.7% (source: chromatin_state_summary.tsv)
- Active_Promoter is 93.0% hypermethylated (4,562 hyper / 4,906 total) (source: chromatin_state_summary.tsv)
- Repressed_Promoter = 1,718 total and is 94.4% hypomethylated (1,621 hypo / 1,718) (source: chromatin_state_summary.tsv)
- Polycomb = 15 DMRs total, 14 hypomethylated (93%); Unmarked = 3,952 total = 36.7% of all DMRs (source: chromatin_state_summary.tsv)
- Mean mC change at Active_Promoter hyper = +3.59% vs hypo = −2.18%; Poised_Enhancer hyper = +4.21% (largest hyper effect) (source: chromatin_state_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] BAP1 loss drives gain of methylation specifically at active regulatory elements (promoters and enhancers), the genomic compartments where 5hmC normally accumulates and TET is most active — consistent with a TET-demethylation block rather than de-novo silencing of already-repressed Polycomb domains. The mirror-image partitioning (hyper→active, hypo→repressed/Polycomb) means a naive pooled enrichment would cancel; the directional split is the signal.
### Plot inventory
- `10a_chromatin_state_distribution` — bar (faceted by direction) + pie of overall state distribution
- `10b_chromatin_by_methylation_direction` — stacked + grouped bars comparing hyper vs hypo state composition
- `10c_chip_mark_overlap_heatmap` — 6-mark overlap % heatmap by direction (All/Hyper/Hypo)
- `10d_coordinated_genes_chromatin` — bar + pie of chromatin state for coordinated mC↑/hmC↓ genes
- `10e_top_genes_chromatin_annotation` — top-20 coordinated genes with state annotation
- `10f_chromatin_stacked_presentation` — presentation-grade stacked bar (state by direction)

## Section 10f: section_10f_expanded_chip_heatmap
### Analysis question
Extend the 6-mark overlap heatmap to 10 chromatin marks (adds H2AK119ub, ATAC, H3K36me2, H3K36me3 as condition-union peak sets) and quantify the % of significant mC DMRs overlapping each mark, split by direction.
### Key results
- H3K4me1 overlap = 91.0% of all significant DMRs (highest); ATAC = 90.3% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K27ac overlap differs sharply by direction: hyper = 64.0% vs hypo = 15.9% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K27me3 overlap is the inverse: hyper = 3.6% vs hypo = 55.4% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H2AK119ub overlap: all = 46.2%, hyper = 38.7%, hypo = 63.4% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K4me3 overlap: hyper = 63.1% vs hypo = 17.9%; CTCF ≈ even (~54% both directions) (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K36me2 = 11.0% and H3K36me3 = 7.3% overall (lowest overlaps) (source: 10f_expanded_chip_overlap_percentages.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The mark profile cleanly separates the two DMR classes: hypermethylated DMRs carry active marks (H3K27ac, H3K4me3, H3K4me1) and less K27me3/K119ub, whereas hypomethylated DMRs carry repressive marks (H3K27me3 55%, K119ub 63%). This is the chromatin signature of a TET block at active loci plus passive demethylation drift at Polycomb-repressed loci, and it foreshadows the directionally opposite ATAC/K119ub couplings quantified in Sections 12-19.
### Plot inventory
- `10f_expanded_chip_heatmap` — 10-mark × 3-row (All/Hyper/Hypo) overlap % heatmap

## Section 11: section_11_mecp2_correlation
### Analysis question
Do mC DMR changes correlate with MeCP2 (CUT&RUN) binding changes? MeCP2 binds 5mCG with high and 5hmCG with low affinity (Mellen 2017), so hypermethylated DMRs are predicted to gain MeCP2. Tested at DMR-overlap, gene-fold, coordinated-gene, binary-gain, stratified, integration, and delta-ratio-regression levels.
### Key results
- DMR-overlap Fisher OR = 5.13, p = 1.27e-33 (hyper DMRs preferentially overlap MeCP2-Up peaks) (source: mecp2_dmr_overlap_summary.tsv)
- Hyper DMRs: 7.16% overlap MeCP2-Up (538/7,513) vs 1.89% MeCP2-Down; hypo DMRs: 5.09% MeCP2-Up vs 6.90% MeCP2-Down (source: mecp2_dmr_overlap_summary.tsv)
- Binary MeCP2 significant-gain enrichment, coordinated vs non-coordinated genes: OR = 1.04, p = 0.795 (NOT significant); coordinated 48/6,069 = 0.79% gain vs non-coordinated 115/15,126 = 0.76% (source: mecp2_binary_gain_coordinated_fisher.tsv)
- Stratified gene-level Spearman (mC% vs MeCP2 nearest fold): Coordinated (n=6,069) rho=+0.024 p=0.065; Discordant (n=1,695) rho=+0.020 p=0.40; mC-only (n=2,060) rho=+0.076 p=5.5e-04 (source: mecp2_stratified_correlations.tsv)
- Delta-ratio linear model (MeCP2 fold ~ delta_ratio + coverage): delta_ratio coefficient = −0.367, p = 1.42e-06 (negative = consistent with prediction); coverage coefficient = +0.0052, p = 1.17e-12 (source: mecp2_delta_ratio_lm_coefficients.tsv)
- Delta-ratio logistic model: delta_ratio OR = 7.20e-07 (≪1, p = 1.97e-05) → strong TET impairment predicts higher P(MeCP2 sig up) (source: mecp2_delta_ratio_glm_odds_ratios.tsv)
- Coordinated mC↑/hmC↓ genes with MeCP2 data = 6,069; of these MeCP2 nearest-fold is up at 1,769 and down at 4,300 (source: mecp2_coordinated_genes.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] MeCP2 shows the predicted enrichment at the level of *which DMRs sit under differential MeCP2 peaks* (OR=5.13), but once aggregated to gene level the directional coupling collapses (stratified rho ≈ 0.02-0.08; binary-gain Fisher NS). The delta-ratio regressions rescue a weak but significant signal: genes with greater TET-efficiency loss have measurably higher MeCP2 occupancy, fitting "more retained 5mC → more MeCP2 recruited." The net reading is that MeCP2 redistribution is a real but quantitatively modest downstream consequence, not a tight 1:1 readout of methylation change.
### Plot inventory
- `11a_mecp2_overlap` — MeCP2-Up/Down overlap % by DMR direction + Fisher
- `11b_mecp2_fold_by_dmr_direction` — MeCP2 fold violin/box by hyper/hypo/NS
- `11c_mc_vs_mecp2_scatter` — gene-level mC% vs MeCP2 fold scatter with quadrant counts
- `11d_mecp2_coordinated_genes` — coordinated-gene MeCP2 box + top-20 bar
- `11e_mecp2_integration_heatmap` — 2×2 mC-direction × MeCP2-direction O/E heatmap (figure-only, no table)
- `11f_mecp2_delta_ratio_lm` — linear-model forest + partial-regression
- `11g_mecp2_delta_ratio_glm` — logistic odds-ratio forest + predicted-probability curve
- `11h_mecp2_binary_gain_fisher` — binary MeCP2-gain enrichment bar
- `11i_mc_vs_mecp2_stratified` — scatter stratified by coordination status

## Section 12: section_12_atac_correlation
### Analysis question
Do mC DMR changes correlate with chromatin-accessibility (ATAC-seq) changes? Prediction (Conway 2021): hyper DMRs lose accessibility (ATAC-down), hypo DMRs gain accessibility (ATAC-up).
### Key results
- DMR-overlap Fisher OR = 0.068, p = 4.44e-178 (strong directional coupling) (source: atac_dmr_overlap_summary.tsv)
- Hyper DMRs: 13.4% overlap ATAC-Down (1,006/7,513) vs 10.2% ATAC-Up; hypo DMRs: 33.4% overlap ATAC-Up (1,088/3,262) vs 2.9% ATAC-Down (source: atac_dmr_overlap_summary.tsv)
- Coordinated mC↑/hmC↓ genes annotated = 6,589; top coordinated genes (e.g. Syt1 net_atac=−10, Gclm net_atac=−5) trend ATAC-down (source: atac_coordinated_genes.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Accessibility is the most cleanly coupled of the four marks at the DMR level: hypermethylated loci close and hypomethylated loci open, exactly as predicted if K119ub-driven compaction blocks TET. The extreme p-value reflects the large directional asymmetry (hypo→ATAC-up is 11× the hypo→ATAC-down rate). This is the chromatin-state correlate of the methylation switch.
### Plot inventory
- `12a_atac_overlap` — ATAC-Up/Down overlap % by DMR direction + Fisher
- `12b_consensus_accessibility` — consensus ctrl-vs-mut accessibility at DMRs
- `12c_atac_coordinated_genes` — ATAC-down box + top-20 coordinated-gene bar
- `12d_mc_vs_atac_scatter` — gene-level mC% vs net-ATAC scatter
- `12e_atac_integration_heatmap` — 2×2 mC × ATAC O/E heatmap

## Section 13: section_13_atac_chromatin_and_loops
### Analysis question
Part A: are differential ATAC peaks enriched in specific chromatin states (especially Polycomb/Repressed_Promoter)? Part B: do Hi-C differential loop anchors show ATAC changes concordant with loop direction (up-in-mutant→ATAC-up, down→ATAC-down)?
### Key results
- ATAC-Up peaks total = 7,620; ATAC-Down total = 3,744 (source: atac_chip_overlap_enrichment.tsv)
- H3K27me3 mark overlap: ATAC-Up = 23.6% vs ATAC-Down = 1.0%, OR (Down/Up) = 0.034, padj < 0.001 *** (source: atac_chip_overlap_enrichment.tsv)
- H3K27ac mark overlap: ATAC-Up = 12.3% vs ATAC-Down = 40.2%, OR = 4.79, padj < 0.001 *** (Down enriched for active enhancers) (source: atac_chip_overlap_enrichment.tsv)
- Chromatin-state distribution: Polycomb is 12.7% of ATAC-Up vs 0.5% of ATAC-Down; Active_Enhancer is 35.8% of ATAC-Down vs 8.6% of ATAC-Up (source: atac_chromatin_state_distribution.tsv)
- Loop-anchor ATAC: 2,910 differential loops (1,723 up-in-mutant, 1,187 down-in-mutant); loops with any-anchor ATAC-Up = 950, any-anchor ATAC-Down = 406 (source: loop_anchor_atac_overlap.tsv)
- Loop-type concordance: Active_Enhancer–Active_Enhancer loops 79.6% concordant (113/142), CTCF_Site–CTCF_Site only 15.4% (52/338) (source: loop_atac_concordance_by_type.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] ATAC-up peaks are unexpectedly enriched at Polycomb/H3K27me3 (12.7% vs 0.5%), suggesting that in the mutant some Polycomb-repressed regions paradoxically open — a candidate for K119ub-driven redistribution of repression rather than uniform compaction. The loop analysis links accessibility to 3D architecture: enhancer-enhancer loops change accessibility concordantly with loop direction (~80%), whereas structural CTCF loops do not (~15%), implying the regulatory (not structural) loop class carries the methylation/accessibility coupling.
### Plot inventory
- `13a_atac_chromatin_state_distribution` — state composition of ATAC-Up vs ATAC-Down
- `13b_atac_chip_overlap_enrichment` — per-mark overlap % + Fisher (Down/Up)
- `13c_atac_chromatin_enrichment_heatmap` — O/E heatmap ATAC-direction × state
- `13d_loop_anchor_atac_overlap` — ATAC overlap at loop anchors by loop direction
- `13e_anchor_atac_by_chromatin_state` — anchor ATAC overlap by 8-category anchor state
- `13f_loop_atac_concordance_by_type` — concordance % per loop type

## Section 14: section_14_h2ak119ub_correlation
### Analysis question
Do mC DMRs correlate with H2AK119ub (the BAP1 PR-DUB substrate)? Prediction (Conway 2021): hyper DMRs gain K119ub (mut), hypo DMRs lose it. Differential peaks are derived by GRanges set operations on ctrl/mut condition-specific peaks.
### Key results
- Condition-overlap (ctrl vs mut) Fisher OR = 0.700, p = 4.55e-16 (source: k119ub_condition_dmr_overlap_summary.tsv)
- Hyper DMRs: 35.4% overlap mutant K119ub vs 27.8% control (gain in mutant); hypo DMRs: 54.1% mut vs 60.6% ctrl (loss in mutant) (source: k119ub_condition_dmr_overlap_summary.tsv)
- Differential (gained/lost) Fisher OR = 4.46, p = 5.80e-97; derived from 6,172 gained and 6,103 lost peaks (source: k119ub_differential_dmr_overlap_summary.tsv)
- Hyper DMRs: 19.6% overlap K119ub-Gained (1,469/7,513) vs 10.4% Lost; hypo DMRs: 12.2% Gained vs 29.0% Lost (946/3,262) (source: k119ub_differential_dmr_overlap_summary.tsv)
- Coordinated mC↑/hmC↓ genes annotated = 6,589; top genes (e.g. Syt1 net_k119ub=+15, Gclm=+3) trend K119ub-gained (source: k119ub_coordinated_genes.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] At the DMR level the predicted convergence holds strongly: hypermethylated sites gain K119ub and hypomethylated sites lose it (differential OR=4.46), the central evidence that K119ub accumulation upstream of TET drives the methylation switch. The condition-overlap OR<1 simply reflects that hypo DMRs start with very high baseline K119ub (60.6% ctrl). Section 17 then re-examines whether this DMR-level signal survives honest gene-level accounting.
### Plot inventory
- `14a_k119ub_condition_overlap` — ctrl-vs-mut K119ub overlap % by DMR direction
- `14b_k119ub_differential_overlap` — gained/lost K119ub overlap % by direction
- `14c_k119ub_coordinated_genes` — K119ub-mut box + top-20 coordinated bar
- `14d_mc_vs_k119ub_scatter` — gene-level mC% vs net-K119ub scatter
- `14e_k119ub_integration_heatmap` — 2×2 mC × K119ub O/E heatmap

## Section 15: section_15_hmc_chromatin_correlations
### Analysis question
Re-run the MeCP2/ATAC/K119ub integration using 5hmC DMR direction (complement to mC). Because mC↑/hmC↓ are 92.3% coordinated, the predicted quadrants flip (hmC-Down ↔ MeCP2-Up/ATAC-Down/K119ub-Gained). Do hmC O/E enrichments match or exceed mC?
### Key results
- hmC × MeCP2 predicted quadrant (hmC Down + MeCP2 Up) O/E = 0.97 (Obs 2,510 / Exp 2,588); whole-table Fisher OR = 0.79, p = 1.49e-05 (source: hmc_mecp2_integration.tsv; hmc_vs_mc_enrichment_comparison.tsv)
- hmC × ATAC predicted quadrants: hmC Down + ATAC Down O/E = 1.37 (961/701); hmC Up + ATAC Up O/E = 1.41 (889/629); Fisher OR = 0.095, p = 1.74e-116 (source: hmc_atac_integration.tsv)
- hmC × K119ub predicted quadrants: hmC Down + K119ub Gained O/E = 1.12 (1,403/1,250); hmC Up + K119ub Lost O/E = 1.41 (520/368); Fisher OR = 2.80, p = 1.12e-35 (source: hmc_k119ub_integration.tsv)
- mC reference O/E for same marks: mC × ATAC = 1.40/1.52, mC × K119ub = 1.30/1.55, mC × MeCP2 = 0.98/0.98 (source: hmc_vs_mc_enrichment_comparison.tsv)
- Across-mark verdict: ATAC and K119ub give strong O/E (1.1-1.6) for both mC and hmC; MeCP2 gives ≈1.0 (no enrichment) for both perspectives (source: hmc_vs_mc_enrichment_comparison.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] The hmC perspective reproduces the mC story for ATAC and K119ub (strong directional O/E), confirming the coupling is real regardless of which modality you anchor on — consistent with mC and hmC being two faces of the same TET-block. MeCP2's O/E hovers at 1.0 in both perspectives, reinforcing that MeCP2 directional coupling is genuinely weak at gene level (its DMR-level OR=5.13 in Section 11 is driven by peak co-location, not net redistribution). The slightly higher hmC-Up + K119ub-Lost O/E (1.41) is a small set but directionally clean.
### Plot inventory
- `15a_hmc_mecp2_heatmap` — 2×2 hmC × MeCP2 O/E
- `15b_hmc_atac_heatmap` — 2×2 hmC × ATAC O/E
- `15c_hmc_k119ub_heatmap` — 2×2 hmC × K119ub O/E
- `15d_mc_vs_hmc_enrichment_comparison` — dot plot of predicted-quadrant O/E, mC vs hmC, all 3 marks

## Section 16: section_16_raw_concordance_barcharts
### Analysis question
Companion to Section 15 O/E heatmaps: show raw concordance rates (% of genes in each methylation group displaying the predicted mark direction) so the dominant-group effect sizes are visible rather than normalized away.
### Key results
- MeCP2 concordance: mC Up = 30.0% MeCP2-Up (2,058/6,862), mC Down = 68.2% MeCP2-Down (2,020/2,962); hmC Down = 29.2% MeCP2-Up, hmC Up = 65.6% MeCP2-Down (source: raw_concordance_all_marks.tsv)
- ATAC concordance: mC Up = 54.9% ATAC-Down (938/1,710), mC Down = 89.7% ATAC-Up (1,119/1,248); hmC Down = 47.9% ATAC-Down, hmC Up = 91.9% ATAC-Up (source: raw_concordance_all_marks.tsv)
- K119ub concordance: mC Up = 68.3% Gained (1,358/1,988), mC Down = 73.7% Lost (889/1,207); hmC Down = 62.3% Gained, hmC Up = 63.0% Lost (source: raw_concordance_all_marks.tsv)
- Dominant-group summary (mC Up vs hmC Down): MeCP2 30.0% vs 29.2%; ATAC 54.9% vs 47.9%; K119ub 68.3% vs 62.3% (source: raw_concordance_summary.tsv)
- The "predicted-positive" group rates (mC Up, hmC Down) are consistently lower than their mirror groups (mC Down, hmC Up), because the down/up mirror groups are small and cleanly concordant (source: raw_concordance_all_marks.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Raw rates reveal that the headline "majority concordant" claim for the dominant mC-Up group is real for K119ub (68%) and ATAC (55%) but barely above chance for MeCP2 (30% predicted-up, i.e. 70% go the other way). The small mirror groups (mC Down→ATAC-Up at 90%) are the most internally consistent, but they involve far fewer genes. This section is the honesty check that motivates Section 17.
### Plot inventory
- `16a_mecp2_concordance_bars` — MeCP2 concordance % across 4 methylation groups
- `16b_atac_concordance_bars` — ATAC concordance % across 4 groups
- `16c_k119ub_concordance_bars` — K119ub concordance % across 4 groups
- `16d_raw_concordance_comparison` — mC-Up vs hmC-Down summary across all 3 marks

## Section 17: section_17_k119ub_honest_assessment
### Analysis question
Honest re-assessment correcting Sections 14/16c framing: include the "No Peaks" category (most genes have no K119ub at all) and report the conditional gain direction *relative to the background gain rate among DMR genes with peaks*. Does the apparent hyper→K119ub-gain enrichment survive?
### Key results
- The majority of genes in every group have NO K119ub peak: mC Up 62.2% no-peaks (4,674/7,513), hmC Down 64.6% (6,150/9,521), mC Down 39.7%, hmC Up 34.5% (source: k119ub_honest_breakdown.tsv)
- Conditional gain among genes WITH peaks: mC Up = 47.8% gained (1,358/2,839), hmC Down = 41.6% (1,403/3,371), mC Down = 16.2%, hmC Up = 23.8% (source: k119ub_honest_breakdown.tsv)
- Background gain rate across all DMR genes with peaks = 33.6% (1,954/5,811) (source: k119ub_honest_breakdown.tsv)
- mC Up gain enrichment above background = +14.2 pp (47.8% vs 33.6%); hmC Down = +8.0 pp (41.6% vs 33.6%) (source: k119ub_honest_breakdown.tsv, derived)
- Background no-peak rate = 58.2% (8,077/13,888 DMR genes have no K119ub peak in either condition) (source: k119ub_honest_breakdown.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] This corrects the Section 14/16c picture: while hypermethylated DMRs are genuinely enriched for K119ub gain (+14 pp over background among peak-bearing genes), the effect is modest and applies to a minority of genes — most DMR genes carry no K119ub at all. The 68.3% "mC Up → K119ub Gained" figure cited from Section 16c is misleading because it excludes the no-peak majority and omits the ~34% background gain rate. The corrected reading: K119ub accumulation is a real but partial contributor to hypermethylation, not a universal gatekeeper.
### Plot inventory
- `17a_k119ub_full_breakdown` — 100% stacked bars including No-Peaks/Equal/Gained/Lost
- `17b_k119ub_conditional_direction` — gain/equal/lost among peak-bearing genes with background reference lines + Fisher vs background
- `17c_k119ub_effect_sizes` — per-gene net-peak strip/box (most changes ±1 peak)

## Section 18: section_18_k119ub_bigwig_signal
### Analysis question
Replace the binary peak-count approach (Section 17) with continuous K119ub BigWig signal: per-gene mean signal, log2(mut/ctrl), and correlation with methylation change. Does continuous signal confirm a genome-wide K119ub increase and is the gene-specific effect strong?
### Key results
- mC Up median K119ub log2FC = +0.0586 (n=6,993), p(vs 0) = 4.67e-96; Cliff's delta vs background = +0.099 (negligible) (source: k119ub_bigwig_signal_summary.tsv)
- mC Down median log2FC = −0.0797 (n=2,974), p(vs 0) = 2.12e-41; Cliff's delta = −0.181 (small) (source: k119ub_bigwig_signal_summary.tsv)
- hmC Down median log2FC = +0.0141 (n=8,757), p(vs 0) = 1.42e-19, but p(vs background) = 0.299 (NS) and Cliff's delta = +0.008 (negligible) (source: k119ub_bigwig_signal_summary.tsv)
- Background (All DMR genes) median log2FC = +0.0070 (n=12,721), p(vs 0) = 1.81e-20 — confirms a global K119ub increase (source: k119ub_bigwig_signal_summary.tsv)
- hmC Up median log2FC = +0.0093 (n=1,842), p(vs background) = 0.770 (NS); Cliff's delta = −0.004 (source: k119ub_bigwig_signal_summary.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] Continuous signal confirms the BAP1-KO expectation of a genome-wide K119ub increase (background median log2FC > 0, p<1e-19) and that hypermethylated (mC Up) genes carry a directionally-correct positive shift while hypomethylated genes shift negative. But the gene-specific effect sizes are tiny (Cliff's delta 0.10-0.18, mostly "negligible"), and hmC-Down genes are statistically indistinguishable from background (p=0.30). The mechanism is a diffuse global K119ub rise, not sharp locus-specific redistribution at DMR genes — the most quantitatively cautious statement of the K119ub-upstream model.
### Plot inventory
- `18a_k119ub_signal_distribution` — ctrl-vs-mut signal violin/box per group (paired Wilcoxon)
- `18b_k119ub_log2fc_distributions` — log2FC violin/box per group vs background
- `18c_methylation_vs_k119ub_scatter` — mod_difference vs K119ub log2FC (Spearman, 5mC + 5hmC facets)
- `18d_k119ub_waterfall` — all genes ranked by K119ub log2FC, DMR genes colored

## Section 19: section_19_h3k27ac_peak_analysis
### Analysis question
Comprehensive H3K27ac condition-specific peak analysis at DMR genes (parallel to K119ub Sections 16-18), 8 panels: status breakdown, conditional direction, effect size, methylation scatter, waterfall, DMR overlap, 2×2 O/E heatmaps, and a 4-mark O/E comparison adding H3K27ac to MeCP2/ATAC/K119ub.
### Key results
- No-peak rates: mC Up 44.2% no-peaks (3,319/7,513), hmC Down 51.7%, mC Down 80.3%, hmC Up 80.9% (source: h3k27ac_peak_analysis_summary.tsv)
- Conditional gain among peak-bearing genes: mC Up = 31.5% gained (1,320/4,194), hmC Down = 34.5%, mC Down = 52.2% (335/642), hmC Up = 52.4%; background gain rate = 35.4% (source: h3k27ac_peak_analysis_summary.tsv)
- Fisher gain-vs-background: mC Down p = 3.20e-16 (OR>1, enriched gain), mC Up p = 4.58e-05, hmC Down p = 0.317 (NS), hmC Up p = 9.03e-11 (source: h3k27ac_peak_analysis_summary.tsv)
- Cliff's delta vs background: mC Down = +0.214 (small, enriched H3K27ac gain), mC Up = −0.052, hmC Up = +0.164, hmC Down = −0.009 (source: h3k27ac_peak_analysis_summary.tsv)
- 2×2 O/E (19g): mC × H3K27ac predicted-enriched cells mC Down + Gained O/E = 1.34 (335/251), mC Up + Lost O/E = 1.11; whole-table Fisher OR = 0.254, p = 2.03e-24 (source: h3k27ac_all_marks_oe_comparison.tsv)
- 4-mark O/E comparison (19h): mC-perspective enriched-quadrant O/E — ATAC 1.40/1.52, K119ub 1.30/1.55, H3K27ac 1.34/1.11, MeCP2 1.04/1.01 (source: h3k27ac_all_marks_oe_comparison.tsv)
### [INTERPRETATION] Biological meaning
[INTERPRETATION] H3K27ac behaves opposite to the simple "active mark lost at hyper" expectation: it is *hypomethylated* (mC Down) genes that significantly gain H3K27ac (OR>1, Cliff's delta +0.21), consistent with hypomethylated loci becoming more active enhancers. Hypermethylated genes show no significant H3K27ac change above background (Cliff's delta −0.05), so loss of activation is not a strong correlate of hypermethylation in this dataset. In the 4-mark ranking, ATAC and K119ub remain the strongest directional couplings; H3K27ac is intermediate and MeCP2 weakest — placing accessibility and the BAP1 substrate at the center of the mechanism.
### Plot inventory
- `19a_h3k27ac_status_breakdown` — 100% stacked status bars
- `19b_h3k27ac_conditional_direction` — gain/equal/lost vs background + Fisher
- `19c_h3k27ac_effect_sizes` — per-gene net-peak strip/box
- `19d_methylation_vs_h3k27ac_scatter` — mod_difference vs net-H3K27ac (Spearman)
- `19e_h3k27ac_waterfall` — all genes ranked by net H3K27ac change
- `19f_h3k27ac_condition_overlap` — ctrl-vs-mut H3K27ac overlap at DMRs + Fisher
- `19g_h3k27ac_oe_heatmaps` — 2×2 mC × H3K27ac and hmC × H3K27ac O/E
- `19h_chromatin_mark_oe_comparison` — 4-mark (MeCP2/ATAC/K119ub/H3K27ac) O/E dot plot

## Cross-section synthesis
These ten sections build the central chromatin-mechanism argument of the paper. Section 10/10f establish *where* the methylation switch occurs — hypermethylation at active promoters/enhancers (carrying H3K27ac, H3K4me3), hypomethylation at Polycomb/Repressed promoters (carrying H3K27me3, high baseline K119ub). Sections 11-14 then test the four candidate chromatin couplings at DMR resolution, finding that ATAC accessibility (OR=0.068) and differential K119ub (OR=4.46) are the strongest and most directionally-correct, MeCP2 is significant only by peak co-location (OR=5.13 at DMRs but ≈1.0 at gene level), and the marks separate cleanly by direction. Section 15/16 confirm the same couplings hold from the hmC perspective (ATAC and K119ub O/E 1.1-1.6; MeCP2 ≈1.0). Crucially, Sections 17-18 are the "honest" correction layer for K119ub: when you include the no-peak majority and a continuous-signal background, the K119ub effect shrinks to a modest +14 pp gain enrichment and a negligible-to-small Cliff's delta atop a genuine *global* K119ub increase — i.e., K119ub is upstream and elevated everywhere, but its gene-specific DMR coupling is diffuse rather than sharp. Section 19 extends the framework to H3K27ac and ranks all four marks, leaving ATAC + K119ub as the load-bearing couplings. Together they support the thesis that BAP1 loss raises H2AK119ub, compacts active chromatin (ATAC-down), and blocks TET-mediated demethylation (mC↑/hmC↓) preferentially at active regulatory elements — with MeCP2 redistribution and H3K27ac change as weaker, partly opposite, downstream effects.

## Tables used
- `chromatin_state_summary.tsv` — per-state hyper/hypo counts and mean/median mC change (Section 10)
- `dmr_chromatin_state_annotation.tsv` — full per-DMR chromatin-state annotation (10,775 sig rows; Section 10)
- `10f_expanded_chip_overlap_percentages.tsv` — 10-mark overlap % by direction (Section 10f)
- `mecp2_dmr_overlap_summary.tsv` — MeCP2-Up/Down overlap counts + Fisher (Section 11a)
- `mecp2_gene_level_correlation.tsv` — per-gene mC vs MeCP2 nearest-fold (9,824 sig genes; Section 11c)
- `mecp2_coordinated_genes.tsv` — coordinated genes with MeCP2 annotation (6,069; Section 11d)
- `mecp2_binary_gain_coordinated_fisher.tsv` — binary MeCP2-gain Fisher (Section 11h)
- `mecp2_stratified_correlations.tsv` — per-coordination-group Spearman (Section 11i)
- `mecp2_delta_ratio_lm_coefficients.tsv` — linear model coefficients (Section 11f)
- `mecp2_delta_ratio_glm_odds_ratios.tsv` — logistic odds ratios (Section 11g)
- `atac_dmr_overlap_summary.tsv` — ATAC-Up/Down overlap + Fisher (Section 12a)
- `atac_coordinated_genes.tsv` — coordinated genes with net-ATAC (Section 12c)
- `atac_chromatin_state_distribution.tsv` — ATAC-direction × chromatin-state counts (Section 13a)
- `atac_chip_overlap_enrichment.tsv` — per-mark ATAC-Up/Down overlap + Fisher (Section 13b)
- `loop_anchor_atac_overlap.tsv` — per-loop anchor ATAC overlap (2,910 loops; Section 13d)
- `loop_atac_concordance_by_type.tsv` — concordance % per loop type (Section 13f)
- `k119ub_condition_dmr_overlap_summary.tsv` — ctrl/mut K119ub overlap + Fisher (Section 14a)
- `k119ub_differential_dmr_overlap_summary.tsv` — gained/lost K119ub overlap + Fisher (Section 14b)
- `k119ub_coordinated_genes.tsv` — coordinated genes with net-K119ub (Section 14c)
- `hmc_mecp2_integration.tsv` / `hmc_atac_integration.tsv` / `hmc_k119ub_integration.tsv` — hmC 2×2 O/E (Section 15a-c)
- `hmc_vs_mc_enrichment_comparison.tsv` — predicted-quadrant O/E, mC vs hmC, 3 marks (Section 15d)
- `raw_concordance_all_marks.tsv` / `raw_concordance_summary.tsv` — raw concordance rates (Section 16)
- `k119ub_honest_breakdown.tsv` — no-peak/gained/equal/lost breakdown + background (Section 17)
- `k119ub_bigwig_signal_summary.tsv` — continuous log2FC, p-values, Cliff's delta (Section 18)
- `h3k27ac_peak_analysis_summary.tsv` — H3K27ac status/gain/Cliff's delta/Fisher (Section 19a-c)
- `h3k27ac_all_marks_oe_comparison.tsv` — 4-mark × 2-perspective O/E (Section 19g-h)

## Data-quality flags
- **TODO.md row 10 says "49.9% DMRs at Active_Promoter."** The table gives 45.5% (4,906/10,775 total) or 60.7% among hyper (4,562/7,513). The 49.9% in TODO.md is not reproduced by `chromatin_state_summary.tsv`. [DISCREPANCY: TODO.md 49.9% vs table 45.5%.]
- **FIGURES.md Figure 15 (11e mC × MeCP2 integration) reports Obs counts 651/1,370/1,830/4,240 with OR=0.91, p=8.42e-02.** Section 11e has no exported table (figure-only), but the run-5 `h3k27ac_all_marks_oe_comparison.tsv` mC × MeCP2 cells are 942 (mC Down + MeCP2 Up) and 4,804 (mC Up + MeCP2 Down) with OR=0.9186, p=0.0773. The OR/p are consistent (~0.92, ~0.08) but the per-cell Obs counts in FIGURES.md do not match the current run-5 tables — likely a prior-run figure description. [UNVERIFIED: 651/1,830 Obs per FIGURES.md, not confirmed in tables.]
- **`atac_chromatin_state_distribution.tsv` uses "Other" not "Unmarked"** for the no-mark category and contains no CTCF_Site row, whereas the DMR chromatin-state system (Section 10) uses "Unmarked." Cosmetic label divergence between the ATAC-peak classification (Section 13a) and the DMR classification (Section 10); both are 7-category but with different residual-bin names.
- **Loop counts: `loop_anchor_atac_overlap.tsv` has 2,910 loops (1,723 up / 1,187 down).** Section 13 console/TODO sometimes cite the same 2,910 "differential loop anchors" — consistent. The per-loop any_atac_up=950 / any_atac_down=406 totals pool both directions; direction-specific concordant counts (up-loop any_atac_up=779; down-loop any_atac_down=370) are recomputed here from the table and are internally consistent.
- **K119ub peak-vs-gene totals differ by construction.** Section 14 derives 6,172 gained / 6,103 lost *peaks*; Section 17 reports 1,954 gained / 1,830 lost *genes* with a 33.6% background gain rate. These are not contradictory — one counts peaks, the other counts genes after ChIPseeker annotation. Worth stating explicitly in the manuscript to avoid apparent conflict.
- **MeCP2 is CUT&RUN, not ChIP** — all MeCP2 references above use "signal"/"mark"/"binding"; the column `mecp2_nearest_fold` is the CUT&RUN differential-binding log2FC.
- **Section 18 requires `data/k119ub_gene_signal.tsv`** (23,150 genes), present and dated 2026-02-09 — confirms BigWig preprocessing ran. No empty/missing tables encountered in this group.
