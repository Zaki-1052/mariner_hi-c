## Section 35: ATAC x Chromatin Features + Loop Anchor Interval Permutation (Confirmatory)
**Key numbers:**
- 38 pairwise regioneReloaded tests at ntimes=100,000 (35A ATAC x 6 ChIP marks=12; 35B ATAC x 7 states=14; 35C loop anchors x 6 features=12) (source: logs/5/biomodal_perm_47851858.out, lines 125-129)
- Loop input: 2,910 loops → 2,327 gained-loop anchors (1,723 up_in_mutant loops) and 1,520 lost-loop anchors (1,187 down_in_mutant loops) (source: logs/5/biomodal_perm_47851858.out, lines 95-97)
- Feature sets: MeCP2 Up 7,686 / Down 1,200; CTCF 32,487; H3K4me1 113,781; Bivalent 318 peaks (source: logs/5/biomodal_perm_47851858.out, lines 79-88)
- All three crosswise objects completed and were cached (source: logs/5/biomodal_perm_47851858.out, line 162)

**What this shows:** Permutation framework to validate the ATAC-vs-mark, ATAC-vs-chromatin-state, and loop-anchor-vs-feature Fisher tests (sections 13b/13c/13d/27/31) against genomic non-independence. The asymmetric anchor counts (more gained than lost) restate the BAP1 loop-gain architectural finding. All three permutation computations ran to completion, but the script crashed at figure 35d (extract_perm_results error), so heatmaps 35a/b/c exist while `permutation_35_comparison.tsv` and figure 35e were not produced. Quantitative confirmation of these relationships is provided by the gene-level section 37 (tests 37-09, 37-11 both Confirmed).

**Figures:**
- 35a_crosswise_atac_x_chip — 2x6 heatmap: ATAC direction x 6 ChIP marks
- 35b_crosswise_atac_x_chromatin_state — 2x7 heatmap: ATAC direction x chromatin states
- 35c_crosswise_loops_x_features — 2x6 heatmap: loop-anchor direction x chromatin features. Panels 35d/35e not generated (run crashed).
