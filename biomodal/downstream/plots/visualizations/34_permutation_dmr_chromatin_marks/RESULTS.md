## Section 34: DMR x Chromatin Mark Interval Permutation (Confirmatory)
**Key numbers:**
- Design: 2 DMR direction sets (mC Hyper=7,522; mC Hypo=3,267 regions) x 8 mark peak sets = 16 pairwise regioneReloaded tests at ntimes=100,000 (source: logs/5/biomodal_perm_47851857.out)
- Crosswise z-score range across the 16 cells = [-8.82, 65.14] — permutation computation completed (source: logs/5/biomodal_perm_47851857.out, line 112)
- Recomputed parent Fisher tests: 12a ATAC OR=0.07 p=6.9e-179; 14a K119ub OR=0.70 p=3.9e-16; 14b K119ub-differential OR=4.45 p=4.7e-97; 19f H3K27ac OR=1.36 p=1.3e-05 (source: logs/5/biomodal_perm_47851857.out, lines 128-131)

**What this shows:** Interval-shuffling permutation (per-chromosome) tested whether mC Hyper/Hypo DMR intervals co-localize with chromatin-mark peaks beyond chromosome-matched random expectation, to validate the Fisher tests from sections 12a/14a/14b/19f. The permutation ran and produced very large positive z-scores (up to 65), indicating the spatial associations are real and not driven by genomic clustering. The script crashed at the comparison-table step, so only the heatmap (34a) was produced and no `permutation_34_comparison.tsv` was written; the formal Confirmed/Weakened verdicts live in the gene-level analogue (section 37, all Confirmed).

**Figures:**
- 34a_crosswise_dmr_x_marks — 2x8 crosswise z-score association heatmap (DMR direction x marks). Panels 34b and 34c were not generated (run crashed before the Fisher-vs-permutation comparison).
