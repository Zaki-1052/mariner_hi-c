## Section 15: 5hmC DMR Correlations with Chromatin Marks
**Key numbers:**
- hmC × MeCP2 predicted quadrant (hmC Down + MeCP2 Up) O/E = 0.97 (2,510/2,588); Fisher OR = 0.79, p = 1.49e-05 (source: hmc_mecp2_integration.tsv)
- hmC × ATAC: hmC Down + ATAC Down O/E = 1.37 (961/701), hmC Up + ATAC Up O/E = 1.41; Fisher OR = 0.095, p = 1.74e-116 (source: hmc_atac_integration.tsv)
- hmC × K119ub: hmC Down + K119ub Gained O/E = 1.12 (1,403/1,250), hmC Up + K119ub Lost O/E = 1.41; Fisher OR = 2.80, p = 1.12e-35 (source: hmc_k119ub_integration.tsv)
- mC reference O/E: ATAC 1.40/1.52, K119ub 1.30/1.55, MeCP2 0.98/0.98 (source: hmc_vs_mc_enrichment_comparison.tsv)

**What this shows:** The 5hmC perspective reproduces the mC couplings for ATAC and K119ub (strong directional O/E 1.1-1.6), confirming the coupling is real regardless of modality. MeCP2 O/E ≈ 1.0 in both perspectives, reinforcing that MeCP2 directional coupling is genuinely weak at gene level.

**Figures:**
- `15a_hmc_mecp2_heatmap` — 2×2 hmC × MeCP2 O/E
- `15b_hmc_atac_heatmap` — 2×2 hmC × ATAC O/E
- `15c_hmc_k119ub_heatmap` — 2×2 hmC × K119ub O/E
- `15d_mc_vs_hmc_enrichment_comparison` — predicted-quadrant O/E dot plot, mC vs hmC
