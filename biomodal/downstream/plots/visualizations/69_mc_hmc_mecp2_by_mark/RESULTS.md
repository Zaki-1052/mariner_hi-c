## Section 69: Does MeCP2 Distinguish 5mC from 5hmC? (by mark)

**Key numbers:**
- Active (K27ac): 5mC rho=0.160 vs 5hmC rho=-0.060; interaction p=7.28e-38 (significant); Fisher z p=1.17e-33 (69_mc_hmc_mecp2_by_mark_stats.tsv)
- Bivalent (K27me3+K27ac): 5mC rho=0.305 vs 5hmC rho=0.032 (NS); interaction p=4.06e-05 (significant) (69_mc_hmc_mecp2_by_mark_stats.tsv)
- Fac. Het (K27me3 only): 5mC rho=0.144 vs 5hmC rho=0.057; interaction p=0.493 — NOT significant (69_mc_hmc_mecp2_by_mark_stats.tsv)
- Neither: 5mC rho=-0.009 (NS) vs 5hmC rho=0.082; interaction p=1.77e-06 (significant) (69_mc_hmc_mecp2_by_mark_stats.tsv)

**What this shows:** In active and bivalent chromatin MeCP2 tracks 5mC positively but 5hmC weakly/negatively, and interaction tests confirm MeCP2 responds differently to the two modifications — consistent with canonical 5mC-reading. The exception is facultative heterochromatin (K27me3 only), where the interaction is non-significant: in this Polycomb context MeCP2 does not cleanly discriminate 5mC from 5hmC, pointing to a non-methylation (K119ub/Polycomb) recruitment mode, consistent with Section 67.

**Figures:**
- 69a_mc_hmc_mecp2_k__ac_only/ — Active context (5mC red + 5hmC blue)
- 69b_mc_hmc_mecp2_k__me____k__ac/ — Bivalent context
- 69c_mc_hmc_mecp2_k__me__only/ — Facultative het context
- 69d_mc_hmc_mecp2_neither/ — Unmarked context
