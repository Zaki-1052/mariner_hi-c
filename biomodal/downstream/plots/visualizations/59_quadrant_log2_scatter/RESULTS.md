## Section 59: Log2 Quadrant Scatters (K119ub x MeCP2 x 5mC x 5hmC)
**Key numbers:**
- Master gene-level table = 23,150 genes; 21,604 K119ub-quantifiable (source: 59_quadrant_master.tsv)
- MeCP2-Up = 79, MeCP2-Down = 34 genes (FDR < 0.05; source: 59_quadrant_master.tsv)
- Master columns integrate K119ub gb_log2fc, MeCP2 fold/FDR, 5mC/5hmC diff, K27ac/K27me3/K119ub DiffBind folds, chromatin state

**What this shows:** Builds the master gene-body integration table joining H2AK119ub BigWig signal, concentration-weighted MeCP2 DiffBind fold, CG 5mC/5hmC, and histone-mark folds, then renders pairwise quadrant scatters with Spearman correlations to encode the BAP1-KO cascade (K119ub up -> 5mC up / 5hmC down / K27ac down / K27me3 up). This table feeds sections 60, 61, 61h, and 61i.

**Figures:**
- 59a_k119ub_vs_mecp2 / 59a2_k119ub_vs_mecp2_peaklevel — K119ub vs MeCP2 (gene-body / peak-level)
- 59b_mecp2_vs_k27ac, 59c_mecp2_vs_k27me3, 59bc_chromatin_composite — euchromatin vs heterochromatin
- 59d_k119ub_vs_5mc, 59e_k119ub_vs_5hmc — K119ub vs methylation change
- 59f_mecp2_vs_5mc, 59g_mecp2_vs_5hmc, 59h_composite — MeCP2 vs methylation + full composite
