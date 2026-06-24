## Section 68: modC vs MeCP2 by Histone-Mark Context

**Key numbers:**
- K27ac only (Active): n=5,929, Spearman rho=0.183, p=1.02e-45 (68_modc_mecp2_by_mark_stats.tsv)
- K27me3+K27ac (Bivalent): n=477, rho=0.461, p~0 — strongest coupling (68_modc_mecp2_by_mark_stats.tsv)
- K27me3 only (Fac. Het): n=3,073, rho=0.189, p=4.10e-26 (68_modc_mecp2_by_mark_stats.tsv)
- Neither (unmarked): n=7,430, rho=0.070, p=1.41e-09 — weakest (68_modc_mecp2_by_mark_stats.tsv)

**What this shows:** Total modC change correlates with MeCP2 binding change only modestly in active and facultative-het contexts (rho ~0.18-0.19) and barely at all at unmarked genes (rho 0.07). The strongest coupling is at bivalent (K27me3+K27ac) promoters (rho 0.46), where the BAP1/Polycomb axis and methylation converge. The context-dependence of correlation strength reinforces that MeCP2's relationship to methylation is gated by chromatin state, not stoichiometric.

**Figures:**
- 68a_modc_mecp2_k__ac_only/ — Active context quadrant scatter
- 68b_modc_mecp2_k__me____k__ac/ — Bivalent context
- 68c_modc_mecp2_k__me__only/ — Facultative het context
- 68d_modc_mecp2_neither/ — Unmarked context
