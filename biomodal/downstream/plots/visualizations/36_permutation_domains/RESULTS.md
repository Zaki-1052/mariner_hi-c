## Section 36: Domain-Level Interval Permutation — Compartments + Polycomb (Confirmatory)
**Key numbers:**
- 24 pairwise regioneReloaded tests at ntimes=100,000 (36A: 4 DMR sets x 4 compartment sets=16; 36B: 4 DMR sets x 2 Polycomb sets=8) (source: logs/5/biomodal_perm_47851859.out, lines 97-100)
- DMR direction sets: mC Hyper 7,522; mC Hypo 3,267; hmC Hypo 9,532; hmC Hyper 1,967 regions (source: logs/5/biomodal_perm_47851859.out, lines 75-78)
- Compartment bins: A 47,241; B 54,480; B→A shift 5,430; A→B shift 2,641 (of 101,721 filtered) (source: logs/5/biomodal_perm_47851859.out, lines 85-89)
- Polycomb: H3K27me3 15,809 peaks; Bivalent 318 peaks (source: logs/5/biomodal_perm_47851859.out, lines 92-93)

**What this shows:** Megabase-scale permutation validating the compartment (section 29) and Polycomb (section 30) Fisher tests. This section ran furthest of the three interval-permutation sections — it produced the compartment heatmap (36a), Polycomb heatmap (36b), and the local z-score curve for mC-hyper DMRs at the A compartment (36c), failing only at the final comparison dot-plot (36d), so no `permutation_36_comparison.tsv` was written. The gene-level companions in section 37 confirm the coupling: B→A shifts anti-associate with hypermethylation (37-10a z=-25.8) and A→B shifts associate with it (37-10b z=+5.1).

**Figures:**
- 36a_crosswise_dmr_x_compartment — 4x4 heatmap: DMR direction x A/B compartment + shifts
- 36b_crosswise_dmr_x_polycomb — 4x2 heatmap: DMR direction x H3K27me3 / Bivalent
- 36c_local_zscore_compartment — local z-score curve, mC Hyper DMRs at A compartment (±50kb). Panel 36d not generated (run crashed).
