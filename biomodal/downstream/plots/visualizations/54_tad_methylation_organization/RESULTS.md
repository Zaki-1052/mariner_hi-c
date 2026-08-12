## Section 54: Non-CG vs CG Methylation TAD Organization
**Key numbers:**
- Inter/intra-TAD variance ratio (Control): CG 5mC=1.018, CG 5hmC=2.078, CHG 5mC=0.0250, CHH 5mC=0.0248, MeCP2=0.0680 (54_tad_methylation_organization_summary.tsv)
- CG 5hmC is the most TAD-organized track (~2x CG 5mC), in both ctrl and mut (54_tad_methylation_organization_summary.tsv)
- CHG/CHH median per-TAD signal = 0 with intra-TAD CV pinned at 4 (empty-signal signature) (54_tad_methylation_organization_summary.tsv)
- MeCP2 median per-TAD signal falls ctrl 0.976 → mut 0.815; variance ratio stays ~0.068-0.070 (not TAD-organized) (54_tad_methylation_organization_summary.tsv)

**What this shows:** Non-CG methylation shows no TAD organization (variance ratio ~0.025, far below 1) because its per-TAD signal is essentially zero. The standout is CG 5hmC, roughly twice as TAD-organized as CG 5mC, reinforcing that the TET-product 5hmC carries the spatially-structured epigenetic information. MeCP2 raw signal is not TAD-organized and drops globally in mutant.

**Figures:**
- `54a_metatad_signal/` — per-TAD mean signal by context × condition
- `54b_boundary_insulation/` — boundary insulation by context
- `54c_intratad_cv/` — intra-TAD CV by context
- `54d_variance_ratio/` — inter/intra-TAD variance ratio (key panel)
- `54e_boundary_correlation/` — TAD score vs methylation insulation
- `54f_differential_overlay/` — log2FC at differential vs non-differential TADs
- `54g_boundary_type/` — insulation by TADCompare boundary type
- `54h_cross_context_fc/` — CG vs non-CG/MeCP2 per-TAD log2FC correlation
- `54j_composite/` — signal + CV + variance-ratio composite
