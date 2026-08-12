## Section 61i: Quadrant Gene Characterization
**Key numbers:**
- 72 quadrant genes: 15 neuronal (21%), 57 non-neuronal (79%) (source: 61i_quadrant_genes_full.tsv; 61i_neuronal_vs_nonneuronal_comparison.tsv)
- Chromatin states: Active_Promoter = 35, Other = 30, Bivalent_Promoter = 4, Repressed_Promoter = 3 (source: 61i_quadrant_genes_full.tsv)
- Only K119ub log2FC differs neuronal vs non-neuronal: 0.411 vs 0.772, q = 0.011 (*) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)
- MeCP2 fold not different: 0.862 vs 0.684, q = 0.497 (ns) (source: 61i_neuronal_vs_nonneuronal_comparison.tsv)

**What this shows:** Within the 72-gene MeCP2-Up + K119ub-Up quadrant, neuronal and non-neuronal members are statistically indistinguishable across MeCP2, 5mC, 5hmC, K27me3, K27ac, and ATAC; the only separator is modestly higher K119ub gain in non-neuronal genes. The set is dominated by Active_Promoter and "Other"/Unmarked states rather than canonical Polycomb domains.

**Figures:**
- 61i_chromatin_state_bar — chromatin-state bar of 72 genes
- 61i_neuronal_vs_nonneuronal — 4-panel boxplots by gene class
- 61i_pergene_strip — per-gene 6-mark strip
- 61i_composite — combined panel
