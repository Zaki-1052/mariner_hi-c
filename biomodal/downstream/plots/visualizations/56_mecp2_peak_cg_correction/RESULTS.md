## Section 56: Peak-Level MeCP2 vs CG Methylation (Non-CG Binding Test)
**Key numbers:**
- Peak regression MeCP2 Fold ~ CG mC log2FC: R²=0.0128, slope=0.1909 (p≈0), n=202,574 (56_mecp2_peak_regression.tsv)
- MeCP2-Up peaks (FDR<0.05, Fold>0): 7,686 — CG-Concordant 4,960 (64.5%), Non-CG Candidate 2,726 (35.5%) (56_mecp2_peak_classification_summary.tsv)
- Non-CG candidates 77.2% Unmarked chromatin vs 62.0% for CG-concordant; concordant peaks far more enhancer-marked (35.7% Active+Poised vs 19.2%) (56_mecp2_peak_chromatin_summary.tsv)
- MeCP2-Up DMR overlap: 2,687 Hyper-DMR vs 4,999 No-CG-DMR (56_mecp2_peak_dmr_crosstab.tsv)

**What this shows:** Peak-level CG methylation explains only ~1.3% of MeCP2 fold-change variance, and 35.5% of MeCP2-Up peaks gain MeCP2 with no CG increase — the "non-CG candidate" set. These candidates are depleted of active/poised enhancer marks and enriched in Unmarked, low-CpG chromatin (77%), matching the mCA-reader model of MeCP2 spreading into gene-poor space (top recurrent locus: Etv1). This 2,726-peak set feeds sections 57-58.

**Figures:**
- `56a_mecp2_vs_cg_scatter/` — MeCP2 Fold vs CG log2FC (all peaks)
- `56b_peak_classification/` — concordant vs non-CG candidate (key panel)
- `56c_cg_levels_by_mecp2/` — CG 5mC at MeCP2-Up vs Down, ctrl/mut
- `56d_residual_distributions/` — residual densities Up vs Down
- `56e_dmr_crosstab/` — MeCP2 direction × CG DMR overlap
- `56f_chromatin_context/` — 7-state composition of the two peak classes
- `56g_composite/` — scatter + classification + residual composite
