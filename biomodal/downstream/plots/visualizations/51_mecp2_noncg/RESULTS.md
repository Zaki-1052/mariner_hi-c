## Section 51: Non-CG Methylation at MeCP2 Binding Sites (gene-body resolution)
**Key numbers:**
- CHG at all MeCP2-bound genes: ctrl 8.08e-06 vs mut 1.10e-05, Wilcoxon p = 0.4857, n = 1,211 genes (mecp2_noncg_summary.tsv)
- CHG at CG-Low MeCP2 genes: ctrl 9.12e-06 vs mut 1.33e-05, p = 0.6612 (mecp2_noncg_summary.tsv)
- CHG at CG-High MeCP2 genes: ctrl 7.04e-06 vs mut 8.69e-06, p = 0.8845 (mecp2_noncg_summary.tsv)
- All effects ~1e-6 to ~4e-6 fractional methylation (detection floor); no test significant (mecp2_noncg_summary.tsv)

**What this shows:** At gene-body resolution there is no detectable non-CG (CHG/CHH) methylation at MeCP2-bound genes and no mutant-vs-control difference, in either CG-low or CG-high strata or genome-wide. Absolute levels sit at the assay noise floor, consistent with global <1% non-CG methylation in evoC. This is a clean negative result for a gene-body-scale non-CG phenotype.

**Figures:**
- `51c_noncg_at_mecp2_by_cg_stratum/` — CHG/CHH at MeCP2 genes, CG-Low vs CG-High × condition
- `51d_mecp2_bound_vs_background_noncg/` — non-CG, MeCP2-bound vs non-MeCP2 genes
- `51a_*`, `51b_*`, `51c_cg_stratified_mecp2_fold/`, `51d_persample_*`, `51e_*` — orphan panels from an earlier script version (not produced by the current code; do not cite)
