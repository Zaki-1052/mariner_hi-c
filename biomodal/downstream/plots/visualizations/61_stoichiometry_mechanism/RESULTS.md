## Section 61: Stoichiometry Test & Mechanism Differentiation
**Key numbers:**
- Global Delta5mC~Delta5hmC slope = -0.959 [-0.978, -0.940], R2 = 0.322, differs from -1 = TRUE (source: 61_stoichiometry_slopes.tsv)
- Total methylation NOT conserved: MeCP2-Up mean Delta = +0.032 (q = 3.2e-11), Coordinated = +0.007 (q = 2.1e-178) (source: 61_stoichiometry_summary.tsv)
- delta_ratio chromatin-model R2: ATAC = 0.095 -> +K27me3 = 0.189 -> +K119ub = 0.229 (n = 3,434) (source: 61_delta_ratio_chromatin_models.tsv)
- Mediation K27me3->ATAC->delta_ratio: indirect = -0.00279 [-0.00337, -0.00226], 18.0% mediated (source: 61_mediation_results.tsv)
- BAP1-KO vs TET-KO delta_ratio: all-genes rho = 0.220 (p = 1.2e-226, n = 20,851) (source: 61_tetko_comparison.tsv)

**What this shows:** The near -1 slope superficially resembles stoichiometric 5hmC->5mC conversion, but total methylation rises at the relevant loci, favoring TET inhibition + de novo DNMT3A over pure dehydroxymethylase. K27me3-driven chromatin compaction (ATAC loss) mediates only ~18% of the TET-efficiency drop, and BAP1-KO partially correlates with TET-KO (shared pathway, not full phenocopy). Defines the 1,149-gene neuronal set used downstream (NARROW/DMR-enriched; bias flagged in section 78).

**Figures:**
- 61a_total_methylation_conservation, 61b_stoichiometry_scatter — conservation + slope test
- 61c_chromatin_model_r2, 61d_mediation_panels — chromatin models + mediation
- 61e_neuronal_venn, 61f_bap1_vs_tetko_scatter, 61g_composite
