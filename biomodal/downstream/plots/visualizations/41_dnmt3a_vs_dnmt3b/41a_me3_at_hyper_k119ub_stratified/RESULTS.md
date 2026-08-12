## Section 41: DNMT3A vs DNMT3B Mechanistic Discrimination
**Key numbers:**
- 7,513 hypermethylated genes classified. Pathway attribution: Unknown 5,498 (73.2%); DNMT3A via H2AK119ub 1,478 (19.7%); Convergent K119ub+me2 268 (3.6%); DNMT3A via H3K36me2 152 (2.0%); DNMT3B via H3K36me3 84 (1.1%); Convergent K119ub+me3 33 (0.4%) (dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- Collapsed: DNMT3A axis (K119ub+me2) = 1,630 (21.7%) vs DNMT3B axis (me3) = 84 (1.1%) → ~19:1 (dnmt3a_vs_dnmt3b_pathway_attribution.tsv)
- 6-mark logistic ORs for hypermethylation: H2AK119ub OR=10.29 (p=1.04e-10); H3K27me3 OR=2.11 (p=8.5e-05); H3K27ac OR=0.44 (p=2.3e-04); H3K36me2 OR=0.62 (p=0.022); ATAC OR=0.50 (ns); H3K36me3 OR=1.00 (p=0.996, ns) (extended_logistic_regression_6mark.tsv)
- Among hyper genes: 1,779 K119ub-gained, 6,589 coordinated mC↑/hmC↓ (dnmt3a_vs_dnmt3b_pathway_attribution.tsv)

**What this shows:** The DNMT3B/SETD2/H3K36me3 mechanism is excluded — H3K36me3 attributable genes are the smallest category (84, 1.1%) and H3K36me3 has zero multivariate predictive value (OR=1.00). Where recruitment evidence exists it is overwhelmingly DNMT3A via H2AK119ub (the dominant positive predictor, OR=10.29), but 73% of hyper genes carry no recruiting-mark gain, consistent with most hypermethylation arising from an upstream TET-demethylation block (6,589 coordinated mC↑/hmC↓ genes) rather than active de novo DNMT recruitment.

**Figures:**
- 41a_me3_at_hyper_k119ub_stratified — me3 fold at hyper genes by K119ub-gain status
- 41b_me2_me3_independent_pathway — me2/me3 at hyper genes lacking K119ub gain
- 41c_logistic_regression_forest — 6-mark logistic OR forest plot (AUCs console-only)
- 41d_pathway_hypermethylation_rate — % hyper by pathway group (chi-square)
- 41e_decision_matrix_heatmap — me3 status × K119ub status O/E decision matrix
- 41f_k119ub_vs_me3_scatter — H2AK119ub vs H3K36me3 fold by mC status
- 41g_pathway_attribution_pie — headline 6-category DNMT pathway pie
- 41h_top50_summary_heatmap — top-50 hyper genes × multi-mark z-scored profile
