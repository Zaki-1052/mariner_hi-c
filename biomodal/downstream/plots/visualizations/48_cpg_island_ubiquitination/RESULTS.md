## Section 48: CpG Island Ubiquitination & De-Novo Methylation
**Key numbers:**
- CpG island universe = 8,910; mC sig = 442 (122 hyper, 320 hypo); hmC sig = 172 (26 hyper, 146 hypo) (48_cpg_island_ubiquitination_summary.tsv)
- De-novo among 122 hyper islands: ctrl mC <0.05 = 7 (5.7%), <0.10 = 12 (9.8%), <0.20 = 18 (14.8%) → 85%+ pre-existing; mean baseline mC hyper = 0.566 vs non-sig 0.127 (48_cpg_island_ubiquitination_summary.tsv)
- K119ub overlap DEPLETED at hyper islands: mut-peak 17.2% (hyper) vs 30.2% (non-sig); gained-K119ub 0.0% vs 1.9% → Fisher OR < 1 (48_cpg_island_ubiquitination_summary.tsv)
- Co-significant (both mC & hmC) = 112; dominant quadrant mC+/hmC- = 51 (48_cpg_island_ubiquitination_summary.tsv)

**What this shows:** Two refining (largely negative) results. CpG-island hypermethylation in BAP1-KO is overwhelmingly amplification of PRE-EXISTING methylation (baseline 0.566; only ~15% de novo), not de-novo promoter silencing. And contrary to a naive "BAP1 loss → K119ub gain → methylation" hypothesis at islands, K119ub is DEPLETED at hypermethylated islands (gained-K119ub overlap 0/122). This pushes the K119ub→methylation coupling to gene bodies / dynamic regions (Sections 47, 11–19) rather than islands; the dominant co-significant quadrant (mC+/hmC-) again reflects the TET-block signature.

**Figures:**
- 48a_k119ub_enrichment, 48b_k119ub_peak_count_violin, 48c_k119ub_diffbind_fold — K119ub at island DMRs
- 48d_baseline_mc_density, 48e_ctrl_vs_mut_scatter — baseline/observed methylation
- 48f_de_novo_classification, 48g_gain_magnitude — de-novo vs pre-existing
- 48h_mc_hmc_scatter, 48i_cosig_heatmap — coordinated mC + hmC
- 48j_k119ub_denovo_forest, 48k_k119ub_denovo_tile — K119ub × baseline Fisher
- 48l_chromatin_overlap_bars, 48m_chromatin_or_heatmap — chromatin-mark context
