## Section 22: Demethylation Efficiency Ratio (TET Conversion)
**Key numbers:**
- WT median ratio = 0.1275, KO median = 0.1193, median delta_ratio = -0.0087 (demethylation_ratio_all_genes.tsv)
- 71.36% of 20,915 genes have a decreased ratio in BAP1-KO (impaired TET) (demethylation_ratio_all_genes.tsv)
- Strongest per-state decrease at Active_Promoter median delta = -0.0281 (***); Repressed_Promoter is the only positive state, +0.0051 (demethylation_ratio_by_chromatin_state.tsv)
- Per-sample medians: ctrl-M 0.1292, ctrl-F 0.1278, mut-M 0.1202, mut-F 0.1168 (demethylation_ratio_per_sample.tsv)

**What this shows:** Per gene, the TET conversion-efficiency ratio 5hmC/(5mC+5hmC) drops in BAP1-KO across most of the genome (71% of genes), with both group-level and per-replicate medians falling and the largest decreases at active chromatin. The shift is partial (0.128 -> 0.118, not to zero), consistent with an indirect rather than ablative TET defect.

**Figures:**
- 22a_wt_ratio_density — WT ratio density distribution.
- 22b_wt_vs_ko_ratio — paired control vs mutant violin+box.
- 22c_delta_ratio_histogram — delta_ratio distribution by sign.
- 22d_delta_ratio_vs_diff — delta_ratio vs mc_diff and hmc_diff scatter.
- 22e_delta_ratio_by_dmr_status — delta_ratio across 4 DMR-status groups.
- 22f_delta_ratio_by_chromatin_state — delta_ratio across 7 chromatin states.
- 22g_top30_decreased_ratio — 30 most-decreased genes (lollipop).
- 22h_per_sample_ratio — per-replicate ratio distributions.
