## Section 38: H3K36me3 Gene Body Analysis (SETD2/DNMT3B axis)
**Key numbers:**
- me3 × mC direction concordance: Fisher OR = 4.30, p = 1.24e-04, n = 207 genes (h3k36me3_direction_concordance.tsv)
- me3 × hmC concordance: Fisher OR = 0.163, p = 1.60e-05, n = 184 (h3k36me3_direction_concordance.tsv)
- me3 fold vs mC diff: Spearman rho = +0.018, p = 0.067 (ns); vs hmC diff: rho = -0.034, p = 5.65e-04 (h3k36me3_direction_concordance.tsv)
- Significant H3K36me3 gene-level peaks: 333 genes (285 gained / 48 lost) of 11,952 (h3k36me3_gene_level_summary.tsv)

**What this shows:** H3K36me3 — the SETD2-deposited mark that recruits DNMT3B — is only a weak, inconsistent correlate of methylation change. The me3 × mC Fisher enrichment (OR=4.30) rests on 207 genes, but the continuous gene-level me3-vs-mC correlation is essentially zero, and me3 tracks hmC (active/hydroxymethylated state) more than mC. DNMT3B does not behave like the primary hypermethylation effector.

**Figures:**
- 38a_h3k36me3_volcano_annotation — DiffBind volcano + significant-peak annotation pie
- 38b_me3_x_mc_oe_heatmap — me3 direction × 5mC DMR direction O/E (Fisher)
- 38c_me3_x_hmc_oe_heatmap — me3 direction × 5hmC DMR direction O/E
- 38d_me3_vs_methylation_scatter — gene-level me3 fold vs 5mC and 5hmC change
- 38e_me3_fold_coordinated_violin — me3 fold at coordinated (mC↑/hmC↓) vs other genes
- 38f_dmr_me3_bed_overlap — mC DMR overlap with me3 gained/lost peak BEDs
- 38g_chromosome_distribution — chromosome distribution of sig me3 peaks (chr13 histone cluster)
- 38h_me3_fold_by_chromatin_state — me3 fold by 7-category chromatin state
