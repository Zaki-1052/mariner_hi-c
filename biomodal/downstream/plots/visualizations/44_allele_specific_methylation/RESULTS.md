## Section 44: Allele-Specific Methylation (ASM)
**Key numbers:**
- Per-sample significant mC ASM: Control mean = 11,295.5, Mutant mean = 22,080.2 → ratio 1.95× (asm_mc_significant_per_sample.tsv)
- Unique significant mC ASM loci = 58,955; mutant-only = 28,641 vs control-only = 9,768; shared = 20,546 (asm_mc_significant_loci.tsv)
- Top locus chr17:39848255 (G→A): mean_meth_diff = 0.359, min_corr_p = 2.06e-139 (asm_mc_significant_loci.tsv)
- Gene-body DMR genes containing ≥1 ASM site = 1,408 (asm_dmr_overlap_summary.tsv)

**What this shows:** BAP1-KO roughly doubles the number of allele-imbalanced mC sites at heterozygous SNVs, and mutant-only loci far outnumber control-only — indicating BAP1 loss induces widespread NEW allelic methylation asymmetry rather than erasing existing imprints (the hypermethylation arm acting unevenly between alleles). 1,408 gene-body DMRs harbor ASM sites, linking the bulk DMR signal to allele-resolved events. Per-sample ASM testing is uncorrected, so absolute counts are exploratory.

**Figures:**
- 44a_filter_distribution_per_sample — ASM filter categories per sample (mC + hmC)
- 44b_pass_site_counts_per_sample — PASS site counts by condition
- 44c_asm_mc_chromosome_distribution — significant loci per chromosome (Fisher)
- 44d_methylation_diff_distributions — allelic methylation-difference density/boxplots
- 44e_mutant_vs_control_asm_comparison — per-sample sig counts, mutant vs control (Wilcoxon)
- 44f_pvalue_distributions — raw and corrected p-value distributions
- 44g_asm_volcano_plot — methylation difference vs significance by condition
- 44h_shared_asm_sites_venn — control vs mutant locus overlap
- 44i_asm_dmr_overlap — ASM sites within gene-body DMRs + DMR-vs-ASM scatter
- 44j_mc_vs_hmc_asm_comparison — mC vs hmC ASM counts + locus Venn
