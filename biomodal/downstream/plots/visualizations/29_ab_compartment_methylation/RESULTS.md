## Section 29: A/B Compartment Methylation Mapping
**Key numbers:**
- mC hyper → A compartment enriched: OR = 13.642 (95% CI 11.92-15.67), p = 0; 48.76% of A genes hyper vs 6.52% of B (compartment_fisher_tests.tsv)
- hmC loss → A enriched: OR = 9.774 (8.83-10.83), p = 0; mC hypo → B: OR = 1.726, p = 8.33e-32 (compartment_fisher_tests.tsv)
- A = 12,351 genes (75.7% hyper among mC-sig, median mc_diff +0.0107), B = 3,743 genes (21.2%, median -0.0042) (compartment_methylation_summary.tsv)
- B→A shifted bins DEPLETED for mC hyper / enriched for hypo: OR = 0.0625, p = 1.37e-66 (n=427: 17 hyper vs 239 hypo) (compartment_fisher_tests.tsv, compartment_methylation_summary.tsv)

**What this shows:** BAP1-KO methylation changes match the Lopez-Moyado TET-KO signature — hypermethylation and hmC loss both concentrate in the A (euchromatin) compartment, consistent with DNMT3A redistribution into open chromatin when TET is impaired. The static A/B asymmetry supports convergent mechanisms with TET-KO. The dynamic compartment-shift result is opposite in sign: bins actively switching B→A in the mutant lose heterochromatic character and gain hypomethylation, not de novo methylation — so the methylation targets are stable-A euchromatin, not newly opened regions. 16,094 genes were assigned to compartments and matched to DMRs.

**Figures:**
- 29a_mc_violin_by_compartment — mC difference, A vs B
- 29b_hmc_violin_by_compartment — hmC difference, A vs B
- 29c_mc_violin_by_shift — mC difference by compartment shift
- 29d_hmc_violin_by_shift — hmC difference by shift
- 29e_dmr_direction_stacked_bar — mC DMR direction proportions by compartment and shift
- 29f_pc1_vs_mc_scatter — control PC1 vs mC difference with loess + Spearman
- 29g_composite_compartment_summary — composite (29a+29b+29e)
