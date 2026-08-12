## Section 66: Methylation by CALDER2 Subcompartment

**Key numbers:**
- A.1 (strong active): 64.3% sig (6,932/10,783), 5,716 hyper vs 1,216 hypo — most affected (66_subcompartment_dmr_summary.tsv)
- B.2 (constitutive het): 21.8% sig (910/4,172), 249 hyper vs 661 hypo — least affected, hypo-dominated (66_subcompartment_dmr_summary.tsv)
- K27ac-only (Active) genes: 77.9% sig, 4,505 hyper / 325 hypo; K27me3-only: 52.8% sig, 92 hyper / 1,581 hypo (66_histone_mark_dmr_summary.tsv)
- All ctrl-vs-mut Wilcoxon tests significant (A.1 5mC & 5hmC p~0; B.2 5mC p=2.14e-62) (66_subcompartment_wilcoxon.tsv)

**What this shows:** BAP1-induced hypermethylation is compartment-specific: active chromatin (A.1/A.2, K27ac-marked) shows the most change and is hypermethylation-dominated, while constitutive heterochromatin (B.2, K27me3-marked) is largely stable and what change occurs is hypomethylation. A purely peak-based H3K27me3/H3K27ac classification cleanly mirrors the Hi-C subcompartment pattern, validating that the effect is driven by underlying active/repressive state. Uses adult/late (250402) CALDER2 labels only.

**Figures:**
- 66a_subcompartment_dmr_fraction/ — % sig DMRs by subcompartment (Chi-squared)
- 66b_direction_by_subcompartment/ — hyper/hypo split per subcompartment
- 66c_methylation_levels_by_subcompartment/ — 5mC & 5hmC violins ctrl vs mut
- 66d_histone_mark_dmr_fraction/ — same analysis by K27me3/K27ac mark
