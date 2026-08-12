## Section 77: MeCP2 developmental aging trajectory (young vs adult)
**Key numbers:**
- Aging-UP peaks: Control 10,930 vs Mutant 23,117 (2.1×); aging-DOWN 2,822 vs 10,646 (3.8×) (77_aging_peak_summary.tsv)
- 7,305 ctrl aging-UP peaks overlap mut (66.8% of ctrl); mut-unique = 15,812 peaks / 1,654 genes (4,274−2,620) (77_aging_overlap.tsv)
- Mut-unique aging genes: 435 sig GO BP terms, 49 neuronal (77_mut_specific_aging_go.tsv)
- Shared loci fold: median ctrl 1.829 vs mut 2.241 (mut +0.41 log2) (77_shared_peak_fold_comparison.tsv)

**What this shows:** MeCP2 normally accumulates with neuronal maturation; the mutant amplifies that program rather than rewiring it. Two-thirds of the control aging-UP peaks are preserved in the mutant, but the mutant adds 1,654 unique aging genes (GO-enriched for neuronal/synaptic terms) and gains more MeCP2 fold even at shared loci. Inputs from collaborator Jai are present, so this section ran fully (it is no longer blocked).

**Figures:**
- `77a_aging_overview/` — UP/DOWN/NS aging peaks by genotype
- `77b_aging_overlap_venn/` — gene overlap (2,620 shared / 288 ctrl / 1,654 mut)
- `77c_mut_specific_go_enrichment/` — GO dotplot of 1,654 mut-unique genes
- `77d_shared_peak_fold_comparison/` — ctrl vs mut fold at 7,305 shared peaks
- `77_composite/` — combined panel
