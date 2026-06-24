## Section 65: MeCP2 vs Methylation Scale (Context-Dependent Reader)

**Key numbers:**
- 20,969 gene bodies tested; 10,775 (51.4%) sig 5mC change; 7,513 mC-hyper (65_methylation_mecp2_cascade.tsv)
- Sig MeCP2 up = 2,052 genes (9.8%) — ~5x fewer than sig-5mC genes (65_methylation_mecp2_cascade.tsv)
- mC-hyper + sig MeCP2 up = only 947 (4.5% of tested; 12.6% of mC-hyper) (65_methylation_mecp2_cascade.tsv)
- Fisher sig-5mC x sig-MeCP2-up: OR=2.63, p=2.23e-61; mecp2_only cell = 359 genes (65_fisher_context_dependence.tsv)

**What this shows:** Methylation change is roughly five times more pervasive than MeCP2 binding change at the gene level, and only ~13% of hypermethylated genes recruit more MeCP2. The two events are coupled above chance (Fisher OR 2.63) but the coupling is weak and context-gated — methylation is necessary but not sufficient for a MeCP2 response. The 359 genes that gain MeCP2 WITHOUT a methylation change (the Fisher mecp2_only cell) are carried into Section 67.

**Figures:**
- 65a_cascade_hierarchy/ — methylation vs MeCP2 funnel bars
- 65b_gene_level_venn/ — mC-hyper / MeCP2-up / coordinated Venn
- 65c_allgene_quadrant_scatter/ — 5mC change vs MeCP2 fold, all genes
- 65d_proportional_funnel/ — tested -> sig -> hyper -> +MeCP2 funnel
- 65e_mecp2_signal_by_meth_status/ — MeCP2 signal violins by methylation status
