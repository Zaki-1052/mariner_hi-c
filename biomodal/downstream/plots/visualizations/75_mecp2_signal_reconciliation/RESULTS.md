## Section 75: MeCP2 signal-direction reconciliation (UP peaks vs gene-body drop)
**Key numbers:**
- 7,686 MeCP2 UP peaks concentrate on only 2,052 genes; 518 genes DOWN-only; 18,625 with no sig peaks (75_mecp2_peak_gene_summary.tsv, 75_peak_annotation_distribution.tsv)
- UP peaks are 51.7% distal-intergenic, 41.0% intron, 2.2% promoter; DOWN peaks 61.8% intron, 8.0% promoter (75_peak_annotation_distribution.tsv)
- K119ub GSEA = 115 sig GO BP terms (3 neuronal); MeCP2 GSEA = 1 sig term (1 neuronal, "synapse assembly") (75_go_term_comparison.tsv)

**What this shows:** The apparent paradox (DiffBind shows ~7,686 UP peaks yet gene-body MeCP2 signal falls genome-wide) is resolved by redistribution: UP peaks pile onto ~2,000 genes (≈10% of MeCP2-bound genes) and are mostly distal-intergenic, while DOWN peaks are genic — MeCP2 vacates gene bodies and concentrates at distal sites. K119ub's functional footprint is broad (115 programs) whereas MeCP2's is narrowly synaptic (a single significant term).

**Figures:**
- `75a_peak_distribution_per_gene/` — genes by MeCP2 peak status
- `75b_genebody_signal_by_peak_status/` — UP-peak genes carry positive gene-body fold
- `75c_gsea_term_comparison/` — K119ub (115) vs MeCP2 (1) sig terms
- `75d_peak_annotation_distribution/` — UP vs DOWN annotation split
- `75_composite/` — combined panel
