## Section 61j: Peak-Level MeCP2-Up GO Enrichment
**Key numbers:**
- MeCP2-Up peaks = 7,686; MeCP2-Down = 1,200; MeCP2-Up + K119ub-Up = 6,849 (source: 62_mecp2_peak_chromatin_signal.tsv)
- MeCP2-Up peaks -> 2,107 genes; Up+K119ub -> 1,897 genes (source: 61j_mecp2up_peak_genes.tsv; 61j_mecp2up_k119up_peak_genes.tsv)
- MeCP2-Up GO BP sig terms: 129 (custom bg) / 818 (genome-wide bg); top = synapse assembly q = 5.3e-22 (source: 61j_mecp2up_custom_go_bp.tsv; 61j_mecp2up_genomewide_go_bp.tsv)
- MeCP2-Down genome-wide: 517 sig terms, top = synapse structure/activity q = 9.7e-11 (source: 61j_mecp2down_genomewide_go_bp.tsv)

**What this shows:** At peak resolution the neuronal/synaptic enrichment is robust and threshold-independent — both MeCP2-Up and MeCP2-Down peaks map to synaptic genes, indicating MeCP2 redistribution concentrates on the neuronal regulatory landscape. Adding the K119ub-Up filter barely changes the gene set (1,897 vs 2,107), consistent with most MeCP2-Up peaks already coinciding with K119ub gain.

**Figures:**
- 61j_mecp2up_custom_dotplot / 61j_mecp2up_genomewide_dotplot — MeCP2-Up GO dotplots
- 61j_up_k119up_custom_dotplot / 61j_up_k119up_genomewide_dotplot — Up+K119ub GO dotplots
- 61j_mecp2down_genomewide_dotplot — MeCP2-Down GO dotplot
- 61j_composite — combined ORA bar composite (from section_61jk_composites.R)
