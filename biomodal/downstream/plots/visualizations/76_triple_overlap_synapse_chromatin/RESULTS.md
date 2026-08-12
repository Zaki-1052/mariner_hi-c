## Section 76: Triple-overlap & synapse/axon-specific chromatin remodeling
**Key numbers:**
- Synapse/axon vs broader-neuronal K27me3: Δ median −0.0444, p 2.95e-3, adj p 6.65e-3 (significant); ATAC p 0.654 NS, K27ac p 0.757 NS (76_synapse_vs_neuronal_stats.tsv)
- Broader-neuronal vs non-neuronal K27me3 is NS (p 0.289) — K27me3-loss specificity sits in the synapse/axon subset (76_synapse_vs_neuronal_stats.tsv)
- K119ub D10 enrichment: Neuronal OR 2.34 (677/4,030, p 1.34e-54), Synapse/axon OR 1.68 (306/2,069, p 2.12e-13), Triple 0/16 OR 0 NS (76_top_decile_fisher.tsv)
- 16 triple-overlap genes incl. Snca, Hcn1, Lgi1, Astn2, Cntn6 (76_triple_overlap_genes.tsv)

**What this shows:** Synapse/axon genes are chromatin-special specifically for Polycomb de-repression — they lose more K27me3 than broader neuronal genes but show no extra accessibility or enhancer gain. The 16 triple-overlap genes are too few for robust statistics and are not at the K119ub extreme (0/16 in D10); instead the broad neuronal and synapse sets are disproportionately K119ub-high, marking them as natural BAP1 substrates.

**Figures:**
- `76a_neuronal_chromatin_with_pvalues/` — neuronal vs non-neuronal violins with p-values
- `76b_four_group_chromatin/` — triple / neuronal-only / coordinated-only / rest
- `76c_synapse_vs_neuronal_chromatin/` — synapse vs broader-neuronal vs non-neuronal
- `76d_k119ub_decile_enrichment/` — D10 odds-ratio forest
- `76_composite/` — combined panel
