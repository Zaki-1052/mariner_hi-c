## Section 78: Stoichiometry & mechanism with the unbiased (broad) neuronal set
**Key numbers:**
- Direction self-correction: narrow neuronal (1,149) mean Δtotal +0.01200 (UP, biased) vs broad neuronal (4,118) −0.00222 (DOWN); broad-only (2,969) −0.00772 (recomputed from diffbind_gene_level_all_marks.tsv + 61/72 gene lists; also 78_total_methylation_summary.tsv for broad)
- Stoichiometry slopes: Neuronal (broad) −0.995 [−1.039, −0.952] and Synapse/axon −1.020 both consistent with −1; All genes −0.959 and Coordinated −1.288 differ (78_stoichiometry_slopes.tsv)
- BAP1-KO vs TET-KO Spearman rho: All 0.220, Non-neuronal 0.246, Neuronal 0.137, Synapse 0.121 — weakest at neuronal (78_tetko_comparison.tsv)
- Quadrant (72 MeCP2-Up + K119ub-Up): broad neuronal 25/72 (35%), synapse/axon 17/72 (24%) (78_quadrant_neuronal_comparison.tsv)

**What this shows:** Section 61's "total methylation increases at neuronal genes" was a DMR-selection artifact: with the unbiased GO set the sign flips to a decrease (−0.0022), and total methylation rises only at coordinated (+0.007) and MeCP2-Up (+0.031) genes. Neuronal and synapse/axon genes change stoichiometrically (slope −1.0, dehydroxymethylase-like) and resemble TET-KO least, distinguishing them from the TET-inhibition pattern at non-neuronal loci.

**Figures:**
- `78a_total_methylation_forest/` — mean Δtotal by gene group
- `78b_stoichiometry_scatter/` — δ-5mC vs δ-5hmC by neuronal class
- `78c_stoichiometry_slope_forest/` — slopes vs −1 reference
- `78d_bap1_vs_tetko/` — BAP1-KO vs TET-KO δ-ratio
- `78e_narrow_vs_broad_bias/` — narrow/broad Δtotal direction flip
- `78g_quadrant_characterization/` — 72 quadrant genes by class
- `78_composite/` — combined panel
