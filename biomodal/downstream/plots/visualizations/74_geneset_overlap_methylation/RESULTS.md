## Section 74: Gene-set overlap (Neuronal × Coordinated × MeCP2-Up) & neuronal methylation levels
**Key numbers:**
- Coordinated × MeCP2-Up: overlap 51, OR 5.16 [3.19, 8.51], adj p 2.82e-12 (74_pairwise_fisher.tsv)
- Neuronal × MeCP2-Up: OR 1.73, adj p 0.034; Neuronal × Coordinated: OR 1.05, adj p 0.141 NS (74_pairwise_fisher.tsv)
- Triple overlap = 16 genes (0.07%); universe 23,150 (74_geneset_overlap_counts.tsv)
- Neuronal genes with methylation data = 4,110; median 5mC 0.6782→0.6858 (↑), 5hmC 0.1272→0.1168 (↓) ctrl→mut (74_neuronal_methylation_levels.tsv)

**What this shows:** Across 23,150 quadrant-master genes, MeCP2-gain co-localizes far more strongly with coordinated mC↑/hmC↓ genes (OR=5.16) than with broad neuronal genes (OR=1.73), and neuronal genes are not preferentially coordinated at all (OR=1.05, NS). At the 4,110 neuronal genes with data, the canonical mC-up/hmC-down signature holds, with the hmC loss slightly exceeding the mC gain so total methylation edges down.

**Figures:**
- `74a_geneset_venn/` — 3-way Venn with Fisher ORs
- `74b_neuronal_total_methylation/` — total (mC+hmC) violins ctrl vs mut
- `74c_neuronal_5mc_levels/` — 5mC violins
- `74d_neuronal_5hmc_levels/` — 5hmC violins
- `74_composite/` — combined panel
