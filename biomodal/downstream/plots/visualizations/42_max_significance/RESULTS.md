## Section 42: Maximum-Significance Gene List
**Key numbers:**
- 5mC genes at q-floor (q < 1e-300) = 100; 5hmC = 67; merged unique = 124 rows / 123 names (max_significance_genes_merged.tsv, max_significance_gene_names.txt)
- Coordinated mC-up/hmC-down = 42 of 124 loci; mC-down/hmC-up = 1; single-modality = 81 (max_significance_genes_merged.tsv)
- Top mC: Syt1 +17.31% (ctrl 0.588 → mut 0.762); top hmC: Syt1 -14.97% (ctrl 0.266 → mut 0.117) (max_significance_genes_mc.tsv, max_significance_genes_hmc.tsv)

**What this shows:** The genes with the most extreme statistical significance are dominated by the coordinated mC-up/hmC-down population (42/124), with synaptic/adhesion neuronal genes (Syt1, Cntnap2, Cadm2, Cdh8) as the canonical floor set. Syt1 is the archetype, gaining ~17% 5mC while losing ~15% 5hmC. PCA panels confirm these genes separate mutant from control per sample.

**Figures:**
- 42_pca_max_significance_genes — PCA (5mC | 5hmC) of the q-floor gene set, per-sample, colored by condition/batch
- 42_pca_all_genes — genome-wide PCA of all gene bodies for comparison
