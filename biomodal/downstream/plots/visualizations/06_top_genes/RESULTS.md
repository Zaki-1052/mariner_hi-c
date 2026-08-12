## Section 06: Top Differentially Methylated Genes

**Key numbers:**
- Top-20 5mC DMRs: 19 hypermethylated / 1 hypomethylated (source: `top_mc_dmrs.tsv`)
- Top-20 5hmC DMRs: 16 hypomethylated / 4 hypermethylated (source: `top_hmc_dmrs.tsv`)
- Total significant: 5mC = 10,775; 5hmC = 11,484; overlap (both) = 8,371 (source: run-5 BEDs / `coordinated_changes.tsv`)
- 5mC-only = 2,404; 5hmC-only = 3,113 (source: derived from run-5 significant sets)

**What this shows:** The most significant 5mC DMRs are uniformly hypermethylated and the most significant 5hmC DMRs predominantly hypomethylated — each ranked list independently recapitulates the global asymmetry. Large overlap (8,371 shared genes) confirms the two marks change together at the same loci.

**Figures:**
- `06a_top_dmrs/` — side-by-side top-20 5mC and top-20 5hmC bar charts.
- `06b_venn_overlap/` — Venn of significant 5mC vs 5hmC gene sets.
- `06_top_genes/` — composite bar charts + Venn.

> Note: `FIGURES.md` Venn ("8836 | 9930 | Both 6722") is OLD-run; current overlap = 8,371.
