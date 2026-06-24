## Section 10: 7-Category Chromatin State Classification of mC DMRs
**Key numbers:**
- 10,775 significant mC DMRs classified (7,513 hyper + 3,262 hypo) (source: chromatin_state_summary.tsv)
- Active_Promoter = 4,906 DMRs (45.5% of all) and is 93.0% hypermethylated (4,562/4,906) (source: chromatin_state_summary.tsv)
- Repressed_Promoter = 1,718 DMRs, 94.4% hypomethylated (1,621/1,718) (source: chromatin_state_summary.tsv)
- Unmarked = 3,952 (36.7%); Polycomb = 15 (14 hypomethylated) (source: chromatin_state_summary.tsv)

**What this shows:** Significant CpG DMRs were assigned to a 7-category chromatin-state system from histone-mark peak overlaps plus TSS distance. Hypermethylation concentrates at active promoters and enhancers; hypomethylation concentrates at Polycomb/Repressed promoters — a mirror-image partition by direction.

**Figures:**
- `10a_chromatin_state_distribution` — bar (faceted by direction) + overall pie
- `10b_chromatin_by_methylation_direction` — stacked + grouped state comparison
- `10c_chip_mark_overlap_heatmap` — 6-mark overlap % by direction
- `10d_coordinated_genes_chromatin` — chromatin state of coordinated mC↑/hmC↓ genes
- `10e_top_genes_chromatin_annotation` — top-20 coordinated genes annotated
- `10f_chromatin_stacked_presentation` — presentation stacked bar
