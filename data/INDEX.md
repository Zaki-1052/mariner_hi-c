# Paper Data & Figure Index

Consolidated ground-truth data and publication figures for the BAP1 Hi-C paper. All files are **copies** from their original analysis locations. This directory serves as a single reference point for reviewers and cross-codebase lookup.

**Generated:** 2026-03-21
**Source mapping:** `docs/Hi-C-paper-annotated.md` (complete provenance for each panel)

---

## Directory Structure

```
data/
├── tsvs/                              # Tabular data files (TSV, BED, TXT)
│   ├── figure_1_tads_boundaries_compartments/   (16 files)
│   ├── figure_2_loop_rewiring/                  (38 files)
│   ├── figure_3_epigenetic_integration/         (15 files)
│   ├── figure_4_abc_analysis/                   (19 files)
│   ├── figure_5_model_functional/               (18 files)
│   └── supplemental/                            (25 files)
├── plots/                             # Publication figures (SVG preferred, JPG for heavy plots)
│   ├── figure_1_tads_boundaries_compartments/   (10 files)
│   ├── figure_2_loop_rewiring/                  (64 files)
│   ├── figure_3_epigenetic_integration/         (28 files)
│   ├── figure_4_abc_analysis/                   (24 files)
│   ├── figure_5_model_functional/               (15 files)
│   └── supplemental/                            (17 files)
├── upstream/                          # Reference input data
│   ├── rna_seq/                                 (1 file)
│   ├── chip_peaks/                              (2 files)
│   └── loop_calls/                              (3 files)
└── INDEX.md                           # This file
```

**Totals:** 137 data files + 158 plot files + 6 upstream = **295 files**

## Naming Convention

Files are named `{panel}_{descriptive_name}.{ext}`:
- **Panel prefix:** `1B`, `2A`, `3C`, etc. maps to the paper figure panel
- **Timepoint prefix:** `early` = P13, `late` = adult
- **Format:** SVG for most plots; JPG for heavy volcano/scatter plots with many points; PDF when SVG/JPG unavailable

---

## Figure 1: Progressive Intensification of 3D Genome Phenotype

TADs, boundaries, and compartment changes between BAP1-KO and wildtype.

### Data (`tsvs/figure_1_tads_boundaries_compartments/`)

| File | Panel | Description | Source |
|------|-------|-------------|--------|
| `1B_tad_all_annotated.tsv` | 1B | All TADs with differential statistics | `tads/tad-pc-analysis/output/tad_analysis/` |
| `1B_tad_significant_relaxed.tsv` | 1B | Significant TADs (FDR<0.15, \|Diff\|>0.15) | same |
| `1B_tad_significant_standard.tsv` | 1B | Significant TADs (FDR<0.05, \|Diff\|>0.30) | same |
| `1B_early_tadcompare_differential.tsv` | 1B | TADCompare results, P13 timepoint | `tads/results/early/tadcompare/` |
| `1B_late_tadcompare_differential.tsv` | 1B | TADCompare results, adult timepoint | `tads/results/late/tadcompare/` |
| `1C_early_differential_boundaries.bed` | 1C | Differential boundaries (P13) | `tads/results/early/final/` |
| `1C_late_differential_boundaries.bed` | 1C | Differential boundaries (adult) | `tads/results/late/final/` |
| `1C_late_shifted_boundaries.bed` | 1C | Shifted boundaries (min 2 bins) | `tads/results/late/final/` |
| `1D_compartment_all_annotated.tsv` | 1D | All compartments with PC1 differential | `tads/tad-pc-analysis/output/compartment_analysis/` |
| `1D_compartment_significant_relaxed.tsv` | 1D | Significant compartment changes (relaxed) | same |
| `1D_compartment_significant_standard.tsv` | 1D | Significant compartment changes (standard) | same |
| `1D_compartment_genome_pct_summary.txt` | 1D | % genome with differential PC1 | same |
| `1E_syt1_nav3_statistics.tsv` | 1E | Syt1/Nav3 locus statistics | `tads/results/visualizations/late/syt1_nav3_focus/` |
| `1F_permutation_test_results.tsv` | 1F | Boundary-loop permutation test | `tads/results/late/boundary_loop_analysis/` |
| `1F_enrichment_statistics.tsv` | 1F | ChIP enrichment at boundaries | same |
| `1F_boundaries_with_chromatin_state.tsv` | 1F | Boundary chromatin state annotations | `tads/results/visualizations/chip/late/` |

### Plots (`plots/figure_1_tads_boundaries_compartments/`)

| File | Panel | Description |
|------|-------|-------------|
| `1B_tad_volcano_relaxed.pdf` | 1B | TAD volcano plot (relaxed thresholds) |
| `1B_tad_volcano_standard.pdf` | 1B | TAD volcano plot (standard thresholds) |
| `1D_compartment_volcano_relaxed.jpg` | 1D | Compartment volcano (JPG - heavy plot) |
| `1D_compartment_volcano_standard.jpg` | 1D | Compartment volcano (JPG - heavy plot) |
| `1D_compartment_genome_pct_pie.svg` | 1D | % genome pie chart by compartment shift |
| `1D_compartment_genome_pct_bar.svg` | 1D | % genome bar chart |
| `1E_syt1_nav3_regional_overview.svg` | 1E | Syt1/Nav3 regional focus |
| `1F_permutation_test.svg` | 1F | Boundary-loop permutation test |
| `1F_chip_mark_overlap_heatmap.svg` | 1F | ChIP mark overlap at boundaries |
| `1F_chromatin_state_by_differential.svg` | 1F | Chromatin state by differential status |

---

## Figure 2: Chromatin Loop Rewiring

Differential loop analysis, APA, anchor annotation.

### Data (`tsvs/figure_2_loop_rewiring/`)

| File | Panel | Description |
|------|-------|-------------|
| `2A_{early,late}_edgeR_{5,10,25}kb.tsv` | 2A | edgeR results per resolution (6 files) |
| `2B_late_distance_shift_summary.tsv` | 2B | Loop distance shift summary |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_enrichment.tsv` | 2C | APA enrichment scores (12 files) |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_stats.tsv` | 2C | APA statistical tests (12 files) |
| `2F_{early,late}_anchor_type_summary.tsv` | 2F | Anchor category counts (2 files) |
| `2G_{early,late}_loop_type_summary.tsv` | 2G | Loop type pair counts (2 files) |
| `2G_{early,late}_extended_characterized_loops.tsv` | 2G | 7-category annotated loops (2 files) |
| `2H_late_characterized_loops.tsv` | 2H | Merged characterized loops (adult) |

### Plots (`plots/figure_2_loop_rewiring/`)

| File | Panel | Description |
|------|-------|-------------|
| `2A_{early,late}_volcano_merged.jpg` | 2A | Multi-resolution merged volcano (JPG - heavy) |
| `2B_{early,late}_distance_cdf.svg` | 2B | Lost vs gained loop eCDF |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_{heatmap_type}.svg` | 2C | APA heatmaps (48 files) |
| `2E_late_cdf_k27me3_anchored.svg` | 2E | K27me3-anchored loop eCDF |
| `2E_late_cdf_H3K27ac_one_anchor.svg` | 2E | H3K27ac-anchored loop eCDF |
| `2E_late_cdf_H3K27me3_one_anchor.svg` | 2E | H3K27me3-anchored loop eCDF |
| `2F_{early,late}_anchor_type_distribution.svg` | 2F | Anchor type bar plots |
| `2G_late_loop_type_by_direction.svg` | 2G | Loop categories by gain/loss |
| `2G_late_loop_type_piechart.svg` | 2G | Loop type pie chart |
| `2H_late_looptype_distance_heatmap.svg` | 2H | Distance x loop type heatmap |
| `2H_late_deg_violin_by_loop_type.svg` | 2H | DEG expression by loop type |
| `2I_late_cdf_{H3K27ac,H3K27me3,CTCF}.svg` | 2I | Mark-specific eCDFs (3 files) |

---

## Figure 3: Integration with Epigenetic Changes

H2AK119ub, ChIP-seq at boundaries, permutation tests, timepoint comparison.

### Data (`tsvs/figure_3_epigenetic_integration/`)

| File | Panel | Description |
|------|-------|-------------|
| `3A_k119ub_anchor_signal.tsv` | 3A | K119ub signal at loop anchors |
| `3A_k119ub_enrichment_global.tsv` | 3A | Global K119ub enrichment |
| `3A_k119ub_enrichment_by_chromstate.tsv` | 3A | K119ub enrichment by chromatin state |
| `3A_k119ub_correlation_results.tsv` | 3A | K119ub-loop correlation |
| `3A_k119ub_logistic_regression.tsv` | 3A | Logistic regression for K119ub |
| `3A_loops_with_k119ub.tsv` | 3A | Loops annotated with K119ub |
| `3B_boundaries_chromatin_state.tsv` | 3B | TAD boundaries with chromatin states |
| `3C_boundary_permutation_results.tsv` | 3C | Boundary permutation test |
| `3C_enrichment_tests_all_loops.tsv` | 3C | Diff ChIP enrichment (all loops) |
| `3C_enrichment_tests_polycomb.tsv` | 3C | Diff ChIP enrichment (polycomb) |
| `3C_overlap_summary_all_loops.tsv` | 3C | ChIP overlap summary (all) |
| `3C_overlap_summary_polycomb.tsv` | 3C | ChIP overlap summary (polycomb) |
| `3C_loop_compartment_annotated.tsv` | 3C | Loops with compartment annotation |
| `3D_timepoint_comparison_stats.tsv` | 3D | Early vs late comparison statistics |
| `3F_syt1_nav3_statistics.tsv` | 3F | Syt1/Nav3 locus stats |

### Plots (`plots/figure_3_epigenetic_integration/`)

| File | Panel | Description |
|------|-------|-------------|
| `3A_cdf_k119ub_up_one_anchor.svg` | 3A | K119ub CDF at gained loop anchors |
| `3A_enrichment_dotplot_by_chromstate.svg` | 3A | K119ub enrichment by chromatin state |
| `3A_scatter_loopFC_vs_k119ub.svg` | 3A | Loop logFC vs K119ub logFC |
| `3A_boxplot_k119ub_by_direction.svg` | 3A | K119ub by loop direction |
| `3A_correlation_summary_heatmap.svg` | 3A | Correlation summary |
| `3B_01-06_*.svg` | 3B | ChIP boundary classification (6 panels) |
| `3C_boundary_permutation_test.svg` | 3C | Boundary permutation test |
| `3C_enrichment_dotplot_*.svg` | 3C | Diff ChIP enrichment dotplots (2 files) |
| `3C_loop_cmpt_q*.svg` | 3C | Loop-compartment crossref (7 files) |
| `3D_*.svg` | 3D | Timepoint comparison plots (6 files) |
| `3F_syt1_nav3_regional_overview.svg` | 3F | Syt1/Nav3 integrated locus |

---

## Figure 4: Enhancer ABC Analysis

Activity-by-Contact analysis, concordance, K119ub integration.

### Data (`tsvs/figure_4_abc_analysis/`)

| File | Panel | Description |
|------|-------|-------------|
| `4A_delta_abc_all_pairs.tsv` | 4A | All enhancer-gene DABC pairs (~180K) |
| `4A_delta_abc_with_rnaseq.tsv` | 4A | DABC merged with RNA-seq (~113K) |
| `4A_activity_contact_summary.tsv` | 4A | Activity vs contact decomposition |
| `4B_discordant_gene_characteristics.tsv` | 4B | Discordant gene features |
| `4B_concordance_by_class.tsv` | 4B | Concordance by enhancer class |
| `4B_discordant_go_bp.tsv` | 4B | GO enrichment for discordant genes |
| `4B_discordant_kegg.tsv` | 4B | KEGG enrichment for discordant genes |
| `4D_class_level_summary.tsv` | 4D | Enhancer class summary stats |
| `4D_enhancer_classes_*.tsv` | 4D | Enhancer class assignments (4 files) |
| `4E_enhancer_class_abc_metrics.tsv` | 4E | ABC metrics by class |
| `4E_enhancer_class_loop_metrics.tsv` | 4E | Loop metrics by class |
| `4E_contact_decay_by_class.tsv` | 4E | Contact decay curves |
| `4F_k119ub_tertile_assignments.tsv` | 4F | K119ub tertile classification |
| `4F_k119ub_tertile_loop_summary.tsv` | 4F | Tertile loop summary |
| `4F_k119ub_abc_correlation_summary.tsv` | 4F | K119ub-ABC correlation |
| `4F_k119ub_abc_enhancer_merged.tsv` | 4F | Merged enhancer K119ub + ABC |

### Plots (`plots/figure_4_abc_analysis/`)

| File | Panel | Description |
|------|-------|-------------|
| `4A_raw_delta_all_pairs.svg` | 4A | Raw delta activity vs contact |
| `4A_log2fc_all_pairs.jpg` | 4A | log2FC scatter (JPG - heavy) |
| `4A_raw_delta_classified.svg` | 4A | Classified by concordance category |
| `4A_raw_delta_promoter_distal.svg` | 4A | Promoter vs distal enhancers |
| `4B_concordance_pie_combined.svg` | 4B | Combined concordance pie chart |
| `4B_concordance_pie_4cat.svg` | 4B | 4-category concordance |
| `4B_discordant_composite.svg` | 4B | Discordant analysis composite |
| `4B_discordant_go_bp.svg` | 4B | GO enrichment for discordant |
| `4B_k119ub_by_concordance.svg` | 4B | K119ub by concordance status |
| `4B_k119ub_significance_rate.svg` | 4B | K119ub significance rate |
| `4C_k119ub_volcano_at_enhancers.jpg` | 4C | K119ub volcano at enhancers (JPG) |
| `4D_summary_patchwork.svg` | 4D | Enhancer model summary |
| `4D_class_composition_bar.svg` | 4D | Class composition |
| `4E_loop_logfc_violin.svg` | 4E | Loop logFC by class |
| `4E_delta_abc_boxplot.svg` | 4E | Delta ABC by class |
| `4E_gene_logfc_by_class.svg` | 4E | Gene logFC by class |
| `4E_contact_decay_{wt,delta}.svg` | 4E | Contact decay curves (2 files) |
| `4E_logFC_vs_deltaABC.svg` | 4E | Loop logFC vs DABC |
| `4F_k119ub_tertile_loop_logfc.svg` | 4F | K119ub tertile loop logFC |
| `4F_k119ub_tertile_delta_abc.svg` | 4F | K119ub tertile DABC |
| `4F_k119ub_vs_loop_logfc_scatter.svg` | 4F | K119ub vs loop logFC |
| `4F_contingency_heatmap.svg` | 4F | K119ub x ABC contingency |
| `4F_boxplot_k119ub_by_abc_category.svg` | 4F | K119ub by ABC category |

---

## Figure 5: Model and Functional Implications

GO analysis, DEG integration, network analysis, top gene heatmaps.

### Data (`tsvs/figure_5_model_functional/`)

| File | Panel | Description |
|------|-------|-------------|
| `5A_boundary_genes.tsv` | 5A | Genes at differential boundaries |
| `5A_go_long_lost_loops.tsv` | 5A | GO for long-range lost loops |
| `5A_go_short_gained_loops.tsv` | 5A | GO for short-range gained loops |
| `5A_abc_go_enrichment.tsv` | 5A | GO at paired ABC anchors |
| `5A_abc_kegg_enrichment.tsv` | 5A | KEGG at paired ABC anchors |
| `5B_deg_boundary_genes.tsv` | 5B | DEGs at TAD boundaries |
| `5B_deg_anchor_genes.tsv` | 5B | DEGs at loop anchors |
| `5B_deg_longrange_vs_shortrange_genes.tsv` | 5B | Long-range lost vs short-range gained genes |
| `5B_gene_level_summary.tsv` | 5B | Gene-level ABC summary (13.6K genes) |
| `5C_gene_structural_profile_filtered.tsv` | 5C | Network node data (filtered) |
| `5C_gene_structural_profile_all.tsv` | 5C | Network node data (all) |
| `5C_edge_list.tsv` | 5C | Network edge list |
| `5C_network_centrality_metrics.tsv` | 5C | Network centrality scores |
| `5C_go_enrichment_results.tsv` | 5C | GO enrichment for network genes |
| `5C_loops_with_gene_assignments.tsv` | 5C | Loop-gene mapping |
| `5C_paired_anchor_matches.tsv` | 5C | Paired anchor matches |
| `5D_top50_combined_score.tsv` | 5D | Top 50 genes by combined score |
| `5D_top50_abc_only.tsv` | 5D | Top 50 genes by ABC score only |

### Plots (`plots/figure_5_model_functional/`)

| File | Panel | Description |
|------|-------|-------------|
| `5A_go_bp_dotplot_boundaries.svg` | 5A | GO at differential boundaries |
| `5A_kegg_dotplot_boundaries.svg` | 5A | KEGG at boundaries |
| `5A_go_comparison_long_vs_short.svg` | 5A | GO: long-range lost vs short gained |
| `5A_go_bp_dotplot_abc.svg` | 5A | GO at paired ABC anchors |
| `5A_kegg_dotplot_abc.svg` | 5A | KEGG at ABC anchors |
| `5B_deg_tad_violin.svg` | 5B | DEG expression at TAD boundaries |
| `5B_deg_loop_anchor_violin.svg` | 5B | DEG expression at loop anchors |
| `5B_deg_loop_permutation_test.svg` | 5B | DEG-loop permutation test |
| `5B_deg_longrange_vs_shortrange.svg` | 5B | Long-range lost vs short-range gained |
| `5C_network_figure.svg` | 5C | Network visualization (Figure 5C) |
| `5C_layer_distribution.svg` | 5C | Layer distribution |
| `5C_edge_type_breakdown.svg` | 5C | Edge type breakdown |
| `5C_go_enrichment_dotplot.svg` | 5C | Network GO enrichment |
| `5D_combined_score_heatmap.svg` | 5D | Top 50 combined structural score |
| `5D_abc_only_heatmap.svg` | 5D | Top 50 ABC-only heatmap |

---

## Supplemental Analyses

### Data (`tsvs/supplemental/`)

| File | Analysis | Description |
|------|----------|-------------|
| `shared_anchor_loops.tsv` | Shared Anchors | Loops sharing anchors with opposing direction |
| `anchor_characterization.tsv` | Shared Anchors | Anchor ChIP-seq characterization |
| `shared_anchor_genes.tsv` | Shared Anchors | Genes at shared anchors |
| `shared_anchor_chip_enrichment.tsv` | Shared Anchors | ChIP enrichment at shared anchors |
| `shared_anchor_paired_distance.tsv` | Shared Anchors | Paired distance statistics |
| `shared_boundary_*.tsv` | Shared Boundaries | Boundary proximity analysis (11 files) |
| `loop_compartment_annotated.tsv` | Loop-Compartment | Loop x compartment crossref |
| `ctcf_stripe_*.tsv` | CTCF Stripes | CTCF stripe x loop crossref (6 files) |
| `timepoint_comparison_stats.tsv` | Timepoint | Early vs late comparison |
| `homer_motif_summary_stats.tsv` | HOMER | Motif enrichment at enhancer subsets |

### Plots (`plots/supplemental/`)

| File | Analysis | Description |
|------|----------|-------------|
| `shared_anchor_distance_violin.svg` | Shared Anchors | Distance distribution |
| `shared_anchor_chip_enrichment.svg` | Shared Anchors | ChIP enrichment dotplot |
| `loop_rewriting_summary.svg` | Loop Rewriting | Rewriting summary figure |
| `ctcf_stripe_*.svg` | CTCF Stripes | Stripe analysis plots (4 files) |
| `paired_anchor_panel.svg` | Paired Anchors | ABC-loop concordance panel |
| `apa_shared_*.svg` | APA Shared | APA for shared anchor loops (8 files) |
| `homer_narrative_composite.svg` | HOMER | Motif narrative composite |

---

## Upstream Reference Data (`upstream/`)

| File | Category | Description |
|------|----------|-------------|
| `rna_seq/adult_rnaseq_results.xlsx` | RNA-seq | BAP1 WT vs KO differential expression |
| `chip_peaks/k119ub_anchor_signal.tsv` | ChIP-seq | H2AK119ub signal at loop anchors |
| `chip_peaks/k119ub_enhancer_signal.tsv` | ChIP-seq | H2AK119ub signal at enhancers |
| `loop_calls/late_characterized_loops.tsv` | Loops | Merged characterized loops (adult) |
| `loop_calls/late_merged_final.bedpe` | Loops | Final loop calls BEDPE (adult) |
| `loop_calls/early_merged_final.bedpe` | Loops | Final loop calls BEDPE (P13) |

---

## Notes

- **Panels not included:** 1A, 1E (Cntnap5a), 2D, 3A (deepTools ARA), 3E (H2Az), 4C (K27ac/ATAC volcanos), 5E (schematic) require HPC-generated screenshots, deepTools heatmaps, or BioRender schematics not stored in this repository.
- **Plot formats:** SVG preferred for editability. JPG used for compartment/loop volcanos and heavy scatter plots. PDF used for TAD volcanos (only format available).
- **Duplicate data:** Some files appear in multiple figure folders (e.g., boundary chromatin state in both Fig 1F and 3B) for self-contained reference.
- **Original locations:** See `docs/Hi-C-paper-annotated.md` for complete source paths and generation scripts.
