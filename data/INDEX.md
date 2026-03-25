# Paper Data & Figure Index

Consolidated ground-truth data and publication figures for the BAP1 Hi-C paper. All files are **copies** from their original analysis locations. This directory serves as a single reference point for reviewers and cross-codebase lookup.

**Generated:** 2026-03-21 | **Updated:** 2026-03-24 (added script provenance)
**Source mapping:** `docs/Hi-C-paper-annotated.md` (complete provenance for each panel)

---

## Directory Structure

```
data/
├── tsvs/                              # Tabular data files (TSV, BED, TXT)
│   ├── figure_1_tads_boundaries_compartments/   (19 files)
│   ├── figure_2_loop_rewiring/                  (38 files)
│   ├── figure_3_epigenetic_integration/         (15 files)
│   ├── figure_4_abc_analysis/                   (19 files)
│   ├── figure_5_model_functional/               (18 files)
│   └── supplemental/                            (25 files)
├── plots/                             # Publication figures (SVG preferred, JPG for heavy plots)
│   ├── figure_1_tads_boundaries_compartments/   (12 files)
│   ├── figure_2_loop_rewiring/                  (64 files)
│   ├── figure_3_epigenetic_integration/         (28 files)
│   ├── figure_4_abc_analysis/                   (24 files)
│   ├── figure_5_model_functional/               (15 files)
│   └── supplemental/                            (17 files)
├── scripts/                           # Copies of figure-generating scripts (by figure)
│   ├── _shared/                                 (1 file: multi_format_output.R)
│   ├── figure_1_tads_boundaries_compartments/   (6 scripts)
│   ├── figure_2_loop_rewiring/                  (9 scripts)
│   ├── figure_3_epigenetic_integration/         (5 scripts)
│   ├── figure_4_abc_analysis/                   (8 scripts)
│   ├── figure_5_model_functional/               (4 scripts)
│   └── supplemental/                            (6 scripts)
├── upstream/                          # Reference input data
│   ├── rna_seq/                                 (1 file)
│   ├── chip_peaks/                              (2 files)
│   ├── diffbind/                                (4 files)
│   └── loop_calls/                              (3 files)
└── INDEX.md                           # This file
```

**Totals:** 140 data files + 160 plot files + 39 scripts + 10 upstream = **349 files**

## Naming Convention

Files are named `{panel}_{descriptive_name}.{ext}`:
- **Panel prefix:** `1B`, `2A`, `3C`, etc. maps to the paper figure panel
- **Timepoint prefix:** `early` = P13, `late` = adult
- **Format:** SVG for most plots; JPG for heavy volcano/scatter plots with many points; PDF when SVG/JPG unavailable

---

## Figure 1: Progressive Intensification of 3D Genome Phenotype

TADs, boundaries, and compartment changes between BAP1-KO and wildtype.

### Data (`tsvs/figure_1_tads_boundaries_compartments/`)

| File | Panel | Description | Source | Script |
|------|-------|-------------|--------|--------|
| `1B_early_tad_all_annotated.tsv` | 1B | All TADs with differential statistics (P13, 5074 TADs) | `tads/tad-pc-analysis/output/tad_analysis/early/` | `scripts/tad_volcano_plot.R` |
| `1B_early_tad_significant_relaxed.tsv` | 1B | Significant TADs, P13 (FDR<0.15, \|Diff\|>0.15; 45 TADs) | same | `scripts/tad_volcano_plot.R` |
| `1B_early_tad_significant_standard.tsv` | 1B | Significant TADs, P13 (FDR<0.05, \|Diff\|>0.30; 0 TADs) | same | `scripts/tad_volcano_plot.R` |
| `1B_late_tad_all_annotated.tsv` | 1B | All TADs with differential statistics (Adult, 4366 TADs) | `tads/tad-pc-analysis/output/tad_analysis/late/` | `scripts/tad_volcano_plot.R` |
| `1B_late_tad_significant_relaxed.tsv` | 1B | Significant TADs, Adult (FDR<0.15, \|Diff\|>0.15; 1283 TADs) | same | `scripts/tad_volcano_plot.R` |
| `1B_late_tad_significant_standard.tsv` | 1B | Significant TADs, Adult (FDR<0.05, \|Diff\|>0.30; 355 TADs) | same | `scripts/tad_volcano_plot.R` |
| `1B_early_tadcompare_differential.tsv` | 1B | TADCompare results, P13 timepoint | `tads/results/early/tadcompare/` | `tads/scripts/02_run_tadcompare.R` |
| `1B_late_tadcompare_differential.tsv` | 1B | TADCompare results, adult timepoint | `tads/results/late/tadcompare/` | `tads/scripts/02_run_tadcompare.R` |
| `1C_early_differential_boundaries.bed` | 1C | Differential boundaries (P13) | `tads/results/early/final/` | `tads/scripts/tad_visualizations.R` |
| `1C_late_differential_boundaries.bed` | 1C | Differential boundaries (adult) | `tads/results/late/final/` | `tads/scripts/tad_visualizations.R` |
| `1C_late_shifted_boundaries.bed` | 1C | Shifted boundaries (min 2 bins) | `tads/results/late/final/` | `tads/scripts/shift_analysis_filtered.R` |
| `1D_compartment_all_annotated.tsv` | 1D | All compartments with PC1 differential | `tads/tad-pc-analysis/output/compartment_analysis/` | `scripts/compartment_volcano_plot.R` |
| `1D_compartment_significant_relaxed.tsv` | 1D | Significant compartment changes (relaxed) | same | `scripts/compartment_volcano_plot.R` |
| `1D_compartment_significant_standard.tsv` | 1D | Significant compartment changes (standard) | same | `scripts/compartment_volcano_plot.R` |
| `1D_compartment_genome_pct_summary.txt` | 1D | % genome with differential PC1 | same | `scripts/compartment_genome_percentage.R` |
| `1E_syt1_nav3_statistics.tsv` | 1E | Syt1/Nav3 locus statistics | `tads/results/visualizations/late/syt1_nav3_focus/` | `tads/scripts/tad_visualizations.R` |
| `1F_permutation_test_results.tsv` | 1F | Boundary-loop permutation test | `tads/results/late/boundary_loop_analysis/` | `tads/scripts/boundary_loop_crossref.R` |
| `1F_enrichment_statistics.tsv` | 1F | ChIP enrichment at boundaries | same | `tads/scripts/boundary_loop_crossref.R` |
| `1F_boundaries_with_chromatin_state.tsv` | 1F | Boundary chromatin state annotations | `tads/results/visualizations/chip/late/` | `tads/scripts/tad_chip_classification.R` |

### Plots (`plots/figure_1_tads_boundaries_compartments/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `1B_early_tad_volcano_relaxed.svg` | 1B | TAD volcano, P13 (relaxed thresholds) | `scripts/tad_volcano_plot.R` |
| `1B_early_tad_volcano_standard.svg` | 1B | TAD volcano, P13 (standard thresholds) | `scripts/tad_volcano_plot.R` |
| `1B_late_tad_volcano_relaxed.svg` | 1B | TAD volcano, Adult (relaxed thresholds) | `scripts/tad_volcano_plot.R` |
| `1B_late_tad_volcano_standard.svg` | 1B | TAD volcano, Adult (standard thresholds) | `scripts/tad_volcano_plot.R` |
| `1D_compartment_volcano_relaxed.jpg` | 1D | Compartment volcano (JPG - heavy plot) | `scripts/compartment_volcano_plot.R` |
| `1D_compartment_volcano_standard.jpg` | 1D | Compartment volcano (JPG - heavy plot) | `scripts/compartment_volcano_plot.R` |
| `1D_compartment_genome_pct_pie.svg` | 1D | % genome pie chart by compartment shift | `scripts/compartment_genome_percentage.R` |
| `1D_compartment_genome_pct_bar.svg` | 1D | % genome bar chart | `scripts/compartment_genome_percentage.R` |
| `1E_syt1_nav3_regional_overview.svg` | 1E | Syt1/Nav3 regional focus | `tads/scripts/tad_visualizations.R` |
| `1F_permutation_test.svg` | 1F | Boundary-loop permutation test | `tads/scripts/boundary_loop_crossref.R` |
| `1F_chip_mark_overlap_heatmap.svg` | 1F | ChIP mark overlap at boundaries | `tads/scripts/tad_chip_classification.R` |
| `1F_chromatin_state_by_differential.svg` | 1F | Chromatin state by differential status | `tads/scripts/tad_chip_classification.R` |

---

## Figure 2: Chromatin Loop Rewiring

Differential loop analysis, APA, anchor annotation.

### Data (`tsvs/figure_2_loop_rewiring/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `2A_{early,late}_edgeR_{5,10,25}kb.tsv` | 2A | edgeR results per resolution (6 files) | `scripts/edgeR.R` |
| `2B_late_distance_shift_summary.tsv` | 2B | Loop distance shift summary | `scripts/loop_distance_analysis.R` |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_enrichment.tsv` | 2C | APA enrichment scores (12 files) | `scripts/apa_analysis.R` |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_stats.tsv` | 2C | APA statistical tests (12 files) | `scripts/apa_analysis.R` |
| `2F_{early,late}_anchor_type_summary.tsv` | 2F | Anchor category counts (2 files) | `peaks/scripts/annotate_loops_extended.R` |
| `2G_{early,late}_loop_type_summary.tsv` | 2G | Loop type pair counts (2 files) | `peaks/scripts/annotate_loops_extended.R` |
| `2G_{early,late}_extended_characterized_loops.tsv` | 2G | 7-category annotated loops (2 files) | `peaks/scripts/annotate_loops_extended.R` |
| `2H_late_characterized_loops.tsv` | 2H | Merged characterized loops (adult) | `scripts/downstream_analysis.R` |

### Plots (`plots/figure_2_loop_rewiring/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `2A_{early,late}_volcano_merged.jpg` | 2A | Multi-resolution merged volcano (JPG - heavy) | `scripts/visualizations.R` |
| `2B_{early,late}_distance_cdf.svg` | 2B | Lost vs gained loop eCDF | `scripts/loop_distance_analysis.R` |
| `2C_{early,late}_apa_{5,10,25}kb_{up,down}_{heatmap_type}.svg` | 2C | APA heatmaps (48 files in `apa/` subdir) | `scripts/apa_analysis.R` |
| `2E_late_cdf_k27me3_anchored.svg` | 2E | K27me3-anchored loop eCDF | `scripts/loop_distance_k27me3_filtered.R` |
| `2E_late_cdf_H3K27ac_one_anchor.svg` | 2E | H3K27ac-anchored loop eCDF | `scripts/loop_distance_mark_filtered.R` |
| `2E_late_cdf_H3K27me3_one_anchor.svg` | 2E | H3K27me3-anchored loop eCDF | `scripts/loop_distance_mark_filtered.R` |
| `2F_{early,late}_anchor_type_distribution.svg` | 2F | Anchor type bar plots | `peaks/scripts/annotate_loops_extended.R` |
| `2G_late_loop_type_by_direction.svg` | 2G | Loop categories by gain/loss | `peaks/scripts/annotate_loops_extended.R` |
| `2G_late_loop_type_piechart.svg` | 2G | Loop type pie chart | `peaks/scripts/annotate_loops_extended.R` |
| `2H_late_looptype_distance_heatmap.svg` | 2H | Distance x loop type heatmap | `scripts/loop_distance_analysis.R` |
| `2H_late_deg_violin_by_loop_type.svg` | 2H | DEG expression by loop type | `scripts/deg_loop_anchor_violin.R` |
| `2I_late_cdf_{H3K27ac,H3K27me3,CTCF}.svg` | 2I | Mark-specific eCDFs (3 files) | `scripts/loop_distance_mark_filtered.R` |

---

## Figure 3: Integration with Epigenetic Changes

H2AK119ub, ChIP-seq at boundaries, permutation tests, timepoint comparison.

### Data (`tsvs/figure_3_epigenetic_integration/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `3A_k119ub_anchor_signal.tsv` | 3A | K119ub signal at loop anchors | `scripts/h2ak119ub_loop_integration.R` |
| `3A_k119ub_enrichment_global.tsv` | 3A | Global K119ub enrichment | `scripts/h2ak119ub_loop_integration.R` |
| `3A_k119ub_enrichment_by_chromstate.tsv` | 3A | K119ub enrichment by chromatin state | `scripts/h2ak119ub_loop_integration.R` |
| `3A_k119ub_correlation_results.tsv` | 3A | K119ub-loop correlation | `scripts/h2ak119ub_loop_integration.R` |
| `3A_k119ub_logistic_regression.tsv` | 3A | Logistic regression for K119ub | `scripts/h2ak119ub_loop_integration.R` |
| `3A_loops_with_k119ub.tsv` | 3A | Loops annotated with K119ub | `scripts/h2ak119ub_loop_integration.R` |
| `3B_boundaries_chromatin_state.tsv` | 3B | TAD boundaries with chromatin states | `tads/scripts/tad_chip_classification.R` |
| `3C_boundary_permutation_results.tsv` | 3C | Boundary permutation test | `tads/scripts/boundary_loop_crossref.R` |
| `3C_enrichment_tests_all_loops.tsv` | 3C | Diff ChIP enrichment (all loops) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_enrichment_tests_polycomb.tsv` | 3C | Diff ChIP enrichment (polycomb) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_overlap_summary_all_loops.tsv` | 3C | ChIP overlap summary (all) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_overlap_summary_polycomb.tsv` | 3C | ChIP overlap summary (polycomb) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_loop_compartment_annotated.tsv` | 3C | Loops with compartment annotation | `scripts/loop_compartment_crossref.R` |
| `3D_timepoint_comparison_stats.tsv` | 3D | Early vs late comparison statistics | `tads/scripts/timepoint_comparison.R` |
| `3F_syt1_nav3_statistics.tsv` | 3F | Syt1/Nav3 locus stats | `tads/scripts/tad_visualizations.R` |

### Plots (`plots/figure_3_epigenetic_integration/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `3A_cdf_k119ub_up_one_anchor.svg` | 3A | K119ub CDF at gained loop anchors | `scripts/h2ak119ub_loop_integration.R` |
| `3A_enrichment_dotplot_by_chromstate.svg` | 3A | K119ub enrichment by chromatin state | `scripts/h2ak119ub_loop_integration.R` |
| `3A_scatter_loopFC_vs_k119ub.svg` | 3A | Loop logFC vs K119ub logFC | `scripts/h2ak119ub_loop_integration.R` |
| `3A_boxplot_k119ub_by_direction.svg` | 3A | K119ub by loop direction | `scripts/h2ak119ub_loop_integration.R` |
| `3A_correlation_summary_heatmap.svg` | 3A | Correlation summary | `scripts/h2ak119ub_loop_integration.R` |
| `3B_01_boundary_chromatin_state_distribution.svg` | 3B | Boundary chromatin state distribution | `tads/scripts/tad_chip_classification.R` |
| `3B_02_chromatin_state_by_differential.svg` | 3B | Chromatin state by differential status | `tads/scripts/tad_chip_classification.R` |
| `3B_03_chromatin_state_by_boundary_type.svg` | 3B | Chromatin state by boundary type | `tads/scripts/tad_chip_classification.R` |
| `3B_04_chromatin_state_by_enrichment.svg` | 3B | Chromatin state by enrichment direction | `tads/scripts/tad_chip_classification.R` |
| `3B_05_chip_mark_overlap_heatmap.svg` | 3B | ChIP mark overlap heatmap | `tads/scripts/tad_chip_classification.R` |
| `3B_06_summary_comparison.svg` | 3B | Summary comparison | `tads/scripts/tad_chip_classification.R` |
| `3C_boundary_permutation_test.svg` | 3C | Boundary permutation test | `tads/scripts/boundary_loop_crossref.R` |
| `3C_enrichment_dotplot_all_loops.svg` | 3C | Diff ChIP enrichment (all loops) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_enrichment_dotplot_polycomb.svg` | 3C | Diff ChIP enrichment (polycomb) | `scripts/diff_chip_polycomb_enrichment.R` |
| `3C_loop_cmpt_q1_shift_overlap_barplot.svg` | 3C | Compartment shift overlap | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q2_ctrl_pc1_violin.svg` | 3C | Control PC1 violin | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q2_distance_vs_pc1diff_scatter.svg` | 3C | Distance vs PC1 diff scatter | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q2_pc1_diff_violin.svg` | 3C | PC1 diff violin | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q3_tad_ir_boxplot.svg` | 3C | TAD inclusion ratio boxplot | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q4_distance_pc1_violin.svg` | 3C | Distance-PC1 violin | `scripts/loop_compartment_crossref.R` |
| `3C_loop_cmpt_q5_polycomb_pc1_violin.svg` | 3C | Polycomb PC1 violin | `scripts/loop_compartment_crossref.R` |
| `3D_boundary_type_comparison.svg` | 3D | Boundary type comparison | `tads/scripts/timepoint_comparison.R` |
| `3D_combined_comparison_summary.svg` | 3D | Combined comparison summary | `tads/scripts/timepoint_comparison.R` |
| `3D_enrichment_direction_comparison.svg` | 3D | Enrichment direction comparison | `tads/scripts/timepoint_comparison.R` |
| `3D_enrichment_direction_counts.svg` | 3D | Enrichment direction counts | `tads/scripts/timepoint_comparison.R` |
| `3D_enrichment_direction_flip.svg` | 3D | Enrichment direction flip | `tads/scripts/timepoint_comparison.R` |
| `3D_net_direction_diverging.svg` | 3D | Net direction diverging bar | `tads/scripts/timepoint_comparison.R` |
| `3F_syt1_nav3_regional_overview.svg` | 3F | Syt1/Nav3 integrated locus | `tads/scripts/tad_visualizations.R` |

---

## Figure 4: Enhancer ABC Analysis

Activity-by-Contact analysis, concordance, K119ub integration.

### Data (`tsvs/figure_4_abc_analysis/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `4A_delta_abc_all_pairs.tsv` | 4A | All enhancer-gene DABC pairs (~180K) | `abc/scripts/step12_activity_contact_scatter.R` |
| `4A_delta_abc_with_rnaseq.tsv` | 4A | DABC merged with RNA-seq (~113K) | `abc/scripts/step12_activity_contact_scatter.R` |
| `4A_activity_contact_summary.tsv` | 4A | Activity vs contact decomposition | `abc/scripts/step12_activity_contact_scatter.R` |
| `4B_discordant_gene_characteristics.tsv` | 4B | Discordant gene features | `abc/scripts/step13_discordant_gene_analysis.R` |
| `4B_concordance_by_class.tsv` | 4B | Concordance by enhancer class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4B_discordant_go_bp.tsv` | 4B | GO enrichment for discordant genes | `abc/scripts/step13b_go_enrichment.R` |
| `4B_discordant_kegg.tsv` | 4B | KEGG enrichment for discordant genes | `abc/scripts/step13b_go_enrichment.R` |
| `4D_class_level_summary.tsv` | 4D | Enhancer class summary stats | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4D_enhancer_classes_*.tsv` | 4D | Enhancer class assignments (4 files) | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_enhancer_class_abc_metrics.tsv` | 4E | ABC metrics by class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_enhancer_class_loop_metrics.tsv` | 4E | Loop metrics by class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_contact_decay_by_class.tsv` | 4E | Contact decay curves | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_k119ub_tertile_assignments.tsv` | 4F | K119ub tertile classification | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_k119ub_tertile_loop_summary.tsv` | 4F | Tertile loop summary | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_k119ub_abc_correlation_summary.tsv` | 4F | K119ub-ABC correlation | `abc/scripts/step10_k119ub_abc_correlation.R` |
| `4F_k119ub_abc_enhancer_merged.tsv` | 4F | Merged enhancer K119ub + ABC | `abc/scripts/step10_k119ub_abc_correlation.R` |

### Plots (`plots/figure_4_abc_analysis/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `4A_raw_delta_all_pairs.svg` | 4A | Raw delta activity vs contact | `abc/scripts/step12_activity_contact_scatter.R` |
| `4A_log2fc_all_pairs.jpg` | 4A | log2FC scatter (JPG - heavy) | `abc/scripts/step12_activity_contact_scatter.R` |
| `4A_raw_delta_classified.svg` | 4A | Classified by concordance category | `abc/scripts/step12_activity_contact_scatter.R` |
| `4A_raw_delta_promoter_distal.svg` | 4A | Promoter vs distal enhancers | `abc/scripts/step12b_promoter_distal_scatter.R` |
| `4B_concordance_pie_combined.svg` | 4B | Combined concordance pie chart | `scripts/concordance_pie_chart.R` |
| `4B_concordance_pie_4cat.svg` | 4B | 4-category concordance | `scripts/concordance_pie_chart.R` |
| `4B_discordant_composite.svg` | 4B | Discordant analysis composite | `abc/scripts/step13_discordant_gene_analysis.R` |
| `4B_discordant_go_bp.svg` | 4B | GO enrichment for discordant | `abc/scripts/step13b_go_enrichment.R` |
| `4B_k119ub_by_concordance.svg` | 4B | K119ub by concordance status | `abc/scripts/step13c_k119ub_concordance.R` |
| `4B_k119ub_significance_rate.svg` | 4B | K119ub significance rate | `abc/scripts/step13c_k119ub_concordance.R` |
| `4C_k119ub_volcano_at_enhancers.jpg` | 4C | K119ub volcano at enhancers (JPG) | `abc/scripts/step10_k119ub_abc_correlation.R` |
| `4D_summary_patchwork.svg` | 4D | Enhancer model summary | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4D_class_composition_bar.svg` | 4D | Class composition | `abc/scripts/step12_activity_contact_scatter.R` |
| `4E_loop_logfc_violin.svg` | 4E | Loop logFC by class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_delta_abc_boxplot.svg` | 4E | Delta ABC by class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_gene_logfc_by_class.svg` | 4E | Gene logFC by class | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_contact_decay_{wt,delta}.svg` | 4E | Contact decay curves (2 files) | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4E_logFC_vs_deltaABC.svg` | 4E | Loop logFC vs DABC | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `4F_k119ub_tertile_loop_logfc.svg` | 4F | K119ub tertile loop logFC | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_k119ub_tertile_delta_abc.svg` | 4F | K119ub tertile DABC | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_k119ub_vs_loop_logfc_scatter.svg` | 4F | K119ub vs loop logFC | `abc/scripts/step11_enhancer_subset_analysis.R` |
| `4F_contingency_heatmap.svg` | 4F | K119ub x ABC contingency | `abc/scripts/step10_k119ub_abc_correlation.R` |
| `4F_boxplot_k119ub_by_abc_category.svg` | 4F | K119ub by ABC category | `abc/scripts/step10_k119ub_abc_correlation.R` |

---

## Figure 5: Model and Functional Implications

GO analysis, DEG integration, network analysis, top gene heatmaps.

### Data (`tsvs/figure_5_model_functional/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `5A_boundary_genes.tsv` | 5A | Genes at differential boundaries | `tads/scripts/tad_visualizations.R` |
| `5A_go_long_lost_loops.tsv` | 5A | GO for long-range lost loops | `scripts/loop_distance_analysis.R` |
| `5A_go_short_gained_loops.tsv` | 5A | GO for short-range gained loops | `scripts/loop_distance_analysis.R` |
| `5A_abc_go_enrichment.tsv` | 5A | GO at paired ABC anchors | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5A_abc_kegg_enrichment.tsv` | 5A | KEGG at paired ABC anchors | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5B_deg_boundary_genes.tsv` | 5B | DEGs at TAD boundaries | `tads/scripts/deg_tad_violin.R` |
| `5B_deg_anchor_genes.tsv` | 5B | DEGs at loop anchors | `scripts/deg_loop_anchor_violin.R` |
| `5B_deg_longrange_vs_shortrange_genes.tsv` | 5B | Long-range lost vs short-range gained genes | `scripts/deg_loop_anchor_violin.R` |
| `5B_gene_level_summary.tsv` | 5B | Gene-level ABC summary (13.6K genes) | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5C_gene_structural_profile_filtered.tsv` | 5C | Network node data (filtered) | `scripts/network_analysis.R` |
| `5C_gene_structural_profile_all.tsv` | 5C | Network node data (all) | `scripts/network_analysis.R` |
| `5C_edge_list.tsv` | 5C | Network edge list | `scripts/network_analysis.R` |
| `5C_network_centrality_metrics.tsv` | 5C | Network centrality scores | `scripts/network_analysis.R` |
| `5C_go_enrichment_results.tsv` | 5C | GO enrichment for network genes | `scripts/network_analysis.R` |
| `5C_loops_with_gene_assignments.tsv` | 5C | Loop-gene mapping | `scripts/network_analysis.R` |
| `5C_paired_anchor_matches.tsv` | 5C | Paired anchor matches | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5D_top50_combined_score.tsv` | 5D | Top 50 genes by combined score | `scripts/structural_heatmap.R` |
| `5D_top50_abc_only.tsv` | 5D | Top 50 genes by ABC score only | `scripts/structural_heatmap.R` |

### Plots (`plots/figure_5_model_functional/`)

| File | Panel | Description | Script |
|------|-------|-------------|--------|
| `5A_go_bp_dotplot_boundaries.svg` | 5A | GO at differential boundaries | `tads/scripts/tad_visualizations.R` |
| `5A_kegg_dotplot_boundaries.svg` | 5A | KEGG at boundaries | `tads/scripts/tad_visualizations.R` |
| `5A_go_comparison_long_vs_short.svg` | 5A | GO: long-range lost vs short gained | `scripts/loop_distance_analysis.R` |
| `5A_go_bp_dotplot_abc.svg` | 5A | GO at paired ABC anchors | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5A_kegg_dotplot_abc.svg` | 5A | KEGG at ABC anchors | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `5B_deg_tad_violin.svg` | 5B | DEG expression at TAD boundaries | `tads/scripts/deg_tad_violin.R` |
| `5B_deg_loop_anchor_violin.svg` | 5B | DEG expression at loop anchors | `scripts/deg_loop_anchor_violin.R` |
| `5B_deg_loop_permutation_test.svg` | 5B | DEG-loop permutation test | `scripts/deg_loop_anchor_violin.R` |
| `5B_deg_longrange_vs_shortrange.svg` | 5B | Long-range lost vs short-range gained | `scripts/deg_loop_anchor_violin.R` |
| `5C_network_figure.svg` | 5C | Network visualization (Figure 5C) | `scripts/network_analysis.R` |
| `5C_layer_distribution.svg` | 5C | Layer distribution | `scripts/network_analysis.R` |
| `5C_edge_type_breakdown.svg` | 5C | Edge type breakdown | `scripts/network_analysis.R` |
| `5C_go_enrichment_dotplot.svg` | 5C | Network GO enrichment | `scripts/network_analysis.R` |
| `5D_combined_score_heatmap.svg` | 5D | Top 50 combined structural score | `scripts/structural_heatmap.R` |
| `5D_abc_only_heatmap.svg` | 5D | Top 50 ABC-only heatmap | `scripts/structural_heatmap.R` |

---

## Supplemental Analyses

### Data (`tsvs/supplemental/`)

| File | Analysis | Description | Script |
|------|----------|-------------|--------|
| `shared_anchor_loops.tsv` | Shared Anchors | Loops sharing anchors with opposing direction | `scripts/shared_anchor_analysis.R` |
| `anchor_characterization.tsv` | Shared Anchors | Anchor ChIP-seq characterization | `scripts/shared_anchor_analysis.R` |
| `shared_anchor_genes.tsv` | Shared Anchors | Genes at shared anchors | `scripts/shared_anchor_analysis.R` |
| `shared_anchor_chip_enrichment.tsv` | Shared Anchors | ChIP enrichment at shared anchors | `scripts/shared_anchor_analysis.R` |
| `shared_anchor_paired_distance.tsv` | Shared Anchors | Paired distance statistics | `scripts/shared_anchor_analysis.R` |
| `shared_boundary_*.tsv` | Shared Boundaries | Boundary proximity analysis (11 files) | `scripts/shared_anchor_boundary_analysis.R` |
| `loop_compartment_annotated.tsv` | Loop-Compartment | Loop x compartment crossref | `scripts/loop_compartment_crossref.R` |
| `ctcf_stripe_*.tsv` | CTCF Stripes | CTCF stripe x loop crossref (6 files) | `scripts/ctcf_stripe_crossref.R` |
| `timepoint_comparison_stats.tsv` | Timepoint | Early vs late comparison | `tads/scripts/timepoint_comparison.R` |
| `homer_motif_summary_stats.tsv` | HOMER | Motif enrichment at enhancer subsets | `abc/scripts/step11b_homer_motif_visualization.R` |

### Plots (`plots/supplemental/`)

| File | Analysis | Description | Script |
|------|----------|-------------|--------|
| `shared_anchor_distance_violin.svg` | Shared Anchors | Distance distribution | `scripts/shared_anchor_analysis.R` |
| `shared_anchor_chip_enrichment.svg` | Shared Anchors | ChIP enrichment dotplot | `scripts/shared_anchor_analysis.R` |
| `loop_rewriting_summary.svg` | Loop Rewriting | Rewriting summary figure | `scripts/loop_distance_analysis.R` |
| `ctcf_stripe_*.svg` | CTCF Stripes | Stripe analysis plots (4 files) | `scripts/ctcf_stripe_crossref.R` |
| `paired_anchor_panel.svg` | Paired Anchors | ABC-loop concordance panel | `abc/scripts/step9b_paired_anchor_analysis.R` |
| `apa_shared_*.svg` | APA Shared | APA for shared anchor loops (8 files) | `scripts/apa_shared_anchors.R` |
| `homer_narrative_composite.svg` | HOMER | Motif narrative composite | `abc/scripts/step11b_homer_motif_visualization.R` |

---

## Upstream Reference Data (`upstream/`)

| File | Category | Description |
|------|----------|-------------|
| `rna_seq/adult_rnaseq_results.xlsx` | RNA-seq | BAP1 WT vs KO differential expression (adult) |
| `rna_seq/young_rnaseq_results.xlsx` | RNA-seq | BAP1 WT vs KO differential expression (P13) |
| `chip_peaks/k119ub_anchor_signal.tsv` | ChIP-seq | H2AK119ub signal at loop anchors |
| `chip_peaks/k119ub_enhancer_signal.tsv` | ChIP-seq | H2AK119ub signal at enhancers |
| `loop_calls/late_characterized_loops.tsv` | Loops | Merged characterized loops (adult) |
| `loop_calls/early_characterized_loops.tsv` | Loops | Merged characterized loops (P13) |
| `loop_calls/late_merged_final.bedpe` | Loops | Final loop calls BEDPE (adult) |
| `loop_calls/early_merged_final.bedpe` | Loops | Final loop calls BEDPE (P13) |
| `homer/Bap1.diff.tad.txt` | HOMER | TAD differential expression input |
| `diffbind/ATAC_allATAC_diffbind_results_summit_appended_ap.txt` | DiffBind | ATAC-seq differential accessibility (all peaks) |
| `diffbind/K119ub_diffbind_results_summit_appended_ap.txt` | DiffBind | H2AK119ub differential binding |
| `diffbind/K27ac_diffbind_results_summit_appended_ap.txt` | DiffBind | H3K27ac differential binding |
| `diffbind/K27me3_diffbind_results_summit_appended_ap.txt` | DiffBind | H3K27me3 differential binding |
| `Go_term_selction.xlsx` | Reference | GO term selection for network analysis |

---

## Script Copies (`scripts/`)

Local copies of the figure-generating scripts, organized by figure. The **Script** column in tables above references the original repo path. The copies below are for self-contained reference; see `REFACTOR_PROMPT.md` for instructions on re-pathing them.

| Directory | Scripts | Original Locations |
|-----------|---------|-------------------|
| `_shared/` | `multi_format_output.R` | `scripts/utils/` |
| `figure_1_.../` | `tad_volcano_plot.R`, `tad_visualizations.R`, `compartment_volcano_plot.R`, `compartment_genome_percentage.R`, `tad_chip_classification.R`, `boundary_loop_crossref.R` | `scripts/`, `tads/scripts/` |
| `figure_2_.../` | `edgeR.R`, `visualizations.R`, `loop_distance_analysis.R`, `apa_analysis.R`, `loop_distance_k27me3_filtered.R`, `loop_distance_mark_filtered.R`, `annotate_loops_extended.R`, `annotate_loops_extended_peaks.R`, `deg_loop_anchor_violin.R` | `scripts/`, `peaks/scripts/` |
| `figure_3_.../` | `h2ak119ub_loop_integration.R`, `preprocess_k119ub_anchor_signal.R`, `diff_chip_polycomb_enrichment.R`, `loop_compartment_crossref.R`, `timepoint_comparison.R` | `scripts/`, `tads/scripts/` |
| `figure_4_.../` | `step12_activity_contact_scatter.R`, `step12b_promoter_distal_scatter.R`, `concordance_pie_chart.R`, `step13_discordant_gene_analysis.R`, `step13b_go_enrichment.R`, `step13c_k119ub_concordance.R`, `step10_k119ub_abc_correlation.R`, `step11_enhancer_subset_analysis.R` | `scripts/`, `abc/scripts/` |
| `figure_5_.../` | `deg_tad_violin.R`, `step9b_paired_anchor_analysis.R`, `network_analysis.R`, `structural_heatmap.R` | `tads/scripts/`, `abc/scripts/`, `scripts/` |
| `supplemental/` | `shared_anchor_analysis.R`, `polycomb_shared_anchor_analysis.R`, `shared_anchor_boundary_analysis.R`, `ctcf_stripe_crossref.R`, `apa_shared_anchors.R`, `step11b_homer_motif_visualization.R` | `scripts/`, `abc/scripts/` |

---

## Notes

- **Panels not included:** 1A, 1E (Cntnap5a), 2D, 3A (deepTools ARA), 3E (H2Az), 4C (K27ac/ATAC volcanos), 5E (schematic) require HPC-generated screenshots, deepTools heatmaps, or BioRender schematics not stored in this repository.
- **Plot formats:** SVG preferred for editability. JPG used for compartment/loop volcanos and heavy scatter plots.
- **Duplicate data:** Some files appear in multiple figure folders (e.g., boundary chromatin state in both Fig 1F and 3B) for self-contained reference.
- **Multi-figure scripts:** Some scripts produce outputs for multiple figures (e.g., `tad_visualizations.R` serves 1C, 1E, 3F, 5A). They are copied once to their primary figure directory.
- **Original locations:** See `docs/Hi-C-paper-annotated.md` for complete source paths and generation scripts.
