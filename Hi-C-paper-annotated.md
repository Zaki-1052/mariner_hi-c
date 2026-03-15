# BAP1 Hi-C Paper — Annotated Figure & Data Inventory

> **Purpose:** Maps every panel in the paper outline to existing figures, scripts, and data tables in this repository.
> Legend: ✅ = exists locally | 🖥️ = on HPC instance | ❓ = needs to be generated | 📐 = conceptual/schematic (Illustrator/BioRender)

---

## Figure 1: Progressive intensification of 3D genome phenotype — TADs, boundaries, representative compartments

### Panel 1A — Contact map of similar site in P13 vs adult where changes are present in both, but more subtle at P13, differential insulation track below each map

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 🖥️ | Juicebox screenshot from .hic files on HPC; insulation track from GENOVA or cooltools |
| **Data** | 🖥️ | `.hic` files at `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/{sample}.hic` |
| **Insulation scores** | 🖥️ | HOMER TAD scores or cooltools insulation; raw scores in `tads/tad-pc-analysis/inputs/late/` |

### Panel 1B — Volcano plot: Differential TAD IDs at P13 and adult

| Item | Status | Path |
|------|--------|------|
| **Figure (Adult/Late)** | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_volcano_relaxed.pdf` |
| | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_volcano_standard.pdf` |
| **Figure (P13/Early)** | ❓ | Needs equivalent early-timepoint TAD volcano (early TADCompare data exists but volcano not generated for HOMER-style differential) |
| **Data (Adult)** | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_significant_relaxed.tsv` |
| | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_significant_standard.tsv` |
| | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_all_annotated.tsv` |
| **Summary** | ✅ | `tads/tad-pc-analysis/output/tad_analysis/tad_volcano_summary.txt` |
| **Script** | ✅ | `scripts/tad_volcano_plot.R` |
| **TADCompare Results (Early)** | ✅ | `tads/results/early/tadcompare/tadcompare_differential_only.tsv` |
| **TADCompare Results (Late)** | ✅ | `tads/results/late/tadcompare/tadcompare_differential_only.tsv` |

### Panel 1C — GENOVA insulation plot: Differential TAD boundaries

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 🖥️ | GENOVA `insulation_plot()` — requires Hi-C matrices on HPC |
| **Boundary data (Early)** | ✅ | `tads/results/early/final/differential_boundaries_final.bed` |
| **Boundary data (Late)** | ✅ | `tads/results/late/final/differential_boundaries_final.bed` |
| **Shifted boundaries** | ✅ | `tads/results/late/final/shifted_boundaries_min2bin.bed` |
| **Visualization scripts** | ✅ | `tads/scripts/tad_visualizations.R` |

### Panel 1D — % of genome with differential PC1 + smaller contact map of obvious differential region

| Item | Status | Path |
|------|--------|------|
| **Compartment volcano (Adult)** | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_volcano_relaxed.pdf` |
| | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_volcano_standard.pdf` |
| **Compartment data** | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_significant_relaxed.tsv` |
| | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_significant_standard.tsv` |
| | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv` |
| **Summary** | ✅ | `tads/tad-pc-analysis/output/compartment_analysis/compartment_volcano_summary.txt` |
| **Script** | ✅ | `scripts/compartment_volcano_plot.R` |
| **Contact map** | 🖥️ | Juicebox screenshot at differential compartment region |
| **% genome calculation** | ✅ | `compartment_genome_pct_bar/`, `compartment_genome_pct_pie/`, `compartment_genome_pct_by_chr/`, `compartment_genome_percentage_summary.txt` |

### Panel 1E — Representative locus for differential TAD, triangular (probably Cntnap5a)

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 🖥️ | Juicebox triangular heatmap at Cntnap5a locus |
| **Syt1/Nav3 regional focus** | ✅ | `tads/results/visualizations/late/syt1_nav3_focus/syt1_nav3_regional_overview/syt1_nav3_regional_overview.pdf` |
| **Syt1/Nav3 stats** | ✅ | `tads/results/visualizations/late/syt1_nav3_statistics.tsv` |

### Panel 1F — Permutation test of differential TADs and differential boundaries in K27ac/K27me3/K119ub

| Item | Status | Path |
|------|--------|------|
| **Boundary-loop permutation (Late)** | ✅ | `tads/results/late/boundary_loop_analysis/plots/permutation_test/permutation_test.pdf` |
| **Boundary-loop permutation data** | ✅ | `tads/results/late/boundary_loop_analysis/permutation_test_results.tsv` |
| **ChIP at boundaries (Late)** | ✅ | `tads/results/visualizations/chip/late/05_chip_mark_overlap_heatmap/05_chip_mark_overlap_heatmap.pdf` |
| | ✅ | `tads/results/visualizations/chip/late/02_chromatin_state_by_differential/02_chromatin_state_by_differential.pdf` |
| **ChIP at boundaries (Early)** | ✅ | `tads/results/visualizations/chip/early/05_chip_mark_overlap_heatmap/05_chip_mark_overlap_heatmap.pdf` |
| **Enrichment stats (Late)** | ✅ | `tads/results/late/boundary_loop_analysis/enrichment_statistics.tsv` |
| **Boundary chromatin state data** | ✅ | `tads/results/visualizations/chip/late/boundaries_with_chromatin_state.tsv` |
| **Diff ChIP-polycomb enrichment** | ✅ | `output/diff_chip_polycomb_enrichment/` (enrichment dotplots, heatmaps) |
| **Scripts** | ✅ | `tads/scripts/tad_chip_classification.R`, `scripts/diff_chip_polycomb_enrichment.R` |

---

## Figure 2: Chromatin loop rewiring — differential loop analysis, APA, annotation w/ ctrl peaks

### Panel 2A — Volcano plots: P13 and adult differential loop

| Item | Status | Path |
|------|--------|------|
| **Adult volcano (per-resolution)** | ✅ | `outputs/250402-late_outputs/edgeR_results_res_{5,10,25}kb/plots/volcano_plot_primary.pdf` |
| **Adult volcano (merged)** | ✅ | `outputs/250402-late_outputs/visualizations/volcano/volcano_merged_multiresolution.pdf` |
| **P13 volcano (per-resolution)** | ✅ | `outputs/250831-early_outputs/edgeR_results_res_{5,10,25}kb/plots/volcano_plot_primary.pdf` |
| **P13 volcano (merged)** | ✅ | `outputs/250831-early_outputs/visualizations/volcano/volcano_merged_multiresolution.pdf` |
| **Data (Adult)** | ✅ | `outputs/250402-late_outputs/edgeR_results_res_{5,10,25}kb/primary_analysis/all_results_primary.tsv` |
| **Data (P13)** | ✅ | `outputs/250831-early_outputs/edgeR_results_res_{5,10,25}kb/primary_analysis/all_results_primary.tsv` |
| **Script** | ✅ | `scripts/edgeR.R`, `scripts/visualizations.R` |

### Panel 2B — Lost vs gained loop eCDF at both P13 and adult

| Item | Status | Path |
|------|--------|------|
| **Figure (Adult)** | ✅ | `output/loops_visualization_extended/late/01_distance_cdf_by_direction/01_distance_cdf_by_direction.pdf` |
| **Figure (P13)** | ✅ | `output/loops_visualization_extended/early/01_distance_cdf_by_direction/01_distance_cdf_by_direction.pdf` |
| **Distance density (Adult)** | ✅ | `output/loops_visualization_extended/late/03_distance_density_overlay/03_distance_density_overlay.pdf` |
| **Distance density (P13)** | ✅ | `output/loops_visualization_extended/early/03_distance_density_overlay/03_distance_density_overlay.pdf` |
| **Summary data** | ✅ | `output/loops_visualization_extended/late/distance_shift_summary.tsv` |
| **Statistics** | ✅ | `output/loops_visualization_extended/late/distance_shift_statistics.txt` |
| **Script** | ✅ | `scripts/loop_distance_analysis.R` |

### Panel 2C — APA analysis at both timepoints

| Item | Status | Path |
|------|--------|------|
| **APA heatmaps (Adult, up loops)** | ✅ | `outputs/250402-late_outputs/apa_results/res_{5,10,25}kb/merged_loops/up_loops/` |
| **APA heatmaps (Adult, down loops)** | ✅ | `outputs/250402-late_outputs/apa_results/res_{5,10,25}kb/merged_loops/down_loops/` |
| **APA heatmaps (P13, up loops)** | ✅ | `outputs/250831-early_outputs/apa_results/res_{5,10,25}kb/merged_loops/up_loops/` |
| **APA heatmaps (P13, down loops)** | ✅ | `outputs/250831-early_outputs/apa_results/res_{5,10,25}kb/merged_loops/down_loops/` |
| **Per-heatmap files** | ✅ | `aggregate_heatmap_ctrl.pdf`, `aggregate_heatmap_mut.pdf`, `difference_heatmap.pdf`, `enrichment_boxplot.pdf` |
| **Enrichment scores** | ✅ | `enrichment_scores.tsv` (per subset) |
| **Statistical tests** | ✅ | `statistical_tests.tsv` (per subset) |
| **Shared anchor APA (Adult)** | ✅ | `output/apa_shared_anchors/late/` (gained short-range + lost long-range) |
| **Script** | ✅ | `scripts/apa_analysis.R`, `scripts/apa_shared_anchors.R` |

### Panel 2D — Example loci of gained and lost loops

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 🖥️ | Juicebox screenshots with BEDPE overlays |
| **BEDPE (Adult)** | ✅ | `outputs/250402-late_outputs/bedpe_final/merged_final.bedpe` |
| **BEDPE (P13)** | ✅ | `outputs/250831-early_outputs/bedpe_final/merged_final.bedpe` |

### Panel 2E — Lost loop length distribution annotation for K27me3 & K27ac → lost loops generally divided b/tw active enhancer loops and long-range polycomb loops

| Item | Status | Path |
|------|--------|------|
| **K27me3-anchored eCDF (Adult)** | ✅ | `output/loops_k27me3_filtered/late/01_cdf_k27me3_anchored/01_cdf_k27me3_anchored.pdf` |
| **K27me3-both eCDF (Adult)** | ✅ | `output/loops_k27me3_filtered/late/01_cdf_k27me3_both/01_cdf_k27me3_both.pdf` |
| **Bivalent eCDF (Adult)** | ✅ | `output/loops_k27me3_filtered/late/01_cdf_bivalent/01_cdf_bivalent.pdf` |
| **K27ac mark-filtered (Adult)** | ✅ | `output/loops_mark_filtered/late/H3K27ac/` (one_anchor + both_anchors CDF/density) |
| **K27me3 mark-filtered (Adult)** | ✅ | `output/loops_mark_filtered/late/H3K27me3/` |
| **K27me3 global (all loops, Adult)** | ✅ | `output/loops_k27me3_global/late/01_cdf_k27me3_anchored/01_cdf_k27me3_anchored.pdf` |
| **K27me3-anchored eCDF (P13)** | ✅ | `output/loops_k27me3_filtered/early/01_cdf_k27me3_anchored/01_cdf_k27me3_anchored.pdf` |
| **Scripts** | ✅ | `scripts/loop_distance_k27me3_filtered.R`, `scripts/loop_distance_k27me3_global.R`, `scripts/loop_distance_mark_filtered.R` |

### Panel 2F — Genomic anchor distribution (combined anchor 1&2)

| Item | Status | Path |
|------|--------|------|
| **Anchor type distribution (Adult)** | ✅ | `peaks/loop_annotation_extended/late/plots/anchor_type_distribution/anchor_type_distribution.pdf` |
| **Anchor type distribution (P13)** | ✅ | `peaks/loop_annotation_extended/early/plots/anchor_type_distribution/anchor_type_distribution.pdf` |
| **Anchor type summary (Adult)** | ✅ | `peaks/loop_annotation_extended/late/anchor_type_summary.tsv` |
| **Anchor type summary (P13)** | ✅ | `peaks/loop_annotation_extended/early/anchor_type_summary.tsv` |
| **Script** | ✅ | `peaks/scripts/annotate_loops_extended.R`, `scripts/annotate_loops_extended.R` |

### Panel 2G — Bar plot: Number of loops by category

| Item | Status | Path |
|------|--------|------|
| **Loop type by direction (Adult)** | ✅ | `peaks/loop_annotation_extended/late/plots/loop_type_by_direction/loop_type_by_direction.pdf` |
| **Loop type by direction (P13)** | ✅ | `peaks/loop_annotation_extended/early/plots/loop_type_by_direction/loop_type_by_direction.pdf` |
| **Pie chart comparison (Adult)** | ✅ | `peaks/loop_annotation_extended/late/plots/loop_type_piechart_comparison/loop_type_piechart_comparison.pdf` |
| **Loop type summary data (Adult)** | ✅ | `peaks/loop_annotation_extended/late/loop_type_summary.tsv` |
| **Loop type summary data (P13)** | ✅ | `peaks/loop_annotation_extended/early/loop_type_summary.tsv` |

### Panel 2H — Loop strength changes by category

| Item | Status | Path |
|------|--------|------|
| **Loop type distance heatmap (Adult)** | ✅ | `output/loops_visualization_extended/late/06_looptype_distance_heatmap/06_looptype_distance_heatmap.pdf` |
| **DEG loop violin by loop type (Adult)** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_violin_by_loop_type/deg_loop_violin_by_loop_type.pdf` |
| **DEG loop violin by chromatin state** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_anchor_violin_by_chromatin_state/deg_loop_anchor_violin_by_chromatin_state.pdf` |
| **ChIP distance trends (Adult)** | ✅ | `output/chip_distance/late/02_mark_trend_lineplot/02_mark_trend_lineplot.pdf` |
| **Script** | ✅ | `scripts/loop_distance_analysis.R`, `scripts/deg_loop_anchor_violin.R` |

### Panel 2I — Lost vs gained loop eCDF for K27ac-anchored, K27me3-anchored, CTCF-CTCF loops

| Item | Status | Path |
|------|--------|------|
| **K27ac-anchored CDF (Adult)** | ✅ | `output/loops_mark_filtered/late/H3K27ac/01_cdf_one_anchor/01_cdf_one_anchor.pdf` |
| **K27me3-anchored CDF (Adult)** | ✅ | `output/loops_mark_filtered/late/H3K27me3/01_cdf_one_anchor/01_cdf_one_anchor.pdf` |
| **CTCF CDF (Adult)** | ✅ | `output/loops_mark_filtered/late/CTCF/01_cdf_one_anchor/01_cdf_one_anchor.pdf` |
| **CTCF stripe cross-ref (Adult)** | ✅ | `output/ctcf_stripe_crossref/late/plots/` (enrichment, distance, permutation, venn) |
| **CTCF stripe statistics** | ✅ | `output/ctcf_stripe_crossref/late/tables/ctcf_stripe_statistics.txt` |
| **Script** | ✅ | `scripts/loop_distance_mark_filtered.R`, `scripts/ctcf_stripe_crossref.R` |

---

## Figure 3: Integration with epigenetic changes

### Panel 3A — ARA for ATAC, H3K27me3, H3K27ac, K119ub

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 🖥️ | Anchor Region Analysis — deepTools `computeMatrix` + `plotHeatmap` on HPC bigwigs |
| **H2AK119ub integration** | ✅ | `output/h2ak119ub_loop_integration/late/` (19 plots + 3 TSVs) |
| **H2AK119ub CDF plots** | ✅ | `output/h2ak119ub_loop_integration/late/01_cdf_k119ub_up_one_anchor/` |
| **H2AK119ub enrichment** | ✅ | `output/h2ak119ub_loop_integration/late/09_enrichment_dotplot_by_chromatin_state/` |
| **H2AK119ub scatter** | ✅ | `output/h2ak119ub_loop_integration/late/12_scatter_loopFC_vs_k119ub_FC/` |
| **K119ub enrichment data** | ✅ | `output/h2ak119ub_loop_integration/late/k119ub_enrichment_summary.tsv` |
| **K119ub signal data** | ✅ | `data/k119ub_anchor_signal.tsv` |
| **Scripts** | ✅ | `scripts/h2ak119ub_loop_integration.R`, `scripts/preprocess_k119ub_anchor_signal.R` |

### Panel 3B — Histone marks at differential TAD boundaries

| Item | Status | Path |
|------|--------|------|
| **ChIP state distribution (Adult)** | ✅ | `tads/results/visualizations/chip/late/01_boundary_chromatin_state_distribution/01_boundary_chromatin_state_distribution.pdf` |
| **ChIP by differential status** | ✅ | `tads/results/visualizations/chip/late/02_chromatin_state_by_differential/02_chromatin_state_by_differential.pdf` |
| **ChIP by boundary type** | ✅ | `tads/results/visualizations/chip/late/03_chromatin_state_by_boundary_type/03_chromatin_state_by_boundary_type.pdf` |
| **ChIP by enrichment direction** | ✅ | `tads/results/visualizations/chip/late/04_chromatin_state_by_enrichment/04_chromatin_state_by_enrichment.pdf` |
| **Overlap heatmap** | ✅ | `tads/results/visualizations/chip/late/05_chip_mark_overlap_heatmap/05_chip_mark_overlap_heatmap.pdf` |
| **Summary comparison** | ✅ | `tads/results/visualizations/chip/late/06_summary_comparison/06_summary_comparison.pdf` |
| **Data** | ✅ | `tads/results/visualizations/chip/late/boundaries_with_chromatin_state.tsv` |
| **Script** | ✅ | `tads/scripts/tad_chip_classification.R` |

### Panel 3C — Permutation tests: adult differential peaks in adult structures → B compartment changes largely indirect of epigenetic changes

| Item | Status | Path |
|------|--------|------|
| **Boundary permutation (Adult)** | ✅ | `tads/results/late/boundary_loop_analysis/plots/permutation_test/permutation_test.pdf` |
| **Boundary permutation data** | ✅ | `tads/results/late/boundary_loop_analysis/permutation_test_results.tsv` |
| **Diff ChIP enrichment (all loops)** | ✅ | `output/diff_chip_polycomb_enrichment/plots/all_loops/02_enrichment_dotplot_all_loops/02_enrichment_dotplot_all_loops.pdf` |
| **Diff ChIP enrichment (polycomb)** | ✅ | `output/diff_chip_polycomb_enrichment/plots/polycomb_shared/02_enrichment_dotplot_polycomb_shared/02_enrichment_dotplot_polycomb_shared.pdf` |
| **Enrichment data** | ✅ | `output/diff_chip_polycomb_enrichment/tables/diff_chip_enrichment_all_loops.tsv` |
| **Loop-compartment cross-ref** | ✅ | `output/loop_compartment_crossref/` (14 plots: PC1 violin, shift overlap, TAD IR boxplot, etc.) |
| **Script** | ✅ | `scripts/diff_chip_polycomb_enrichment.R`, `scripts/loop_compartment_crossref.R` |

### Panel 3D — Permutation tests: P12 marks predict adult structures (P12 and adult loops against adult loops/boundaries/TADs/compartments)

| Item | Status | Notes |
|------|--------|-------|
| **Cross-timepoint permutation** | ❓ | Needs cross-timepoint comparison: P12 ChIP peaks → adult structural features. Partial data in `tads/results/visualizations/comparison/` |
| **Timepoint comparison plots** | ✅ | `tads/results/visualizations/comparison/boundary_type_comparison/`, `enrichment_direction_comparison/`, etc. |
| **Timepoint comparison data** | ✅ | `tads/results/visualizations/comparison/timepoint_comparison_stats.tsv` |
| **Script** | ✅ | `tads/scripts/timepoint_comparison.R` |

### Panel 3E — H2Az variant dynamics at structural features (conditional on publication)

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | ❓ | Conditional — depends on H2Az data availability by summer |
| **Data** | ❓ | Not yet in repository |

### Panel 3F — One integrated locus, probably Syt1

| Item | Status | Path |
|------|--------|------|
| **Syt1/Nav3 regional overview** | ✅ | `tads/results/visualizations/late/syt1_nav3_focus/syt1_nav3_regional_overview/syt1_nav3_regional_overview.pdf` |
| **Syt1/Nav3 statistics** | ✅ | `tads/results/visualizations/late/syt1_nav3_statistics.tsv` |
| **Contact map** | 🖥️ | Juicebox screenshot at Syt1 locus with loop/boundary overlays |

---

## Figure 4: Enhancer ABC analysis

### Panel 4A — ABC-RNA correlation using unnormalized score

| Item | Status | Path |
|------|--------|------|
| **ΔABC vs log2FC scatter** | ✅ | `abc/results/figures/delta_abc_vs_log2fc.pdf` |
| **Sum ΔABC vs log2FC** | ✅ | `abc/results/figures/sum_delta_abc_vs_log2fc.pdf` |
| **Activity-contact scatter (raw delta, all)** | ✅ | `abc/results/figures/activity_contact_scatter/02_raw_delta_all_pairs/02_raw_delta_all_pairs.pdf` |
| **Activity-contact scatter (log2FC, all)** | ✅ | `abc/results/figures/activity_contact_scatter/04_log2fc_all_pairs/04_log2fc_all_pairs.pdf` |
| **Activity-contact scatter (classified)** | ✅ | `abc/results/figures/activity_contact_scatter/01_raw_delta_classified/01_raw_delta_classified.pdf` |
| **Promoter vs distal** | ✅ | `abc/results/figures/activity_contact_scatter/08_raw_delta_promoter_distal/08_raw_delta_promoter_distal.pdf` |
| **Summary data** | ✅ | `abc/results/figures/activity_contact_scatter/activity_contact_summary.tsv` |
| **Data** | ✅ | `abc/results/delta_abc_all_pairs.tsv` (180K pairs) |
| | ✅ | `abc/results/delta_abc_with_rnaseq.tsv` (113K pairs w/ RNA-seq) |
| **Script** | ✅ | `abc/scripts/step12_activity_contact_scatter.R`, `abc/scripts/step12b_promoter_distal_scatter.R` |

### Panel 4B — Concordance analysis for 957 DEGs (pie chart of concordant vs discordant, 4 categories; characteristics of discordant in supplemental)

| Item | Status | Path |
|------|--------|------|
| **Concordance pie chart (binary)** | ✅ | `abc/results/figures/concordance_pie/concordance_pie_binary/concordance_pie_binary.pdf` |
| **Concordance pie chart (4-cat)** | ✅ | `abc/results/figures/concordance_pie/concordance_pie_4cat/concordance_pie_4cat.pdf` |
| **Concordance pie chart (combined)** | ✅ | `abc/results/figures/concordance_pie/concordance_pie_combined/concordance_pie_combined.pdf` |
| **Concordance bar plot** | ✅ | `abc/results/enhancer_subset_analysis/16_abc_rnaseq_concordance_bar/16_abc_rnaseq_concordance_bar.pdf` |
| **Discordant composite figure** | ✅ | `abc/results/figures/discordant_analysis/01_discordant_composite/01_discordant_composite.pdf` |
| **Enhancer agreement** | ✅ | `abc/results/figures/discordant_analysis/02_enhancer_agreement/02_enhancer_agreement.pdf` |
| **ΔABC magnitude** | ✅ | `abc/results/figures/discordant_analysis/03_dabc_magnitude/03_dabc_magnitude.pdf` |
| **log2FC magnitude** | ✅ | `abc/results/figures/discordant_analysis/04_log2fc_magnitude/04_log2fc_magnitude.pdf` |
| **Class enrichment** | ✅ | `abc/results/figures/discordant_analysis/05_class_enrichment/05_class_enrichment.pdf` |
| **Distance analysis** | ✅ | `abc/results/figures/discordant_analysis/06_distance/06_distance.pdf` |
| **ΔABC vs log2FC scatter** | ✅ | `abc/results/figures/discordant_analysis/08_dabc_vs_log2fc_scatter/08_dabc_vs_log2fc_scatter.pdf` |
| **GO enrichment** | ✅ | `abc/results/figures/discordant_analysis/09_go_bp_compareCluster/09_go_bp_compareCluster.pdf` |
| **KEGG enrichment** | ✅ | `abc/results/figures/discordant_analysis/10_kegg_compareCluster/10_kegg_compareCluster.pdf` |
| **K119ub by concordance** | ✅ | `abc/results/figures/discordant_analysis/11_k119ub_by_concordance_boxplot/11_k119ub_by_concordance_boxplot.pdf` |
| **K119ub significance rate** | ✅ | `abc/results/figures/discordant_analysis/12_k119ub_significance_rate/12_k119ub_significance_rate.pdf` |
| **Data** | ✅ | `abc/results/discordant_gene_characteristics.tsv` |
| | ✅ | `abc/results/enhancer_subset_analysis/abc_rnaseq_concordance_by_class.tsv` |
| | ✅ | `abc/results/figures/discordant_analysis/go_bp_enrichment_results.tsv` |
| | ✅ | `abc/results/figures/discordant_analysis/kegg_enrichment_results.tsv` |
| **Scripts** | ✅ | `abc/scripts/step13_discordant_gene_analysis.R`, `abc/scripts/step13b_go_enrichment.R`, `abc/scripts/step13c_k119ub_concordance.R` |

### Panel 4C — Volcano plots for K27ac, ATAC, K119ub within combined ATAC peaks

| Item | Status | Notes |
|------|--------|-------|
| **K27ac volcano** | 🖥️ | DiffBind results — likely generated on HPC or in separate epigenomics analysis |
| **ATAC volcano** | 🖥️ | DiffBind ATAC differential — separate analysis pipeline |
| **K119ub volcano at enhancers** | ✅ | `abc/results/figures/k119ub_correlation/09_k119ub_volcano_at_enhancers/09_k119ub_volcano_at_enhancers.pdf` |
| **Differential peaks** | ✅ | `peaks/new/H2AK119ub_up.bed`, `H2AK119ub_down.bed` |
| | ✅ | `peaks/atac_seq/ATAC_up.bed`, `ATAC_down.bed` |

### Panel 4D — Model for subsetting enhancers

| Item | Status | Path |
|------|--------|------|
| **Summary patchwork** | ✅ | `abc/results/enhancer_subset_analysis/12_summary_patchwork/12_summary_patchwork.pdf` |
| **Class composition bar** | ✅ | `abc/results/figures/activity_contact_scatter/10_class_composition_bar/10_class_composition_bar.pdf` |
| **Class-level summary** | ✅ | `abc/results/enhancer_subset_analysis/class_level_summary.tsv` |
| **Enhancer class files** | ✅ | `abc/enhancer_subsets/enhancer_classes_activity_lost.tsv` |
| | ✅ | `abc/enhancer_subsets/enhancer_classes_k119ub_only.tsv` |
| | ✅ | `abc/enhancer_subsets/enhancer_classes_activity_gain.tsv` |
| | ✅ | `abc/enhancer_subsets/enhancer_classes_stable.tsv` |
| **Script** | ✅ | `abc/scripts/step11_enhancer_subset_analysis.R` |

### Panel 4E — Loop logFC, delta activity, delta contacts, gene logFC

| Item | Status | Path |
|------|--------|------|
| **Loop logFC violin** | ✅ | `abc/results/enhancer_subset_analysis/05_loop_logfc_violin/05_loop_logfc_violin.pdf` |
| **Loop logFC density** | ✅ | `abc/results/enhancer_subset_analysis/06_loop_logfc_density/06_loop_logfc_density.pdf` |
| **Delta ABC boxplot** | ✅ | `abc/results/enhancer_subset_analysis/07_delta_abc_boxplot/07_delta_abc_boxplot.pdf` |
| **Delta unnorm boxplot** | ✅ | `abc/results/enhancer_subset_analysis/08_delta_unnorm_boxplot/08_delta_unnorm_boxplot.pdf` |
| **Gene logFC by class** | ✅ | `abc/results/enhancer_subset_analysis/04_gene_logfc_by_class/04_gene_logfc_by_class.pdf` |
| **Contact decay WT** | ✅ | `abc/results/enhancer_subset_analysis/10_contact_decay_wt/10_contact_decay_wt.pdf` |
| **Contact decay delta** | ✅ | `abc/results/enhancer_subset_analysis/11_contact_decay_delta/11_contact_decay_delta.pdf` |
| **ABC direction stacked bar** | ✅ | `abc/results/enhancer_subset_analysis/09_abc_direction_stacked_bar/09_abc_direction_stacked_bar.pdf` |
| **ΔABC vs gene logFC scatter** | ✅ | `abc/results/enhancer_subset_analysis/17_delta_abc_vs_gene_logfc_scatter/17_delta_abc_vs_gene_logfc_scatter.pdf` |
| **Paired anchor logFC vs ΔABC** | ✅ | `abc/results/paired_anchor_plots/logFC_vs_deltaABC/logFC_vs_deltaABC.pdf` |
| **Paired anchor logFC vs delta unnorm** | ✅ | `abc/results/paired_anchor_plots/logFC_vs_delta_unnorm/logFC_vs_delta_unnorm.pdf` |
| **Data** | ✅ | `abc/results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv` |
| | ✅ | `abc/results/enhancer_subset_analysis/enhancer_class_loop_metrics.tsv` |
| | ✅ | `abc/results/enhancer_subset_analysis/contact_decay_by_class.tsv` |
| **Script** | ✅ | `abc/scripts/step11_enhancer_subset_analysis.R` |

### Panel 4F — ABC category vs K119ub significance

| Item | Status | Path |
|------|--------|------|
| **K119ub tertile loop logFC** | ✅ | `abc/results/enhancer_subset_analysis/13_k119ub_tertile_loop_logfc/13_k119ub_tertile_loop_logfc.pdf` |
| **K119ub tertile delta ABC** | ✅ | `abc/results/enhancer_subset_analysis/14_k119ub_tertile_delta_abc/14_k119ub_tertile_delta_abc.pdf` |
| **K119ub vs loop logFC scatter** | ✅ | `abc/results/enhancer_subset_analysis/15_k119ub_vs_loop_logfc_scatter/15_k119ub_vs_loop_logfc_scatter.pdf` |
| **K119ub target gene logFC** | ✅ | `abc/results/enhancer_subset_analysis/18_k119ub_target_gene_logfc_hist/18_k119ub_target_gene_logfc_hist.pdf` |
| **Contingency heatmap** | ✅ | `abc/results/figures/k119ub_correlation/10_contingency_heatmap/10_contingency_heatmap.pdf` |
| **K119ub by ABC category boxplot** | ✅ | `abc/results/figures/k119ub_correlation/04_boxplot_k119ub_by_abc_category/04_boxplot_k119ub_by_abc_category.pdf` |
| **K119ub violin by category** | ✅ | `abc/results/figures/k119ub_correlation/06_violin_k119ub_by_category/06_violin_k119ub_by_category.pdf` |
| **Delta activity by K119ub sig** | ✅ | `abc/results/figures/k119ub_correlation/11_delta_activity_by_k119ub_sig/11_delta_activity_by_k119ub_sig.pdf` |
| **Scatter colored by significance** | ✅ | `abc/results/figures/k119ub_correlation/12_scatter_colored_by_significance/12_scatter_colored_by_significance.pdf` |
| **Data** | ✅ | `abc/results/enhancer_subset_analysis/k119ub_tertile_assignments.tsv` |
| | ✅ | `abc/results/enhancer_subset_analysis/k119ub_tertile_loop_summary.tsv` |
| | ✅ | `abc/results/k119ub_abc_correlation_summary.tsv` |
| | ✅ | `abc/results/k119ub_abc_enhancer_merged.tsv` |
| **Scripts** | ✅ | `abc/scripts/step10_k119ub_abc_correlation.R`, `abc/scripts/step11_enhancer_subset_analysis.R` |

---

## Figure 5: Model and functional implications

### Panel 5A — GO analysis: diff boundaries, loops, AxC connections

| Item | Status | Path |
|------|--------|------|
| **GO at differential boundaries (Adult)** | ✅ | `tads/results/visualizations/late/go_bp_dotplot/go_bp_dotplot.pdf` |
| **GO MF at boundaries** | ✅ | `tads/results/visualizations/late/go_mf_dotplot/go_mf_dotplot.pdf` |
| **KEGG at boundaries** | ✅ | `tads/results/visualizations/late/kegg_dotplot/kegg_dotplot.pdf` |
| **Boundary genes** | ✅ | `tads/results/visualizations/late/boundary_genes.tsv` |
| **GO at loops (long vs short)** | ✅ | `output/loops_visualization_extended/late/08_go_comparison_long_vs_short/08_go_comparison_long_vs_short.pdf` |
| **GO long-range lost data** | ✅ | `output/loops_visualization_extended/late/go_long_lost_loops.tsv` |
| **GO short-range gained data** | ✅ | `output/loops_visualization_extended/late/go_short_gained_loops.tsv` |
| **GO at paired ABC anchors** | ✅ | `abc/results/paired_anchor_plots/go_bp_dotplot/go_bp_dotplot.pdf` |
| **KEGG at paired anchors** | ✅ | `abc/results/paired_anchor_plots/kegg_dotplot/kegg_dotplot.pdf` |
| **GO enrichment data** | ✅ | `abc/results/paired_anchor_go_enrichment.tsv` |
| **KEGG enrichment data** | ✅ | `abc/results/paired_anchor_kegg_enrichment.tsv` |
| **Loop enrichment** | ✅ | `outputs/250402-late_outputs/visualizations/enrichment/` (go_bp_dotplot, go_mf_dotplot, kegg_dotplot) |

### Panel 5B — DEG logFC for genes near diff boundaries, at diff loops, with lost AxC connections

| Item | Status | Path |
|------|--------|------|
| **DEG-TAD violin (Adult)** | ✅ | `tads/results/visualizations/late/deg_violin/deg_tad_violin_late.pdf` |
| **DEG-TAD boundary genes** | ✅ | `tads/results/visualizations/late/deg_violin/deg_boundary_genes_late.tsv` |
| **DEG-TAD statistics** | ✅ | `tads/results/visualizations/late/deg_violin/deg_tad_statistics_late.txt` |
| **DEG at loop anchors (Adult)** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_anchor_violin/deg_loop_anchor_violin.pdf` |
| **DEG expression correlation** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_expression_correlation/deg_loop_expression_correlation.pdf` |
| **DEG by distance stratified** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_violin_distance_2x2_stratified/deg_loop_violin_distance_2x2_stratified.pdf` |
| **DEG long-range lost vs short-range gained** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_violin_longrange_lost_vs_shortrange_gained/deg_loop_violin_longrange_lost_vs_shortrange_gained.pdf` |
| **DEG polycomb** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_anchor_violin_polycomb/deg_loop_anchor_violin_polycomb.pdf` |
| **DEG permutation test** | ✅ | `output/deg_loop_violin/late/plots/deg_loop_permutation_test/deg_loop_permutation_test.pdf` |
| **DEG correlation data** | ✅ | `output/deg_loop_violin/late/tables/correlation_summary.tsv` |
| **DEG permutation data** | ✅ | `output/deg_loop_violin/late/tables/permutation_results.tsv` |
| **ABC gene-level summary** | ✅ | `abc/results/gene_level_summary.tsv` (13,588 genes) |
| **Scripts** | ✅ | `scripts/deg_loop_anchor_violin.R`, `tads/scripts/deg_tad_violin.R` |

### Panel 5C — Network analysis (Nodes = genes with combined structural changes; Node size = AxC score change; Node color = Gene expression logFC)

| Item | Status | Notes |
|------|--------|-------|
| **Network figure** | ❓ | Needs to be generated — combine boundary, loop, and ABC structural scores per gene |
| **Input: boundary genes** | ✅ | `tads/results/visualizations/late/boundary_genes.tsv` |
| **Input: loop-gene assignments** | ✅ | `abc/results/loops_with_gene_assignments.tsv` |
| **Input: ABC gene summary** | ✅ | `abc/results/gene_level_summary.tsv` |
| **Input: characterized loops** | ✅ | `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv` |
| **Input: RNA-seq** | ✅ | `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` |

### Panel 5D — Heatmap of top 50 genes by combined structural score (Columns: logFC genes, logFC AxC)

| Item | Status | Notes |
|------|--------|-------|
| **Heatmap** | ❓ | Needs to be generated — rank genes by combined structural disruption score |
| **Input: gene-level ABC** | ✅ | `abc/results/gene_level_summary.tsv` |
| **Input: loop data** | ✅ | `abc/results/paired_anchor_matches.tsv` |
| **Input: boundary data** | ✅ | `tads/results/late/final/tadcompare_final_annotated.tsv` |

### Panel 5E — 3×2 model of epigenetic layer, structural layer, and functional (DEG) layer at P13 and adult

| Item | Status | Notes |
|------|--------|-------|
| **Figure** | 📐 | Conceptual/schematic — Illustrator or BioRender |
| **Supporting quantitative data** | ✅ | All data exists to populate the model numerically |

---

## Supplemental Analyses (Available but not in main figures)

### Shared Anchor / Loop Rewriting Analysis

| Item | Path |
|------|------|
| **Shared anchor loops (Adult)** | `output/shared_anchor_analysis/late/shared_anchor_loops.tsv` |
| **Anchor characterization** | `output/shared_anchor_analysis/late/anchor_characterization.tsv` |
| **Distance violin** | `output/shared_anchor_analysis/late/plots/distance_violin_shared/distance_violin_shared.pdf` |
| **ChIP enrichment** | `output/shared_anchor_analysis/late/plots/chip_enrichment_dotplot/chip_enrichment_dotplot.pdf` |
| **Polycomb-specific** | `output/shared_anchor_analysis/late/polycomb_specific/` (5 plots + 5 tables) |
| **Shared anchor boundaries** | `output/shared_anchor_boundary_analysis/late/` (7 plots + 10 tables) |
| **Loop rewriting summary** | `output/loops_visualization_extended/late/09_loop_rewriting_summary/09_loop_rewriting_summary.pdf` |
| **Scripts** | `scripts/shared_anchor_analysis.R`, `scripts/polycomb_shared_anchor_analysis.R`, `scripts/shared_anchor_boundary_analysis.R` |

### HOMER Motif Enrichment (ABC enhancer subsets)

| Item | Path |
|------|------|
| **Top motifs by comparison** | `abc/results/enhancer_subset_analysis/homer_motif_viz/01-06_top_motifs_*.pdf` |
| **Cross-comparison heatmap** | `abc/results/enhancer_subset_analysis/homer_motif_viz/07_cross_comparison_heatmap/` |
| **Dose-response** | `abc/results/enhancer_subset_analysis/homer_motif_viz/08_dose_comparison_dotplot/` |
| **Narrative composite** | `abc/results/enhancer_subset_analysis/homer_motif_viz/12_narrative_composite/` |
| **Summary stats** | `abc/results/enhancer_subset_analysis/homer_motif_summary_stats.tsv` |
| **Script** | `abc/scripts/step11_homer_motif_enrichment.sb`, `abc/scripts/step11b_homer_motif_visualization.R` |

### Paired Anchor ABC-Loop Concordance

| Item | Path |
|------|------|
| **Concordance panel** | `abc/results/paired_anchor_plots/paired_anchor_panel/paired_anchor_panel.pdf` |
| **FDR-stratified** | `abc/results/paired_anchor_plots/fdr_stratified_concordance/fdr_stratified_concordance.pdf` |
| **RNA-seq concordance** | `abc/results/paired_anchor_plots/rnaseq_concordance/rnaseq_concordance.pdf` |
| **K119ub at paired enhancers** | `abc/results/paired_anchor_plots/k119ub_at_paired_enhancers/k119ub_at_paired_enhancers.pdf` |
| **Distance concordance** | `abc/results/paired_anchor_plots/distance_concordance/distance_concordance.pdf` |
| **Summary** | `abc/results/paired_anchor_summary.txt` |
| **Script** | `abc/scripts/step9b_paired_anchor_analysis.R` |

### Loop-Compartment Cross-Reference

| Item | Path |
|------|------|
| **Shift overlap barplot** | `output/loop_compartment_crossref/q1_shift_overlap_barplot/q1_shift_overlap_barplot.pdf` |
| **PC1 diff violin** | `output/loop_compartment_crossref/q2_pc1_diff_violin/q2_pc1_diff_violin.pdf` |
| **Ctrl PC1 violin** | `output/loop_compartment_crossref/q2_ctrl_pc1_violin/q2_ctrl_pc1_violin.pdf` |
| **TAD IR boxplot** | `output/loop_compartment_crossref/q3_tad_ir_boxplot/q3_tad_ir_boxplot.pdf` |
| **Polycomb PC1** | `output/loop_compartment_crossref/q5_polycomb_pc1_violin/q5_polycomb_pc1_violin.pdf` |
| **Script** | `scripts/loop_compartment_crossref.R` |

### Timepoint Comparison (Early vs Late)

| Item | Path |
|------|------|
| **Boundary type comparison** | `tads/results/visualizations/comparison/boundary_type_comparison/boundary_type_comparison.pdf` |
| **Enrichment direction comparison** | `tads/results/visualizations/comparison/enrichment_direction_comparison/enrichment_direction_comparison.pdf` |
| **Combined summary** | `tads/results/visualizations/comparison/combined_comparison_summary/combined_comparison_summary.pdf` |
| **Net direction diverging** | `tads/results/visualizations/comparison/net_direction_diverging/net_direction_diverging.pdf` |
| **Data** | `tads/results/visualizations/comparison/timepoint_comparison_stats.tsv` |
| **Script** | `tads/scripts/timepoint_comparison.R` |

### Methylation / Biomodal Analysis

| Item | Path |
|------|------|
| **EM-seq analysis pipeline** | `biomodal/downstream/scripts/viz_sections/` (30 analysis sections) |
| **Plots** | `biomodal/downstream/plots/` |
| **QC** | `biomodal/downstream/qc/` |

### Stripes Analysis

| Item | Path |
|------|------|
| **Stripe outputs** | `stripes/outputs/` |
| **Stripe scripts** | `stripes/scripts/` |

---

## Summary: Items Needing Generation

| Panel | What's Missing | Data Available? |
|-------|---------------|-----------------|
| **1A** Contact maps at P13 vs adult | Juicebox screenshots with insulation tracks | 🖥️ .hic files on HPC |
| **1B** P13 TAD volcano | Early HOMER-style TAD volcano (TADCompare data exists) | ✅ `tads/results/early/tadcompare/` |
| **1D** % genome calculation | ✅ Done — `compartment_genome_pct_bar/`, summary stats | ✅ `compartment_significant_*.tsv` |
| **1E** Cntnap5a locus | Juicebox triangular heatmap | 🖥️ .hic files |
| **2D** Example loop loci | Juicebox screenshots with BEDPE overlays | ✅ BEDPE files exist |
| **3A** ARA heatmaps | deepTools heatmaps at anchors for 4 marks | 🖥️ bigwig files on HPC |
| **3D** Cross-timepoint permutation | P12 peaks → adult structure permutation test | ✅ All input data exists |
| **3E** H2Az dynamics | Conditional on publication | ❓ No data yet |
| **4C** K27ac/ATAC volcanoes in ATAC peaks | DiffBind volcano within combined ATAC regions | 🖥️ DiffBind results |
| **5C** Network analysis | Integrate boundary + loop + ABC per gene → network | ✅ All component data exists |
| **5D** Top 50 genes heatmap | Rank genes by combined structural score → heatmap | ✅ All component data exists |
| **5E** 3×2 conceptual model | Schematic figure | 📐 Illustrator/BioRender |

---

## Key Data Tables Reference

| File | Location | Rows | Description |
|------|----------|------|-------------|
| `all_results_primary.tsv` | `outputs/250402-late_outputs/edgeR_results_res_10kb/primary_analysis/` | ~39K | All tested loops (adult) |
| `characterized_loops.tsv` | `outputs/250402-late_outputs/merged_loops/` | ~2.9K | Annotated differential loops (adult) |
| `extended_characterized_loops.tsv` | `peaks/loop_annotation_extended/late/` | ~2.9K | 7-category annotation (adult) |
| `delta_abc_all_pairs.tsv` | `abc/results/` | 180K | All enhancer-gene pairs with ΔABC |
| `delta_abc_with_rnaseq.tsv` | `abc/results/` | 113K | ΔABC merged with RNA-seq |
| `gene_level_summary.tsv` | `abc/results/` | 13.6K | Gene-level strongest enhancer + DE |
| `tadcompare_final_annotated.tsv` | `tads/results/late/final/` | — | Final TAD boundaries with annotations |
| `boundaries_with_chromatin_state.tsv` | `tads/results/visualizations/chip/late/` | — | TAD boundaries + ChIP-seq states |
| `k119ub_anchor_signal.tsv` | `data/` | — | K119ub signal at loop anchors |
| `k119ub_abc_enhancer_merged.tsv` | `abc/results/` | ~9K | Enhancer K119ub + ABC scores |
| `anchor_type_summary.tsv` | `peaks/loop_annotation_extended/late/` | 7 | Anchor category counts |
| `loop_type_summary.tsv` | `peaks/loop_annotation_extended/late/` | ~28 | Loop type pair counts |

---

## Key Scripts Reference

| Script | Language | Purpose |
|--------|----------|---------|
| `scripts/edgeR.R` | R | Differential loop analysis (quasi-likelihood GLM) |
| `scripts/apa_analysis.R` | R | Aggregate Peak Analysis heatmaps |
| `scripts/loop_distance_analysis.R` | R | Loop distance eCDF/density by direction |
| `scripts/loop_distance_mark_filtered.R` | R | Mark-specific loop distance analysis |
| `scripts/annotate_loops_extended.R` | R | 7-category anchor classification |
| `scripts/h2ak119ub_loop_integration.R` | R | H2AK119ub integration with loops |
| `scripts/deg_loop_anchor_violin.R` | R | DEG expression at loop anchors |
| `scripts/shared_anchor_analysis.R` | R | Loop rewriting / shared anchor analysis |
| `scripts/diff_chip_polycomb_enrichment.R` | R | Differential ChIP enrichment at loops |
| `scripts/loop_compartment_crossref.R` | R | Loop × compartment cross-reference |
| `scripts/ctcf_stripe_crossref.R` | R | CTCF stripe × loop cross-reference |
| `scripts/tad_volcano_plot.R` | R | TAD volcano plots |
| `scripts/compartment_volcano_plot.R` | R | Compartment volcano plots |
| `tads/scripts/tad_chip_classification.R` | R | ChIP-seq at TAD boundaries |
| `tads/scripts/deg_tad_violin.R` | R | DEGs at TAD boundaries |
| `tads/scripts/timepoint_comparison.R` | R | Early vs late comparison |
| `abc/scripts/step6_delta_abc.py` | Python | Compute ΔABC between conditions |
| `abc/scripts/step11_enhancer_subset_analysis.R` | R | Stratified enhancer class analysis |
| `abc/scripts/step10_k119ub_abc_correlation.R` | R | K119ub-ABC correlation |
| `abc/scripts/step12_activity_contact_scatter.R` | R | Decompose ΔABC into activity vs contact |
| `abc/scripts/step13_discordant_gene_analysis.R` | R | Concordance/discordance analysis |

---

*Generated: 2026-03-12 | All paths relative to `/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/`*
