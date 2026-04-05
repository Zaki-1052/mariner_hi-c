#!/bin/bash
# biomodal/downstream/scripts/copy_presentation_plots.sh
# Copies all plots used in the presentation into a numbered presentation/ folder
# Plots are numbered to match presentation slide order

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(dirname "$SCRIPT_DIR")"
VIZ="$BASE/plots/visualizations"
COMP="$VIZ/comparison_shallow_vs_deep"
DMR_PLOTS="$BASE/modality/DMR_processed/summary/plots"
DEST="$BASE/presentation"

rm -rf "$DEST"
mkdir -p "$DEST"

N=0
copy() {
  N=$((N + 1))
  local src="$1"
  local label="$2"
  local num=$(printf "%03d" $N)
  local ext="${src##*.}"
  if [[ -f "$src" ]]; then
    cp "$src" "$DEST/${num}_${label}.${ext}"
  else
    echo "WARNING: missing $src"
  fi
}

# =============================================================================
# COMPARISON: Run-4 vs Run-5 (slides 1-5)
# =============================================================================
copy "$COMP/20a_dmr_counts_by_region/20a_dmr_counts_by_region.jpg"         "dmr_counts_comparison"
copy "$COMP/20b_gene_status_summary/20b_gene_status_summary.jpg"           "gene_status_comparison"
copy "$COMP/20d_effect_size_scatter/20d_effect_size_scatter.jpg"            "effect_size_stability"
copy "$COMP/20e_qvalue_improvement/20e_qvalue_improvement.jpg"             "qvalue_improvement"
copy "$COMP/20g_coordinated_pattern_stability/20g_coordinated_pattern_stability.jpg" "coordinated_stability"

# =============================================================================
# QC & OVERVIEW (slides 6-7)
# =============================================================================
copy "$VIZ/01_qc_overview/01_qc_overview.jpg"                              "qc_overview"
copy "$VIZ/42_max_significance/42_pca_all_genes/42_pca_all_genes.jpg"      "pca_all_genes"

# =============================================================================
# DMR STATISTICS (slides 8-10)
# =============================================================================
copy "$VIZ/03_dmr_region_statistics/03_dmr_region_statistics.jpg"           "dmr_region_statistics"
copy "$VIZ/03a_dmr_by_region/03a_dmr_by_region.jpg"                        "dmr_by_region"
copy "$VIZ/03b_direction_comparison/03b_direction_comparison.jpg"           "direction_comparison"
copy "$DMR_PLOTS/01_summary_significant_dmrs.png"                          "summary_significant_dmrs_by_context"
copy "$DMR_PLOTS/05_heatmap_significant_dmrs.png"                          "heatmap_significant_dmrs"

# =============================================================================
# EFFECT SIZE (slide 11)
# =============================================================================
copy "$VIZ/07_effect_size_distributions/07_effect_size_distributions.jpg"   "effect_size_distributions"

# =============================================================================
# KEY GENES & COORDINATED CHANGES (slides 12-14)
# =============================================================================
copy "$VIZ/05a_mc_hmc_scatter/05a_mc_hmc_scatter.jpg"                      "mc_hmc_scatter"
copy "$VIZ/28b_mc_hmc_concordance_scatter/28b_mc_hmc_concordance_scatter.jpg" "mc_hmc_concordance_scatter"
copy "$VIZ/05_coordinated_changes/05_coordinated_changes.jpg"              "coordinated_changes"

# =============================================================================
# VOLCANO PLOTS (slide 15)
# =============================================================================
copy "$VIZ/04_volcano_plots/04_volcano_plots.jpg"                           "volcano_plots"

# =============================================================================
# TOP GENES & VENN (slides 16-17)
# =============================================================================
copy "$VIZ/05b_top_coordinated_genes/05b_top_coordinated_genes.jpg"        "top_coordinated_genes"
copy "$VIZ/06_top_genes/06_top_genes.jpg"                                  "top_genes"
copy "$VIZ/06a_top_dmrs/06a_top_dmrs.jpg"                                  "top_dmrs"
copy "$VIZ/06b_venn_overlap/06b_venn_overlap.jpg"                          "venn_overlap"

# =============================================================================
# CHROMATIN STATE (slides 18-20)
# =============================================================================
copy "$VIZ/10a_chromatin_state_distribution/10a_chromatin_state_distribution.jpg"   "chromatin_state_distribution"
copy "$VIZ/10d_coordinated_genes_chromatin/10d_coordinated_genes_chromatin.jpg"     "coordinated_genes_chromatin"
copy "$VIZ/10c_chip_mark_overlap_heatmap/10c_chip_mark_overlap_heatmap.jpg"        "chip_mark_overlap_heatmap"
copy "$VIZ/10b_chromatin_by_methylation_direction/10b_chromatin_by_methylation_direction.jpg" "chromatin_by_direction"
copy "$VIZ/10e_top_genes_chromatin_annotation/10e_top_genes_chromatin_annotation.jpg" "top_genes_chromatin"

# =============================================================================
# MeCP2 INTEGRATION (slides 21-22)
# =============================================================================
copy "$VIZ/11a_mecp2_overlap/11a_mecp2_overlap.jpg"                        "mecp2_overlap"
copy "$VIZ/11h_mecp2_binary_gain_fisher/11h_mecp2_binary_gain_fisher.jpg"  "mecp2_fisher"
copy "$VIZ/11c_mc_vs_mecp2_scatter/11c_mc_vs_mecp2_scatter.jpg"           "mc_vs_mecp2_scatter"
copy "$VIZ/11e_mecp2_integration_heatmap/11e_mecp2_integration_heatmap.jpg" "mecp2_integration_heatmap"
copy "$VIZ/11b_mecp2_fold_by_dmr_direction/11b_mecp2_fold_by_dmr_direction.jpg" "mecp2_fold_by_direction"

# =============================================================================
# ATAC-seq (slides 23-26)
# =============================================================================
copy "$VIZ/12a_atac_overlap/12a_atac_overlap.jpg"                          "atac_overlap"
copy "$VIZ/12e_atac_integration_heatmap/12e_atac_integration_heatmap.jpg"  "atac_integration_heatmap"
copy "$VIZ/12c_atac_coordinated_genes/12c_atac_coordinated_genes.jpg"      "atac_coordinated_genes"
copy "$VIZ/12d_mc_vs_atac_scatter/12d_mc_vs_atac_scatter.jpg"             "mc_vs_atac_scatter"
copy "$VIZ/13c_atac_chromatin_enrichment_heatmap/13c_atac_chromatin_enrichment_heatmap.jpg" "atac_chromatin_enrichment"

# =============================================================================
# H2AK119ub (slides 27-28)
# =============================================================================
copy "$VIZ/14b_k119ub_differential_overlap/14b_k119ub_differential_overlap.jpg"   "k119ub_differential_overlap"
copy "$VIZ/14e_k119ub_integration_heatmap/14e_k119ub_integration_heatmap.jpg"     "k119ub_integration_heatmap"

# =============================================================================
# 5hmC INTEGRATION HEATMAPS (slide 29)
# =============================================================================
copy "$VIZ/15a_hmc_mecp2_heatmap/15a_hmc_mecp2_heatmap.jpg"               "hmc_mecp2_heatmap"
copy "$VIZ/15c_hmc_k119ub_heatmap/15c_hmc_k119ub_heatmap.jpg"            "hmc_k119ub_heatmap"
copy "$VIZ/15b_hmc_atac_heatmap/15b_hmc_atac_heatmap.jpg"                 "hmc_atac_heatmap"
copy "$VIZ/15d_mc_vs_hmc_enrichment_comparison/15d_mc_vs_hmc_enrichment_comparison.jpg" "mc_vs_hmc_enrichment"

# =============================================================================
# K119ub COORDINATED & SCATTER (slides 30-31)
# =============================================================================
copy "$VIZ/14c_k119ub_coordinated_genes/14c_k119ub_coordinated_genes.jpg"  "k119ub_coordinated_genes"
copy "$VIZ/14d_mc_vs_k119ub_scatter/14d_mc_vs_k119ub_scatter.jpg"         "mc_vs_k119ub_scatter"

# =============================================================================
# H3K27ac (slides 32-33)
# =============================================================================
copy "$VIZ/19f_h3k27ac_condition_overlap/19f_h3k27ac_condition_overlap.jpg" "h3k27ac_condition_overlap"
copy "$VIZ/19g_h3k27ac_oe_heatmaps/19g_h3k27ac_oe_heatmaps.jpg"          "h3k27ac_oe_heatmaps"
copy "$VIZ/19d_methylation_vs_h3k27ac_scatter/19d_methylation_vs_h3k27ac_scatter.jpg" "h3k27ac_scatter"
copy "$VIZ/19h_chromatin_mark_oe_comparison/19h_chromatin_mark_oe_comparison.jpg"     "chromatin_mark_oe_comparison"

# =============================================================================
# CONCORDANCE (slides 34-35)
# =============================================================================
copy "$VIZ/16d_raw_concordance_comparison/16d_raw_concordance_comparison.jpg" "raw_concordance_comparison"
copy "$VIZ/16a_mecp2_concordance_bars/16a_mecp2_concordance_bars.jpg"      "mecp2_concordance"
copy "$VIZ/16b_atac_concordance_bars/16b_atac_concordance_bars.jpg"        "atac_concordance"
copy "$VIZ/16c_k119ub_concordance_bars/16c_k119ub_concordance_bars.jpg"    "k119ub_concordance"

# =============================================================================
# K119ub DETAILED (slides 36-38)
# =============================================================================
copy "$VIZ/18c_methylation_vs_k119ub_scatter/18c_methylation_vs_k119ub_scatter.jpg" "methylation_vs_k119ub"
copy "$VIZ/18b_k119ub_log2fc_distributions/18b_k119ub_log2fc_distributions.jpg"     "k119ub_log2fc"
copy "$VIZ/17a_k119ub_full_breakdown/17a_k119ub_full_breakdown.jpg"        "k119ub_full_breakdown"
copy "$VIZ/17b_k119ub_conditional_direction/17b_k119ub_conditional_direction.jpg" "k119ub_conditional"
copy "$VIZ/17c_k119ub_effect_sizes/17c_k119ub_effect_sizes.jpg"           "k119ub_effect_sizes"

# =============================================================================
# GO/KEGG ENRICHMENT (slides 39-40)
# =============================================================================
copy "$VIZ/08a_enrichment_go_bp/08a_enrichment_go_bp.jpg"                  "go_bp"
copy "$VIZ/08b_enrichment_go_cc/08b_enrichment_go_cc.jpg"                  "go_cc"
copy "$VIZ/08c_enrichment_go_mf/08c_enrichment_go_mf.jpg"                  "go_mf"
copy "$VIZ/08d_enrichment_kegg/08d_enrichment_kegg.jpg"                    "kegg"
copy "$VIZ/08e_enrichment_delta_ratio_compare_go_bp/08e_enrichment_delta_ratio_compare_go_bp.jpg" "delta_ratio_go_bp"
copy "$VIZ/08f_enrichment_delta_ratio_compare_kegg/08f_enrichment_delta_ratio_compare_kegg.jpg"   "delta_ratio_kegg"

# =============================================================================
# EXPRESSION (slides 41-43)
# =============================================================================
copy "$VIZ/20d_mc_expression_heatmap/20d_mc_expression_heatmap.jpg"        "mc_expression_heatmap"
copy "$VIZ/20a_coordinated_expression_breakdown/20a_coordinated_expression_breakdown.jpg" "expression_breakdown"
copy "$VIZ/20c_log2fc_violin_comparison/20c_log2fc_violin_comparison.jpg"  "log2fc_violin"
copy "$VIZ/20b_methylation_vs_expression_scatter/20b_methylation_vs_expression_scatter.jpg" "methylation_vs_expression"

# =============================================================================
# DISCORDANT GENES (slides 44-45)
# =============================================================================
copy "$VIZ/21a_discordant_composite/21a_discordant_composite.jpg"          "discordant_composite"
copy "$VIZ/21b_mc_hmc_concordance_scatter/21b_mc_hmc_concordance_scatter.jpg" "discordant_concordance_scatter"

# =============================================================================
# CONCORDANT / ALL QUADRANTS (slides 46-48)
# =============================================================================
copy "$VIZ/21e_all_quadrants_comprehensive/21e_all_quadrants_comprehensive.jpg" "all_quadrants_comprehensive"
copy "$VIZ/28a_coordinated_composite/28a_coordinated_composite.jpg"        "coordinated_composite"
copy "$VIZ/28e_all_quadrants_comprehensive/28e_all_quadrants_comprehensive.jpg" "coordinated_all_quadrants"

# =============================================================================
# LOOP INTEGRATION (slides 49-52)
# =============================================================================
copy "$VIZ/27_methylation_loop_integration/27b_mc_direction_by_loop_direction/27b_mc_direction_by_loop_direction.jpg" "hypermeth_rate_by_loop"
copy "$VIZ/27_methylation_loop_integration/27b_delta_ratio_violin_by_loop_direction/27b_delta_ratio_violin_by_loop_direction.jpg" "delta_ratio_by_loop"
copy "$VIZ/27_methylation_loop_integration/27b_mc_diff_violin_by_loop_direction/27b_mc_diff_violin_by_loop_direction.jpg" "mc_diff_by_loop"
copy "$VIZ/27_methylation_loop_integration/27c_k119ub_methylation_loop_convergence_heatmap/27c_k119ub_methylation_loop_convergence_heatmap.jpg" "k119ub_loop_convergence"
copy "$VIZ/27_methylation_loop_integration/27d_logistic_regression_forest_plot/27d_logistic_regression_forest_plot.jpg" "loop_logistic_regression"

# =============================================================================
# DEMETHYLATION RATIO (slides 53-55)
# =============================================================================
copy "$VIZ/22c_delta_ratio_histogram/22c_delta_ratio_histogram.jpg"        "delta_ratio_histogram"
copy "$VIZ/22f_delta_ratio_by_chromatin_state/22f_delta_ratio_by_chromatin_state.jpg" "delta_ratio_by_chromatin"
copy "$VIZ/23d_dose_response_scatter/23d_dose_response_scatter.jpg"        "dose_response_scatter"

# =============================================================================
# PREDICTIVE MODELS (slides 56-60)
# =============================================================================
copy "$VIZ/23a_roc_curves/23a_roc_curves.jpg"                              "roc_curves_s23"
copy "$VIZ/25b_binary_vs_continuous_s23/25b_binary_vs_continuous_s23.jpg"   "binary_vs_continuous_s23"
copy "$VIZ/25a_dose_response_delta_ratio/25a_dose_response_delta_ratio.jpg" "dose_response_delta_ratio"
copy "$VIZ/25f_predicted_vs_observed/25f_predicted_vs_observed.jpg"         "predicted_vs_observed"
copy "$VIZ/24b_roc_curves/24b_roc_curves.jpg"                              "roc_curves_s24"
copy "$VIZ/24c_feature_importance/24c_feature_importance.jpg"               "feature_importance"
copy "$VIZ/24i_interaction_k119ub_hmc/24i_interaction_k119ub_hmc.jpg"      "interaction_k119ub_hmc"

# =============================================================================
# TET-KO COMPARISON (slides 61-63)
# =============================================================================
copy "$VIZ/26i_response_decomposition/26i_response_decomposition.jpg"      "response_decomposition"
copy "$VIZ/26c_qq_plot/26c_qq_plot.jpg"                                    "qq_plot_tet"
copy "$VIZ/26b_delta_ratio_density/26b_delta_ratio_density.jpg"            "delta_ratio_density_tet"

# =============================================================================
# A/B COMPARTMENT (slides 64-66)
# =============================================================================
copy "$VIZ/29_ab_compartment_methylation/29g_composite_compartment_summary/29g_composite_compartment_summary.jpg" "compartment_composite"
copy "$VIZ/29_ab_compartment_methylation/29f_pc1_vs_mc_scatter/29f_pc1_vs_mc_scatter.jpg"                 "pc1_vs_mc_scatter"
copy "$VIZ/29_ab_compartment_methylation/29c_mc_violin_by_shift/29c_mc_violin_by_shift.jpg"               "mc_by_compartment_shift"

# =============================================================================
# POLYCOMB (slides 67-68)
# =============================================================================
copy "$VIZ/30_polycomb_enrichment/30a_polycomb_vs_non_polycomb_stacked_bar/30a_polycomb_vs_non_polycomb_stacked_bar.jpg" "polycomb_stacked_bar"
copy "$VIZ/30_polycomb_enrichment/30d_per_state_hypermethylation_rate/30d_per_state_hypermethylation_rate.jpg" "hypermeth_rate_by_state"

# =============================================================================
# MeCP2 AT LOOP ANCHORS (slide 69)
# =============================================================================
copy "$VIZ/31_mecp2_loop_anchor_integration/31a_mecp2_peak_overlap_at_loop_anchors/31a_mecp2_peak_overlap_at_loop_anchors.jpg" "mecp2_at_loop_anchors"

# =============================================================================
# DIFFBIND 4-MARK (slides 70-72)
# =============================================================================
copy "$VIZ/33_multi_mark_diffbind/33a_diffbind_volcano_plots/33a_diffbind_volcano_plots.jpg"       "diffbind_volcanos"
copy "$VIZ/33_multi_mark_diffbind/33b_cross_mark_correlation_heatmap/33b_cross_mark_correlation_heatmap.jpg" "cross_mark_correlation"
copy "$VIZ/33_multi_mark_diffbind/33c_quantitative_oe_dotplot/33c_quantitative_oe_dotplot.jpg"     "quantitative_oe_dotplot"
copy "$VIZ/33_multi_mark_diffbind/33d_methylation_vs_mark_scatters/33d_methylation_vs_mark_scatters.jpg" "methylation_vs_marks"
copy "$VIZ/33_multi_mark_diffbind/33e_logistic_regression_forest/33e_logistic_regression_forest.jpg" "diffbind_logistic_forest"

# =============================================================================
# H3K36me3 (slides 73-76)
# =============================================================================
copy "$VIZ/38_h3k36me3_gene_body/38a_h3k36me3_volcano_annotation/38a_h3k36me3_volcano_annotation.jpg" "h3k36me3_volcano"
copy "$VIZ/38_h3k36me3_gene_body/38b_me3_x_mc_oe_heatmap/38b_me3_x_mc_oe_heatmap.jpg"                "me3_x_mc_oe"
copy "$VIZ/38_h3k36me3_gene_body/38c_me3_x_hmc_oe_heatmap/38c_me3_x_hmc_oe_heatmap.jpg"              "me3_x_hmc_oe"
copy "$VIZ/38_h3k36me3_gene_body/38d_me3_vs_methylation_scatter/38d_me3_vs_methylation_scatter.jpg"    "me3_vs_methylation"
copy "$VIZ/38_h3k36me3_gene_body/38e_me3_fold_coordinated_violin/38e_me3_fold_coordinated_violin.jpg"  "me3_fold_coordinated"
copy "$VIZ/38_h3k36me3_gene_body/38f_dmr_me3_bed_overlap/38f_dmr_me3_bed_overlap.jpg"                  "dmr_me3_overlap"
copy "$VIZ/38_h3k36me3_gene_body/38h_me3_fold_by_chromatin_state/38h_me3_fold_by_chromatin_state.jpg"  "me3_by_chromatin_state"

# =============================================================================
# H3K36me2 (slides 77-79)
# =============================================================================
copy "$VIZ/39_h3k36me2_boundary/39a_h3k36me2_volcano_annotation/39a_h3k36me2_volcano_annotation.jpg" "h3k36me2_volcano"
copy "$VIZ/39_h3k36me2_boundary/39b_me2_x_mc_oe_heatmap/39b_me2_x_mc_oe_heatmap.jpg"                "me2_x_mc_oe"
copy "$VIZ/39_h3k36me2_boundary/39c_me2_at_k27me3_boundary/39c_me2_at_k27me3_boundary.jpg"          "me2_k27me3_boundary"
copy "$VIZ/39_h3k36me2_boundary/39d_me2_x_k27me3_oe_heatmap/39d_me2_x_k27me3_oe_heatmap.jpg"        "me2_x_k27me3_oe"
copy "$VIZ/39_h3k36me2_boundary/39e_me2_vs_methylation_scatter/39e_me2_vs_methylation_scatter.jpg"    "me2_vs_methylation"
copy "$VIZ/39_h3k36me2_boundary/39f_me2_by_chromatin_state/39f_me2_by_chromatin_state.jpg"            "me2_by_chromatin_state"

# =============================================================================
# H3K36me2/me3 COMBINED (slides 80-82)
# =============================================================================
copy "$VIZ/40_h3k36me2_me3_combined/40a_expanded_correlation_heatmap/40a_expanded_correlation_heatmap.jpg" "expanded_correlation_heatmap"
copy "$VIZ/40_h3k36me2_me3_combined/40b_me2_vs_me3_scatter/40b_me2_vs_me3_scatter.jpg"                    "me2_vs_me3_scatter"
copy "$VIZ/40_h3k36me2_me3_combined/40c_me2_me3_ratio_delta/40c_me2_me3_ratio_delta.jpg"                  "me2_me3_ratio_delta"
copy "$VIZ/40_h3k36me2_me3_combined/40d_three_way_venn/40d_three_way_venn.jpg"                              "me2_me3_mc_venn"
copy "$VIZ/40_h3k36me2_me3_combined/40e_direction_flow/40e_direction_flow.jpg"                              "me2_me3_direction_flow"
copy "$VIZ/40_h3k36me2_me3_combined/40f_go_comparison/40f_go_comparison.jpg"                                "me2_me3_go_comparison"

# =============================================================================
# DNMT3A vs DNMT3B (slides 83-86)
# =============================================================================
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41a_me3_at_hyper_k119ub_stratified/41a_me3_at_hyper_k119ub_stratified.jpg" "me3_hyper_k119ub_stratified"
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41b_me2_me3_independent_pathway/41b_me2_me3_independent_pathway.jpg"       "me2_me3_independent"
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41c_logistic_regression_forest/41c_logistic_regression_forest.jpg"          "dnmt_logistic_forest"
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41d_pathway_hypermethylation_rate/41d_pathway_hypermethylation_rate.jpg"    "pathway_hypermeth_rate"
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41e_decision_matrix_heatmap/41e_decision_matrix_heatmap.jpg"                "decision_matrix"
copy "$VIZ/41_dnmt3a_vs_dnmt3b/41h_top50_summary_heatmap/41h_top50_summary_heatmap.jpg"                    "top50_multimark_heatmap"

# =============================================================================
# CHROMOSOME ANALYSIS (slides 87-92)
# =============================================================================
copy "$VIZ/43_cg_exploratory/43a_cg_direction_breakdown/43a_cg_direction_breakdown.jpg"        "cg_direction_breakdown"
copy "$VIZ/43_cg_exploratory/43b_cg_mc_hmc_venn/43b_cg_mc_hmc_venn.jpg"                        "cg_mc_hmc_venn"
copy "$VIZ/43_cg_exploratory/43c_cg_mc_chromosome_distribution/43c_cg_mc_chromosome_distribution.jpg" "mc_chromosome_distribution"
copy "$VIZ/43_cg_exploratory/43d_cg_hmc_chromosome_distribution/43d_cg_hmc_chromosome_distribution.jpg" "hmc_chromosome_distribution"
copy "$VIZ/43_cg_exploratory/43e_cg_mc_vs_hmc_chr_comparison/43e_cg_mc_vs_hmc_chr_comparison.jpg" "mc_vs_hmc_chr_comparison"
copy "$VIZ/43_cg_exploratory/43g_cg_mc_chr_direction/43g_cg_mc_chr_direction.jpg"              "mc_chr_direction"
copy "$VIZ/43_cg_exploratory/43h_cg_hmc_chr_direction/43h_cg_hmc_chr_direction.jpg"            "hmc_chr_direction"
copy "$VIZ/43_cg_exploratory/43k_cg_mc_sig_rate_vs_gene_density/43k_cg_mc_sig_rate_vs_gene_density.jpg" "mc_sig_rate_vs_density"
copy "$VIZ/43_cg_exploratory/43l_cg_mc_effect_size_by_chr/43l_cg_mc_effect_size_by_chr.jpg"    "mc_effect_by_chr"
copy "$VIZ/43_cg_exploratory/43m_chrX_direction_comparison/43m_chrX_direction_comparison.jpg"   "chrX_direction"

# =============================================================================
# ALLELE-SPECIFIC METHYLATION (slides 93-95)
# =============================================================================
copy "$VIZ/44_allele_specific_methylation/44c_asm_mc_chromosome_distribution/44c_asm_mc_chromosome_distribution.jpg" "asm_chr_distribution"
copy "$VIZ/44_allele_specific_methylation/44e_mutant_vs_control_asm_comparison/44e_mutant_vs_control_asm_comparison.jpg" "asm_mut_vs_ctrl"
copy "$VIZ/44_allele_specific_methylation/44h_shared_asm_sites_venn/44h_shared_asm_sites_venn.jpg" "asm_venn"
copy "$VIZ/44_allele_specific_methylation/44i_asm_dmr_overlap/44i_asm_dmr_overlap.jpg"         "asm_dmr_overlap"
copy "$VIZ/44_allele_specific_methylation/44j_mc_vs_hmc_asm_comparison/44j_mc_vs_hmc_asm_comparison.jpg" "asm_mc_vs_hmc"

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "Copied $N plots to $DEST/"
echo "Files are numbered 001-$(printf '%03d' $N) in presentation order."
