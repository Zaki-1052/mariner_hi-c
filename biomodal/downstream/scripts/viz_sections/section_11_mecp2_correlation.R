# biomodal/downstream/scripts/viz_sections/section_11_mecp2_correlation.R
# Section 11: MeCP2 Correlation at DMRs
# Standalone script - sources shared config for all dependencies and data
# Tests whether DMR methylation changes correlate with MeCP2 binding changes.
#
# Biological prediction (Mellen et al. 2017):
#   MeCP2 binds 5mCG with HIGH affinity, 5hmCG with LOW affinity.
#   BAP1 loss -> H2AK119ub accumulation -> TET blocked -> 5mC up / 5hmC down
#   Therefore: hypermethylated DMRs should show INCREASED MeCP2 occupancy.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_11_mecp2_correlation.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# MeCP2-SPECIFIC HELPER FUNCTIONS
# =============================================================================

#' Load MeCP2 annotated differential binding data
#' @param filepath Path to MeCP2_annotated.txt
#' @return data.frame with MeCP2 peak-level data
load_mecp2_annotated <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("MeCP2 annotated file not found: ", filepath)
    return(NULL)
  }

  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "")

  # Ensure numeric columns are properly typed
  numeric_cols <- c("start", "end", "width", "Conc", "Conc_mut", "Conc_ctrl",
                    "Fold", "p.value", "FDR", "geneStart", "geneEnd",
                    "geneLength", "distanceToTSS")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  cat(sprintf("  Loaded %d MeCP2 peaks\n", nrow(df)))
  cat(sprintf("  Significant (FDR < 0.05): %d up, %d down\n",
              sum(df$FDR < 0.05 & df$Fold > 0, na.rm = TRUE),
              sum(df$FDR < 0.05 & df$Fold < 0, na.rm = TRUE)))

  return(df)
}

#' Aggregate MeCP2 peaks per gene
#' @param mecp2_df MeCP2 annotated data.frame
#' @return data.frame with gene-level MeCP2 summary
aggregate_mecp2_by_gene <- function(mecp2_df) {
  gene_summary <- mecp2_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      mean_fold = mean(Fold, na.rm = TRUE),
      max_fold = max(Fold, na.rm = TRUE),
      min_fold = min(Fold, na.rm = TRUE),
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr = FDR[which.min(abs(distanceToTSS))],
      min_distance_tss = min(abs(distanceToTSS), na.rm = TRUE),
      n_peaks = n(),
      any_sig = any(FDR < 0.05),
      .groups = "drop"
    )

  return(gene_summary)
}

# =============================================================================
# Re-compute coordinated changes (from Section 5) for independence
# =============================================================================

mc_sig_coord <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_ctrl = mean_mod_group1, mc_mut = mean_mod_group2,
                chr, start, end)

hmc_sig_coord <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2)

coordinated <- inner_join(mc_sig_coord, hmc_sig_coord, by = "gene")

coordinated <- coordinated %>%
  mutate(
    coordinated_pattern = (mc_diff > 0 & hmc_diff < 0),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  ) %>%
  arrange(desc(combined_effect))

# =============================================================================
# SECTION 11: MeCP2 CORRELATION AT DMRs
# =============================================================================

cat("================================================================================\n")
cat("SECTION 11: MeCP2 CORRELATION AT DMRs\n")
cat("================================================================================\n\n")

# Load MeCP2 data
cat("Loading MeCP2 differential binding data...\n")
mecp2_annotated <- load_mecp2_annotated(MECP2_FILES$annotated)
mecp2_up_gr <- load_chip_peaks(MECP2_FILES$up, "MeCP2 Up")
mecp2_down_gr <- load_chip_peaks(MECP2_FILES$down, "MeCP2 Down")

if (is.null(mecp2_annotated) || is.null(mecp2_up_gr) || is.null(mecp2_down_gr)) {
  cat("  ERROR: MeCP2 data files not found. Skipping Section 11.\n")
} else {

  # Aggregate MeCP2 data to gene level
  cat("\nAggregating MeCP2 peaks per gene...\n")
  mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated)
  cat(sprintf("  %d unique genes with MeCP2 peaks\n", nrow(mecp2_gene)))

  # -----------------------------------------------------------------------
  # FIGURE 11a: MeCP2 Peak Overlap at DMRs
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 11a: MeCP2 peak overlap at DMRs...\n")

  # Get significant mC DMRs split by direction
  mc_hyper <- mc_dmr %>% dplyr::filter(significant & mod_difference > 0)
  mc_hypo <- mc_dmr %>% dplyr::filter(significant & mod_difference < 0)

  # Convert to GRanges
  hyper_gr <- dmr_to_granges(mc_hyper)
  hypo_gr <- dmr_to_granges(mc_hypo)

  # Compute overlaps with MeCP2 peaks
  hyper_up <- sum(countOverlaps(hyper_gr, mecp2_up_gr) > 0)
  hyper_down <- sum(countOverlaps(hyper_gr, mecp2_down_gr) > 0)
  hypo_up <- sum(countOverlaps(hypo_gr, mecp2_up_gr) > 0)
  hypo_down <- sum(countOverlaps(hypo_gr, mecp2_down_gr) > 0)

  n_hyper <- length(hyper_gr)
  n_hypo <- length(hypo_gr)

  cat(sprintf("  Hypermethylated DMRs: %d/%d (%.1f%%) overlap MeCP2 Up, %d/%d (%.1f%%) overlap MeCP2 Down\n",
              hyper_up, n_hyper, 100 * hyper_up / n_hyper,
              hyper_down, n_hyper, 100 * hyper_down / n_hyper))
  cat(sprintf("  Hypomethylated DMRs: %d/%d (%.1f%%) overlap MeCP2 Up, %d/%d (%.1f%%) overlap MeCP2 Down\n",
              hypo_up, n_hypo, 100 * hypo_up / n_hypo,
              hypo_down, n_hypo, 100 * hypo_down / n_hypo))

  # Fisher's exact test: hyper/hypo x up/down
  fisher_mat <- matrix(c(hyper_up, hyper_down, hypo_up, hypo_down), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Hypermethylated", "Hypomethylated"),
                                       c("MeCP2 Up", "MeCP2 Down")))
  fisher_result <- fisher.test(fisher_mat)

  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_result$estimate, fisher_result$p.value))

  # Build plot data
  overlap_df <- data.frame(
    DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
    MeCP2_Direction = rep(c("MeCP2 Up", "MeCP2 Down"), 2),
    Count = c(hyper_up, hyper_down, hypo_up, hypo_down),
    Total = rep(c(n_hyper, n_hypo), each = 2)
  ) %>%
    mutate(Percentage = 100 * Count / Total)

  overlap_df$DMR_Direction <- factor(overlap_df$DMR_Direction,
                                     levels = c("Hypermethylated", "Hypomethylated"))
  overlap_df$MeCP2_Direction <- factor(overlap_df$MeCP2_Direction,
                                       levels = c("MeCP2 Up", "MeCP2 Down"))

  p_11a <- ggplot(overlap_df, aes(x = DMR_Direction, y = Percentage, fill = MeCP2_Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = COLORS$mecp2[c("MeCP2 Up", "MeCP2 Down")],
                      name = "MeCP2 Binding") +
    scale_y_continuous(limits = c(0, max(overlap_df$Percentage) * 1.3), expand = c(0, 0)) +
    labs(
      title = "MeCP2 Peak Overlap at Differentially Methylated Regions",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e",
                         fisher_result$estimate, fisher_result$p.value),
      x = "DMR Direction (5mC)", y = "% of DMRs Overlapping MeCP2 Peaks"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_11a, file.path(OUTPUT_DIR, "11a_mecp2_overlap"),
                          width = 8, height = 7)

  # Save overlap summary table
  overlap_summary <- data.frame(
    DMR_Direction = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    MeCP2_Direction = c("MeCP2 Up", "MeCP2 Down", "MeCP2 Up", "MeCP2 Down"),
    Overlap_Count = c(hyper_up, hyper_down, hypo_up, hypo_down),
    Total_DMRs = c(n_hyper, n_hyper, n_hypo, n_hypo),
    Percentage = c(100 * hyper_up / n_hyper, 100 * hyper_down / n_hyper,
                   100 * hypo_up / n_hypo, 100 * hypo_down / n_hypo),
    Fisher_OR = fisher_result$estimate,
    Fisher_pvalue = fisher_result$p.value
  )
  write.table(overlap_summary, file.path(TABLES_DIR, "mecp2_dmr_overlap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: mecp2_dmr_overlap_summary.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 11b: MeCP2 Fold Change Distribution by DMR Direction
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 11b: MeCP2 fold change by DMR direction...\n")

  # Join MeCP2 gene-level data to DMR data
  mc_dmr_mecp2 <- mc_dmr %>%
    left_join(mecp2_gene %>% dplyr::select(SYMBOL, nearest_fold, nearest_fdr, mean_fold, n_peaks),
              by = c("gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(nearest_fold))

  cat(sprintf("  %d DMR genes matched to MeCP2 data\n", nrow(mc_dmr_mecp2)))

  # Categorize DMR genes
  mc_dmr_mecp2 <- mc_dmr_mecp2 %>%
    mutate(dmr_category = case_when(
      significant & mod_difference > 0 ~ "Hypermethylated",
      significant & mod_difference < 0 ~ "Hypomethylated",
      TRUE ~ "Not Significant"
    ))

  mc_dmr_mecp2$dmr_category <- factor(mc_dmr_mecp2$dmr_category,
                                      levels = c("Hypermethylated", "Hypomethylated", "Not Significant"))

  # Wilcoxon test between groups
  hyper_fold <- mc_dmr_mecp2$nearest_fold[mc_dmr_mecp2$dmr_category == "Hypermethylated"]
  hypo_fold <- mc_dmr_mecp2$nearest_fold[mc_dmr_mecp2$dmr_category == "Hypomethylated"]
  ns_fold <- mc_dmr_mecp2$nearest_fold[mc_dmr_mecp2$dmr_category == "Not Significant"]

  wilcox_hyper_hypo <- wilcox.test(hyper_fold, hypo_fold)
  wilcox_hyper_ns <- wilcox.test(hyper_fold, ns_fold)

  cat(sprintf("  Median MeCP2 fold change - Hyper: %.3f, Hypo: %.3f, NS: %.3f\n",
              median(hyper_fold), median(hypo_fold), median(ns_fold)))
  cat(sprintf("  Wilcoxon Hyper vs Hypo: p = %.2e\n", wilcox_hyper_hypo$p.value))
  cat(sprintf("  Wilcoxon Hyper vs NS: p = %.2e\n", wilcox_hyper_ns$p.value))

  # Build significance annotation
  max_fold_val <- max(abs(mc_dmr_mecp2$nearest_fold), na.rm = TRUE)
  y_upper <- max_fold_val * 1.1

  p_11b <- ggplot(mc_dmr_mecp2, aes(x = dmr_category, y = nearest_fold, fill = dmr_category)) +
    geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Hypermethylated" = COLORS$direction["Hypermethylated"],
                                 "Hypomethylated" = COLORS$direction["Hypomethylated"],
                                 "Not Significant" = "grey70")) +
    annotate("text", x = 1.5, y = y_upper,
             label = sprintf("Hyper vs Hypo: p = %.1e", wilcox_hyper_hypo$p.value),
             size = 3.5, fontface = "italic") +
    annotate("text", x = 2, y = y_upper * 0.9,
             label = sprintf("Hyper vs NS: p = %.1e", wilcox_hyper_ns$p.value),
             size = 3.5, fontface = "italic") +
    labs(
      title = "MeCP2 Fold Change Distribution by DMR Direction",
      subtitle = "Prediction: Hypermethylated DMR genes show higher MeCP2 binding (positive fold change)",
      x = "DMR Category (5mC)", y = "MeCP2 log2 Fold Change (Mutant/Control)"
    ) +
    theme_biomodal() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p_11b, file.path(OUTPUT_DIR, "11b_mecp2_fold_by_dmr_direction"),
                          width = 9, height = 7)

  # -----------------------------------------------------------------------
  # FIGURE 11c: Scatter - mC Change vs MeCP2 Fold Change
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 11c: mC change vs MeCP2 fold change scatter...\n")

  scatter_data <- mc_dmr_mecp2 %>%
    dplyr::filter(significant) %>%
    mutate(
      mc_pct = mod_difference * 100,
      label_gene = ifelse(gene %in% KEY_GENES, gene, "")
    )

  # Spearman correlation
  cor_test <- cor.test(scatter_data$mc_pct, scatter_data$nearest_fold,
                       method = "spearman")
  cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))

  # Count genes per quadrant
  q1 <- sum(scatter_data$mc_pct > 0 & scatter_data$nearest_fold > 0)  # mC up, MeCP2 up (predicted)
  q2 <- sum(scatter_data$mc_pct < 0 & scatter_data$nearest_fold > 0)  # mC down, MeCP2 up
  q3 <- sum(scatter_data$mc_pct < 0 & scatter_data$nearest_fold < 0)  # mC down, MeCP2 down (predicted)
  q4 <- sum(scatter_data$mc_pct > 0 & scatter_data$nearest_fold < 0)  # mC up, MeCP2 down

  p_11c <- ggplot(scatter_data, aes(x = mc_pct, y = nearest_fold)) +
    geom_point(aes(color = dmr_category), alpha = 0.4, size = 1.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 15,
                    fontface = "italic", color = "grey20",
                    segment.color = "grey60", segment.size = 0.3) +
    scale_color_manual(values = c("Hypermethylated" = COLORS$direction["Hypermethylated"],
                                  "Hypomethylated" = COLORS$direction["Hypomethylated"]),
                       name = "DMR Direction") +
    # Quadrant annotations
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = max(scatter_data$nearest_fold) * 0.9,
             label = sprintf("mC\u2191 MeCP2\u2191\nn=%d (predicted)", q1),
             size = 3, color = "#D95F02", fontface = "bold") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = min(scatter_data$nearest_fold) * 0.9,
             label = sprintf("mC\u2193 MeCP2\u2193\nn=%d (predicted)", q3),
             size = 3, color = "#7570B3", fontface = "bold") +
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = min(scatter_data$nearest_fold) * 0.9,
             label = sprintf("mC\u2191 MeCP2\u2193\nn=%d", q4),
             size = 3, color = "grey50") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = max(scatter_data$nearest_fold) * 0.9,
             label = sprintf("mC\u2193 MeCP2\u2191\nn=%d", q2),
             size = 3, color = "grey50") +
    labs(
      title = "5mC Change vs MeCP2 Fold Change at Significant DMRs",
      subtitle = sprintf("Spearman \u03C1 = %.3f, p = %.2e | Gene-level nearest peak",
                         cor_test$estimate, cor_test$p.value),
      x = "5mC Change (Mutant - Control, %)",
      y = "MeCP2 log2 Fold Change (Mutant/Control)"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_11c, file.path(OUTPUT_DIR, "11c_mc_vs_mecp2_scatter"),
                          width = 10, height = 9)

  # Save gene-level correlation table
  gene_corr_table <- mc_dmr_mecp2 %>%
    dplyr::filter(significant) %>%
    dplyr::select(gene, chr, start, end, mod_difference, dmr_qvalue,
                  dmr_category, mecp2_nearest_fold = nearest_fold,
                  mecp2_nearest_fdr = nearest_fdr, mecp2_mean_fold = mean_fold,
                  mecp2_n_peaks = n_peaks) %>%
    arrange(dmr_qvalue)

  write.table(gene_corr_table, file.path(TABLES_DIR, "mecp2_gene_level_correlation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: mecp2_gene_level_correlation.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 11d: MeCP2 Signal at Coordinated mC up / hmC down Genes
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 11d: MeCP2 at coordinated genes...\n")

  # Get coordinated genes (mC up, hmC down) with MeCP2 data
  coordinated_mecp2 <- coordinated %>%
    dplyr::filter(coordinated_pattern) %>%
    left_join(mecp2_gene %>% dplyr::select(SYMBOL, nearest_fold, nearest_fdr, mean_fold, n_peaks),
              by = c("gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(nearest_fold))

  cat(sprintf("  %d coordinated genes matched to MeCP2 data\n", nrow(coordinated_mecp2)))

  # Compare coordinated vs all other genes
  coord_gene_list <- coordinated$gene[coordinated$coordinated_pattern]
  mecp2_gene_classified <- mecp2_gene %>%
    mutate(category = ifelse(SYMBOL %in% coord_gene_list,
                             "Coordinated\n(mC\u2191/hmC\u2193)",
                             "All Other Genes"))

  wilcox_coord <- wilcox.test(
    mecp2_gene_classified$nearest_fold[mecp2_gene_classified$category == "Coordinated\n(mC\u2191/hmC\u2193)"],
    mecp2_gene_classified$nearest_fold[mecp2_gene_classified$category == "All Other Genes"]
  )
  cat(sprintf("  Wilcoxon coordinated vs other: p = %.2e\n", wilcox_coord$p.value))

  # Boxplot comparison
  p_11d_box <- ggplot(mecp2_gene_classified, aes(x = category, y = nearest_fold, fill = category)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Coordinated\n(mC\u2191/hmC\u2193)" = "#D95F02",
                                 "All Other Genes" = "grey70")) +
    annotate("text", x = 1.5, y = max(mecp2_gene_classified$nearest_fold, na.rm = TRUE) * 0.95,
             label = sprintf("Wilcoxon p = %.2e", wilcox_coord$p.value),
             size = 3.5, fontface = "italic") +
    labs(
      title = "MeCP2 Binding at Coordinated Genes",
      subtitle = "Genes with mC increase + hmC decrease vs all other genes",
      x = "", y = "MeCP2 log2 Fold Change\n(Mutant/Control)"
    ) +
    theme_biomodal() +
    theme(legend.position = "none")

  # Top coordinated genes bar chart with MeCP2 annotation
  top_coord_mecp2 <- coordinated_mecp2 %>%
    arrange(desc(combined_effect)) %>%
    head(20) %>%
    mutate(
      gene = factor(gene, levels = rev(unique(gene))),
      mecp2_label = sprintf("MeCP2: %+.2f", nearest_fold),
      mecp2_direction = ifelse(nearest_fold > 0, "Up", "Down")
    )

  p_11d_bar <- ggplot(top_coord_mecp2, aes(x = gene, y = combined_effect * 100)) +
    geom_bar(aes(fill = mecp2_direction), stat = "identity", width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_text(aes(label = mecp2_label), hjust = -0.1, size = 2.8, color = "grey30") +
    scale_fill_manual(values = c("Up" = "#D95F02", "Down" = "#7570B3"),
                      name = "MeCP2 Direction") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "Top 20 Coordinated Genes",
      subtitle = "mC\u2191/hmC\u2193 genes with MeCP2 fold change annotation",
      x = "", y = "Combined Effect (|mC| + |hmC| change, %)"
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  p_11d <- p_11d_box + p_11d_bar +
    plot_layout(widths = c(1, 2)) +
    plot_annotation(
      title = "MeCP2 Signal at Coordinated mC\u2191/hmC\u2193 Genes",
      subtitle = "Do genes with impaired TET-mediated demethylation show increased MeCP2 occupancy?",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_11d, file.path(OUTPUT_DIR, "11d_mecp2_coordinated_genes"),
                          width = 15, height = 9)

  # Save coordinated genes with MeCP2 annotations
  coord_export <- coordinated_mecp2 %>%
    dplyr::select(gene, mc_diff, hmc_diff, mc_q, hmc_q, combined_effect,
                  mecp2_nearest_fold = nearest_fold, mecp2_nearest_fdr = nearest_fdr,
                  mecp2_mean_fold = mean_fold, mecp2_n_peaks = n_peaks) %>%
    arrange(desc(combined_effect))

  write.table(coord_export, file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: mecp2_coordinated_genes.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 11e: Summary Heatmap - mC Direction x MeCP2 Direction
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 11e: Integration heatmap...\n")

  # Stratify genes into 4 quadrants by DMR direction x MeCP2 direction
  sig_with_mecp2 <- mc_dmr_mecp2 %>%
    dplyr::filter(significant) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      mecp2_direction = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
    )

  # Contingency table
  quad_table <- table(sig_with_mecp2$mc_direction, sig_with_mecp2$mecp2_direction)
  cat("  Contingency table:\n")
  print(quad_table)

  # Fisher's exact test on the 2x2 table
  fisher_quad <- fisher.test(quad_table)
  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_quad$estimate, fisher_quad$p.value))

  # Compute expected counts and enrichment (observed/expected)
  total_n <- sum(quad_table)
  row_totals <- rowSums(quad_table)
  col_totals <- colSums(quad_table)
  expected <- outer(row_totals, col_totals) / total_n

  enrichment <- quad_table / expected

  # Build heatmap data
  heatmap_df <- as.data.frame(as.table(quad_table))
  colnames(heatmap_df) <- c("mc_direction", "mecp2_direction", "Observed")

  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("mc_direction", "mecp2_direction", "Expected")

  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("mc_direction", "mecp2_direction", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("mc_direction", "mecp2_direction")) %>%
    left_join(enrichment_df, by = c("mc_direction", "mecp2_direction")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nOR: %.2f", Observed, Expected, Enrichment),
      mc_direction = factor(mc_direction, levels = c("mC Up", "mC Down")),
      mecp2_direction = factor(mecp2_direction, levels = c("MeCP2 Up", "MeCP2 Down")),
      is_predicted = (mc_direction == "mC Up" & mecp2_direction == "MeCP2 Up") |
                     (mc_direction == "mC Down" & mecp2_direction == "MeCP2 Down")
    )

  p_11e <- ggplot(heatmap_data, aes(x = mecp2_direction, y = mc_direction, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2) +
    # Highlight predicted quadrants
    geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
              aes(x = mecp2_direction, y = mc_direction),
              fill = NA, color = "black", linewidth = 1.5, linetype = "solid") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(
      title = "Integration: 5mC Direction x MeCP2 Binding Direction",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Black borders = predicted quadrants",
                         fisher_quad$estimate, fisher_quad$p.value),
      x = "MeCP2 Binding Direction (Mutant vs Control)",
      y = "5mC DMR Direction"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  save_multiformat_ggplot(p_11e, file.path(OUTPUT_DIR, "11e_mecp2_integration_heatmap"),
                          width = 9, height = 7)

  # -----------------------------------------------------------------------
  # Print summary
  # -----------------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 11 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Total MeCP2 peaks analyzed: %d\n", nrow(mecp2_annotated)))
  cat(sprintf("Genes with MeCP2 peaks: %d\n", nrow(mecp2_gene)))
  cat(sprintf("DMR genes matched to MeCP2: %d\n", nrow(mc_dmr_mecp2 %>% dplyr::filter(significant))))
  cat(sprintf("Coordinated genes with MeCP2 data: %d\n", nrow(coordinated_mecp2)))
  cat(sprintf("\nOverlap analysis: Fisher OR = %.2f, p = %.2e\n",
              fisher_result$estimate, fisher_result$p.value))
  cat(sprintf("Fold change correlation: Spearman rho = %.3f, p = %.2e\n",
              cor_test$estimate, cor_test$p.value))
  cat(sprintf("Coordinated vs other: Wilcoxon p = %.2e\n", wilcox_coord$p.value))
  cat(sprintf("Integration: Fisher OR = %.2f, p = %.2e\n",
              fisher_quad$estimate, fisher_quad$p.value))
  cat(sprintf("\nMedian MeCP2 fold change:\n"))
  cat(sprintf("  Hypermethylated DMRs: %.3f\n", median(hyper_fold)))
  cat(sprintf("  Hypomethylated DMRs:  %.3f\n", median(hypo_fold)))
  cat(sprintf("  Not significant:      %.3f\n", median(ns_fold)))
}

cat("\nSection 11 complete.\n\n")
