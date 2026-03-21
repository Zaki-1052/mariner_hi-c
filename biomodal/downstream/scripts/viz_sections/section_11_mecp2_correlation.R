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
  # FIGURE 11h: Fisher's Exact on Binary MeCP2 Gain (Coordinated vs Other)
  # -----------------------------------------------------------------------
  # Tests whether coordinated genes (mC up/hmC down) are enriched for
  # significant MeCP2 binding increase compared to non-coordinated genes.
  cat("\nCreating Figure 11h: Fisher's exact on binary MeCP2 gain...\n")

  # Define binary MeCP2 significant gain per gene
  mecp2_gene_binary <- mecp2_gene %>%
    mutate(mecp2_sig_gain = nearest_fdr < 0.05 & nearest_fold > 0)

  # Classify genes: coordinated vs non-coordinated
  mecp2_gene_binary <- mecp2_gene_binary %>%
    mutate(gene_category = ifelse(SYMBOL %in% coord_gene_list,
                                  "Coordinated\n(mC\u2191/hmC\u2193)",
                                  "Non-coordinated"))

  # Counts for 2x2 table
  coord_gain    <- sum(mecp2_gene_binary$gene_category == "Coordinated\n(mC\u2191/hmC\u2193)" &
                       mecp2_gene_binary$mecp2_sig_gain)
  coord_no_gain <- sum(mecp2_gene_binary$gene_category == "Coordinated\n(mC\u2191/hmC\u2193)" &
                       !mecp2_gene_binary$mecp2_sig_gain)
  other_gain    <- sum(mecp2_gene_binary$gene_category == "Non-coordinated" &
                       mecp2_gene_binary$mecp2_sig_gain)
  other_no_gain <- sum(mecp2_gene_binary$gene_category == "Non-coordinated" &
                       !mecp2_gene_binary$mecp2_sig_gain)

  fisher_binary_mat <- matrix(c(coord_gain, other_gain, coord_no_gain, other_no_gain),
                              nrow = 2,
                              dimnames = list(c("Coordinated", "Non-coordinated"),
                                              c("MeCP2 Gain", "No MeCP2 Gain")))
  fisher_binary <- fisher.test(fisher_binary_mat)

  cat(sprintf("  Coordinated: %d/%d (%.1f%%) with MeCP2 sig gain\n",
              coord_gain, coord_gain + coord_no_gain,
              100 * coord_gain / (coord_gain + coord_no_gain)))
  cat(sprintf("  Non-coordinated: %d/%d (%.1f%%) with MeCP2 sig gain\n",
              other_gain, other_gain + other_no_gain,
              100 * other_gain / (other_gain + other_no_gain)))
  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_binary$estimate, fisher_binary$p.value))

  # Build plot data
  binary_plot_df <- data.frame(
    gene_category = c("Coordinated\n(mC\u2191/hmC\u2193)", "Coordinated\n(mC\u2191/hmC\u2193)",
                       "Non-coordinated", "Non-coordinated"),
    mecp2_status = rep(c("MeCP2 Sig Gain", "No MeCP2 Sig Gain"), 2),
    count = c(coord_gain, coord_no_gain, other_gain, other_no_gain),
    total = c(coord_gain + coord_no_gain, coord_gain + coord_no_gain,
              other_gain + other_no_gain, other_gain + other_no_gain),
    stringsAsFactors = FALSE
  ) %>%
    mutate(percentage = 100 * count / total)

  binary_plot_df$gene_category <- factor(binary_plot_df$gene_category,
                                          levels = c("Coordinated\n(mC\u2191/hmC\u2193)",
                                                     "Non-coordinated"))
  binary_plot_df$mecp2_status <- factor(binary_plot_df$mecp2_status,
                                         levels = c("MeCP2 Sig Gain", "No MeCP2 Sig Gain"))

  p_11h <- ggplot(binary_plot_df %>% dplyr::filter(mecp2_status == "MeCP2 Sig Gain"),
                  aes(x = gene_category, y = percentage, fill = gene_category)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", percentage, count, total)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Coordinated\n(mC\u2191/hmC\u2193)" = "#D95F02",
                                  "Non-coordinated" = "grey70")) +
    scale_y_continuous(limits = c(0, max(binary_plot_df$percentage[binary_plot_df$mecp2_status == "MeCP2 Sig Gain"]) * 1.4),
                       expand = c(0, 0)) +
    labs(
      title = "Binary MeCP2 Gain Enrichment at Coordinated Genes",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | MeCP2 sig gain = FDR < 0.05 & Fold > 0",
                         fisher_binary$estimate, fisher_binary$p.value),
      x = "Gene Category", y = "% of Genes with MeCP2 Significant Gain"
    ) +
    theme_biomodal() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p_11h, file.path(OUTPUT_DIR, "11h_mecp2_binary_gain_fisher"),
                          width = 8, height = 7)

  # Save Fisher's exact table
  fisher_binary_export <- data.frame(
    gene_category = c("Coordinated", "Coordinated", "Non-coordinated", "Non-coordinated"),
    mecp2_gain = c(TRUE, FALSE, TRUE, FALSE),
    count = c(coord_gain, coord_no_gain, other_gain, other_no_gain),
    total = c(coord_gain + coord_no_gain, coord_gain + coord_no_gain,
              other_gain + other_no_gain, other_gain + other_no_gain),
    percentage = c(100 * coord_gain / (coord_gain + coord_no_gain),
                   100 * coord_no_gain / (coord_gain + coord_no_gain),
                   100 * other_gain / (other_gain + other_no_gain),
                   100 * other_no_gain / (other_gain + other_no_gain)),
    fisher_OR = fisher_binary$estimate,
    fisher_pvalue = fisher_binary$p.value,
    stringsAsFactors = FALSE
  )
  write.table(fisher_binary_export, file.path(TABLES_DIR, "mecp2_binary_gain_coordinated_fisher.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: mecp2_binary_gain_coordinated_fisher.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 11i: Scatter Stratified by Coordinated vs Discordant
  # -----------------------------------------------------------------------
  # Tests whether the mC-to-MeCP2 relationship is stronger among
  # coordinated genes (mC up/hmC down) vs discordant or mC-only genes.
  cat("\nCreating Figure 11i: mC vs MeCP2 scatter stratified by coordination status...\n")

  # Build hmC-significant gene lookup for stratification
  hmc_sig_genes <- hmc_sig_coord$gene

  # For each significant mC DMR gene, determine coordination status
  stratified_scatter <- mc_dmr_mecp2 %>%
    dplyr::filter(significant) %>%
    mutate(mc_pct = mod_difference * 100) %>%
    left_join(
      coordinated %>% dplyr::select(gene, hmc_diff, coordinated_pattern),
      by = "gene"
    ) %>%
    mutate(
      coord_status = case_when(
        !is.na(coordinated_pattern) & coordinated_pattern ~ "Coordinated\n(mC\u2191/hmC\u2193)",
        !is.na(coordinated_pattern) & !coordinated_pattern ~ "Discordant",
        TRUE ~ "mC-only"
      )
    )

  stratified_scatter$coord_status <- factor(stratified_scatter$coord_status,
                                             levels = c("Coordinated\n(mC\u2191/hmC\u2193)",
                                                        "Discordant", "mC-only"))

  # Per-group Spearman correlations
  strat_groups <- levels(stratified_scatter$coord_status)
  strat_cor_results <- data.frame(
    group = character(), n_genes = integer(),
    spearman_rho = numeric(), spearman_p = numeric(),
    median_mc_diff = numeric(), median_mecp2_fold = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp in strat_groups) {
    grp_data <- stratified_scatter %>% dplyr::filter(coord_status == grp)
    if (nrow(grp_data) >= 5) {
      ct <- cor.test(grp_data$mc_pct, grp_data$nearest_fold, method = "spearman")
      strat_cor_results <- rbind(strat_cor_results, data.frame(
        group = grp,
        n_genes = nrow(grp_data),
        spearman_rho = ct$estimate,
        spearman_p = ct$p.value,
        median_mc_diff = median(grp_data$mc_pct, na.rm = TRUE),
        median_mecp2_fold = median(grp_data$nearest_fold, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  %s (n=%d): rho = %.3f, p = %.2e\n",
                  grp, nrow(grp_data), ct$estimate, ct$p.value))
    } else {
      cat(sprintf("  %s (n=%d): too few genes for correlation\n", grp, nrow(grp_data)))
    }
  }

  # Color palette for coordination status
  coord_colors <- c(
    "Coordinated\n(mC\u2191/hmC\u2193)" = "#4DAF4A",
    "Discordant" = "#984EA3",
    "mC-only" = "grey60"
  )

  # Build annotation label for rho values
  rho_labels <- paste(
    sapply(seq_len(nrow(strat_cor_results)), function(i) {
      sprintf("%s: \u03C1=%.3f, p=%.1e (n=%d)",
              gsub("\n", " ", strat_cor_results$group[i]),
              strat_cor_results$spearman_rho[i],
              strat_cor_results$spearman_p[i],
              strat_cor_results$n_genes[i])
    }),
    collapse = "\n"
  )

  p_11i <- ggplot(stratified_scatter, aes(x = mc_pct, y = nearest_fold, color = coord_status)) +
    geom_point(alpha = 0.3, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_manual(values = coord_colors, name = "Coordination Status") +
    annotate("text", x = min(stratified_scatter$mc_pct, na.rm = TRUE) * 0.5,
             y = max(stratified_scatter$nearest_fold, na.rm = TRUE) * 0.95,
             label = rho_labels, size = 3, hjust = 0, vjust = 1, lineheight = 1.2) +
    labs(
      title = "5mC Change vs MeCP2 Fold Change: Stratified by Coordination Status",
      subtitle = "Prediction: Coordinated genes (mC\u2191/hmC\u2193) show stronger mC\u2192MeCP2 relationship",
      x = "5mC Change (Mutant - Control, %)",
      y = "MeCP2 log2 Fold Change (Mutant/Control)"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_11i, file.path(OUTPUT_DIR, "11i_mc_vs_mecp2_stratified"),
                          width = 12, height = 9)

  # Save per-group correlation table
  write.table(strat_cor_results, file.path(TABLES_DIR, "mecp2_stratified_correlations.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: mecp2_stratified_correlations.tsv\n")

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
  # FIGURES 11f-11g: DELTA-RATIO REGRESSION MODELS
  # -----------------------------------------------------------------------
  # Test whether delta_ratio (TET efficiency change) predicts MeCP2 binding.
  # Linear: does delta_ratio predict magnitude of MeCP2 fold change?
  # Logistic: does delta_ratio predict probability of significant MeCP2 increase?
  #
  # Expected: negative coefficient (more TET impairment -> more MeCP2 binding)
  # because impaired TET -> more mC retained -> more MeCP2 recruited.

  cat("\n")
  cat("================================================================================\n")
  cat("DELTA-RATIO REGRESSION MODELS (Figures 11f-11g)\n")
  cat("================================================================================\n\n")

  stopifnot("broom package required for regression models" =
              requireNamespace("broom", quietly = TRUE))
  suppressPackageStartupMessages(library(broom))

  delta_ratio_path <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
  stopifnot("demethylation_ratio_all_genes.tsv not found — run Section 22 first" =
              file.exists(delta_ratio_path))

  delta_ratio_df <- read.table(delta_ratio_path, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
  cat(sprintf("Loaded %d genes with delta_ratio values\n", nrow(delta_ratio_df)))

  # Join delta_ratio to mc_dmr_mecp2 (already built at line 213)
  reg_data <- mc_dmr_mecp2 %>%
    dplyr::filter(significant) %>%
    left_join(delta_ratio_df %>% dplyr::select(gene, delta_ratio),
              by = "gene") %>%
    dplyr::filter(!is.na(delta_ratio) & !is.na(nearest_fold) & !is.na(mean_coverage))

  cat(sprintf("Regression dataset: %d significant DMR genes with delta_ratio + MeCP2 data\n",
              nrow(reg_data)))

  if (nrow(reg_data) >= 20) {

    # -------------------------------------------------------------------
    # FIGURE 11f: Linear Model — MeCP2 fold ~ delta_ratio + coverage
    # -------------------------------------------------------------------
    cat("\nFitting linear model: nearest_fold ~ delta_ratio + mean_coverage...\n")

    lm_fit <- lm(nearest_fold ~ delta_ratio + mean_coverage, data = reg_data)
    lm_tidy <- tidy(lm_fit, conf.int = TRUE)
    lm_glance <- glance(lm_fit)

    cat(sprintf("  R-squared: %.4f, Adj R-squared: %.4f\n",
                lm_glance$r.squared, lm_glance$adj.r.squared))
    cat(sprintf("  F-statistic: %.2f, p = %.2e\n", lm_glance$statistic, lm_glance$p.value))
    cat("  Coefficients:\n")
    for (i in seq_len(nrow(lm_tidy))) {
      cat(sprintf("    %s: %.4f (SE=%.4f, p=%.2e)\n",
                  lm_tidy$term[i], lm_tidy$estimate[i],
                  lm_tidy$std.error[i], lm_tidy$p.value[i]))
    }

    # Save coefficients table
    write.table(lm_tidy, file.path(TABLES_DIR, "mecp2_delta_ratio_lm_coefficients.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Saved: mecp2_delta_ratio_lm_coefficients.tsv\n")

    # Left panel: Coefficient forest plot (exclude intercept)
    lm_coef_df <- lm_tidy %>%
      dplyr::filter(term != "(Intercept)") %>%
      mutate(
        term_label = case_when(
          term == "delta_ratio" ~ "Delta-Ratio\n(TET efficiency change)",
          term == "mean_coverage" ~ "Mean Coverage\n(technical covariate)"
        ),
        sig_label = ifelse(p.value < 0.05, "p < 0.05", "p >= 0.05")
      )

    p_11f_forest <- ggplot(lm_coef_df, aes(x = estimate, y = term_label, color = sig_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 1, linewidth = 0.8) +
      scale_color_manual(values = c("p < 0.05" = "#E41A1C", "p >= 0.05" = "grey50"),
                         name = "Significance") +
      labs(
        title = "Linear Model Coefficients",
        subtitle = sprintf("R\u00B2 = %.3f, F(%.0f,%.0f) = %.1f, p = %.2e",
                           lm_glance$r.squared, lm_glance$df, lm_glance$df.residual,
                           lm_glance$statistic, lm_glance$p.value),
        x = "Coefficient Estimate (95% CI)", y = ""
      ) +
      theme_biomodal() +
      theme(legend.position = "bottom")

    # Right panel: Partial regression plot (delta_ratio effect)
    # Residualize both Y and X on coverage to isolate delta_ratio effect
    reg_data$resid_fold <- residuals(lm(nearest_fold ~ mean_coverage, data = reg_data))
    reg_data$resid_delta <- residuals(lm(delta_ratio ~ mean_coverage, data = reg_data))

    p_11f_partial <- ggplot(reg_data, aes(x = resid_delta, y = resid_fold)) +
      geom_point(alpha = 0.3, size = 1.5, color = "grey40") +
      geom_smooth(method = "lm", color = "#E41A1C", linewidth = 0.8, se = TRUE, alpha = 0.15) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      labs(
        title = "Partial Regression: Delta-Ratio Effect",
        subtitle = "After residualizing out mean coverage",
        x = "Delta-Ratio (residualized)",
        y = "MeCP2 Fold Change (residualized)"
      ) +
      theme_biomodal()

    p_11f <- p_11f_forest + p_11f_partial +
      plot_layout(widths = c(1, 1.3)) +
      plot_annotation(
        title = "Linear Model: MeCP2 Fold Change ~ Delta-Ratio + Coverage",
        subtitle = "Does TET efficiency change predict MeCP2 binding magnitude at significant DMRs?",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
        )
      )

    save_multiformat_ggplot(p_11f, file.path(OUTPUT_DIR, "11f_mecp2_delta_ratio_lm"),
                            width = 14, height = 7)

    # -------------------------------------------------------------------
    # FIGURE 11g: Logistic Model — MeCP2 sig up ~ delta_ratio + coverage
    # -------------------------------------------------------------------
    cat("\nFitting logistic model: mecp2_up ~ delta_ratio + mean_coverage...\n")

    reg_data$mecp2_up <- as.integer(reg_data$nearest_fold > 0 & reg_data$nearest_fdr < 0.05)

    n_up <- sum(reg_data$mecp2_up == 1)
    n_down <- sum(reg_data$mecp2_up == 0)
    cat(sprintf("  MeCP2 significant up: %d, other: %d\n", n_up, n_down))

    if (n_up >= 5 && n_down >= 5) {
      glm_fit <- glm(mecp2_up ~ delta_ratio + mean_coverage,
                      data = reg_data, family = binomial)
      glm_tidy <- tidy(glm_fit, conf.int = TRUE, exponentiate = TRUE)
      glm_glance <- glance(glm_fit)

      cat(sprintf("  AIC: %.1f, Deviance: %.1f\n", glm_glance$AIC, glm_glance$deviance))
      cat("  Odds ratios:\n")
      for (i in seq_len(nrow(glm_tidy))) {
        cat(sprintf("    %s: OR = %.4f (95%% CI: %.4f-%.4f, p = %.2e)\n",
                    glm_tidy$term[i], glm_tidy$estimate[i],
                    glm_tidy$conf.low[i], glm_tidy$conf.high[i],
                    glm_tidy$p.value[i]))
      }

      # Save odds ratios table
      write.table(glm_tidy, file.path(TABLES_DIR, "mecp2_delta_ratio_glm_odds_ratios.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved: mecp2_delta_ratio_glm_odds_ratios.tsv\n")

      # Left panel: Odds ratio forest plot (exclude intercept)
      glm_coef_df <- glm_tidy %>%
        dplyr::filter(term != "(Intercept)") %>%
        mutate(
          term_label = case_when(
            term == "delta_ratio" ~ "Delta-Ratio\n(TET efficiency change)",
            term == "mean_coverage" ~ "Mean Coverage\n(technical covariate)"
          ),
          sig_label = ifelse(p.value < 0.05, "p < 0.05", "p >= 0.05")
        )

      p_11g_forest <- ggplot(glm_coef_df, aes(x = estimate, y = term_label, color = sig_label)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 1, linewidth = 0.8) +
        scale_color_manual(values = c("p < 0.05" = "#E41A1C", "p >= 0.05" = "grey50"),
                           name = "Significance") +
        scale_x_log10() +
        labs(
          title = "Logistic Model: Odds Ratios",
          subtitle = sprintf("AIC = %.1f | Reference line at OR = 1", glm_glance$AIC),
          x = "Odds Ratio (95% CI, log scale)", y = ""
        ) +
        theme_biomodal() +
        theme(legend.position = "bottom")

      # Right panel: Predicted probability curve vs delta_ratio
      # Hold mean_coverage at median
      median_coverage <- median(reg_data$mean_coverage, na.rm = TRUE)
      delta_seq <- seq(min(reg_data$delta_ratio), max(reg_data$delta_ratio), length.out = 200)
      pred_df <- data.frame(delta_ratio = delta_seq, mean_coverage = median_coverage)
      pred_df$prob <- predict(glm_fit, newdata = pred_df, type = "response")

      # Get raw coefficients for annotation
      glm_tidy_raw <- tidy(glm_fit, conf.int = TRUE, exponentiate = FALSE)
      delta_coef <- glm_tidy_raw$estimate[glm_tidy_raw$term == "delta_ratio"]

      p_11g_prob <- ggplot() +
        geom_rug(data = reg_data, aes(x = delta_ratio, color = factor(mecp2_up)),
                 sides = "b", alpha = 0.3) +
        geom_line(data = pred_df, aes(x = delta_ratio, y = prob),
                  color = "#E41A1C", linewidth = 1) +
        geom_ribbon(data = pred_df, aes(x = delta_ratio, ymin = 0, ymax = prob),
                    fill = "#E41A1C", alpha = 0.08) +
        scale_color_manual(values = c("0" = "grey60", "1" = "#D95F02"),
                           labels = c("Not sig up", "MeCP2 sig up"),
                           name = "") +
        labs(
          title = "Predicted Probability of MeCP2 Sig. Increase",
          subtitle = sprintf("Coverage held at median (%.0f) | \u03B2(delta_ratio) = %.3f",
                             median_coverage, delta_coef),
          x = "Delta-Ratio (TET efficiency change)",
          y = "P(MeCP2 significantly up)"
        ) +
        theme_biomodal() +
        theme(legend.position = "bottom")

      p_11g <- p_11g_forest + p_11g_prob +
        plot_layout(widths = c(1, 1.3)) +
        plot_annotation(
          title = "Logistic Model: MeCP2 Significant Increase ~ Delta-Ratio + Coverage",
          subtitle = "Does TET efficiency change predict probability of significant MeCP2 binding increase?",
          theme = theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
          )
        )

      save_multiformat_ggplot(p_11g, file.path(OUTPUT_DIR, "11g_mecp2_delta_ratio_glm"),
                              width = 14, height = 7)

    } else {
      cat("  Not enough events in both classes for logistic regression\n")
      glm_fit <- NULL
    }

  } else {
    cat("  Not enough genes for regression models (need >= 20)\n")
    lm_fit <- NULL
    glm_fit <- NULL
  }

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
  cat(sprintf("Binary MeCP2 gain (coordinated vs other): Fisher OR = %.2f, p = %.2e\n",
              fisher_binary$estimate, fisher_binary$p.value))
  cat(sprintf("Integration: Fisher OR = %.2f, p = %.2e\n",
              fisher_quad$estimate, fisher_quad$p.value))
  cat(sprintf("\nMedian MeCP2 fold change:\n"))
  cat(sprintf("  Hypermethylated DMRs: %.3f\n", median(hyper_fold)))
  cat(sprintf("  Hypomethylated DMRs:  %.3f\n", median(hypo_fold)))
  cat(sprintf("  Not significant:      %.3f\n", median(ns_fold)))

  # Delta-ratio regression summary
  if (exists("lm_fit") && !is.null(lm_fit)) {
    lm_g <- glance(lm_fit)
    lm_t <- tidy(lm_fit)
    delta_row <- lm_t[lm_t$term == "delta_ratio", ]
    cat(sprintf("\nDelta-ratio linear model:\n"))
    cat(sprintf("  R-squared: %.4f, F-statistic p = %.2e\n",
                lm_g$r.squared, lm_g$p.value))
    cat(sprintf("  delta_ratio coefficient: %.4f (p = %.2e)\n",
                delta_row$estimate, delta_row$p.value))
    if (delta_row$estimate < 0) {
      cat("  Direction: NEGATIVE (consistent with prediction — more TET impairment -> more MeCP2)\n")
    } else {
      cat("  Direction: POSITIVE (unexpected — opposite to prediction)\n")
    }
  }

  if (exists("glm_fit") && !is.null(glm_fit)) {
    glm_t <- tidy(glm_fit, exponentiate = TRUE)
    delta_or <- glm_t[glm_t$term == "delta_ratio", ]
    cat(sprintf("\nDelta-ratio logistic model:\n"))
    cat(sprintf("  delta_ratio OR: %.4f (95%% CI: %.4f-%.4f, p = %.2e)\n",
                delta_or$estimate, delta_or$conf.low, delta_or$conf.high,
                delta_or$p.value))
    if (delta_or$estimate < 1) {
      cat("  Direction: OR < 1 (consistent with prediction — more TET impairment -> higher MeCP2 up probability)\n")
    } else {
      cat("  Direction: OR > 1 (unexpected — opposite to prediction)\n")
    }
  }
}

cat("\nSection 11 complete.\n\n")
