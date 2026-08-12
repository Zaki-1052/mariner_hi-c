# biomodal/downstream/scripts/viz_sections/section_12_atac_correlation.R
# Section 12: ATAC-seq Accessibility Correlation at DMRs
# Standalone script - sources shared config for all dependencies and data
# Tests whether DMR methylation changes correlate with chromatin accessibility changes.
#
# Biological prediction (Conway et al. 2021):
#   BAP1 loss -> H2AK119ub accumulation -> chromatin compaction
#   Therefore:
#     - Hypermethylated DMRs (mC up / hmC down) -> DECREASED accessibility (ATAC down)
#     - Hypomethylated DMRs (mC down)           -> INCREASED accessibility (ATAC up)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_12_atac_correlation.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# Load ChIPseeker for gene annotation of ATAC peaks (not in shared config)
suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# ATAC-SPECIFIC HELPER FUNCTIONS
# =============================================================================

#' Annotate ATAC peaks to nearest genes using ChIPseeker
#' @param atac_gr GRanges of ATAC peaks
#' @param txdb TxDb object for gene annotation
#' @return data.frame with peak-to-gene mapping
annotate_atac_to_genes <- function(atac_gr, txdb) {
  anno <- suppressMessages(annotatePeak(
    atac_gr,
    TxDb = txdb,
    annoDb = "org.Mm.eg.db",
    level = "gene"
  ))
  as.data.frame(anno)
}

#' Aggregate ATAC peaks per gene
#' @param anno_df Annotated ATAC peak data.frame (from ChIPseeker)
#' @param direction Character label for ATAC direction ("up" or "down")
#' @return data.frame with gene-level ATAC peak counts
aggregate_atac_by_gene <- function(anno_df, direction) {
  anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      !!paste0("n_atac_", direction) := n(),
      !!paste0("min_dist_tss_", direction) := min(abs(distanceToTSS), na.rm = TRUE),
      .groups = "drop"
    )
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
# SECTION 12: ATAC-seq CORRELATION AT DMRs
# =============================================================================

cat("================================================================================\n")
cat("SECTION 12: ATAC-seq ACCESSIBILITY CORRELATION AT DMRs\n")
cat("================================================================================\n\n")

# Load ATAC differential peak data
cat("Loading ATAC-seq differential peaks...\n")
atac_up_gr <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

if (is.null(atac_up_gr) || is.null(atac_down_gr)) {
  cat("  ERROR: ATAC differential peak files not found. Skipping Section 12.\n")
} else {

  # Get significant mC DMRs split by direction
  mc_hyper <- mc_dmr %>% dplyr::filter(significant & mod_difference > 0)
  mc_hypo <- mc_dmr %>% dplyr::filter(significant & mod_difference < 0)

  # Convert to GRanges
  hyper_gr <- dmr_to_granges(mc_hyper)
  hypo_gr <- dmr_to_granges(mc_hypo)

  # -----------------------------------------------------------------------
  # FIGURE 12a: ATAC Differential Peak Overlap at DMRs
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 12a: ATAC differential peak overlap at DMRs...\n")

  # Compute overlaps
  hyper_up <- sum(countOverlaps(hyper_gr, atac_up_gr) > 0)
  hyper_down <- sum(countOverlaps(hyper_gr, atac_down_gr) > 0)
  hypo_up <- sum(countOverlaps(hypo_gr, atac_up_gr) > 0)
  hypo_down <- sum(countOverlaps(hypo_gr, atac_down_gr) > 0)

  n_hyper <- length(hyper_gr)
  n_hypo <- length(hypo_gr)

  cat(sprintf("  Hypermethylated DMRs: %d/%d (%.1f%%) overlap ATAC Up, %d/%d (%.1f%%) overlap ATAC Down\n",
              hyper_up, n_hyper, 100 * hyper_up / n_hyper,
              hyper_down, n_hyper, 100 * hyper_down / n_hyper))
  cat(sprintf("  Hypomethylated DMRs: %d/%d (%.1f%%) overlap ATAC Up, %d/%d (%.1f%%) overlap ATAC Down\n",
              hypo_up, n_hypo, 100 * hypo_up / n_hypo,
              hypo_down, n_hypo, 100 * hypo_down / n_hypo))

  # Fisher's exact test: hyper/hypo x up/down
  fisher_mat <- matrix(c(hyper_up, hyper_down, hypo_up, hypo_down), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Hypermethylated", "Hypomethylated"),
                                       c("ATAC Up", "ATAC Down")))
  fisher_result <- fisher.test(fisher_mat)

  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_result$estimate, fisher_result$p.value))

  # Build plot data
  overlap_df <- data.frame(
    DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
    ATAC_Direction = rep(c("ATAC Up", "ATAC Down"), 2),
    Count = c(hyper_up, hyper_down, hypo_up, hypo_down),
    Total = rep(c(n_hyper, n_hypo), each = 2)
  ) %>%
    mutate(Percentage = 100 * Count / Total)

  overlap_df$DMR_Direction <- factor(overlap_df$DMR_Direction,
                                     levels = c("Hypermethylated", "Hypomethylated"))
  overlap_df$ATAC_Direction <- factor(overlap_df$ATAC_Direction,
                                      levels = c("ATAC Up", "ATAC Down"))

  p_12a <- ggplot(overlap_df, aes(x = DMR_Direction, y = Percentage, fill = ATAC_Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = COLORS$atac[c("ATAC Up", "ATAC Down")],
                      name = "ATAC-seq Change") +
    scale_y_continuous(limits = c(0, max(overlap_df$Percentage) * 1.3), expand = c(0, 0)) +
    labs(
      title = "ATAC-seq Peak Overlap at Differentially Methylated Regions",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Prediction: Hyper DMRs \u2192 ATAC Down",
                         fisher_result$estimate, fisher_result$p.value),
      x = "DMR Direction (5mC)", y = "% of DMRs Overlapping ATAC Peaks"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_12a, file.path(OUTPUT_DIR, "12a_atac_overlap"),
                          width = 8, height = 7)

  # Save overlap summary table
  overlap_summary <- data.frame(
    DMR_Direction = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    ATAC_Direction = c("ATAC Up", "ATAC Down", "ATAC Up", "ATAC Down"),
    Overlap_Count = c(hyper_up, hyper_down, hypo_up, hypo_down),
    Total_DMRs = c(n_hyper, n_hyper, n_hypo, n_hypo),
    Percentage = c(100 * hyper_up / n_hyper, 100 * hyper_down / n_hyper,
                   100 * hypo_up / n_hypo, 100 * hypo_down / n_hypo),
    Fisher_OR = fisher_result$estimate,
    Fisher_pvalue = fisher_result$p.value
  )
  write.table(overlap_summary, file.path(TABLES_DIR, "atac_dmr_overlap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: atac_dmr_overlap_summary.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 12b: Consensus Accessibility at DMRs
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 12b: Consensus accessibility at DMRs...\n")

  # Load consensus peaks (skip with warning if not yet generated)
  consensus_ctrl_gr <- NULL
  consensus_mut_gr <- NULL

  if (file.exists(ATAC_FILES$consensus_ctrl) && file.exists(ATAC_FILES$consensus_mut)) {
    consensus_ctrl_gr <- load_chip_peaks(ATAC_FILES$consensus_ctrl, "Consensus Control")
    consensus_mut_gr <- load_chip_peaks(ATAC_FILES$consensus_mut, "Consensus Mutant")
  } else {
    cat("  WARNING: Consensus peak files not found. Run generate_consensus.sh first.\n")
    cat("  Skipping Figure 12b.\n")
  }

  if (!is.null(consensus_ctrl_gr) && !is.null(consensus_mut_gr)) {
    # Compute fraction of hyper/hypo DMRs overlapping each consensus set
    hyper_ctrl_overlap <- sum(countOverlaps(hyper_gr, consensus_ctrl_gr) > 0)
    hyper_mut_overlap <- sum(countOverlaps(hyper_gr, consensus_mut_gr) > 0)
    hypo_ctrl_overlap <- sum(countOverlaps(hypo_gr, consensus_ctrl_gr) > 0)
    hypo_mut_overlap <- sum(countOverlaps(hypo_gr, consensus_mut_gr) > 0)

    cat(sprintf("  Hypermethylated DMRs: %.1f%% ctrl accessible, %.1f%% mut accessible\n",
                100 * hyper_ctrl_overlap / n_hyper, 100 * hyper_mut_overlap / n_hyper))
    cat(sprintf("  Hypomethylated DMRs: %.1f%% ctrl accessible, %.1f%% mut accessible\n",
                100 * hypo_ctrl_overlap / n_hypo, 100 * hypo_mut_overlap / n_hypo))

    consensus_df <- data.frame(
      DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
      Genotype = rep(c("Control", "Mutant"), 2),
      Count = c(hyper_ctrl_overlap, hyper_mut_overlap, hypo_ctrl_overlap, hypo_mut_overlap),
      Total = rep(c(n_hyper, n_hypo), each = 2)
    ) %>%
      mutate(Percentage = 100 * Count / Total)

    consensus_df$DMR_Direction <- factor(consensus_df$DMR_Direction,
                                         levels = c("Hypermethylated", "Hypomethylated"))
    consensus_df$Genotype <- factor(consensus_df$Genotype,
                                     levels = c("Control", "Mutant"))

    p_12b <- ggplot(consensus_df, aes(x = DMR_Direction, y = Percentage, fill = Genotype)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.6, color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
                position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$condition, name = "Consensus Peaks") +
      scale_y_continuous(limits = c(0, max(consensus_df$Percentage) * 1.3), expand = c(0, 0)) +
      labs(
        title = "Consensus Chromatin Accessibility at DMRs",
        subtitle = "Prediction: Hypermethylated DMRs lose accessibility (higher ctrl than mut overlap)",
        x = "DMR Direction (5mC)", y = "% of DMRs in Accessible Regions"
      ) +
      theme_biomodal() +
      theme(legend.position = "top")

    save_multiformat_ggplot(p_12b, file.path(OUTPUT_DIR, "12b_consensus_accessibility"),
                            width = 8, height = 7)
  }

  # -----------------------------------------------------------------------
  # FIGURE 12c: ATAC at Coordinated mC up / hmC down Genes
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 12c: ATAC at coordinated genes...\n")

  # Annotate ATAC peaks to genes
  cat("  Annotating ATAC Up peaks to genes...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  atac_up_anno <- annotate_atac_to_genes(atac_up_gr, txdb)
  cat("  Annotating ATAC Down peaks to genes...\n")
  atac_down_anno <- annotate_atac_to_genes(atac_down_gr, txdb)

  # Aggregate per gene
  atac_up_gene <- aggregate_atac_by_gene(atac_up_anno, "up")
  atac_down_gene <- aggregate_atac_by_gene(atac_down_anno, "down")

  cat(sprintf("  ATAC Up peaks mapped to %d unique genes\n", nrow(atac_up_gene)))
  cat(sprintf("  ATAC Down peaks mapped to %d unique genes\n", nrow(atac_down_gene)))

  # Join ATAC gene data
  atac_gene <- full_join(atac_up_gene, atac_down_gene, by = "SYMBOL") %>%
    mutate(
      n_atac_up = replace_na(n_atac_up, 0),
      n_atac_down = replace_na(n_atac_down, 0),
      net_atac = n_atac_up - n_atac_down
    )

  # Classify coordinated genes
  coord_gene_list <- coordinated$gene[coordinated$coordinated_pattern]
  atac_gene_classified <- atac_gene %>%
    mutate(category = ifelse(SYMBOL %in% coord_gene_list,
                             "Coordinated\n(mC\u2191/hmC\u2193)",
                             "All Other Genes"))

  # Count ATAC_down peaks per gene for coordinated vs all other
  wilcox_coord <- wilcox.test(
    atac_gene_classified$n_atac_down[atac_gene_classified$category == "Coordinated\n(mC\u2191/hmC\u2193)"],
    atac_gene_classified$n_atac_down[atac_gene_classified$category == "All Other Genes"]
  )
  cat(sprintf("  Wilcoxon ATAC Down at coordinated vs other: p = %.2e\n", wilcox_coord$p.value))

  # Boxplot comparison
  p_12c_box <- ggplot(atac_gene_classified, aes(x = category, y = n_atac_down, fill = category)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = c("Coordinated\n(mC\u2191/hmC\u2193)" = "#E6AB02",
                                 "All Other Genes" = "grey70")) +
    annotate("text", x = 1.5, y = max(atac_gene_classified$n_atac_down, na.rm = TRUE) * 0.95,
             label = sprintf("Wilcoxon p = %.2e", wilcox_coord$p.value),
             size = 3.5, fontface = "italic") +
    labs(
      title = "ATAC Down Peaks at Coordinated Genes",
      subtitle = "Genes with mC\u2191 + hmC\u2193 vs all other genes",
      x = "", y = "Number of ATAC Down Peaks per Gene"
    ) +
    theme_biomodal() +
    theme(legend.position = "none")

  # Top coordinated genes bar chart with ATAC annotation
  coordinated_atac <- coordinated %>%
    dplyr::filter(coordinated_pattern) %>%
    left_join(atac_gene %>% dplyr::select(SYMBOL, n_atac_up, n_atac_down, net_atac),
              by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_atac_up = replace_na(n_atac_up, 0),
      n_atac_down = replace_na(n_atac_down, 0),
      net_atac = replace_na(net_atac, 0)
    )

  top_coord_atac <- coordinated_atac %>%
    arrange(desc(combined_effect)) %>%
    head(20) %>%
    mutate(
      gene = factor(gene, levels = rev(unique(gene))),
      atac_label = sprintf("ATAC: %+d", net_atac),
      atac_direction = case_when(
        net_atac < 0 ~ "Net Down",
        net_atac > 0 ~ "Net Up",
        TRUE ~ "No Change"
      )
    )

  p_12c_bar <- ggplot(top_coord_atac, aes(x = gene, y = combined_effect * 100)) +
    geom_bar(aes(fill = atac_direction), stat = "identity", width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_text(aes(label = atac_label), hjust = -0.1, size = 2.8, color = "grey30") +
    scale_fill_manual(values = c("Net Down" = "#66A61E", "Net Up" = "#E6AB02",
                                 "No Change" = "grey70"),
                      name = "Net ATAC Direction") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "Top 20 Coordinated Genes",
      subtitle = "mC\u2191/hmC\u2193 genes with net ATAC change annotation",
      x = "", y = "Combined Effect (|mC| + |hmC| change, %)"
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  p_12c <- p_12c_box + p_12c_bar +
    plot_layout(widths = c(1, 2)) +
    plot_annotation(
      title = "ATAC-seq Signal at Coordinated mC\u2191/hmC\u2193 Genes",
      subtitle = "Do genes with impaired TET-mediated demethylation show decreased accessibility?",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_12c, file.path(OUTPUT_DIR, "12c_atac_coordinated_genes"),
                          width = 15, height = 9)

  # Save coordinated genes with ATAC annotations
  coord_export <- coordinated_atac %>%
    dplyr::select(gene, mc_diff, hmc_diff, mc_q, hmc_q, combined_effect,
                  n_atac_up, n_atac_down, net_atac) %>%
    arrange(desc(combined_effect))

  write.table(coord_export, file.path(TABLES_DIR, "atac_coordinated_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: atac_coordinated_genes.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 12d: Gene-Level Scatter - mC Change vs Net ATAC Change
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 12d: mC change vs net ATAC change scatter...\n")

  # Join ATAC gene-level data to significant DMR data
  mc_dmr_atac <- mc_dmr %>%
    dplyr::filter(significant) %>%
    left_join(atac_gene %>% dplyr::select(SYMBOL, n_atac_up, n_atac_down, net_atac),
              by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_atac_up = replace_na(n_atac_up, 0),
      n_atac_down = replace_na(n_atac_down, 0),
      net_atac = replace_na(net_atac, 0)
    )

  # Filter to genes with at least some ATAC signal for scatter
  scatter_data <- mc_dmr_atac %>%
    dplyr::filter(net_atac != 0) %>%
    mutate(
      mc_pct = mod_difference * 100,
      label_gene = ifelse(gene %in% KEY_GENES, gene, ""),
      dmr_category = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
    )

  cat(sprintf("  %d significant DMR genes with non-zero net ATAC change\n", nrow(scatter_data)))

  # Spearman correlation
  cor_test <- cor.test(scatter_data$mc_pct, scatter_data$net_atac,
                       method = "spearman")
  cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))

  # Count genes per quadrant
  q1 <- sum(scatter_data$mc_pct > 0 & scatter_data$net_atac > 0)  # mC up, ATAC up
  q2 <- sum(scatter_data$mc_pct < 0 & scatter_data$net_atac > 0)  # mC down, ATAC up (predicted)
  q3 <- sum(scatter_data$mc_pct < 0 & scatter_data$net_atac < 0)  # mC down, ATAC down
  q4 <- sum(scatter_data$mc_pct > 0 & scatter_data$net_atac < 0)  # mC up, ATAC down (predicted)

  p_12d <- ggplot(scatter_data, aes(x = mc_pct, y = net_atac)) +
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
    # Quadrant annotations - predicted quadrants in bold color
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = min(scatter_data$net_atac) * 0.9,
             label = sprintf("mC\u2191 ATAC\u2193\nn=%d (predicted)", q4),
             size = 3, color = "#66A61E", fontface = "bold") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = max(scatter_data$net_atac) * 0.9,
             label = sprintf("mC\u2193 ATAC\u2191\nn=%d (predicted)", q2),
             size = 3, color = "#E6AB02", fontface = "bold") +
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = max(scatter_data$net_atac) * 0.9,
             label = sprintf("mC\u2191 ATAC\u2191\nn=%d", q1),
             size = 3, color = "grey50") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = min(scatter_data$net_atac) * 0.9,
             label = sprintf("mC\u2193 ATAC\u2193\nn=%d", q3),
             size = 3, color = "grey50") +
    labs(
      title = "5mC Change vs Net ATAC-seq Change at Significant DMRs",
      subtitle = sprintf("Spearman \u03C1 = %.3f, p = %.2e | Prediction: Negative correlation",
                         cor_test$estimate, cor_test$p.value),
      x = "5mC Change (Mutant - Control, %)",
      y = "Net ATAC Change (n_up - n_down peaks per gene)"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_12d, file.path(OUTPUT_DIR, "12d_mc_vs_atac_scatter"),
                          width = 10, height = 9)

  # Save gene-level ATAC annotation table
  gene_atac_table <- mc_dmr_atac %>%
    dplyr::select(gene, chr, start, end, mod_difference, dmr_qvalue,
                  direction, n_atac_up, n_atac_down, net_atac) %>%
    arrange(dmr_qvalue)

  write.table(gene_atac_table, file.path(TABLES_DIR, "atac_gene_level_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: atac_gene_level_annotation.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 12e: 2x2 Enrichment Heatmap (mC Direction x ATAC Direction)
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 12e: Integration heatmap...\n")

  # Stratify genes into 4 quadrants (exclude net_atac == 0)
  sig_with_atac <- mc_dmr_atac %>%
    dplyr::filter(net_atac != 0) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      atac_direction = ifelse(net_atac > 0, "ATAC Up", "ATAC Down")
    )

  # Contingency table
  quad_table <- table(sig_with_atac$mc_direction, sig_with_atac$atac_direction)
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
  colnames(heatmap_df) <- c("mc_direction", "atac_direction", "Observed")

  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("mc_direction", "atac_direction", "Expected")

  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("mc_direction", "atac_direction", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("mc_direction", "atac_direction")) %>%
    left_join(enrichment_df, by = c("mc_direction", "atac_direction")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nOR: %.2f", Observed, Expected, Enrichment),
      mc_direction = factor(mc_direction, levels = c("mC Up", "mC Down")),
      atac_direction = factor(atac_direction, levels = c("ATAC Up", "ATAC Down")),
      is_predicted = (mc_direction == "mC Up" & atac_direction == "ATAC Down") |
                     (mc_direction == "mC Down" & atac_direction == "ATAC Up")
    )

  p_12e <- ggplot(heatmap_data, aes(x = atac_direction, y = mc_direction, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2) +
    # Highlight predicted quadrants with black borders
    geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
              aes(x = atac_direction, y = mc_direction),
              fill = NA, color = "black", linewidth = 1.5, linetype = "solid") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(
      title = "Integration: 5mC Direction x ATAC-seq Direction",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Black borders = predicted quadrants",
                         fisher_quad$estimate, fisher_quad$p.value),
      x = "ATAC-seq Direction (Mutant vs Control)",
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

  save_multiformat_ggplot(p_12e, file.path(OUTPUT_DIR, "12e_atac_integration_heatmap"),
                          width = 9, height = 7)

  # -----------------------------------------------------------------------
  # Print summary
  # -----------------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 12 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Total ATAC Up peaks: %d\n", length(atac_up_gr)))
  cat(sprintf("Total ATAC Down peaks: %d\n", length(atac_down_gr)))
  cat(sprintf("Genes with ATAC Up peaks: %d\n", nrow(atac_up_gene)))
  cat(sprintf("Genes with ATAC Down peaks: %d\n", nrow(atac_down_gene)))
  cat(sprintf("Significant DMR genes with non-zero net ATAC: %d\n", nrow(scatter_data)))
  cat(sprintf("Coordinated genes (mC\u2191/hmC\u2193) with ATAC data: %d\n",
              sum(coordinated_atac$gene %in% atac_gene$SYMBOL)))
  cat(sprintf("\nOverlap analysis: Fisher OR = %.2f, p = %.2e\n",
              fisher_result$estimate, fisher_result$p.value))
  cat(sprintf("mC vs ATAC correlation: Spearman rho = %.3f, p = %.2e\n",
              cor_test$estimate, cor_test$p.value))
  cat(sprintf("ATAC Down at coordinated vs other: Wilcoxon p = %.2e\n", wilcox_coord$p.value))
  cat(sprintf("Integration: Fisher OR = %.2f, p = %.2e\n",
              fisher_quad$estimate, fisher_quad$p.value))
  cat(sprintf("\nATAC overlap at DMRs:\n"))
  cat(sprintf("  Hypermethylated: %.1f%% ATAC Up, %.1f%% ATAC Down\n",
              100 * hyper_up / n_hyper, 100 * hyper_down / n_hyper))
  cat(sprintf("  Hypomethylated:  %.1f%% ATAC Up, %.1f%% ATAC Down\n",
              100 * hypo_up / n_hypo, 100 * hypo_down / n_hypo))

  if (!is.null(consensus_ctrl_gr) && !is.null(consensus_mut_gr)) {
    cat(sprintf("\nConsensus accessibility at DMRs:\n"))
    cat(sprintf("  Hypermethylated: %.1f%% ctrl, %.1f%% mut\n",
                100 * hyper_ctrl_overlap / n_hyper, 100 * hyper_mut_overlap / n_hyper))
    cat(sprintf("  Hypomethylated:  %.1f%% ctrl, %.1f%% mut\n",
                100 * hypo_ctrl_overlap / n_hypo, 100 * hypo_mut_overlap / n_hypo))
  }
}

cat("\nSection 12 complete.\n\n")
