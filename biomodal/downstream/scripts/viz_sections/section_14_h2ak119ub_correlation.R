# biomodal/downstream/scripts/viz_sections/section_14_h2ak119ub_correlation.R
# Section 14: H2AK119ub Correlation at DMRs
# Standalone script - sources shared config for all dependencies and data
# Tests whether DMR methylation changes correlate with H2AK119ub changes.
#
# Biological prediction (Conway et al. 2021):
#   BAP1 (PR-DUB) removes H2AK119ub. In BAP1-KO, H2AK119ub accumulates,
#   causing chromatin compaction that blocks TET access -> 5mC up / 5hmC down.
#   Therefore:
#     - Hypermethylated DMRs should show INCREASED H2AK119ub (gained in mutant)
#     - Hypomethylated DMRs should show DECREASED H2AK119ub (lost in mutant)
#
# Data: Condition-specific peaks (ctrl and mut BEDs), not pre-computed differential.
#   Gained/lost peaks derived via GenomicRanges set operations.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_14_h2ak119ub_correlation.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# Load ChIPseeker for gene annotation of K119ub peaks (not in shared config)
suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# K119ub-SPECIFIC HELPER FUNCTIONS
# =============================================================================

#' Derive differential peaks from condition-specific peak sets
#' Gained = mut-only (no overlap with ctrl), Lost = ctrl-only (no overlap with mut)
#' @param ctrl_gr GRanges of control peaks
#' @param mut_gr GRanges of mutant peaks
#' @return list with gained, lost, and shared GRanges
derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  # Find overlaps between mut and ctrl
  mut_hits <- countOverlaps(mut_gr, ctrl_gr)
  ctrl_hits <- countOverlaps(ctrl_gr, mut_gr)

  gained <- mut_gr[mut_hits == 0]   # Mut peaks with no ctrl overlap
  lost <- ctrl_gr[ctrl_hits == 0]    # Ctrl peaks with no mut overlap
  shared_mut <- mut_gr[mut_hits > 0] # Mut peaks overlapping ctrl

  list(
    gained = gained,
    lost = lost,
    shared = shared_mut
  )
}

#' Annotate K119ub peaks to nearest genes using ChIPseeker
#' @param k119ub_gr GRanges of K119ub peaks
#' @param txdb TxDb object for gene annotation
#' @return data.frame with peak-to-gene mapping
annotate_k119ub_to_genes <- function(k119ub_gr, txdb) {
  anno <- suppressMessages(annotatePeak(
    k119ub_gr,
    TxDb = txdb,
    annoDb = "org.Mm.eg.db",
    level = "gene"
  ))
  as.data.frame(anno)
}

#' Aggregate K119ub peaks per gene
#' @param anno_df Annotated K119ub peak data.frame (from ChIPseeker)
#' @param condition Character label for condition ("ctrl", "mut", "gained", "lost")
#' @return data.frame with gene-level K119ub peak counts
aggregate_k119ub_by_gene <- function(anno_df, condition) {
  anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      !!paste0("n_k119ub_", condition) := n(),
      !!paste0("min_dist_tss_", condition) := min(abs(distanceToTSS), na.rm = TRUE),
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
# SECTION 14: H2AK119ub CORRELATION AT DMRs
# =============================================================================

cat("================================================================================\n")
cat("SECTION 14: H2AK119ub CORRELATION AT DMRs\n")
cat("================================================================================\n\n")

# Load K119ub condition-specific peaks
cat("Loading H2AK119ub condition-specific peaks...\n")
k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Control")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mutant")

if (is.null(k119ub_ctrl_gr) || is.null(k119ub_mut_gr)) {
  cat("  ERROR: K119ub peak files not found. Skipping Section 14.\n")
} else {

  # Derive gained/lost peaks
  cat("\nDeriving differential K119ub peaks...\n")
  diff_peaks <- derive_differential_peaks(k119ub_ctrl_gr, k119ub_mut_gr)

  cat(sprintf("  Gained (mut-only): %d peaks\n", length(diff_peaks$gained)))
  cat(sprintf("  Lost (ctrl-only):  %d peaks\n", length(diff_peaks$lost)))
  cat(sprintf("  Shared:            %d peaks\n", length(diff_peaks$shared)))

  # Get significant mC DMRs split by direction
  mc_hyper <- mc_dmr %>% dplyr::filter(significant & mod_difference > 0)
  mc_hypo <- mc_dmr %>% dplyr::filter(significant & mod_difference < 0)

  # Convert to GRanges
  hyper_gr <- dmr_to_granges(mc_hyper)
  hypo_gr <- dmr_to_granges(mc_hypo)
  n_hyper <- length(hyper_gr)
  n_hypo <- length(hypo_gr)

  # -----------------------------------------------------------------------
  # FIGURE 14a: Condition-Specific H2AK119ub Overlap at DMRs
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 14a: Condition-specific K119ub overlap at DMRs...\n")

  # Compute overlaps: % of hyper/hypo DMRs overlapping K119ub in ctrl vs mut
  hyper_ctrl <- sum(countOverlaps(hyper_gr, k119ub_ctrl_gr) > 0)
  hyper_mut <- sum(countOverlaps(hyper_gr, k119ub_mut_gr) > 0)
  hypo_ctrl <- sum(countOverlaps(hypo_gr, k119ub_ctrl_gr) > 0)
  hypo_mut <- sum(countOverlaps(hypo_gr, k119ub_mut_gr) > 0)

  cat(sprintf("  Hypermethylated DMRs: %.1f%% ctrl K119ub, %.1f%% mut K119ub\n",
              100 * hyper_ctrl / n_hyper, 100 * hyper_mut / n_hyper))
  cat(sprintf("  Hypomethylated DMRs: %.1f%% ctrl K119ub, %.1f%% mut K119ub\n",
              100 * hypo_ctrl / n_hypo, 100 * hypo_mut / n_hypo))

  # Fisher's exact test: DMR direction x condition overlap
  fisher_14a_mat <- matrix(c(hyper_ctrl, hyper_mut, hypo_ctrl, hypo_mut), nrow = 2, byrow = TRUE,
                           dimnames = list(c("Hypermethylated", "Hypomethylated"),
                                           c("Control K119ub", "Mutant K119ub")))
  fisher_14a <- fisher.test(fisher_14a_mat)

  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_14a$estimate, fisher_14a$p.value))

  # Build plot data
  overlap_14a_df <- data.frame(
    DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
    Genotype = rep(c("Control", "Mutant"), 2),
    Count = c(hyper_ctrl, hyper_mut, hypo_ctrl, hypo_mut),
    Total = rep(c(n_hyper, n_hypo), each = 2)
  ) %>%
    mutate(Percentage = 100 * Count / Total)

  overlap_14a_df$DMR_Direction <- factor(overlap_14a_df$DMR_Direction,
                                          levels = c("Hypermethylated", "Hypomethylated"))
  overlap_14a_df$Genotype <- factor(overlap_14a_df$Genotype,
                                     levels = c("Control", "Mutant"))

  p_14a <- ggplot(overlap_14a_df, aes(x = DMR_Direction, y = Percentage, fill = Genotype)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = COLORS$condition, name = "K119ub Peaks") +
    scale_y_continuous(limits = c(0, max(overlap_14a_df$Percentage) * 1.3), expand = c(0, 0)) +
    labs(
      title = "H2AK119ub Peak Overlap at Differentially Methylated Regions",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Prediction: Hyper DMRs \u2192 higher mut overlap",
                         fisher_14a$estimate, fisher_14a$p.value),
      x = "DMR Direction (5mC)", y = "% of DMRs Overlapping K119ub Peaks"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_14a, file.path(OUTPUT_DIR, "14a_k119ub_condition_overlap"),
                          width = 8, height = 7)

  # Save condition overlap summary table
  overlap_14a_summary <- data.frame(
    DMR_Direction = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    Condition = c("Control", "Mutant", "Control", "Mutant"),
    Overlap_Count = c(hyper_ctrl, hyper_mut, hypo_ctrl, hypo_mut),
    Total_DMRs = c(n_hyper, n_hyper, n_hypo, n_hypo),
    Percentage = c(100 * hyper_ctrl / n_hyper, 100 * hyper_mut / n_hyper,
                   100 * hypo_ctrl / n_hypo, 100 * hypo_mut / n_hypo),
    Fisher_OR = fisher_14a$estimate,
    Fisher_pvalue = fisher_14a$p.value
  )
  write.table(overlap_14a_summary, file.path(TABLES_DIR, "k119ub_condition_dmr_overlap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: k119ub_condition_dmr_overlap_summary.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 14b: Gained/Lost H2AK119ub Overlap at DMRs
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 14b: Gained/lost K119ub overlap at DMRs...\n")

  if (length(diff_peaks$gained) < 10 || length(diff_peaks$lost) < 10) {
    cat("  WARNING: Very few gained or lost peaks (<10). Results may be unreliable.\n")
  }

  # Compute overlaps with gained/lost K119ub peaks
  hyper_gained <- sum(countOverlaps(hyper_gr, diff_peaks$gained) > 0)
  hyper_lost <- sum(countOverlaps(hyper_gr, diff_peaks$lost) > 0)
  hypo_gained <- sum(countOverlaps(hypo_gr, diff_peaks$gained) > 0)
  hypo_lost <- sum(countOverlaps(hypo_gr, diff_peaks$lost) > 0)

  cat(sprintf("  Hypermethylated DMRs: %d/%d (%.1f%%) K119ub Gained, %d/%d (%.1f%%) K119ub Lost\n",
              hyper_gained, n_hyper, 100 * hyper_gained / n_hyper,
              hyper_lost, n_hyper, 100 * hyper_lost / n_hyper))
  cat(sprintf("  Hypomethylated DMRs: %d/%d (%.1f%%) K119ub Gained, %d/%d (%.1f%%) K119ub Lost\n",
              hypo_gained, n_hypo, 100 * hypo_gained / n_hypo,
              hypo_lost, n_hypo, 100 * hypo_lost / n_hypo))

  # Fisher's exact test: DMR direction x K119ub differential
  fisher_14b_mat <- matrix(c(hyper_gained, hyper_lost, hypo_gained, hypo_lost), nrow = 2, byrow = TRUE,
                           dimnames = list(c("Hypermethylated", "Hypomethylated"),
                                           c("K119ub Gained", "K119ub Lost")))
  fisher_14b <- fisher.test(fisher_14b_mat)

  cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
              fisher_14b$estimate, fisher_14b$p.value))

  # Build plot data
  overlap_14b_df <- data.frame(
    DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
    K119ub_Direction = rep(c("K119ub Gained", "K119ub Lost"), 2),
    Count = c(hyper_gained, hyper_lost, hypo_gained, hypo_lost),
    Total = rep(c(n_hyper, n_hypo), each = 2)
  ) %>%
    mutate(Percentage = 100 * Count / Total)

  overlap_14b_df$DMR_Direction <- factor(overlap_14b_df$DMR_Direction,
                                          levels = c("Hypermethylated", "Hypomethylated"))
  overlap_14b_df$K119ub_Direction <- factor(overlap_14b_df$K119ub_Direction,
                                             levels = c("K119ub Gained", "K119ub Lost"))

  p_14b <- ggplot(overlap_14b_df, aes(x = DMR_Direction, y = Percentage, fill = K119ub_Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = COLORS$k119ub[c("K119ub Gained", "K119ub Lost")],
                      name = "K119ub Change") +
    scale_y_continuous(limits = c(0, max(overlap_14b_df$Percentage) * 1.3), expand = c(0, 0)) +
    labs(
      title = "Differential H2AK119ub Peak Overlap at DMRs",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Prediction: Hyper DMRs \u2192 K119ub Gained",
                         fisher_14b$estimate, fisher_14b$p.value),
      x = "DMR Direction (5mC)", y = "% of DMRs Overlapping Differential K119ub Peaks"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_14b, file.path(OUTPUT_DIR, "14b_k119ub_differential_overlap"),
                          width = 8, height = 7)

  # Save differential overlap summary table
  overlap_14b_summary <- data.frame(
    DMR_Direction = c("Hypermethylated", "Hypermethylated", "Hypomethylated", "Hypomethylated"),
    K119ub_Direction = c("K119ub Gained", "K119ub Lost", "K119ub Gained", "K119ub Lost"),
    Overlap_Count = c(hyper_gained, hyper_lost, hypo_gained, hypo_lost),
    Total_DMRs = c(n_hyper, n_hyper, n_hypo, n_hypo),
    Percentage = c(100 * hyper_gained / n_hyper, 100 * hyper_lost / n_hyper,
                   100 * hypo_gained / n_hypo, 100 * hypo_lost / n_hypo),
    N_Gained_Peaks = length(diff_peaks$gained),
    N_Lost_Peaks = length(diff_peaks$lost),
    Fisher_OR = fisher_14b$estimate,
    Fisher_pvalue = fisher_14b$p.value
  )
  write.table(overlap_14b_summary, file.path(TABLES_DIR, "k119ub_differential_dmr_overlap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: k119ub_differential_dmr_overlap_summary.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 14c: K119ub at Coordinated mC up / hmC down Genes
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 14c: K119ub at coordinated genes...\n")

  # Annotate condition-specific peaks to genes
  cat("  Annotating K119ub Control peaks to genes...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  k119ub_ctrl_anno <- annotate_k119ub_to_genes(k119ub_ctrl_gr, txdb)
  cat("  Annotating K119ub Mutant peaks to genes...\n")
  k119ub_mut_anno <- annotate_k119ub_to_genes(k119ub_mut_gr, txdb)

  # Aggregate per gene
  k119ub_ctrl_gene <- aggregate_k119ub_by_gene(k119ub_ctrl_anno, "ctrl")
  k119ub_mut_gene <- aggregate_k119ub_by_gene(k119ub_mut_anno, "mut")

  cat(sprintf("  K119ub Control peaks mapped to %d unique genes\n", nrow(k119ub_ctrl_gene)))
  cat(sprintf("  K119ub Mutant peaks mapped to %d unique genes\n", nrow(k119ub_mut_gene)))

  # Join K119ub gene data
  k119ub_gene <- full_join(k119ub_ctrl_gene, k119ub_mut_gene, by = "SYMBOL") %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = n_k119ub_mut - n_k119ub_ctrl
    )

  # Classify coordinated genes
  coord_gene_list <- coordinated$gene[coordinated$coordinated_pattern]
  k119ub_gene_classified <- k119ub_gene %>%
    mutate(category = ifelse(SYMBOL %in% coord_gene_list,
                             "Coordinated\n(mC\u2191/hmC\u2193)",
                             "All Other Genes"))

  # Wilcoxon: mut peak count at coordinated vs all other genes
  wilcox_coord <- wilcox.test(
    k119ub_gene_classified$n_k119ub_mut[k119ub_gene_classified$category == "Coordinated\n(mC\u2191/hmC\u2193)"],
    k119ub_gene_classified$n_k119ub_mut[k119ub_gene_classified$category == "All Other Genes"]
  )
  cat(sprintf("  Wilcoxon K119ub Mut peaks at coordinated vs other: p = %.2e\n", wilcox_coord$p.value))

  # Violin + boxplot comparison
  p_14c_box <- ggplot(k119ub_gene_classified, aes(x = category, y = n_k119ub_mut, fill = category)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = c("Coordinated\n(mC\u2191/hmC\u2193)" = "#756BB1",
                                 "All Other Genes" = "grey70")) +
    annotate("text", x = 1.5, y = max(k119ub_gene_classified$n_k119ub_mut, na.rm = TRUE) * 0.95,
             label = sprintf("Wilcoxon p = %.2e", wilcox_coord$p.value),
             size = 3.5, fontface = "italic") +
    labs(
      title = "K119ub Mutant Peaks at Coordinated Genes",
      subtitle = "Genes with mC\u2191 + hmC\u2193 vs all other genes",
      x = "", y = "Number of K119ub Mutant Peaks per Gene"
    ) +
    theme_biomodal() +
    theme(legend.position = "none")

  # Top coordinated genes bar chart with K119ub annotation
  coordinated_k119ub <- coordinated %>%
    dplyr::filter(coordinated_pattern) %>%
    left_join(k119ub_gene %>% dplyr::select(SYMBOL, n_k119ub_ctrl, n_k119ub_mut, net_k119ub),
              by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = replace_na(net_k119ub, 0)
    )

  top_coord_k119ub <- coordinated_k119ub %>%
    arrange(desc(combined_effect)) %>%
    head(20) %>%
    mutate(
      gene = factor(gene, levels = rev(unique(gene))),
      k119ub_label = sprintf("K119ub: %+d", net_k119ub),
      k119ub_direction = case_when(
        net_k119ub > 0 ~ "Net Gained",
        net_k119ub < 0 ~ "Net Lost",
        TRUE ~ "No Change"
      )
    )

  p_14c_bar <- ggplot(top_coord_k119ub, aes(x = gene, y = combined_effect * 100)) +
    geom_bar(aes(fill = k119ub_direction), stat = "identity", width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_text(aes(label = k119ub_label), hjust = -0.1, size = 2.8, color = "grey30") +
    scale_fill_manual(values = c("Net Gained" = "#756BB1", "Net Lost" = "#74C476",
                                 "No Change" = "grey70"),
                      name = "Net K119ub Direction") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      title = "Top 20 Coordinated Genes",
      subtitle = "mC\u2191/hmC\u2193 genes with net K119ub change annotation",
      x = "", y = "Combined Effect (|mC| + |hmC| change, %)"
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  p_14c <- p_14c_box + p_14c_bar +
    plot_layout(widths = c(1, 2)) +
    plot_annotation(
      title = "H2AK119ub Signal at Coordinated mC\u2191/hmC\u2193 Genes",
      subtitle = "Do genes with impaired TET-mediated demethylation show increased H2AK119ub?",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_14c, file.path(OUTPUT_DIR, "14c_k119ub_coordinated_genes"),
                          width = 15, height = 9)

  # Save coordinated genes with K119ub annotations
  coord_export <- coordinated_k119ub %>%
    dplyr::select(gene, mc_diff, hmc_diff, mc_q, hmc_q, combined_effect,
                  n_k119ub_ctrl, n_k119ub_mut, net_k119ub) %>%
    arrange(desc(combined_effect))

  write.table(coord_export, file.path(TABLES_DIR, "k119ub_coordinated_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: k119ub_coordinated_genes.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 14d: Scatter - mC Change vs Net K119ub Change
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 14d: mC change vs net K119ub change scatter...\n")

  # Join K119ub gene-level data to significant DMR data
  mc_dmr_k119ub <- mc_dmr %>%
    dplyr::filter(significant) %>%
    left_join(k119ub_gene %>% dplyr::select(SYMBOL, n_k119ub_ctrl, n_k119ub_mut, net_k119ub),
              by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = replace_na(net_k119ub, 0)
    )

  # Filter to genes with non-zero net K119ub for scatter
  scatter_data <- mc_dmr_k119ub %>%
    dplyr::filter(net_k119ub != 0) %>%
    mutate(
      mc_pct = mod_difference * 100,
      label_gene = ifelse(gene %in% KEY_GENES, gene, ""),
      dmr_category = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
    )

  cat(sprintf("  %d significant DMR genes with non-zero net K119ub change\n", nrow(scatter_data)))

  # Spearman correlation
  cor_test <- cor.test(scatter_data$mc_pct, scatter_data$net_k119ub,
                       method = "spearman")
  cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))

  # Count genes per quadrant
  q1 <- sum(scatter_data$mc_pct > 0 & scatter_data$net_k119ub > 0)  # mC up, K119ub gained (predicted)
  q2 <- sum(scatter_data$mc_pct < 0 & scatter_data$net_k119ub > 0)  # mC down, K119ub gained
  q3 <- sum(scatter_data$mc_pct < 0 & scatter_data$net_k119ub < 0)  # mC down, K119ub lost (predicted)
  q4 <- sum(scatter_data$mc_pct > 0 & scatter_data$net_k119ub < 0)  # mC up, K119ub lost

  p_14d <- ggplot(scatter_data, aes(x = mc_pct, y = net_k119ub)) +
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
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = max(scatter_data$net_k119ub) * 0.9,
             label = sprintf("mC\u2191 K119ub\u2191\nn=%d (predicted)", q1),
             size = 3, color = "#756BB1", fontface = "bold") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = min(scatter_data$net_k119ub) * 0.9,
             label = sprintf("mC\u2193 K119ub\u2193\nn=%d (predicted)", q3),
             size = 3, color = "#74C476", fontface = "bold") +
    annotate("text", x = max(scatter_data$mc_pct) * 0.7, y = min(scatter_data$net_k119ub) * 0.9,
             label = sprintf("mC\u2191 K119ub\u2193\nn=%d", q4),
             size = 3, color = "grey50") +
    annotate("text", x = min(scatter_data$mc_pct) * 0.7, y = max(scatter_data$net_k119ub) * 0.9,
             label = sprintf("mC\u2193 K119ub\u2191\nn=%d", q2),
             size = 3, color = "grey50") +
    labs(
      title = "5mC Change vs Net H2AK119ub Change at Significant DMRs",
      subtitle = sprintf("Spearman \u03C1 = %.3f, p = %.2e | Prediction: Positive correlation",
                         cor_test$estimate, cor_test$p.value),
      x = "5mC Change (Mutant - Control, %)",
      y = "Net K119ub Change (n_mut - n_ctrl peaks per gene)"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_14d, file.path(OUTPUT_DIR, "14d_mc_vs_k119ub_scatter"),
                          width = 10, height = 9)

  # Save gene-level K119ub annotation table
  gene_k119ub_table <- mc_dmr_k119ub %>%
    dplyr::select(gene, chr, start, end, mod_difference, dmr_qvalue,
                  direction, n_k119ub_ctrl, n_k119ub_mut, net_k119ub) %>%
    arrange(dmr_qvalue)

  write.table(gene_k119ub_table, file.path(TABLES_DIR, "k119ub_gene_level_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: k119ub_gene_level_annotation.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 14e: 2x2 Enrichment Heatmap (mC Direction x K119ub Direction)
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 14e: Integration heatmap...\n")

  # Stratify genes into 4 quadrants (exclude net_k119ub == 0)
  sig_with_k119ub <- mc_dmr_k119ub %>%
    dplyr::filter(net_k119ub != 0) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      k119ub_direction = ifelse(net_k119ub > 0, "K119ub Gained", "K119ub Lost")
    )

  # Contingency table
  quad_table <- table(sig_with_k119ub$mc_direction, sig_with_k119ub$k119ub_direction)
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
  colnames(heatmap_df) <- c("mc_direction", "k119ub_direction", "Observed")

  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("mc_direction", "k119ub_direction", "Expected")

  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("mc_direction", "k119ub_direction", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("mc_direction", "k119ub_direction")) %>%
    left_join(enrichment_df, by = c("mc_direction", "k119ub_direction")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nOR: %.2f", Observed, Expected, Enrichment),
      mc_direction = factor(mc_direction, levels = c("mC Up", "mC Down")),
      k119ub_direction = factor(k119ub_direction, levels = c("K119ub Gained", "K119ub Lost")),
      is_predicted = (mc_direction == "mC Up" & k119ub_direction == "K119ub Gained") |
                     (mc_direction == "mC Down" & k119ub_direction == "K119ub Lost")
    )

  p_14e <- ggplot(heatmap_data, aes(x = k119ub_direction, y = mc_direction, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2) +
    # Highlight predicted quadrants with black borders
    geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
              aes(x = k119ub_direction, y = mc_direction),
              fill = NA, color = "black", linewidth = 1.5, linetype = "solid") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(
      title = "Integration: 5mC Direction x H2AK119ub Direction",
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Black borders = predicted quadrants",
                         fisher_quad$estimate, fisher_quad$p.value),
      x = "H2AK119ub Direction (Mutant vs Control)",
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

  save_multiformat_ggplot(p_14e, file.path(OUTPUT_DIR, "14e_k119ub_integration_heatmap"),
                          width = 9, height = 7)

  # -----------------------------------------------------------------------
  # Print summary
  # -----------------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 14 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("Total K119ub Control peaks: %d\n", length(k119ub_ctrl_gr)))
  cat(sprintf("Total K119ub Mutant peaks: %d\n", length(k119ub_mut_gr)))
  cat(sprintf("Derived differential peaks:\n"))
  cat(sprintf("  Gained (mut-only): %d\n", length(diff_peaks$gained)))
  cat(sprintf("  Lost (ctrl-only):  %d\n", length(diff_peaks$lost)))
  cat(sprintf("  Shared:            %d\n", length(diff_peaks$shared)))
  cat(sprintf("\nGenes with K119ub Control peaks: %d\n", nrow(k119ub_ctrl_gene)))
  cat(sprintf("Genes with K119ub Mutant peaks: %d\n", nrow(k119ub_mut_gene)))
  cat(sprintf("Significant DMR genes with non-zero net K119ub: %d\n", nrow(scatter_data)))
  cat(sprintf("Coordinated genes (mC\u2191/hmC\u2193) with K119ub data: %d\n",
              sum(coordinated_k119ub$gene %in% k119ub_gene$SYMBOL)))
  cat(sprintf("\nCondition overlap (14a): Fisher OR = %.2f, p = %.2e\n",
              fisher_14a$estimate, fisher_14a$p.value))
  cat(sprintf("Differential overlap (14b): Fisher OR = %.2f, p = %.2e\n",
              fisher_14b$estimate, fisher_14b$p.value))
  cat(sprintf("mC vs K119ub correlation (14d): Spearman rho = %.3f, p = %.2e\n",
              cor_test$estimate, cor_test$p.value))
  cat(sprintf("K119ub Mut at coordinated vs other (14c): Wilcoxon p = %.2e\n", wilcox_coord$p.value))
  cat(sprintf("Integration (14e): Fisher OR = %.2f, p = %.2e\n",
              fisher_quad$estimate, fisher_quad$p.value))
  cat(sprintf("\nCondition-specific K119ub overlap at DMRs:\n"))
  cat(sprintf("  Hypermethylated: %.1f%% ctrl, %.1f%% mut\n",
              100 * hyper_ctrl / n_hyper, 100 * hyper_mut / n_hyper))
  cat(sprintf("  Hypomethylated:  %.1f%% ctrl, %.1f%% mut\n",
              100 * hypo_ctrl / n_hypo, 100 * hypo_mut / n_hypo))
  cat(sprintf("\nDifferential K119ub overlap at DMRs:\n"))
  cat(sprintf("  Hypermethylated: %.1f%% gained, %.1f%% lost\n",
              100 * hyper_gained / n_hyper, 100 * hyper_lost / n_hyper))
  cat(sprintf("  Hypomethylated:  %.1f%% gained, %.1f%% lost\n",
              100 * hypo_gained / n_hypo, 100 * hypo_lost / n_hypo))
}

cat("\nSection 14 complete.\n\n")
