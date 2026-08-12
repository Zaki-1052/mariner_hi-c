# biomodal/downstream/scripts/viz_sections/section_16_raw_concordance_barcharts.R
# Section 16: Raw Concordance Bar Charts for Predicted Chromatin-Methylation Trends
# Companion to Section 15 (O/E heatmaps) — presents the same data as direct
# concordance rates (%) so the actual effect sizes in the dominant groups
# (mC Up, hmC Down) are visible rather than obscured by O/E normalization.
#
# For each mark, asks: "What % of genes in each methylation group show the
# predicted chromatin mark direction?"
#
# Predictions:
#   mC Up    → MeCP2 Up, ATAC Down, K119ub Gained
#   mC Down  → MeCP2 Down, ATAC Up, K119ub Lost
#   hmC Down → MeCP2 Up, ATAC Down, K119ub Gained
#   hmC Up   → MeCP2 Down, ATAC Up, K119ub Lost
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_16_raw_concordance_barcharts.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# HELPER FUNCTIONS (same as section 15)
# =============================================================================

load_mecp2_annotated <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, quote = "",
                   comment.char = "")
  for (col in c("Fold", "FDR", "p.value", "distanceToTSS")) {
    if (col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
  }
  df
}

aggregate_mecp2_by_gene <- function(mecp2_df) {
  mecp2_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      n_peaks = n(),
      .groups = "drop"
    )
}

annotate_peaks_to_genes <- function(gr, txdb) {
  anno <- suppressMessages(annotatePeak(
    gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
  ))
  as.data.frame(anno)
}

aggregate_peaks_by_gene <- function(anno_df, label) {
  anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(!!paste0("n_", label) := n(), .groups = "drop")
}

derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  list(
    gained = mut_gr[countOverlaps(mut_gr, ctrl_gr) == 0],
    lost = ctrl_gr[countOverlaps(ctrl_gr, mut_gr) == 0]
  )
}

# =============================================================================
# SECTION 16: RAW CONCORDANCE BAR CHARTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 16: RAW CONCORDANCE BAR CHARTS\n")
cat("================================================================================\n\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# -----------------------------------------------------------------------
# Prepare methylation gene groups
# -----------------------------------------------------------------------
cat("Preparing methylation gene groups...\n")

mc_sig <- mc_dmr %>% dplyr::filter(significant)
hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

mc_up_genes <- mc_sig$gene[mc_sig$mod_difference > 0]
mc_down_genes <- mc_sig$gene[mc_sig$mod_difference < 0]
hmc_down_genes <- hmc_sig$gene[hmc_sig$mod_difference < 0]
hmc_up_genes <- hmc_sig$gene[hmc_sig$mod_difference > 0]

cat(sprintf("  mC Up:    %d genes\n", length(mc_up_genes)))
cat(sprintf("  mC Down:  %d genes\n", length(mc_down_genes)))
cat(sprintf("  hmC Down: %d genes\n", length(hmc_down_genes)))
cat(sprintf("  hmC Up:   %d genes\n", length(hmc_up_genes)))

# -----------------------------------------------------------------------
# Compute concordance rate for a given mark
# For each methylation group, what % of genes show the predicted mark direction?
# Returns a data.frame with one row per methylation group
# -----------------------------------------------------------------------
compute_concordance <- function(mark_gene_df, mark_col, mark_name,
                                predicted_positive, predicted_negative) {
  # mark_gene_df: data.frame with SYMBOL and mark_col columns
  # predicted_positive: direction label when mark value > 0 (e.g., "MeCP2 Up")
  # predicted_negative: direction label when mark value < 0 (e.g., "MeCP2 Down")

  mark_vals <- setNames(mark_gene_df[[mark_col]], mark_gene_df$SYMBOL)

  compute_group <- function(gene_list, group_name, predicted_dir) {
    matched <- gene_list[gene_list %in% names(mark_vals)]
    vals <- mark_vals[matched]

    if (predicted_dir == predicted_positive) {
      concordant <- sum(vals > 0)
    } else {
      concordant <- sum(vals < 0)
    }

    data.frame(
      group = group_name,
      mark = mark_name,
      predicted_direction = predicted_dir,
      n_concordant = concordant,
      n_total = length(vals),
      pct_concordant = 100 * concordant / length(vals),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    compute_group(mc_up_genes, "mC Up", predicted_positive),
    compute_group(mc_down_genes, "mC Down", predicted_negative),
    compute_group(hmc_down_genes, "hmC Down", predicted_positive),
    compute_group(hmc_up_genes, "hmC Up", predicted_negative)
  )
}

# -----------------------------------------------------------------------
# Build a per-mark concordance bar chart
# -----------------------------------------------------------------------
build_concordance_barchart <- function(conc_df, mark_name, mark_subtitle) {
  # Order: mC Up, mC Down, hmC Down, hmC Up
  conc_df$group <- factor(conc_df$group,
                          levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))

  # Color by methylation type
  conc_df$met_type <- ifelse(grepl("^mC ", conc_df$group), "5mC", "5hmC")

  # Fisher test: predicted-up group vs predicted-down group within each perspective
  mc_up_row <- conc_df[conc_df$group == "mC Up", ]
  mc_down_row <- conc_df[conc_df$group == "mC Down", ]
  hmc_down_row <- conc_df[conc_df$group == "hmC Down", ]
  hmc_up_row <- conc_df[conc_df$group == "hmC Up", ]

  # mC perspective: mC Up concordant rate vs mC Down concordant rate
  fisher_mc <- fisher.test(matrix(c(
    mc_up_row$n_concordant, mc_up_row$n_total - mc_up_row$n_concordant,
    mc_down_row$n_concordant, mc_down_row$n_total - mc_down_row$n_concordant
  ), nrow = 2, byrow = TRUE))

  # hmC perspective: hmC Down concordant rate vs hmC Up concordant rate
  fisher_hmc <- fisher.test(matrix(c(
    hmc_down_row$n_concordant, hmc_down_row$n_total - hmc_down_row$n_concordant,
    hmc_up_row$n_concordant, hmc_up_row$n_total - hmc_up_row$n_concordant
  ), nrow = 2, byrow = TRUE))

  p <- ggplot(conc_df, aes(x = group, y = pct_concordant, fill = met_type)) +
    geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)",
                                  pct_concordant, n_concordant, n_total)),
              vjust = -0.3, size = 3.2, lineheight = 0.9) +
    geom_text(aes(label = predicted_direction, y = 1),
              size = 2.5, color = "grey30", fontface = "italic") +
    # Significance brackets
    annotate("segment", x = 1, xend = 2,
             y = max(conc_df$pct_concordant) * 1.15,
             yend = max(conc_df$pct_concordant) * 1.15,
             linewidth = 0.4) +
    annotate("text", x = 1.5,
             y = max(conc_df$pct_concordant) * 1.18,
             label = sprintf("p = %.1e", fisher_mc$p.value),
             size = 3, fontface = "italic") +
    annotate("segment", x = 3, xend = 4,
             y = max(conc_df$pct_concordant) * 1.15,
             yend = max(conc_df$pct_concordant) * 1.15,
             linewidth = 0.4) +
    annotate("text", x = 3.5,
             y = max(conc_df$pct_concordant) * 1.18,
             label = sprintf("p = %.1e", fisher_hmc$p.value),
             size = 3, fontface = "italic") +
    scale_fill_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
                      name = "Methylation Type") +
    scale_y_continuous(limits = c(0, max(conc_df$pct_concordant) * 1.35),
                       expand = c(0, 0)) +
    labs(
      title = sprintf("%s: %% of Genes with Predicted Direction", mark_name),
      subtitle = mark_subtitle,
      x = "Methylation Group (Significant DMR Genes)",
      y = "% Showing Predicted Mark Direction"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  list(plot = p, data = conc_df, fisher_mc = fisher_mc, fisher_hmc = fisher_hmc)
}

# =============================================================================
# FIGURE 16a: MeCP2 Concordance
# =============================================================================
cat("\n--- FIGURE 16a: MeCP2 Concordance ---\n\n")

mecp2_annotated <- load_mecp2_annotated(MECP2_FILES$annotated)

if (is.null(mecp2_annotated)) {
  cat("  ERROR: MeCP2 data not found. Skipping 16a.\n")
  result_16a <- NULL
} else {
  mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated)
  cat(sprintf("  %d genes with MeCP2 data\n", nrow(mecp2_gene)))

  conc_mecp2 <- compute_concordance(mecp2_gene, "nearest_fold", "MeCP2",
                                    "MeCP2 Up", "MeCP2 Down")

  result_16a <- build_concordance_barchart(
    conc_mecp2, "MeCP2",
    "mC Up / hmC Down genes predicted to show increased MeCP2 binding"
  )

  save_multiformat_ggplot(result_16a$plot,
                          file.path(OUTPUT_DIR, "16a_mecp2_concordance_bars"),
                          width = 9, height = 7)

  cat("  Results:\n")
  for (i in 1:nrow(conc_mecp2)) {
    cat(sprintf("    %s: %.1f%% %s (%d/%d)\n",
                conc_mecp2$group[i], conc_mecp2$pct_concordant[i],
                conc_mecp2$predicted_direction[i],
                conc_mecp2$n_concordant[i], conc_mecp2$n_total[i]))
  }
  cat(sprintf("  Fisher mC Up vs Down: p = %.2e\n", result_16a$fisher_mc$p.value))
  cat(sprintf("  Fisher hmC Down vs Up: p = %.2e\n", result_16a$fisher_hmc$p.value))
}

# =============================================================================
# FIGURE 16b: ATAC Concordance
# =============================================================================
cat("\n--- FIGURE 16b: ATAC Concordance ---\n\n")

atac_up_gr <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

if (is.null(atac_up_gr) || is.null(atac_down_gr)) {
  cat("  ERROR: ATAC data not found. Skipping 16b.\n")
  result_16b <- NULL
} else {
  cat("  Annotating ATAC peaks to genes...\n")
  atac_up_gene <- aggregate_peaks_by_gene(
    annotate_peaks_to_genes(atac_up_gr, txdb), "atac_up")
  atac_down_gene <- aggregate_peaks_by_gene(
    annotate_peaks_to_genes(atac_down_gr, txdb), "atac_down")

  atac_gene <- full_join(atac_up_gene, atac_down_gene, by = "SYMBOL") %>%
    mutate(
      n_atac_up = replace_na(n_atac_up, 0),
      n_atac_down = replace_na(n_atac_down, 0),
      net_atac = n_atac_up - n_atac_down
    ) %>%
    dplyr::filter(net_atac != 0)

  cat(sprintf("  %d genes with non-zero net ATAC\n", nrow(atac_gene)))

  # For ATAC: predicted direction for mC Up is ATAC Down (negative net_atac)
  conc_atac <- compute_concordance(atac_gene, "net_atac", "ATAC",
                                   "ATAC Down", "ATAC Up")
  # Note: predicted_positive="ATAC Down" means when net_atac > 0 that's NOT concordant
  # We need to flip: for mC Up and hmC Down, concordant = net_atac < 0
  # For mC Down and hmC Up, concordant = net_atac > 0
  # The compute_concordance function checks vals > 0 for predicted_positive
  # Since ATAC Down means net_atac < 0, we need to negate

  # Recompute manually for ATAC since the sign convention is inverted
  atac_vals <- setNames(atac_gene$net_atac, atac_gene$SYMBOL)

  compute_atac_group <- function(gene_list, group_name, want_negative) {
    matched <- gene_list[gene_list %in% names(atac_vals)]
    vals <- atac_vals[matched]
    if (want_negative) {
      concordant <- sum(vals < 0)
      pred_dir <- "ATAC Down"
    } else {
      concordant <- sum(vals > 0)
      pred_dir <- "ATAC Up"
    }
    data.frame(
      group = group_name, mark = "ATAC",
      predicted_direction = pred_dir,
      n_concordant = concordant,
      n_total = length(vals),
      pct_concordant = 100 * concordant / length(vals),
      stringsAsFactors = FALSE
    )
  }

  conc_atac <- rbind(
    compute_atac_group(mc_up_genes, "mC Up", TRUE),
    compute_atac_group(mc_down_genes, "mC Down", FALSE),
    compute_atac_group(hmc_down_genes, "hmC Down", TRUE),
    compute_atac_group(hmc_up_genes, "hmC Up", FALSE)
  )

  result_16b <- build_concordance_barchart(
    conc_atac, "ATAC-seq",
    "mC Up / hmC Down genes predicted to show decreased accessibility"
  )

  save_multiformat_ggplot(result_16b$plot,
                          file.path(OUTPUT_DIR, "16b_atac_concordance_bars"),
                          width = 9, height = 7)

  cat("  Results:\n")
  for (i in 1:nrow(conc_atac)) {
    cat(sprintf("    %s: %.1f%% %s (%d/%d)\n",
                conc_atac$group[i], conc_atac$pct_concordant[i],
                conc_atac$predicted_direction[i],
                conc_atac$n_concordant[i], conc_atac$n_total[i]))
  }
  cat(sprintf("  Fisher mC Up vs Down: p = %.2e\n", result_16b$fisher_mc$p.value))
  cat(sprintf("  Fisher hmC Down vs Up: p = %.2e\n", result_16b$fisher_hmc$p.value))
}

# =============================================================================
# FIGURE 16c: K119ub Concordance
# =============================================================================
cat("\n--- FIGURE 16c: K119ub Concordance ---\n\n")

k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Control")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mutant")

if (is.null(k119ub_ctrl_gr) || is.null(k119ub_mut_gr)) {
  cat("  ERROR: K119ub data not found. Skipping 16c.\n")
  result_16c <- NULL
} else {
  cat("  Annotating K119ub peaks to genes...\n")
  k119ub_ctrl_gene <- aggregate_peaks_by_gene(
    annotate_peaks_to_genes(k119ub_ctrl_gr, txdb), "k119ub_ctrl")
  k119ub_mut_gene <- aggregate_peaks_by_gene(
    annotate_peaks_to_genes(k119ub_mut_gr, txdb), "k119ub_mut")

  k119ub_gene <- full_join(k119ub_ctrl_gene, k119ub_mut_gene, by = "SYMBOL") %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = n_k119ub_mut - n_k119ub_ctrl
    ) %>%
    dplyr::filter(net_k119ub != 0)

  cat(sprintf("  %d genes with non-zero net K119ub\n", nrow(k119ub_gene)))

  conc_k119ub <- compute_concordance(k119ub_gene, "net_k119ub", "K119ub",
                                     "K119ub Gained", "K119ub Lost")

  result_16c <- build_concordance_barchart(
    conc_k119ub, "H2AK119ub",
    "mC Up / hmC Down genes predicted to show K119ub accumulation"
  )

  save_multiformat_ggplot(result_16c$plot,
                          file.path(OUTPUT_DIR, "16c_k119ub_concordance_bars"),
                          width = 9, height = 7)

  cat("  Results:\n")
  for (i in 1:nrow(conc_k119ub)) {
    cat(sprintf("    %s: %.1f%% %s (%d/%d)\n",
                conc_k119ub$group[i], conc_k119ub$pct_concordant[i],
                conc_k119ub$predicted_direction[i],
                conc_k119ub$n_concordant[i], conc_k119ub$n_total[i]))
  }
  cat(sprintf("  Fisher mC Up vs Down: p = %.2e\n", result_16c$fisher_mc$p.value))
  cat(sprintf("  Fisher hmC Down vs Up: p = %.2e\n", result_16c$fisher_hmc$p.value))
}

# =============================================================================
# FIGURE 16d: Summary Comparison — Raw Concordance Rates (mC Up vs hmC Down)
# =============================================================================
cat("\n--- FIGURE 16d: Summary Comparison ---\n\n")

# Collect the dominant predicted groups: mC Up and hmC Down
summary_rows <- list()

collect_main_groups <- function(conc_df) {
  conc_df %>%
    dplyr::select(group, mark, predicted_direction,
                  n_concordant, n_total, pct_concordant) %>%
    dplyr::filter(group %in% c("mC Up", "hmC Down")) %>%
    mutate(perspective = ifelse(grepl("^mC", group), "5mC", "5hmC"))
}

if (!is.null(result_16a)) {
  summary_rows[[1]] <- collect_main_groups(result_16a$data)
}
if (!is.null(result_16b)) {
  summary_rows[[2]] <- collect_main_groups(conc_atac)
}
if (!is.null(result_16c)) {
  summary_rows[[3]] <- collect_main_groups(conc_k119ub)
}

summary_df <- do.call(rbind, summary_rows)

if (is.null(summary_df) || nrow(summary_df) == 0) {
  cat("  ERROR: No data for summary. Skipping 16d.\n")
} else {
  summary_df$mark <- factor(summary_df$mark, levels = c("MeCP2", "ATAC", "K119ub"))
  summary_df$perspective <- factor(summary_df$perspective, levels = c("5mC", "5hmC"))

  p_16d <- ggplot(summary_df,
                  aes(x = mark, y = pct_concordant, fill = perspective)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.6),
             width = 0.5, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)",
                                  pct_concordant, n_concordant, n_total)),
              position = position_dodge(width = 0.6),
              vjust = -0.3, size = 3, lineheight = 0.9) +
    geom_text(aes(label = predicted_direction, y = 2),
              position = position_dodge(width = 0.6),
              size = 2.3, color = "grey40", fontface = "italic") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.55, y = 51.5, label = "50% (chance)",
             size = 2.5, color = "grey50", hjust = 0) +
    scale_fill_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
                      labels = c("5mC" = "mC Up (n=3,616)",
                                 "5hmC" = "hmC Down (n=4,361)"),
                      name = "Dominant Predicted Group") +
    scale_y_continuous(limits = c(0, max(summary_df$pct_concordant) * 1.3),
                       expand = c(0, 0)) +
    labs(
      title = "Raw Concordance: mC Up vs hmC Down Across Chromatin Marks",
      subtitle = "% of genes in each dominant methylation group showing predicted mark direction",
      x = "Chromatin Mark",
      y = "% with Predicted Direction"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_16d,
                          file.path(OUTPUT_DIR, "16d_raw_concordance_comparison"),
                          width = 10, height = 8)

  # Save all concordance data
  base_cols <- c("group", "mark", "predicted_direction",
                 "n_concordant", "n_total", "pct_concordant")
  all_conc <- do.call(rbind, list(
    if (!is.null(result_16a)) result_16a$data[, base_cols] else NULL,
    if (!is.null(result_16b)) conc_atac[, base_cols] else NULL,
    if (!is.null(result_16c)) conc_k119ub[, base_cols] else NULL
  ))

  write.table(all_conc, file.path(TABLES_DIR, "raw_concordance_all_marks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: raw_concordance_all_marks.tsv\n")

  write.table(summary_df, file.path(TABLES_DIR, "raw_concordance_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: raw_concordance_summary.tsv\n")

  # Print comparison
  cat("\n  mC Up vs hmC Down comparison:\n")
  for (m in c("MeCP2", "ATAC", "K119ub")) {
    mc_row <- summary_df[summary_df$perspective == "5mC" & as.character(summary_df$mark) == m, ]
    hmc_row <- summary_df[summary_df$perspective == "5hmC" & as.character(summary_df$mark) == m, ]
    if (nrow(mc_row) > 0 && nrow(hmc_row) > 0) {
      diff <- hmc_row$pct_concordant - mc_row$pct_concordant
      winner <- ifelse(diff > 0, "hmC Down", "mC Up")
      cat(sprintf("    %s: mC Up = %.1f%%, hmC Down = %.1f%% (%+.1f pp, %s stronger)\n",
                  m, mc_row$pct_concordant, hmc_row$pct_concordant, diff, winner))
    }
  }
}

# =============================================================================
# SECTION 16 SUMMARY
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("SECTION 16 SUMMARY\n")
cat("================================================================================\n")
cat("Concordance = % of genes in methylation group showing predicted mark direction\n\n")

if (!is.null(result_16a)) {
  cat("MeCP2:\n")
  for (i in 1:nrow(result_16a$data)) {
    r <- result_16a$data[i, ]
    cat(sprintf("  %-10s %5.1f%% %s (%d/%d)\n",
                r$group, r$pct_concordant, r$predicted_direction,
                r$n_concordant, r$n_total))
  }
  cat(sprintf("  Fisher: mC p=%.1e, hmC p=%.1e\n\n",
              result_16a$fisher_mc$p.value, result_16a$fisher_hmc$p.value))
}

if (!is.null(result_16b)) {
  cat("ATAC:\n")
  for (i in 1:nrow(conc_atac)) {
    r <- conc_atac[i, ]
    cat(sprintf("  %-10s %5.1f%% %s (%d/%d)\n",
                r$group, r$pct_concordant, r$predicted_direction,
                r$n_concordant, r$n_total))
  }
  cat(sprintf("  Fisher: mC p=%.1e, hmC p=%.1e\n\n",
              result_16b$fisher_mc$p.value, result_16b$fisher_hmc$p.value))
}

if (!is.null(result_16c)) {
  cat("K119ub:\n")
  for (i in 1:nrow(conc_k119ub)) {
    r <- conc_k119ub[i, ]
    cat(sprintf("  %-10s %5.1f%% %s (%d/%d)\n",
                r$group, r$pct_concordant, r$predicted_direction,
                r$n_concordant, r$n_total))
  }
  cat(sprintf("  Fisher: mC p=%.1e, hmC p=%.1e\n\n",
              result_16c$fisher_mc$p.value, result_16c$fisher_hmc$p.value))
}

cat("Section 16 complete.\n\n")
