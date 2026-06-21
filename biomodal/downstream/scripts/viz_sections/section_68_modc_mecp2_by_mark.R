# biomodal/downstream/scripts/viz_sections/section_68_modc_mecp2_by_mark.R
# Section 68: modC vs MeCP2 Scatter by Histone Mark Category
# Standalone script - sources shared config for all dependencies and data
#
# Four quadrant scatter plots (following section 59 pattern):
#   68a: K27ac only (Active) — modC change vs MeCP2 fold
#   68b: K27me3 + K27ac (Bivalent)
#   68c: K27me3 only (Fac. Het)
#   68d: Neither
# Dots colored by methylation direction (red = hyper, blue = hypo).
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_68_modc_mecp2_by_mark.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# QUADRANT HELPER FUNCTIONS (from section 59, copied for independence)
# =============================================================================

assign_quadrant <- function(dx, dy) {
  ifelse(dx > 0 & dy > 0, "Q1_up_up",
  ifelse(dx < 0 & dy > 0, "Q2_down_up",
  ifelse(dx < 0 & dy < 0, "Q3_down_down",
  ifelse(dx > 0 & dy < 0, "Q4_up_down", "origin"))))
}

make_quadrant_labels <- function(x, y, x_range, y_range) {
  quads <- assign_quadrant(x, y)
  counts <- table(factor(quads, levels = c("Q1_up_up", "Q2_down_up",
                                           "Q3_down_down", "Q4_up_down")))
  total <- sum(counts)

  pad_x <- diff(x_range) * 0.04
  pad_y <- diff(y_range) * 0.04

  data.frame(
    label = c(
      sprintf("n=%d\n(%.0f%%)", counts["Q2_down_up"],
              100 * counts["Q2_down_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q1_up_up"],
              100 * counts["Q1_up_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q3_down_down"],
              100 * counts["Q3_down_down"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q4_up_down"],
              100 * counts["Q4_up_down"] / total)
    ),
    x = c(x_range[1] + pad_x, x_range[2] - pad_x,
          x_range[1] + pad_x, x_range[2] - pad_x),
    y = c(y_range[2] - pad_y, y_range[2] - pad_y,
          y_range[1] + pad_y, y_range[1] + pad_y),
    hjust = c(0, 1, 0, 1),
    vjust = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
}

clip_symmetric <- function(x, pctl = 0.995) {
  lim <- quantile(abs(x), pctl, na.rm = TRUE)
  c(-lim, lim)
}

scatter_quadrant_plot <- function(df, x_col, y_col,
                                  color_col = NULL, color_values = NULL,
                                  x_lab, y_lab, title,
                                  label_genes = KEY_GENES,
                                  point_alpha = 0.2, point_size = 0.5) {

  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]

  rho <- cor(x_vals, y_vals, method = "spearman", use = "complete.obs")
  rho_p <- cor.test(x_vals, y_vals, method = "spearman")$p.value
  rho_label <- if (rho_p < 2.2e-16) {
    sprintf("rho = %.3f, p < 2.2e-16", rho)
  } else {
    sprintf("rho = %.3f, p = %.2e", rho, rho_p)
  }

  x_lim <- clip_symmetric(x_vals)
  y_lim <- clip_symmetric(y_vals)

  valid <- !is.na(x_vals) & !is.na(y_vals) & x_vals != 0 & y_vals != 0
  q_labels <- make_quadrant_labels(x_vals[valid], y_vals[valid], x_lim, y_lim)

  if (!is.null(color_col)) {
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                        color = .data[[color_col]]))
  } else {
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]]))
  }

  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_text(data = q_labels, aes(x = x, y = y, label = label,
              hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 2.8, color = "grey30") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(
      title = title,
      subtitle = sprintf("Spearman %s | n = %s genes",
                         rho_label, format(sum(valid), big.mark = ",")),
      x = x_lab,
      y = y_lab
    ) +
    theme_biomodal()

  if (!is.null(color_col) && !is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values,
                                name = NULL,
                                na.value = "grey80") +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
  } else if (is.null(color_col)) {
    p <- p + aes(color = NULL)
    p$layers[[3]]$aes_params$colour <- "#404040"
  }

  if (!is.null(label_genes) && length(label_genes) > 0) {
    label_df <- df %>%
      dplyr::filter(gene %in% label_genes) %>%
      dplyr::filter(!is.na(.data[[x_col]]) & !is.na(.data[[y_col]]))
    if (nrow(label_df) > 0) {
      p <- p + geom_text_repel(
        data = label_df,
        aes(x = .data[[x_col]], y = .data[[y_col]], label = gene),
        inherit.aes = FALSE, size = 2.5, max.overlaps = 15,
        segment.color = "grey50", segment.size = 0.3,
        fontface = "italic", color = "black"
      )
    }
  }

  list(plot = p, rho = rho, rho_p = rho_p, n = sum(valid), x_lim = x_lim, y_lim = y_lim)
}

# =============================================================================
# SECTION 68: modC vs MeCP2 BY HISTONE MARK
# =============================================================================

cat("================================================================================\n")
cat("SECTION 68: modC vs MeCP2 SCATTER BY HISTONE MARK CATEGORY\n")
cat("================================================================================\n\n")

# ---- Load data ---------------------------------------------------------------

DIFFBIND_GENE_PATH <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")

stopifnot(
  "diffbind_gene_level_all_marks.tsv not found (run section_33 first)" =
    file.exists(DIFFBIND_GENE_PATH),
  "MeCP2 annotated file not found" =
    file.exists(MECP2_FILES$annotated)
)

cat("Loading gene-level data...\n")

diffbind_genes <- read.table(DIFFBIND_GENE_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind gene table: %d genes\n", nrow(diffbind_genes)))

# Aggregate MeCP2 from full annotated file (21k+ genes, unfiltered)
cat("Loading and aggregating MeCP2 from full annotated file...\n")
mecp2_raw <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, fill = TRUE, quote = "")
mecp2_raw$Fold <- as.numeric(mecp2_raw$Fold)

mecp2_gene <- mecp2_raw %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    mecp2_mean_fold = mean(Fold, na.rm = TRUE),
    .groups = "drop"
  )
cat(sprintf("  MeCP2: %d unique genes with peaks\n", nrow(mecp2_gene)))

# ---- Build master table ------------------------------------------------------

cat("\nBuilding master gene-level table...\n")

master <- diffbind_genes %>%
  dplyr::select(gene, chr, start, end,
                mc_ctrl, mc_mut, mc_diff, mc_sig,
                hmc_ctrl, hmc_mut, hmc_diff, hmc_sig,
                total_ctrl, total_mut) %>%
  left_join(mecp2_gene, by = c("gene" = "SYMBOL")) %>%
  dplyr::mutate(
    modc_diff = total_mut - total_ctrl,
    modc_direction = ifelse(modc_diff > 0, "Hypermethylated", "Hypomethylated")
  ) %>%
  dplyr::filter(is.finite(modc_diff) & is.finite(mecp2_mean_fold))

cat(sprintf("  Master table: %d genes with both modC and MeCP2 data\n", nrow(master)))

# ---- Classify by histone mark ------------------------------------------------

cat("\nClassifying genes by H3K27me3 / H3K27ac overlap...\n")

k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
k27ac_gr  <- load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac")

gene_gr <- GRanges(
  seqnames = master$chr,
  ranges = IRanges(start = master$start, end = master$end)
)

master$k27me3_overlap <- countOverlaps(gene_gr, k27me3_gr) > 0
master$k27ac_overlap  <- countOverlaps(gene_gr, k27ac_gr) > 0

master$mark_category <- case_when(
  master$k27me3_overlap & master$k27ac_overlap ~ "K27me3 + K27ac (Bivalent)",
  master$k27me3_overlap                        ~ "K27me3 only (Fac. Het)",
  master$k27ac_overlap                         ~ "K27ac only (Active)",
  TRUE                                         ~ "Neither"
)

MARK_ORDER <- c("K27ac only (Active)", "K27me3 + K27ac (Bivalent)",
                "K27me3 only (Fac. Het)", "Neither")
master$mark_category <- factor(master$mark_category, levels = MARK_ORDER)

cat("  Gene counts per category:\n")
for (cat_name in MARK_ORDER) {
  n <- sum(master$mark_category == cat_name)
  cat(sprintf("    %s: %d\n", cat_name, n))
}

# Arrange so significant dots render on top
master <- master %>%
  arrange(desc(modc_direction == "Hypermethylated"))

# ---- Color scheme ------------------------------------------------------------

DIR_COLORS <- c("Hypermethylated" = "#D7191C", "Hypomethylated" = "#2C7BB6")

# ---- Generate 4 scatter panels -----------------------------------------------

cat("\nCreating scatter plots...\n")

results_list <- list()

for (cat_name in MARK_ORDER) {
  cat_df <- master %>% dplyr::filter(mark_category == cat_name)

  short_name <- gsub(" \\(.*\\)", "", cat_name)
  panel_id <- switch(cat_name,
    "K27ac only (Active)" = "68a",
    "K27me3 + K27ac (Bivalent)" = "68b",
    "K27me3 only (Fac. Het)" = "68c",
    "Neither" = "68d"
  )

  cat(sprintf("  %s: %s (%d genes)...\n", panel_id, cat_name, nrow(cat_df)))

  result <- scatter_quadrant_plot(
    df = cat_df,
    x_col = "modc_diff",
    y_col = "mecp2_mean_fold",
    color_col = "modc_direction",
    color_values = DIR_COLORS,
    x_lab = "modC Change (Mut - Ctrl)",
    y_lab = "MeCP2 Mean Fold (log2)",
    title = sprintf("%s: modC vs MeCP2 (%s)", panel_id, cat_name),
    point_alpha = 0.3,
    point_size = 0.8
  )

  save_multiformat_ggplot(result$plot,
                          file.path(OUTPUT_DIR, sprintf("%s_modc_mecp2_%s",
                                    panel_id, gsub("[^a-z]", "_", tolower(short_name)))),
                          width = 9, height = 8)

  results_list[[panel_id]] <- list(
    plot = result$plot,
    mark_category = cat_name,
    n_genes = result$n,
    spearman_rho = result$rho,
    spearman_p = result$rho_p
  )

  # Quadrant counts for table
  valid_df <- cat_df %>%
    dplyr::filter(is.finite(modc_diff) & is.finite(mecp2_mean_fold) &
                  modc_diff != 0 & mecp2_mean_fold != 0)
  quads <- assign_quadrant(valid_df$modc_diff, valid_df$mecp2_mean_fold)
  q_counts <- table(factor(quads, levels = c("Q1_up_up", "Q2_down_up",
                                             "Q3_down_down", "Q4_up_down")))
  results_list[[panel_id]]$n_Q1 <- as.integer(q_counts["Q1_up_up"])
  results_list[[panel_id]]$n_Q2 <- as.integer(q_counts["Q2_down_up"])
  results_list[[panel_id]]$n_Q3 <- as.integer(q_counts["Q3_down_down"])
  results_list[[panel_id]]$n_Q4 <- as.integer(q_counts["Q4_up_down"])
}

# ---- Combined figure ---------------------------------------------------------

cat("\nCreating combined Figure 68...\n")

p_68_combined <- (results_list[["68a"]]$plot | results_list[["68b"]]$plot) /
                 (results_list[["68c"]]$plot | results_list[["68d"]]$plot) +
  plot_annotation(
    title = "modC Change vs MeCP2 Binding by Chromatin Context",
    subtitle = "Red = hypermethylated (modC up), Blue = hypomethylated (modC down)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_68_combined,
                        file.path(OUTPUT_DIR, "68_modc_mecp2_by_mark"),
                        width = 16, height = 14)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

# Stats summary
stats_df <- do.call(rbind, lapply(names(results_list), function(pid) {
  r <- results_list[[pid]]
  data.frame(
    panel = pid,
    mark_category = r$mark_category,
    n_genes = r$n_genes,
    spearman_rho = r$spearman_rho,
    spearman_p = r$spearman_p,
    n_Q1_hyper_mecp2up = r$n_Q1,
    n_Q2_hypo_mecp2up = r$n_Q2,
    n_Q3_hypo_mecp2down = r$n_Q3,
    n_Q4_hyper_mecp2down = r$n_Q4,
    stringsAsFactors = FALSE
  )
}))

write.table(stats_df,
            file.path(TABLES_DIR, "68_modc_mecp2_by_mark_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 68_modc_mecp2_by_mark_stats.tsv\n")

# Per-gene table
gene_out <- master %>%
  dplyr::select(gene, chr, start, end, modc_diff, mecp2_mean_fold,
                mark_category, modc_direction,
                mc_diff, hmc_diff, total_ctrl, total_mut,
                k27me3_overlap, k27ac_overlap)

write.table(gene_out,
            file.path(TABLES_DIR, "68_per_gene_modc_mecp2_mark.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 68_per_gene_modc_mecp2_mark.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 68 COMPLETE\n")
cat("================================================================================\n")
