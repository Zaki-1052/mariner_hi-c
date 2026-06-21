# biomodal/downstream/scripts/viz_sections/section_70_mc_hmc_sig_category_by_mark.R
# Section 70: Methylation vs MeCP2 Colored by Significance Category
# Standalone script - sources shared config for all dependencies and data
#
# Same 4-panel layout as section 69, but dots colored by which modification
# changed significantly:
#   Red    = only 5mC significant (1 dot per gene)
#   Blue   = only 5hmC significant (1 dot per gene)
#   Purple = both 5mC and 5hmC significant (2 dots per gene)
#   Grey   = neither significant (excluded from main plot)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_70_mc_hmc_sig_category_by_mark.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# QUADRANT HELPER FUNCTIONS (from section 59)
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
      sprintf("n=%d\n(%.0f%%)", counts["Q2_down_up"],    100 * counts["Q2_down_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q1_up_up"],      100 * counts["Q1_up_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q3_down_down"],  100 * counts["Q3_down_down"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q4_up_down"],    100 * counts["Q4_up_down"] / total)
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

# =============================================================================
# SECTION 70
# =============================================================================

cat("================================================================================\n")
cat("SECTION 70: METHYLATION vs MeCP2 BY SIGNIFICANCE CATEGORY\n")
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

# Aggregate MeCP2 from full annotated file (21k+ genes, not the 9.8k filtered table)
cat("Loading and aggregating MeCP2 from full annotated file...\n")
mecp2_raw <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, fill = TRUE, quote = "")
mecp2_raw$Fold <- as.numeric(mecp2_raw$Fold)
mecp2_raw$FDR <- as.numeric(mecp2_raw$FDR)
mecp2_raw$distanceToTSS <- as.numeric(mecp2_raw$distanceToTSS)

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
                mc_diff, hmc_diff, mc_sig, hmc_sig) %>%
  left_join(mecp2_gene, by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(is.finite(mc_diff) & is.finite(hmc_diff) &
                is.finite(mecp2_mean_fold))

master$mc_sig <- as.logical(master$mc_sig)
master$hmc_sig <- as.logical(master$hmc_sig)

master$sig_category <- case_when(
  master$mc_sig & master$hmc_sig   ~ "Both",
  master$mc_sig                    ~ "5mC only",
  master$hmc_sig                   ~ "5hmC only",
  TRUE                             ~ "Neither"
)

cat(sprintf("  %d genes total\n", nrow(master)))
cat("  Significance categories:\n")
for (cat_name in c("Both", "5mC only", "5hmC only", "Neither")) {
  cat(sprintf("    %s: %d\n", cat_name, sum(master$sig_category == cat_name)))
}

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

# ---- Build long-format plot data ---------------------------------------------
# "Both" genes: 2 purple dots (one for mc_diff, one for hmc_diff)
# "5mC only": 1 red dot at mc_diff
# "5hmC only": 1 blue dot at hmc_diff
# "Neither": excluded

cat("\nBuilding plot data...\n")

both_genes <- master %>% dplyr::filter(sig_category == "Both")
mc_only_genes <- master %>% dplyr::filter(sig_category == "5mC only")
hmc_only_genes <- master %>% dplyr::filter(sig_category == "5hmC only")

make_plot_rows <- function(genes_df, diff_col, plot_cat, color_grp) {
  if (nrow(genes_df) == 0) return(NULL)
  data.frame(gene = genes_df$gene,
             meth_diff = genes_df[[diff_col]],
             mecp2_fold = genes_df$mecp2_mean_fold,
             mark_category = genes_df$mark_category,
             plot_category = plot_cat,
             color_group = color_grp,
             stringsAsFactors = FALSE)
}

plot_data <- do.call(rbind, Filter(Negate(is.null), list(
  make_plot_rows(both_genes, "mc_diff", "Both (5mC dot)", "Both 5mC + 5hmC"),
  make_plot_rows(both_genes, "hmc_diff", "Both (5hmC dot)", "Both 5mC + 5hmC"),
  make_plot_rows(mc_only_genes, "mc_diff", "5mC only", "5mC only"),
  make_plot_rows(hmc_only_genes, "hmc_diff", "5hmC only", "5hmC only")
)))

plot_data$color_group <- factor(plot_data$color_group,
                                levels = c("5mC only", "5hmC only", "Both 5mC + 5hmC"))

cat(sprintf("  Plot data: %d dots (%d genes sig in both x2 + %d 5mC-only + %d 5hmC-only)\n",
            nrow(plot_data), nrow(both_genes), nrow(mc_only_genes), nrow(hmc_only_genes)))

SIG_COLORS <- c("5mC only" = "#E41A1C", "5hmC only" = "#377EB8",
                "Both 5mC + 5hmC" = "#984EA3")

# ---- Generate 4 scatter panels -----------------------------------------------

cat("\nCreating scatter plots...\n")

results_list <- list()

for (mark_name in MARK_ORDER) {
  mark_df <- plot_data %>% dplyr::filter(mark_category == mark_name)

  panel_id <- switch(mark_name,
    "K27ac only (Active)" = "70a",
    "K27me3 + K27ac (Bivalent)" = "70b",
    "K27me3 only (Fac. Het)" = "70c",
    "Neither" = "70d"
  )

  n_genes_unique <- length(unique(mark_df$gene))
  cat(sprintf("  %s: %s (%d unique genes, %d dots)...\n",
              panel_id, mark_name, n_genes_unique, nrow(mark_df)))

  # Per-category counts in this mark
  cat_counts <- mark_df %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    count(color_group)
  for (i in seq_len(nrow(cat_counts))) {
    cat(sprintf("    %s: %d genes\n", cat_counts$color_group[i], cat_counts$n[i]))
  }

  # Spearman on all dots combined
  rho <- cor(mark_df$meth_diff, mark_df$mecp2_fold, method = "spearman", use = "complete.obs")
  rho_p <- cor.test(mark_df$meth_diff, mark_df$mecp2_fold, method = "spearman")$p.value
  rho_label <- ifelse(rho_p < 2.2e-16,
                      sprintf("rho = %.3f, p < 2.2e-16", rho),
                      sprintf("rho = %.3f, p = %.1e", rho, rho_p))

  x_lim <- clip_symmetric(mark_df$meth_diff)
  y_lim <- clip_symmetric(mark_df$mecp2_fold)

  valid <- !is.na(mark_df$meth_diff) & !is.na(mark_df$mecp2_fold) &
           mark_df$meth_diff != 0 & mark_df$mecp2_fold != 0
  q_labels <- make_quadrant_labels(mark_df$meth_diff[valid], mark_df$mecp2_fold[valid],
                                   x_lim, y_lim)

  # Render "Both" behind, then single-sig on top
  mark_df <- mark_df %>% arrange(color_group == "Both 5mC + 5hmC")

  p <- ggplot(mark_df, aes(x = meth_diff, y = mecp2_fold, color = color_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_text(data = q_labels, aes(x = x, y = y, label = label,
              hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 2.8, color = "grey30") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    scale_color_manual(values = SIG_COLORS, name = "Sig Category") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(
      title = sprintf("%s: %s", panel_id, mark_name),
      subtitle = sprintf("Spearman %s | n = %s genes (%d dots)",
                         rho_label, format(n_genes_unique, big.mark = ","), nrow(mark_df)),
      x = "Methylation Change (Mut - Ctrl)",
      y = "MeCP2 Mean Fold (log2)"
    ) +
    theme_biomodal()

  save_multiformat_ggplot(p,
                          file.path(OUTPUT_DIR, sprintf("%s_sig_category_%s",
                                    panel_id, gsub("[^a-z]", "_", tolower(gsub(" \\(.*\\)", "", mark_name))))),
                          width = 9, height = 8)

  results_list[[panel_id]] <- list(
    plot = p,
    mark_category = mark_name,
    n_genes = n_genes_unique,
    n_dots = nrow(mark_df),
    rho = rho,
    rho_p = rho_p
  )
}

# ---- Combined figure ---------------------------------------------------------

cat("\nCreating combined Figure 70...\n")

p_70_combined <- (results_list[["70a"]]$plot | results_list[["70b"]]$plot) /
                 (results_list[["70c"]]$plot | results_list[["70d"]]$plot) +
  plot_annotation(
    title = "Methylation vs MeCP2 by Significance Category",
    subtitle = "Red = 5mC only sig | Blue = 5hmC only sig | Purple = both sig (2 dots per gene)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_70_combined,
                        file.path(OUTPUT_DIR, "70_sig_category_by_mark"),
                        width = 16, height = 14)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

stats_df <- do.call(rbind, lapply(names(results_list), function(pid) {
  r <- results_list[[pid]]
  data.frame(
    panel = pid,
    mark_category = r$mark_category,
    n_genes = r$n_genes,
    n_dots = r$n_dots,
    spearman_rho = r$rho,
    spearman_p = r$rho_p,
    stringsAsFactors = FALSE
  )
}))

# Add per-mark sig category breakdown
for (mark_name in MARK_ORDER) {
  mark_master <- master %>% dplyr::filter(mark_category == mark_name)
  pid <- switch(mark_name,
    "K27ac only (Active)" = "70a",
    "K27me3 + K27ac (Bivalent)" = "70b",
    "K27me3 only (Fac. Het)" = "70c",
    "Neither" = "70d"
  )
  stats_df$n_both[stats_df$panel == pid] <- sum(mark_master$sig_category == "Both")
  stats_df$n_mc_only[stats_df$panel == pid] <- sum(mark_master$sig_category == "5mC only")
  stats_df$n_hmc_only[stats_df$panel == pid] <- sum(mark_master$sig_category == "5hmC only")
  stats_df$n_neither[stats_df$panel == pid] <- sum(mark_master$sig_category == "Neither")
}

write.table(stats_df,
            file.path(TABLES_DIR, "70_sig_category_by_mark_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 70_sig_category_by_mark_stats.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 70 COMPLETE\n")
cat("================================================================================\n")
