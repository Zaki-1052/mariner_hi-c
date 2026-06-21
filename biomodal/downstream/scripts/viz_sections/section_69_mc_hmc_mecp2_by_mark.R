# biomodal/downstream/scripts/viz_sections/section_69_mc_hmc_mecp2_by_mark.R
# Section 69: 5mC vs 5hmC at MeCP2 — Does MeCP2 Distinguish Modification Type?
# Standalone script - sources shared config for all dependencies and data
#
# Same 4-panel layout as section 68 but each gene plotted TWICE:
#   Red dot  = 5mC change (x) vs MeCP2 fold (y)
#   Blue dot = 5hmC change (x) vs MeCP2 fold (y)
# If red and blue are garbled together, MeCP2 doesn't distinguish 5mC from 5hmC.
#
#   Panel 69a: K27ac only (Active)
#   Panel 69b: K27me3 + K27ac (Bivalent)
#   Panel 69c: K27me3 only (Fac. Het)
#   Panel 69d: Neither
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_69_mc_hmc_mecp2_by_mark.R

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
# SECTION 69: 5mC vs 5hmC AT MeCP2 BY HISTONE MARK
# =============================================================================

cat("================================================================================\n")
cat("SECTION 69: 5mC vs 5hmC AT MeCP2 — MODIFICATION TYPE BY MARK CATEGORY\n")
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
                mc_diff, hmc_diff) %>%
  left_join(mecp2_gene, by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(is.finite(mc_diff) & is.finite(hmc_diff) &
                is.finite(mecp2_mean_fold))

cat(sprintf("  Master table: %d genes\n", nrow(master)))

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
  cat(sprintf("    %s: %d\n", cat_name, sum(master$mark_category == cat_name)))
}

# ---- Pivot to long: each gene appears twice (5mC and 5hmC) ------------------

cat("\nPivoting to long format (2 rows per gene: 5mC + 5hmC)...\n")

master_long <- rbind(
  data.frame(gene = master$gene, chr = master$chr, start = master$start, end = master$end,
             meth_diff = master$mc_diff,
             mecp2_fold = master$mecp2_mean_fold,
             mark_category = master$mark_category,
             modality = "5mC",
             stringsAsFactors = FALSE),
  data.frame(gene = master$gene, chr = master$chr, start = master$start, end = master$end,
             meth_diff = master$hmc_diff,
             mecp2_fold = master$mecp2_mean_fold,
             mark_category = master$mark_category,
             modality = "5hmC",
             stringsAsFactors = FALSE)
) %>%
  dplyr::mutate(
    modality = factor(modality, levels = c("5mC", "5hmC"))
  )

cat(sprintf("  Long table: %d rows (%d genes x 2 modalities)\n",
            nrow(master_long), nrow(master)))

MOD_COLORS <- COLORS$methylation  # "5mC" = "#E41A1C", "5hmC" = "#377EB8"

# ---- Generate 4 scatter panels -----------------------------------------------

cat("\nCreating scatter plots...\n")

results_list <- list()

for (cat_name in MARK_ORDER) {
  cat_df <- master_long %>% dplyr::filter(mark_category == cat_name)

  panel_id <- switch(cat_name,
    "K27ac only (Active)" = "69a",
    "K27me3 + K27ac (Bivalent)" = "69b",
    "K27me3 only (Fac. Het)" = "69c",
    "Neither" = "69d"
  )

  n_genes <- nrow(cat_df) / 2
  cat(sprintf("  %s: %s (%d genes, %d dots)...\n",
              panel_id, cat_name, n_genes, nrow(cat_df)))

  # Compute per-modality Spearman
  mc_sub <- cat_df %>% dplyr::filter(modality == "5mC")
  hmc_sub <- cat_df %>% dplyr::filter(modality == "5hmC")

  rho_mc <- cor(mc_sub$meth_diff, mc_sub$mecp2_fold, method = "spearman", use = "complete.obs")
  rho_hmc <- cor(hmc_sub$meth_diff, hmc_sub$mecp2_fold, method = "spearman", use = "complete.obs")
  rho_mc_p <- cor.test(mc_sub$meth_diff, mc_sub$mecp2_fold, method = "spearman")$p.value
  rho_hmc_p <- cor.test(hmc_sub$meth_diff, hmc_sub$mecp2_fold, method = "spearman")$p.value

  rho_mc_label <- ifelse(rho_mc_p < 2.2e-16,
                         sprintf("5mC rho = %.3f, p < 2.2e-16", rho_mc),
                         sprintf("5mC rho = %.3f, p = %.1e", rho_mc, rho_mc_p))
  rho_hmc_label <- ifelse(rho_hmc_p < 2.2e-16,
                          sprintf("5hmC rho = %.3f, p < 2.2e-16", rho_hmc),
                          sprintf("5hmC rho = %.3f, p = %.1e", rho_hmc, rho_hmc_p))

  # Symmetric axis limits (shared across both modalities)
  x_lim <- clip_symmetric(cat_df$meth_diff)
  y_lim <- clip_symmetric(cat_df$mecp2_fold)

  # Quadrant annotations (combined, all dots)
  valid <- !is.na(cat_df$meth_diff) & !is.na(cat_df$mecp2_fold) &
           cat_df$meth_diff != 0 & cat_df$mecp2_fold != 0
  q_labels <- make_quadrant_labels(cat_df$meth_diff[valid], cat_df$mecp2_fold[valid],
                                   x_lim, y_lim)

  # Render 5hmC first (blue behind), 5mC on top (red)
  cat_df <- cat_df %>% arrange(modality == "5mC")

  p <- ggplot(cat_df, aes(x = meth_diff, y = mecp2_fold, color = modality)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_text(data = q_labels, aes(x = x, y = y, label = label,
              hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 2.8, color = "grey30") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    scale_color_manual(values = MOD_COLORS, name = "Modification") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(
      title = sprintf("%s: %s", panel_id, cat_name),
      subtitle = sprintf("%s | %s | n = %s genes",
                         rho_mc_label, rho_hmc_label,
                         format(n_genes, big.mark = ",")),
      x = "Methylation Change (Mut - Ctrl)",
      y = "MeCP2 Mean Fold (log2)"
    ) +
    theme_biomodal()

  # Key gene labels (from 5mC layer only to avoid duplicates)
  mc_label_df <- mc_sub %>%
    dplyr::filter(gene %in% KEY_GENES & is.finite(meth_diff) & is.finite(mecp2_fold))
  if (nrow(mc_label_df) > 0) {
    p <- p + geom_text_repel(
      data = mc_label_df,
      aes(x = meth_diff, y = mecp2_fold, label = gene),
      inherit.aes = FALSE, size = 2.5, max.overlaps = 15,
      segment.color = "grey50", segment.size = 0.3,
      fontface = "italic", color = "black"
    )
  }

  save_multiformat_ggplot(p,
                          file.path(OUTPUT_DIR, sprintf("%s_mc_hmc_mecp2_%s",
                                    panel_id, gsub("[^a-z]", "_", tolower(gsub(" \\(.*\\)", "", cat_name))))),
                          width = 9, height = 8)

  # Interaction test: does MeCP2 respond differently to 5mC vs 5hmC?
  interaction_model <- lm(mecp2_fold ~ meth_diff * modality, data = cat_df)
  interaction_summary <- summary(interaction_model)
  interaction_coefs <- coef(interaction_summary)

  # The interaction term is "meth_diff:modality5hmC"
  interaction_row <- interaction_coefs["meth_diff:modality5hmC", , drop = FALSE]
  interaction_p <- interaction_row[1, "Pr(>|t|)"]
  interaction_est <- interaction_row[1, "Estimate"]

  cat(sprintf("    Interaction test (meth_diff:modality): est=%+.4f, p=%.2e\n",
              interaction_est, interaction_p))
  cat(sprintf("    %s\n", ifelse(interaction_p < 0.05,
              "=> MeCP2 DOES distinguish 5mC from 5hmC in this context",
              "=> MeCP2 does NOT distinguish 5mC from 5hmC (interaction NS)")))

  # Fisher z-test for comparing two Spearman correlations
  z_mc <- atanh(rho_mc)
  z_hmc <- atanh(rho_hmc)
  se_diff <- sqrt(1 / (nrow(mc_sub) - 3) + 1 / (nrow(hmc_sub) - 3))
  z_stat <- (z_mc - z_hmc) / se_diff
  fisher_z_p <- 2 * pnorm(-abs(z_stat))

  cat(sprintf("    Fisher z-test (rho_5mC vs rho_5hmC): z=%.3f, p=%.2e\n",
              z_stat, fisher_z_p))

  results_list[[panel_id]] <- list(
    plot = p,
    mark_category = cat_name,
    n_genes = as.integer(n_genes),
    rho_mc = rho_mc, rho_mc_p = rho_mc_p,
    rho_hmc = rho_hmc, rho_hmc_p = rho_hmc_p,
    interaction_est = interaction_est,
    interaction_p = interaction_p,
    fisher_z = z_stat,
    fisher_z_p = fisher_z_p
  )
}

# ---- Combined figure ---------------------------------------------------------

cat("\nCreating combined Figure 69...\n")

p_69_combined <- (results_list[["69a"]]$plot | results_list[["69b"]]$plot) /
                 (results_list[["69c"]]$plot | results_list[["69d"]]$plot) +
  plot_annotation(
    title = "5mC vs 5hmC at MeCP2 by Chromatin Context",
    subtitle = "Red = 5mC change, Blue = 5hmC change | If garbled: MeCP2 doesn't distinguish modification type",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_69_combined,
                        file.path(OUTPUT_DIR, "69_mc_hmc_mecp2_by_mark"),
                        width = 16, height = 14)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

stats_df <- do.call(rbind, lapply(names(results_list), function(pid) {
  r <- results_list[[pid]]
  data.frame(
    panel = pid,
    mark_category = r$mark_category,
    n_genes = r$n_genes,
    spearman_rho_5mC = r$rho_mc,
    spearman_p_5mC = r$rho_mc_p,
    spearman_rho_5hmC = r$rho_hmc,
    spearman_p_5hmC = r$rho_hmc_p,
    interaction_estimate = r$interaction_est,
    interaction_p = r$interaction_p,
    interaction_significant = r$interaction_p < 0.05,
    fisher_z_rho_comparison = r$fisher_z,
    fisher_z_p = r$fisher_z_p,
    stringsAsFactors = FALSE
  )
}))

write.table(stats_df,
            file.path(TABLES_DIR, "69_mc_hmc_mecp2_by_mark_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 69_mc_hmc_mecp2_by_mark_stats.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 69 COMPLETE\n")
cat("================================================================================\n")
