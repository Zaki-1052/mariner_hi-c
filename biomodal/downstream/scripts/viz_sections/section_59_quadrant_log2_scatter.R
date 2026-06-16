# biomodal/downstream/scripts/viz_sections/section_59_quadrant_log2_scatter.R
# Log2 quadrant scatter plots: K119ub × MeCP2 × 5mC × 5hmC at gene-body resolution.
#
# Five pairwise comparisons testing the BAP1-KO mechanistic model:
#   59a: K119ub log2(mut/ctrl) vs MeCP2 DiffBind fold
#   59b: MeCP2 fold vs H3K27ac fold (euchromatin genes only)
#   59c: MeCP2 fold vs H3K27me3 fold (heterochromatin genes only)
#   59d: K119ub log2(mut/ctrl) vs 5mC mod_difference
#   59e: K119ub log2(mut/ctrl) vs 5hmC mod_difference
#   59f: Composite panel (patchwork)
#
# Prerequisites:
#   data/k119ub_gene_signal.tsv (from preprocess_k119ub_bigwig.R on HPC)
#   tables/diffbind_gene_level_all_marks.tsv (from section_33)
#   tables/mecp2_gene_level_correlation.tsv (from section_11)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_59_quadrant_log2_scatter.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 59: LOG2 QUADRANT SCATTER PLOTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 59: LOG2 QUADRANT SCATTER PLOTS\n")
cat("  K119ub x MeCP2 x 5mC x 5hmC gene-body integration\n")
cat("================================================================================\n\n")

# =============================================================================
# INPUT VALIDATION
# =============================================================================

K119UB_SIGNAL_PATH <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")
MECP2_GENE_PATH    <- file.path(TABLES_DIR, "mecp2_gene_level_correlation.tsv")
DIFFBIND_GENE_PATH <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")

cat("Validating input files...\n")
stopifnot(
  "K119ub gene signal file not found (run preprocess_k119ub_bigwig.R on HPC)" =
    file.exists(K119UB_SIGNAL_PATH)
)
cat(sprintf("  [OK] K119ub signal: %s\n", K119UB_SIGNAL_PATH))

stopifnot(
  "MeCP2 gene-level table not found (run section_11 first)" =
    file.exists(MECP2_GENE_PATH)
)
cat(sprintf("  [OK] MeCP2 gene table: %s\n", MECP2_GENE_PATH))

stopifnot(
  "DiffBind gene-level table not found (run section_33 first)" =
    file.exists(DIFFBIND_GENE_PATH)
)
cat(sprintf("  [OK] DiffBind gene table: %s\n", DIFFBIND_GENE_PATH))

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\nLoading data...\n")

k119ub <- read.table(K119UB_SIGNAL_PATH, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
cat(sprintf("  K119ub signal: %d genes (%d quantifiable)\n",
            nrow(k119ub), sum(k119ub$gb_signal_class == "quantifiable")))

mecp2 <- read.table(MECP2_GENE_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  MeCP2 gene table: %d genes\n", nrow(mecp2)))

diffbind <- read.table(DIFFBIND_GENE_PATH, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind gene table: %d genes\n", nrow(diffbind)))

# =============================================================================
# MASTER JOIN
# =============================================================================

cat("\nBuilding master gene-level data frame...\n")

mecp2_slim <- mecp2 %>%
  dplyr::select(gene, mecp2_nearest_fold, mecp2_nearest_fdr,
                mecp2_mean_fold, mecp2_n_peaks)

diffbind_slim <- diffbind %>%
  dplyr::select(gene, mc_ctrl, mc_mut, mc_diff, mc_q, mc_sig,
                hmc_ctrl, hmc_mut, hmc_diff, hmc_q, hmc_sig,
                chromatin_state, k27ac_fold, k27ac_fdr,
                k27me3_fold, k27me3_fdr, k119ub_fold, k119ub_fdr)

master <- k119ub %>%
  dplyr::rename(gene = symbol) %>%
  dplyr::left_join(mecp2_slim, by = "gene") %>%
  dplyr::left_join(diffbind_slim, by = "gene")

n_with_mecp2 <- sum(!is.na(master$mecp2_mean_fold))
n_with_diffbind <- sum(!is.na(master$mc_diff))
n_quantifiable <- sum(master$gb_signal_class == "quantifiable", na.rm = TRUE)

cat(sprintf("  Master table: %d genes\n", nrow(master)))
cat(sprintf("  With MeCP2 fold:   %d (%.1f%%)\n",
            n_with_mecp2, 100 * n_with_mecp2 / nrow(master)))
cat(sprintf("  With methylation:  %d (%.1f%%)\n",
            n_with_diffbind, 100 * n_with_diffbind / nrow(master)))
cat(sprintf("  K119ub quantifiable: %d\n", n_quantifiable))

# =============================================================================
# QUADRANT HELPER FUNCTIONS
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

# Factory function for quadrant scatter plots
scatter_quadrant_plot <- function(df, x_col, y_col,
                                  color_col = NULL, color_values = NULL,
                                  x_lab, y_lab, title,
                                  label_genes = KEY_GENES,
                                  point_alpha = 0.2, point_size = 0.5) {

  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]

  # Spearman correlation
  rho <- cor(x_vals, y_vals, method = "spearman", use = "complete.obs")
  rho_p <- cor.test(x_vals, y_vals, method = "spearman")$p.value
  rho_label <- if (rho_p < 2.2e-16) {
    sprintf("rho = %.3f, p < 2.2e-16", rho)
  } else {
    sprintf("rho = %.3f, p = %.2e", rho, rho_p)
  }

  # Symmetric axis clipping
 x_lim <- clip_symmetric(x_vals)
  y_lim <- clip_symmetric(y_vals)

  # Quadrant annotations
  valid <- !is.na(x_vals) & !is.na(y_vals) & x_vals != 0 & y_vals != 0
  q_labels <- make_quadrant_labels(x_vals[valid], y_vals[valid], x_lim, y_lim)

  # Base plot
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

  # Color scale
  if (!is.null(color_col) && !is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values,
                                name = NULL,
                                na.value = "grey80") +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
  } else if (is.null(color_col)) {
    p <- p + aes(color = NULL)
    p$layers[[3]]$aes_params$colour <- "#404040"
  }

  # Key gene labels
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

SEC59_DIR <- file.path(OUTPUT_DIR, "59_quadrant_log2_scatter")
dir.create(SEC59_DIR, recursive = TRUE, showWarnings = FALSE)

save_section_plot <- function(p, name, w = 10, h = 9) {
  save_multiformat_ggplot(p, base_path = file.path(SEC59_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

# =============================================================================
# PLOT 59a: K119UB log2(mut/ctrl) vs MeCP2 DiffBind FOLD
# =============================================================================

cat("\n--- Plot 59a: K119ub vs MeCP2 ---\n")

df_59a <- master %>%
  dplyr::filter(gb_signal_class == "quantifiable",
                !is.na(mecp2_mean_fold)) %>%
  dplyr::mutate(
    mecp2_status = dplyr::case_when(
      mecp2_nearest_fdr < Q_THRESHOLD & mecp2_mean_fold > 0 ~ "MeCP2 Up",
      mecp2_nearest_fdr < Q_THRESHOLD & mecp2_mean_fold < 0 ~ "MeCP2 Down",
      TRUE ~ "Not Significant"
    ),
    mecp2_status = factor(mecp2_status,
                          levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant"))
  )

# Plot NS points first, significant on top
df_59a <- df_59a %>% dplyr::arrange(desc(mecp2_status))

cat(sprintf("  Subset: %d genes\n", nrow(df_59a)))
cat(sprintf("  MeCP2 Up:   %d\n", sum(df_59a$mecp2_status == "MeCP2 Up")))
cat(sprintf("  MeCP2 Down: %d\n", sum(df_59a$mecp2_status == "MeCP2 Down")))
cat(sprintf("  NS:         %d\n", sum(df_59a$mecp2_status == "Not Significant")))

res_59a <- scatter_quadrant_plot(
  df_59a,
  x_col = "gb_log2fc", y_col = "mecp2_mean_fold",
  color_col = "mecp2_status", color_values = COLORS$mecp2,
  x_lab = expression("H2AK119ub " * log[2] * "(mut/ctrl)"),
  y_lab = "MeCP2 DiffBind fold change",
  title = "H2AK119ub vs MeCP2 at Gene Bodies"
)

save_section_plot(res_59a$plot, "59a_k119ub_vs_mecp2")
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59a$rho, res_59a$rho_p))

# =============================================================================
# PLOT 59b: MeCP2 vs H3K27ac (EUCHROMATIN)
# =============================================================================

cat("\n--- Plot 59b: MeCP2 vs H3K27ac ---\n")

df_59b <- master %>%
  dplyr::filter(!is.na(mecp2_mean_fold),
                !is.na(k27ac_fold))

cat(sprintf("  Genes with MeCP2 + H3K27ac fold: %d\n", nrow(df_59b)))

res_59b <- scatter_quadrant_plot(
  df_59b,
  x_col = "mecp2_mean_fold", y_col = "k27ac_fold",
  x_lab = "MeCP2 DiffBind fold change",
  y_lab = "H3K27ac DiffBind fold change",
  title = "MeCP2 vs H3K27ac (Euchromatin Mark)"
)

save_section_plot(res_59b$plot, "59b_mecp2_vs_k27ac", w = 9, h = 8)
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59b$rho, res_59b$rho_p))

# =============================================================================
# PLOT 59c: MeCP2 vs H3K27me3 (HETEROCHROMATIN)
# =============================================================================

cat("\n--- Plot 59c: MeCP2 vs H3K27me3 ---\n")

df_59c <- master %>%
  dplyr::filter(!is.na(mecp2_mean_fold),
                !is.na(k27me3_fold))

cat(sprintf("  Genes with MeCP2 + H3K27me3 fold: %d\n", nrow(df_59c)))

res_59c <- scatter_quadrant_plot(
  df_59c,
  x_col = "mecp2_mean_fold", y_col = "k27me3_fold",
  x_lab = "MeCP2 DiffBind fold change",
  y_lab = "H3K27me3 DiffBind fold change",
  title = "MeCP2 vs H3K27me3 (Heterochromatin Mark)"
)

save_section_plot(res_59c$plot, "59c_mecp2_vs_k27me3", w = 9, h = 8)
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59c$rho, res_59c$rho_p))

# =============================================================================
# COMPOSITE 59bc: EUCHROMATIN + HETEROCHROMATIN SIDE-BY-SIDE
# =============================================================================

cat("\n--- Composite 59bc ---\n")

p_59bc <- res_59b$plot + res_59c$plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "MeCP2 vs Euchromatin and Heterochromatin Marks",
    tag_levels = "A"
  )

save_section_plot(p_59bc, "59bc_chromatin_composite", w = 18, h = 8)

# =============================================================================
# PLOT 59d: K119UB log2(mut/ctrl) vs 5mC mod_difference
# =============================================================================

cat("\n--- Plot 59d: K119ub vs 5mC ---\n")

df_59d <- master %>%
  dplyr::filter(gb_signal_class == "quantifiable",
                !is.na(mc_diff)) %>%
  dplyr::mutate(
    mc_status = ifelse(mc_sig, "5mC DMR", "Not DMR"),
    mc_status = factor(mc_status, levels = c("5mC DMR", "Not DMR"))
  )

# NS first, DMR on top
df_59d <- df_59d %>% dplyr::arrange(mc_status)

cat(sprintf("  Subset: %d genes\n", nrow(df_59d)))
cat(sprintf("  5mC DMR: %d | Not DMR: %d\n",
            sum(df_59d$mc_status == "5mC DMR"),
            sum(df_59d$mc_status == "Not DMR")))

MC_STATUS_COLORS <- c("5mC DMR" = "#E41A1C", "Not DMR" = "grey70")

res_59d <- scatter_quadrant_plot(
  df_59d,
  x_col = "gb_log2fc", y_col = "mc_diff",
  color_col = "mc_status", color_values = MC_STATUS_COLORS,
  x_lab = expression("H2AK119ub " * log[2] * "(mut/ctrl)"),
  y_lab = expression(Delta * "5mC (mutant - control)"),
  title = "H2AK119ub vs 5mC Change at Gene Bodies"
)

save_section_plot(res_59d$plot, "59d_k119ub_vs_5mc")
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59d$rho, res_59d$rho_p))

# =============================================================================
# PLOT 59e: K119UB log2(mut/ctrl) vs 5hmC mod_difference
# =============================================================================

cat("\n--- Plot 59e: K119ub vs 5hmC ---\n")

df_59e <- master %>%
  dplyr::filter(gb_signal_class == "quantifiable",
                !is.na(hmc_diff)) %>%
  dplyr::mutate(
    hmc_status = ifelse(hmc_sig, "5hmC DMR", "Not DMR"),
    hmc_status = factor(hmc_status, levels = c("5hmC DMR", "Not DMR"))
  )

# NS first, DMR on top
df_59e <- df_59e %>% dplyr::arrange(hmc_status)

cat(sprintf("  Subset: %d genes\n", nrow(df_59e)))
cat(sprintf("  5hmC DMR: %d | Not DMR: %d\n",
            sum(df_59e$hmc_status == "5hmC DMR"),
            sum(df_59e$hmc_status == "Not DMR")))

HMC_STATUS_COLORS <- c("5hmC DMR" = "#377EB8", "Not DMR" = "grey70")

res_59e <- scatter_quadrant_plot(
  df_59e,
  x_col = "gb_log2fc", y_col = "hmc_diff",
  color_col = "hmc_status", color_values = HMC_STATUS_COLORS,
  x_lab = expression("H2AK119ub " * log[2] * "(mut/ctrl)"),
  y_lab = expression(Delta * "5hmC (mutant - control)"),
  title = "H2AK119ub vs 5hmC Change at Gene Bodies"
)

save_section_plot(res_59e$plot, "59e_k119ub_vs_5hmc")
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59e$rho, res_59e$rho_p))

# =============================================================================
# PLOT 59f: MeCP2 DiffBind FOLD vs 5mC mod_difference
# =============================================================================

cat("\n--- Plot 59f: MeCP2 vs 5mC ---\n")

df_59f <- master %>%
  dplyr::filter(!is.na(mecp2_mean_fold),
                !is.na(mc_diff)) %>%
  dplyr::mutate(
    mc_status = ifelse(mc_sig, "5mC DMR", "Not DMR"),
    mc_status = factor(mc_status, levels = c("5mC DMR", "Not DMR"))
  )

df_59f <- df_59f %>% dplyr::arrange(mc_status)

cat(sprintf("  Genes with MeCP2 + 5mC: %d\n", nrow(df_59f)))
cat(sprintf("  5mC DMR: %d | Not DMR: %d\n",
            sum(df_59f$mc_status == "5mC DMR"),
            sum(df_59f$mc_status == "Not DMR")))

MC_STATUS_COLORS_F <- c("5mC DMR" = "#E41A1C", "Not DMR" = "grey70")

res_59f <- scatter_quadrant_plot(
  df_59f,
  x_col = "mecp2_mean_fold", y_col = "mc_diff",
  color_col = "mc_status", color_values = MC_STATUS_COLORS_F,
  x_lab = "MeCP2 DiffBind fold change",
  y_lab = expression(Delta * "5mC (mutant - control)"),
  title = "MeCP2 vs 5mC Change at Gene Bodies"
)

save_section_plot(res_59f$plot, "59f_mecp2_vs_5mc")
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59f$rho, res_59f$rho_p))

# =============================================================================
# PLOT 59g: MeCP2 DiffBind FOLD vs 5hmC mod_difference
# =============================================================================

cat("\n--- Plot 59g: MeCP2 vs 5hmC ---\n")

df_59g <- master %>%
  dplyr::filter(!is.na(mecp2_mean_fold),
                !is.na(hmc_diff)) %>%
  dplyr::mutate(
    hmc_status = ifelse(hmc_sig, "5hmC DMR", "Not DMR"),
    hmc_status = factor(hmc_status, levels = c("5hmC DMR", "Not DMR"))
  )

df_59g <- df_59g %>% dplyr::arrange(hmc_status)

cat(sprintf("  Genes with MeCP2 + 5hmC: %d\n", nrow(df_59g)))
cat(sprintf("  5hmC DMR: %d | Not DMR: %d\n",
            sum(df_59g$hmc_status == "5hmC DMR"),
            sum(df_59g$hmc_status == "Not DMR")))

HMC_STATUS_COLORS_G <- c("5hmC DMR" = "#377EB8", "Not DMR" = "grey70")

res_59g <- scatter_quadrant_plot(
  df_59g,
  x_col = "mecp2_mean_fold", y_col = "hmc_diff",
  color_col = "hmc_status", color_values = HMC_STATUS_COLORS_G,
  x_lab = "MeCP2 DiffBind fold change",
  y_lab = expression(Delta * "5hmC (mutant - control)"),
  title = "MeCP2 vs 5hmC Change at Gene Bodies"
)

save_section_plot(res_59g$plot, "59g_mecp2_vs_5hmc")
cat(sprintf("  Spearman rho = %.3f (p = %.2e)\n", res_59g$rho, res_59g$rho_p))

# =============================================================================
# COMPOSITE 59h: ALL 7 PLOTS
# =============================================================================

cat("\n--- Composite 59h ---\n")

p_59h <- (res_59a$plot + plot_spacer()) /
         (res_59b$plot + res_59c$plot) /
         (res_59d$plot + res_59e$plot) /
         (res_59f$plot + res_59g$plot) +
  plot_annotation(
    title = "Ubiquitin-MeCP2-Methylation Quadrant Analysis",
    tag_levels = "A"
  )

save_section_plot(p_59h, "59h_composite", w = 20, h = 32)

# =============================================================================
# EXPORT MASTER TABLE
# =============================================================================

cat("\n--- Exporting master table ---\n")

export_cols <- c("gene", "chr", "start", "end",
                 "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc", "gb_signal_class",
                 "mecp2_mean_fold", "mecp2_nearest_fdr", "mecp2_n_peaks",
                 "mc_ctrl", "mc_mut", "mc_diff", "mc_q", "mc_sig",
                 "hmc_ctrl", "hmc_mut", "hmc_diff", "hmc_q", "hmc_sig",
                 "chromatin_state",
                 "k27ac_fold", "k27ac_fdr",
                 "k27me3_fold", "k27me3_fdr",
                 "k119ub_fold", "k119ub_fdr")

existing_cols <- intersect(export_cols, colnames(master))
export_df <- master[, existing_cols]

export_path <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
write.table(export_df, export_path, sep = "\t", row.names = FALSE,
            quote = FALSE)
cat(sprintf("  Saved: %s (%d genes, %d columns)\n",
            export_path, nrow(export_df), ncol(export_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 59 SUMMARY\n")
cat("================================================================================\n\n")

results <- list(
  "59a K119ub vs MeCP2"               = res_59a,
  "59b MeCP2 vs K27ac (euchromatin)"   = res_59b,
  "59c MeCP2 vs K27me3 (heterochrom.)" = res_59c,
  "59d K119ub vs 5mC"                  = res_59d,
  "59e K119ub vs 5hmC"                 = res_59e,
  "59f MeCP2 vs 5mC"                   = res_59f,
  "59g MeCP2 vs 5hmC"                  = res_59g
)

for (name in names(results)) {
  r <- results[[name]]
  cat(sprintf("  %-40s  n=%6s  rho=%+.3f\n",
              name, format(r$n, big.mark = ","), r$rho))
}

cat("\nAll plots saved to:", SEC59_DIR, "\n")
cat("Master table saved to:", export_path, "\n")
cat("Section 59 complete.\n")
