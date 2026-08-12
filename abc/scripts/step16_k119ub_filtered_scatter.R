# abc/scripts/step16_k119ub_filtered_scatter.R
# Filter K119ub-vs-ΔABC scatter to differentially ubiquitinated enhancers.
#
# Jesse Dixon noted that 90% of per-enhancer data clusters near zero, hiding
# the real K119ub-ABC trend. Filtering to FDR < 0.05 K119ub sites (from
# DiffBind) should produce a much cleaner scatter showing the expected negative
# correlation between K119ub gain and ABC score loss.
#
# Inputs (relative to abc/ working directory):
#   ../data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_enhancer_merged.tsv
#
# Outputs:
#   results/figures/k119ub_filtered_scatter/  — 5 plots (PDF+SVG+JPG)
#   results/figures/k119ub_filtered_scatter/k119ub_filter_summary.tsv
#
# Usage:
#   cd abc && Rscript scripts/step16_k119ub_filtered_scatter.R

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE <- "../data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_enhancer_merged.tsv"
OUTPUT_DIR <- "results/figures/k119ub_filtered_scatter"
MULTIFORMAT_UTIL <- "../scripts/utils/multi_format_output.R"

FDR_PRIMARY   <- 0.05
FDR_SECONDARY <- 0.10

SIG_COLORS <- c(Sig_Down = "#2166AC", NS = "grey70", Sig_Up = "#B2182B")

cat("================================================================================\n")
cat("STEP 16: K119UB FILTERED SCATTER (FDR < 0.05)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggpointdensity)
})

stopifnot(file.exists(MULTIFORMAT_UTIL))
source(MULTIFORMAT_UTIL)
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
stopifnot(file.exists(INPUT_FILE))
cat(sprintf("  [OK] %s\n", INPUT_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

run_spearman <- function(x, y, label) {
  valid <- !is.na(x) & !is.na(y)
  n <- sum(valid)
  if (n < 10) {
    cat(sprintf("  %-50s: n=%d (too few)\n", label, n))
    return(data.frame(label = label, spearman_rho = NA, spearman_p = NA,
                      n = n, ci_lo = NA, ci_hi = NA,
                      stringsAsFactors = FALSE))
  }
  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  z <- atanh(test$estimate)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1.96, 1.96) * se)
  cat(sprintf("  %-50s: rho = %+.4f [%+.4f, %+.4f], p = %.2e, n = %d\n",
              label, test$estimate, ci[1], ci[2], test$p.value, n))
  data.frame(
    label = label,
    spearman_rho = as.numeric(test$estimate),
    spearman_p = test$p.value,
    n = n,
    ci_lo = ci[1], ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

format_rho <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  sprintf("rho = %+.3f\np = %.1e\nn = %s",
          test$estimate, test$p.value, format(sum(valid), big.mark = ","))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
enh <- read.table(INPUT_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  Total enhancers: %d\n", nrow(enh)))

# =============================================================================
# DEFINE SUBSETS
# =============================================================================

cat("\nDefining significance subsets...\n")

enh$k119ub_sig <- ifelse(
  enh$FDR < FDR_PRIMARY & enh$Fold > 0, "Sig_Up",
  ifelse(enh$FDR < FDR_PRIMARY & enh$Fold < 0, "Sig_Down", "NS")
)
enh$k119ub_sig <- factor(enh$k119ub_sig, levels = c("Sig_Down", "NS", "Sig_Up"))

sig_05 <- enh[enh$FDR < FDR_PRIMARY, ]
sig_10 <- enh[enh$FDR < FDR_SECONDARY, ]

cat(sprintf("  All enhancers:  %d\n", nrow(enh)))
cat(sprintf("  FDR < %.2f:     %d (%.1f%%)\n",
            FDR_PRIMARY, nrow(sig_05), 100 * nrow(sig_05) / nrow(enh)))
cat(sprintf("  FDR < %.2f:     %d (%.1f%%)\n",
            FDR_SECONDARY, nrow(sig_10), 100 * nrow(sig_10) / nrow(enh)))

sig_counts <- table(enh$k119ub_sig)
for (nm in names(sig_counts)) {
  cat(sprintf("    %-10s: %d\n", nm, sig_counts[nm]))
}

# =============================================================================
# CORRELATIONS
# =============================================================================

cat("\n--- Correlation analysis ---\n")

corr_results <- list()

corr_results[[1]] <- run_spearman(
  enh$Fold, enh$mean_delta_unnorm,
  "K119ub Fold vs mean_delta_unnorm [all]"
)
corr_results[[2]] <- run_spearman(
  sig_05$Fold, sig_05$mean_delta_unnorm,
  sprintf("K119ub Fold vs mean_delta_unnorm [FDR<%.2f]", FDR_PRIMARY)
)
corr_results[[3]] <- run_spearman(
  sig_10$Fold, sig_10$mean_delta_unnorm,
  sprintf("K119ub Fold vs mean_delta_unnorm [FDR<%.2f]", FDR_SECONDARY)
)
corr_results[[4]] <- run_spearman(
  sig_05$Fold, sig_05$delta_activity,
  sprintf("K119ub Fold vs delta_activity [FDR<%.2f]", FDR_PRIMARY)
)
corr_results[[5]] <- run_spearman(
  sig_05$Fold, sig_05$mean_delta_abc,
  sprintf("K119ub Fold vs mean_delta_abc [FDR<%.2f]", FDR_PRIMARY)
)

corr_table <- do.call(rbind, corr_results)

corr_out <- file.path(OUTPUT_DIR, "k119ub_filter_summary.tsv")
write.table(corr_table, corr_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", corr_out))

# =============================================================================
# PLOT 01: ALL ENHANCERS (baseline, the "90% blob")
# =============================================================================

cat("\n--- Plot 01: All enhancers (baseline) ---\n")

p01 <- ggplot(enh, aes(x = Fold, y = mean_delta_unnorm)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh$Fold, enh$mean_delta_unnorm),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    y = expression("Mean " * Delta * "(A " %*% " C) per enhancer"),
    title = "All enhancers (n = 35,777)"
  ) +
  theme_pub

save_plot(p01, "01_all_enhancers_baseline", w = 7, h = 6)

# =============================================================================
# PLOT 02: FDR < 0.05 FILTERED (primary result)
# =============================================================================

cat("\n--- Plot 02: FDR < 0.05 filtered ---\n")

p02 <- ggplot(sig_05, aes(x = Fold, y = mean_delta_unnorm)) +
  geom_pointdensity(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(sig_05$Fold, sig_05$mean_delta_unnorm),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    y = expression("Mean " * Delta * "(A " %*% " C) per enhancer"),
    title = sprintf("Differentially ubiquitinated enhancers (FDR < %.2f, n = %s)",
                    FDR_PRIMARY, format(nrow(sig_05), big.mark = ","))
  ) +
  theme_pub

save_plot(p02, "02_fdr05_filtered", w = 7, h = 6)

# =============================================================================
# PLOT 03: SIDE-BY-SIDE COMPOSITE (the key figure)
# =============================================================================

cat("\n--- Plot 03: Before/after composite ---\n")

rho_all <- corr_table$spearman_rho[grepl("\\[all\\]", corr_table$label)]
rho_filt <- corr_table$spearman_rho[grepl(
  sprintf("\\[FDR<%.2f\\]", FDR_PRIMARY), corr_table$label) &
  grepl("delta_unnorm", corr_table$label)]

p03a <- ggplot(enh, aes(x = Fold, y = mean_delta_unnorm)) +
  geom_pointdensity(size = 0.2, alpha = 0.6) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh$Fold, enh$mean_delta_unnorm),
           hjust = 1.1, vjust = 1.3, size = 3.2, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (log"[2] * "FC)"),
    y = expression("Mean " * Delta * "(A " %*% " C)"),
    title = sprintf("All enhancers (rho = %+.3f)", rho_all)
  ) +
  theme_pub +
  theme(legend.position = "none")

p03b <- ggplot(sig_05, aes(x = Fold, y = mean_delta_unnorm)) +
  geom_pointdensity(size = 0.4, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(sig_05$Fold, sig_05$mean_delta_unnorm),
           hjust = 1.1, vjust = 1.3, size = 3.2, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (log"[2] * "FC)"),
    y = expression("Mean " * Delta * "(A " %*% " C)"),
    title = sprintf("FDR < %.2f only (rho = %+.3f)", FDR_PRIMARY, rho_filt)
  ) +
  theme_pub +
  theme(legend.position = "none")

p03 <- p03a + p03b +
  plot_annotation(
    title = "K119ub fold change vs ABC score change: effect of filtering",
    subtitle = sprintf("Filtering to differentially ubiquitinated enhancers removes the uninformative blob"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(p03, "03_composite_before_after", w = 13, h = 6)

# =============================================================================
# PLOT 04: COLORED BY DIRECTION (Sig_Up / Sig_Down)
# =============================================================================

cat("\n--- Plot 04: FDR < 0.05, colored by direction ---\n")

sig_05_dir <- sig_05
sig_05_dir$direction <- ifelse(sig_05_dir$Fold > 0, "Sig_Up", "Sig_Down")
sig_05_dir$direction <- factor(sig_05_dir$direction,
                               levels = c("Sig_Down", "Sig_Up"))

n_up <- sum(sig_05_dir$direction == "Sig_Up")
n_down <- sum(sig_05_dir$direction == "Sig_Down")

p04 <- ggplot(sig_05_dir, aes(x = Fold, y = mean_delta_unnorm,
                               color = direction)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", linewidth = 0.6, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(Sig_Down = "#2166AC", Sig_Up = "#B2182B"),
    labels = c(sprintf("K119ub Down (n=%s)", format(n_down, big.mark = ",")),
               sprintf("K119ub Up (n=%s)", format(n_up, big.mark = ","))),
    name = NULL
  ) +
  labs(
    x = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    y = expression("Mean " * Delta * "(A " %*% " C) per enhancer"),
    title = sprintf("K119ub direction at differentially ubiquitinated enhancers (FDR < %.2f)",
                    FDR_PRIMARY)
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p04, "04_fdr05_by_direction", w = 7, h = 6)

# =============================================================================
# PLOT 05: ACTIVITY COMPONENT AT SIGNIFICANT SITES
# =============================================================================

cat("\n--- Plot 05: K119ub Fold vs delta_activity (FDR < 0.05) ---\n")

p05 <- ggplot(sig_05, aes(x = Fold, y = delta_activity)) +
  geom_pointdensity(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(sig_05$Fold, sig_05$delta_activity),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    y = expression(Delta * "Activity (KO - WT)"),
    title = sprintf("K119ub vs activity change at sig. enhancers (FDR < %.2f)",
                    FDR_PRIMARY)
  ) +
  theme_pub

save_plot(p05, "05_fdr05_activity_component", w = 7, h = 6)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 16 COMPLETE\n")
cat(sprintf("Outputs in: %s\n", OUTPUT_DIR))
cat("\nKey result — signal improvement after filtering:\n")
cat(sprintf("  All enhancers (n=%d):        rho = %+.4f\n",
            nrow(enh), rho_all))
cat(sprintf("  FDR < %.2f (n=%d):       rho = %+.4f\n",
            FDR_PRIMARY, nrow(sig_05), rho_filt))
cat("================================================================================\n")
