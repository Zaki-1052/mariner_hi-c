# abc/scripts/step15_abc_component_vs_expression.R
# Decompose ΔABC into activity and contact components vs gene expression.
#
# Jesse Dixon (Salk) suggested separating delta_activity and delta_contact to
# show the ΔABC-expression correlation is mostly activity-driven. The existing
# Fig 4A shows summed Δ(A×C) vs log2FC; this script produces two companion
# panels: sum(delta_activity) vs log2FC and sum(delta_contact) vs log2FC.
#
# Inputs (relative to abc/ working directory):
#   ../data/tsvs/figure_4_abc_analysis/4A_delta_abc_with_rnaseq.tsv
#
# Outputs:
#   results/figures/abc_component_expression/  — 5 plots (PDF+SVG+JPG)
#   results/figures/abc_component_expression/component_correlation_summary.tsv
#   results/figures/abc_component_expression/gene_level_components.tsv
#
# Usage:
#   cd abc && Rscript scripts/step15_abc_component_vs_expression.R

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE <- "../data/tsvs/figure_4_abc_analysis/4A_delta_abc_with_rnaseq.tsv"
OUTPUT_DIR <- "results/figures/abc_component_expression"
MULTIFORMAT_UTIL <- "../scripts/utils/multi_format_output.R"

DE_PADJ_CUTOFF <- 0.05

cat("================================================================================\n")
cat("STEP 15: ABC COMPONENT DECOMPOSITION vs GENE EXPRESSION\n")
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
                      pearson_r = NA, pearson_p = NA, n = n,
                      ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE))
  }
  sp_test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  pe_test <- cor.test(x[valid], y[valid], method = "pearson")
  z <- atanh(sp_test$estimate)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1.96, 1.96) * se)
  cat(sprintf("  %-50s: rho = %+.4f [%+.4f, %+.4f], r = %+.4f, n = %d\n",
              label, sp_test$estimate, ci[1], ci[2], pe_test$estimate, n))
  data.frame(
    label = label,
    spearman_rho = as.numeric(sp_test$estimate),
    spearman_p = sp_test$p.value,
    pearson_r = as.numeric(pe_test$estimate),
    pearson_p = pe_test$p.value,
    n = n,
    ci_lo = ci[1], ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

format_stats <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  sp <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  pe <- cor.test(x[valid], y[valid], method = "pearson")
  sprintf("rho = %+.3f\nr = %+.3f\nn = %s",
          sp$estimate, pe$estimate, format(sum(valid), big.mark = ","))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
pairs <- read.table(INPUT_FILE, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
cat(sprintf("  E-G pairs: %d rows\n", nrow(pairs)))

# =============================================================================
# COMPUTE PER-PAIR DELTAS
# =============================================================================

cat("\nComputing per-pair deltas...\n")
pairs$delta_activity <- pairs$activity_base_KO - pairs$activity_base_WT
pairs$delta_contact  <- pairs$hic_contact_pl_scaled_adj_KO -
                        pairs$hic_contact_pl_scaled_adj_WT

cat(sprintf("  delta_activity range: [%.4f, %.4f]\n",
            min(pairs$delta_activity), max(pairs$delta_activity)))
cat(sprintf("  delta_contact  range: [%.6f, %.6f]\n",
            min(pairs$delta_contact), max(pairs$delta_contact)))

# =============================================================================
# GENE-LEVEL AGGREGATION
# =============================================================================

cat("\nAggregating to gene level...\n")

gene_agg <- do.call(rbind, lapply(split(pairs, pairs$TargetGene), function(df) {
  data.frame(
    TargetGene         = df$TargetGene[1],
    log2FC             = df$log2FC[1],
    padj               = df$padj[1],
    baseMean           = df$baseMean[1],
    sum_delta_activity = sum(df$delta_activity),
    sum_delta_contact  = sum(df$delta_contact),
    sum_delta_unnorm   = sum(df$delta_unnorm),
    sum_delta_abc      = sum(df$delta_ABC),
    n_enhancers        = nrow(df),
    stringsAsFactors   = FALSE
  )
}))
rownames(gene_agg) <- NULL

cat(sprintf("  Total genes: %d\n", nrow(gene_agg)))

gene_agg$is_de <- !is.na(gene_agg$padj) & gene_agg$padj < DE_PADJ_CUTOFF
n_de <- sum(gene_agg$is_de)
cat(sprintf("  DE genes (padj < %.2f): %d\n", DE_PADJ_CUTOFF, n_de))

de_genes <- gene_agg[gene_agg$is_de, ]

# Save gene-level table
gene_out <- file.path(OUTPUT_DIR, "gene_level_components.tsv")
write.table(gene_agg, gene_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", gene_out))

# =============================================================================
# CORRELATIONS
# =============================================================================

cat("\n--- Correlation analysis ---\n")

corr_results <- list()

cat("\n  All genes:\n")
corr_results[[1]] <- run_spearman(
  gene_agg$sum_delta_unnorm, gene_agg$log2FC,
  "sum_delta_unnorm vs log2FC [all]"
)
corr_results[[2]] <- run_spearman(
  gene_agg$sum_delta_activity, gene_agg$log2FC,
  "sum_delta_activity vs log2FC [all]"
)
corr_results[[3]] <- run_spearman(
  gene_agg$sum_delta_contact, gene_agg$log2FC,
  "sum_delta_contact vs log2FC [all]"
)

cat("\n  DE genes only:\n")
corr_results[[4]] <- run_spearman(
  de_genes$sum_delta_unnorm, de_genes$log2FC,
  "sum_delta_unnorm vs log2FC [DE]"
)
corr_results[[5]] <- run_spearman(
  de_genes$sum_delta_activity, de_genes$log2FC,
  "sum_delta_activity vs log2FC [DE]"
)
corr_results[[6]] <- run_spearman(
  de_genes$sum_delta_contact, de_genes$log2FC,
  "sum_delta_contact vs log2FC [DE]"
)

corr_table <- do.call(rbind, corr_results)

corr_out <- file.path(OUTPUT_DIR, "component_correlation_summary.tsv")
write.table(corr_table, corr_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", corr_out))

# =============================================================================
# PLOT 01: ACTIVITY COMPONENT vs EXPRESSION (DE genes)
# =============================================================================

cat("\n--- Plot 01: Activity component vs expression ---\n")

p01 <- ggplot(de_genes, aes(x = sum_delta_activity, y = log2FC)) +
  geom_pointdensity(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_stats(de_genes$sum_delta_activity, de_genes$log2FC),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Sigma * Delta * "Activity per gene (KO - WT)"),
    y = expression("Gene expression log"[2] * "FC"),
    title = "Activity component vs gene expression"
  ) +
  theme_pub

save_plot(p01, "01_activity_vs_expression", w = 7, h = 6)

# =============================================================================
# PLOT 02: CONTACT COMPONENT vs EXPRESSION (DE genes)
# =============================================================================

cat("\n--- Plot 02: Contact component vs expression ---\n")

p02 <- ggplot(de_genes, aes(x = sum_delta_contact, y = log2FC)) +
  geom_pointdensity(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_stats(de_genes$sum_delta_contact, de_genes$log2FC),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Sigma * Delta * "Contact per gene (KO - WT)"),
    y = expression("Gene expression log"[2] * "FC"),
    title = "Contact component vs gene expression"
  ) +
  theme_pub

save_plot(p02, "02_contact_vs_expression", w = 7, h = 6)

# =============================================================================
# PLOT 03: SIDE-BY-SIDE COMPOSITE (the key figure)
# =============================================================================

cat("\n--- Plot 03: Side-by-side composite ---\n")

# Extract rho values for panel labels
rho_act <- corr_table$spearman_rho[corr_table$label ==
  "sum_delta_activity vs log2FC [DE]"]
rho_con <- corr_table$spearman_rho[corr_table$label ==
  "sum_delta_contact vs log2FC [DE]"]

p03a <- ggplot(de_genes, aes(x = sum_delta_activity, y = log2FC)) +
  geom_pointdensity(size = 0.4, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_stats(de_genes$sum_delta_activity, de_genes$log2FC),
           hjust = 1.1, vjust = 1.3, size = 3.2, fontface = "bold") +
  labs(
    x = expression(Sigma * Delta * "Activity per gene"),
    y = expression("Gene expression log"[2] * "FC"),
    title = sprintf("Activity (rho = %+.3f)", rho_act)
  ) +
  theme_pub +
  theme(legend.position = "none")

p03b <- ggplot(de_genes, aes(x = sum_delta_contact, y = log2FC)) +
  geom_pointdensity(size = 0.4, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_stats(de_genes$sum_delta_contact, de_genes$log2FC),
           hjust = 1.1, vjust = 1.3, size = 3.2, fontface = "bold") +
  labs(
    x = expression(Sigma * Delta * "Contact per gene"),
    y = expression("Gene expression log"[2] * "FC"),
    title = sprintf("Contact (rho = %+.3f)", rho_con)
  ) +
  theme_pub +
  theme(legend.position = "none")

p03 <- p03a + p03b +
  plot_annotation(
    title = "ABC decomposition: Activity vs Contact components",
    subtitle = sprintf("DE genes (padj < %.2f, n = %s)",
                       DE_PADJ_CUTOFF,
                       format(n_de, big.mark = ",")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(p03, "03_composite_activity_vs_contact", w = 13, h = 6)

# =============================================================================
# PLOT 04: BASELINE — COMBINED Δ(A×C) vs EXPRESSION
# =============================================================================

cat("\n--- Plot 04: Baseline delta_unnorm vs expression ---\n")

p04 <- ggplot(de_genes, aes(x = sum_delta_unnorm, y = log2FC)) +
  geom_pointdensity(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_stats(de_genes$sum_delta_unnorm, de_genes$log2FC),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Sigma * Delta * "(Activity " %*% " Contact) per gene"),
    y = expression("Gene expression log"[2] * "FC"),
    title = expression("Combined " * Delta * "(A " %*% " C) vs gene expression (baseline)")
  ) +
  theme_pub

save_plot(p04, "04_baseline_unnorm_vs_expression", w = 7, h = 6)

# =============================================================================
# PLOT 05: BAR CHART COMPARING CORRELATION STRENGTHS
# =============================================================================

cat("\n--- Plot 05: Correlation comparison bar chart ---\n")

bar_data <- corr_table[grepl("\\[DE\\]", corr_table$label), ]
bar_data$metric <- c("Combined\ndelta(AxC)", "Activity\nonly", "Contact\nonly")
bar_data$metric <- factor(bar_data$metric,
                          levels = c("Activity\nonly", "Contact\nonly",
                                     "Combined\ndelta(AxC)"))

p05 <- ggplot(bar_data, aes(x = metric, y = spearman_rho, fill = metric)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.5) +
  geom_text(aes(label = sprintf("rho = %+.3f", spearman_rho)),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Activity\nonly" = "#D95F02",
                                "Contact\nonly" = "#7570B3",
                                "Combined\ndelta(AxC)" = "#1B9E77")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = NULL,
    y = expression("Spearman " * rho * " with gene expression log"[2] * "FC"),
    title = "Correlation strength: Activity vs Contact vs Combined",
    subtitle = sprintf("DE genes (n = %s) | Error bars = 95%% CI",
                       format(n_de, big.mark = ","))
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 11))

save_plot(p05, "05_correlation_comparison", w = 6, h = 5.5)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 15 COMPLETE\n")
cat(sprintf("Outputs in: %s\n", OUTPUT_DIR))
cat("\nKey result (DE genes):\n")
cat(sprintf("  Activity component:  rho = %+.4f\n", rho_act))
cat(sprintf("  Contact component:   rho = %+.4f\n", rho_con))
rho_comb <- corr_table$spearman_rho[corr_table$label ==
  "sum_delta_unnorm vs log2FC [DE]"]
cat(sprintf("  Combined delta(AxC): rho = %+.4f\n", rho_comb))
cat("================================================================================\n")
