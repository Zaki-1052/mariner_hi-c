# scripts/tad_k119ub_correlation.R
#
# Correlate differential TADs with H2AK119ub status
# (Jesse Dixon suggestion: which TADs gaining density sit in ub-enriched regions?)
#
# Input:
#   - data/tsvs/figure_1_tads_boundaries_compartments/1B_late_tad_all_annotated.tsv
#   - data/upstream/diffbind/K119ub_diffbind_results_summit_appended_ap.txt
#
# Output:
#   - outputs/tad_k119ub_analysis/ (plots + annotated TSV + statistics)
#
# Usage:
#   Rscript scripts/tad_k119ub_correlation.R

# =============================================================================
# 1. SETUP
# =============================================================================

cat("=== TAD–K119ub Correlation Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(ggpubr)
})

source("scripts/utils/multi_format_output.R")

# =============================================================================
# 2. CONFIGURATION
# =============================================================================

TAD_FILE <- "data/tsvs/figure_1_tads_boundaries_compartments/1B_late_tad_all_annotated.tsv"
UB_FILE  <- "data/upstream/diffbind/K119ub_diffbind_results_summit_appended_ap.txt"

OUTPUT_DIR <- "outputs/tad_k119ub_analysis"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Significance thresholds
TAD_FDR_THRESHOLD <- 0.05
UB_FDR_THRESHOLD  <- 0.05

# Plot colors: increased (mutant-enriched) = red, decreased = blue
DIRECTION_COLORS <- c(
  "Increased_in_Mutant" = "#d73027",
  "Decreased_in_Mutant" = "#4575b4"
)

# =============================================================================
# 3. LOAD DATA
# =============================================================================

cat("Loading data...\n")

# TAD differential results
tad_df <- read.table(TAD_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  TADs: %d rows\n", nrow(tad_df)))

# K119ub DiffBind peaks
ub_df <- read.table(UB_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  K119ub peaks: %d rows\n", nrow(ub_df)))

# =============================================================================
# 4. CONVERT TO GRANGES
# =============================================================================

cat("Converting to GRanges...\n")

tad_gr <- GRanges(
  seqnames = tad_df$chr1,
  ranges   = IRanges(start = tad_df$start1, end = tad_df$end1),
  TAD_name   = tad_df$TAD_name,
  Difference = tad_df$Difference,
  pvalue     = tad_df$pvalue,
  adj_pvalue = tad_df$adj_pvalue,
  TAD_size   = tad_df$TAD_size,
  direction  = tad_df$direction
)

ub_gr <- GRanges(
  seqnames = ub_df$Peak_Chr,
  ranges   = IRanges(start = ub_df$Peak_Start, end = ub_df$Peak_End),
  Fold = ub_df$Fold,
  FDR  = ub_df$FDR
)

cat(sprintf("  TAD GRanges: %d\n", length(tad_gr)))
cat(sprintf("  K119ub GRanges: %d\n", length(ub_gr)))

# =============================================================================
# 5. COMPUTE K119UB METRICS PER TAD
# =============================================================================

cat("Computing K119ub overlap metrics per TAD...\n")

# Find all overlaps
hits <- findOverlaps(tad_gr, ub_gr)
cat(sprintf("  Total overlaps: %d\n", length(hits)))

# Build per-TAD summary from overlap hits
ub_metrics <- tibble(
  tad_idx = queryHits(hits),
  ub_fold = ub_gr$Fold[subjectHits(hits)],
  ub_fdr  = ub_gr$FDR[subjectHits(hits)]
) %>%
  group_by(tad_idx) %>%
  summarise(
    peak_count          = n(),
    mean_ub_fold        = mean(ub_fold),
    median_ub_fold      = median(ub_fold),
    n_sig_peaks         = sum(ub_fdr < UB_FDR_THRESHOLD),
    mean_sig_ub_fold    = ifelse(sum(ub_fdr < UB_FDR_THRESHOLD) > 0,
                                 mean(ub_fold[ub_fdr < UB_FDR_THRESHOLD]),
                                 NA_real_),
    .groups = "drop"
  )

# Merge back to full TAD data
tad_annotated <- tad_df %>%
  mutate(tad_idx = row_number()) %>%
  left_join(ub_metrics, by = "tad_idx") %>%
  mutate(
    peak_count     = replace_na(peak_count, 0),
    n_sig_peaks    = replace_na(n_sig_peaks, 0),
    peak_density   = peak_count / (TAD_size / 1e6),
    sig_peak_density = n_sig_peaks / (TAD_size / 1e6)
  ) %>%
  select(-tad_idx)

cat(sprintf("  TADs with >= 1 K119ub peak: %d / %d (%.1f%%)\n",
            sum(tad_annotated$peak_count > 0), nrow(tad_annotated),
            100 * mean(tad_annotated$peak_count > 0)))
cat(sprintf("  TADs with >= 1 sig K119ub peak: %d / %d (%.1f%%)\n",
            sum(tad_annotated$n_sig_peaks > 0), nrow(tad_annotated),
            100 * mean(tad_annotated$n_sig_peaks > 0)))

# =============================================================================
# 6. STATISTICAL TESTS
# =============================================================================

cat("\nStatistical tests...\n")

inc <- tad_annotated %>% filter(direction == "Increased_in_Mutant")
dec <- tad_annotated %>% filter(direction == "Decreased_in_Mutant")

# Peak density
wilcox_density <- wilcox.test(inc$peak_density, dec$peak_density)
cat(sprintf("  Peak density — Increased median: %.2f, Decreased median: %.2f, p = %.2e\n",
            median(inc$peak_density), median(dec$peak_density), wilcox_density$p.value))

# Sig peak density
wilcox_sig_density <- wilcox.test(inc$sig_peak_density, dec$sig_peak_density)
cat(sprintf("  Sig peak density — Increased median: %.2f, Decreased median: %.2f, p = %.2e\n",
            median(inc$sig_peak_density), median(dec$sig_peak_density), wilcox_sig_density$p.value))

# Mean ub fold change (only TADs with peaks)
inc_fold <- inc %>% filter(!is.na(mean_ub_fold))
dec_fold <- dec %>% filter(!is.na(mean_ub_fold))
wilcox_fold <- wilcox.test(inc_fold$mean_ub_fold, dec_fold$mean_ub_fold)
cat(sprintf("  Mean K119ub fold — Increased median: %.3f, Decreased median: %.3f, p = %.2e\n",
            median(inc_fold$mean_ub_fold), median(dec_fold$mean_ub_fold), wilcox_fold$p.value))

# Spearman correlation: Difference vs K119ub metrics
cor_density <- cor.test(tad_annotated$Difference, tad_annotated$peak_density, method = "spearman")
cor_fold <- cor.test(
  tad_annotated$Difference[!is.na(tad_annotated$mean_ub_fold)],
  tad_annotated$mean_ub_fold[!is.na(tad_annotated$mean_ub_fold)],
  method = "spearman"
)
cat(sprintf("  Spearman(Difference, peak_density): rho = %.3f, p = %.2e\n",
            cor_density$estimate, cor_density$p.value))
cat(sprintf("  Spearman(Difference, mean_ub_fold): rho = %.3f, p = %.2e\n",
            cor_fold$estimate, cor_fold$p.value))

# Significant TADs only
sig_tads <- tad_annotated %>% filter(adj_pvalue < TAD_FDR_THRESHOLD)
sig_inc <- sig_tads %>% filter(direction == "Increased_in_Mutant")
sig_dec <- sig_tads %>% filter(direction == "Decreased_in_Mutant")

wilcox_sig_tad_density <- wilcox.test(sig_inc$peak_density, sig_dec$peak_density)
cat(sprintf("\n  Significant TADs (FDR < %.2f): n = %d (Inc: %d, Dec: %d)\n",
            TAD_FDR_THRESHOLD, nrow(sig_tads), nrow(sig_inc), nrow(sig_dec)))
cat(sprintf("  Sig TAD peak density — Increased median: %.2f, Decreased median: %.2f, p = %.2e\n",
            median(sig_inc$peak_density), median(sig_dec$peak_density),
            wilcox_sig_tad_density$p.value))

# =============================================================================
# 7. PLOTS
# =============================================================================

cat("\nGenerating plots...\n")

# --- Helper: format p-value for annotation ---
format_pval <- function(p) {
  if (p < 1e-16) return("p < 1e-16")
  if (p < 0.00001) return(sprintf("p = %.2e", p))
  if (p < 0.001) return(sprintf("p = %.2e", p))
  sprintf("p = %.4f", p)
}

# --- Plot 1: K119ub peak density by TAD direction (all TADs) ---
p1 <- ggplot(tad_annotated, aes(x = direction, y = peak_density, fill = direction)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = DIRECTION_COLORS) +
  scale_x_discrete(labels = c("Increased_in_Mutant" = "Increased", "Decreased_in_Mutant" = "Decreased")) +
  annotate("text", x = 1.5, y = max(tad_annotated$peak_density) * 0.95,
           label = format_pval(wilcox_density$p.value), size = 4) +
  annotate("segment", x = 1, xend = 2,
           y = max(tad_annotated$peak_density) * 0.88,
           yend = max(tad_annotated$peak_density) * 0.88,
           color = "black", linewidth = 0.5) +
  labs(
    title = "K119ub Peak Density by TAD Interaction Change",
    subtitle = sprintf("All TADs (n = %d increased, %d decreased)", nrow(inc), nrow(dec)),
    x = "TAD Interaction Density",
    y = "K119ub peaks per Mb"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "none",
    panel.grid = element_blank()
  )
save_multiformat_ggplot(p1, file.path(OUTPUT_DIR, "tad_k119ub_peak_density_violin"), width = 6, height = 7)

# --- Plot 2: Mean K119ub fold change by TAD direction ---
fold_data <- tad_annotated %>% filter(!is.na(mean_ub_fold))

p2 <- ggplot(fold_data, aes(x = direction, y = mean_ub_fold, fill = direction)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_manual(values = DIRECTION_COLORS) +
  scale_x_discrete(labels = c("Increased_in_Mutant" = "Increased", "Decreased_in_Mutant" = "Decreased")) +
  annotate("text", x = 1.5, y = max(fold_data$mean_ub_fold) * 0.95,
           label = format_pval(wilcox_fold$p.value), size = 4) +
  annotate("segment", x = 1, xend = 2,
           y = max(fold_data$mean_ub_fold) * 0.88,
           yend = max(fold_data$mean_ub_fold) * 0.88,
           color = "black", linewidth = 0.5) +
  labs(
    title = "Mean K119ub Fold Change by TAD Interaction Change",
    subtitle = sprintf("TADs with >= 1 K119ub peak (n = %d increased, %d decreased)",
                       nrow(inc_fold), nrow(dec_fold)),
    x = "TAD Interaction Density",
    y = "Mean K119ub fold change (mut/ctrl)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "none",
    panel.grid = element_blank()
  )
save_multiformat_ggplot(p2, file.path(OUTPUT_DIR, "tad_k119ub_fold_change_violin"), width = 6, height = 7)

# --- Plot 3: Scatter — TAD Difference vs mean K119ub fold change ---
scatter_data <- tad_annotated %>%
  filter(!is.na(mean_ub_fold)) %>%
  mutate(significant = adj_pvalue < TAD_FDR_THRESHOLD)

p3 <- ggplot(scatter_data, aes(x = Difference, y = mean_ub_fold)) +
  geom_point(aes(color = significant), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_color_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "gray60"),
    labels = c("TRUE" = sprintf("FDR < %.2f", TAD_FDR_THRESHOLD), "FALSE" = "NS"),
    name = "TAD significance"
  ) +
  annotate("text", x = max(scatter_data$Difference) * 0.7,
           y = max(scatter_data$mean_ub_fold) * 0.9,
           label = sprintf("rho = %.3f\n%s", cor_fold$estimate, format_pval(cor_fold$p.value)),
           size = 3.5, hjust = 0) +
  labs(
    title = "TAD Interaction Change vs. K119ub Enrichment",
    x = "TAD Interaction Density Difference (mut - ctrl)",
    y = "Mean K119ub fold change (mut/ctrl)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = c(0.15, 0.9),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid = element_blank()
  )
save_multiformat_ggplot(p3, file.path(OUTPUT_DIR, "tad_k119ub_scatter"), width = 8, height = 7)

# --- Plot 4: Significant TADs only — peak density and fold change ---
if (nrow(sig_inc) >= 5 && nrow(sig_dec) >= 5) {
  sig_fold_data <- sig_tads %>% filter(!is.na(mean_ub_fold))
  sig_inc_fold <- sig_fold_data %>% filter(direction == "Increased_in_Mutant")
  sig_dec_fold <- sig_fold_data %>% filter(direction == "Decreased_in_Mutant")

  wilcox_sig_fold <- wilcox.test(sig_inc_fold$mean_ub_fold, sig_dec_fold$mean_ub_fold)

  p4a <- ggplot(sig_tads, aes(x = direction, y = peak_density, fill = direction)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    scale_x_discrete(labels = c("Increased_in_Mutant" = "Increased", "Decreased_in_Mutant" = "Decreased")) +
    annotate("text", x = 1.5, y = max(sig_tads$peak_density) * 0.95,
             label = format_pval(wilcox_sig_tad_density$p.value), size = 4) +
    annotate("segment", x = 1, xend = 2,
             y = max(sig_tads$peak_density) * 0.88,
             yend = max(sig_tads$peak_density) * 0.88,
             color = "black", linewidth = 0.5) +
    labs(
      title = "K119ub Peak Density — Significant TADs",
      subtitle = sprintf("FDR < %.2f (n = %d increased, %d decreased)",
                         TAD_FDR_THRESHOLD, nrow(sig_inc), nrow(sig_dec)),
      x = "TAD Interaction Density",
      y = "K119ub peaks per Mb"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )
  save_multiformat_ggplot(p4a, file.path(OUTPUT_DIR, "tad_k119ub_sig_tads_density_violin"), width = 6, height = 7)

  if (nrow(sig_inc_fold) >= 3 && nrow(sig_dec_fold) >= 3) {
    p4b <- ggplot(sig_fold_data, aes(x = direction, y = mean_ub_fold, fill = direction)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      scale_x_discrete(labels = c("Increased_in_Mutant" = "Increased", "Decreased_in_Mutant" = "Decreased")) +
      annotate("text", x = 1.5, y = max(sig_fold_data$mean_ub_fold) * 0.95,
               label = format_pval(wilcox_sig_fold$p.value), size = 4) +
      annotate("segment", x = 1, xend = 2,
               y = max(sig_fold_data$mean_ub_fold) * 0.88,
               yend = max(sig_fold_data$mean_ub_fold) * 0.88,
               color = "black", linewidth = 0.5) +
      labs(
        title = "Mean K119ub Fold Change — Significant TADs",
        subtitle = sprintf("FDR < %.2f (n = %d increased, %d decreased)",
                           TAD_FDR_THRESHOLD, nrow(sig_inc_fold), nrow(sig_dec_fold)),
        x = "TAD Interaction Density",
        y = "Mean K119ub fold change (mut/ctrl)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.grid = element_blank()
      )
    save_multiformat_ggplot(p4b, file.path(OUTPUT_DIR, "tad_k119ub_sig_tads_fold_violin"), width = 6, height = 7)
  }
} else {
  cat("  Skipping significant-only plots (too few TADs per group)\n")
}

# =============================================================================
# 8. SAVE OUTPUTS
# =============================================================================

cat("\nSaving outputs...\n")

# Annotated TSV
out_tsv <- file.path(OUTPUT_DIR, "tad_k119ub_annotated.tsv")
write.table(tad_annotated, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Wrote: %s (%d rows)\n", out_tsv, nrow(tad_annotated)))

# Statistics summary
stats_file <- file.path(OUTPUT_DIR, "tad_k119ub_statistics.txt")
sink(stats_file)
cat("=== TAD-K119ub Correlation Statistics ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("--- Data Summary ---\n")
cat(sprintf("Total TADs: %d\n", nrow(tad_annotated)))
cat(sprintf("  Increased_in_Mutant: %d\n", nrow(inc)))
cat(sprintf("  Decreased_in_Mutant: %d\n", nrow(dec)))
cat(sprintf("Total K119ub peaks: %d\n", length(ub_gr)))
cat(sprintf("TADs with >= 1 K119ub peak: %d (%.1f%%)\n",
            sum(tad_annotated$peak_count > 0), 100 * mean(tad_annotated$peak_count > 0)))
cat(sprintf("TADs with >= 1 significant K119ub peak: %d (%.1f%%)\n\n",
            sum(tad_annotated$n_sig_peaks > 0), 100 * mean(tad_annotated$n_sig_peaks > 0)))

cat("--- All TADs: K119ub Peak Density (peaks/Mb) ---\n")
cat(sprintf("  Increased median: %.2f\n", median(inc$peak_density)))
cat(sprintf("  Decreased median: %.2f\n", median(dec$peak_density)))
cat(sprintf("  Wilcoxon p = %.2e\n\n", wilcox_density$p.value))

cat("--- All TADs: Sig K119ub Peak Density (peaks/Mb) ---\n")
cat(sprintf("  Increased median: %.2f\n", median(inc$sig_peak_density)))
cat(sprintf("  Decreased median: %.2f\n", median(dec$sig_peak_density)))
cat(sprintf("  Wilcoxon p = %.2e\n\n", wilcox_sig_density$p.value))

cat("--- All TADs: Mean K119ub Fold Change (TADs with peaks) ---\n")
cat(sprintf("  Increased median: %.3f (n = %d)\n", median(inc_fold$mean_ub_fold), nrow(inc_fold)))
cat(sprintf("  Decreased median: %.3f (n = %d)\n", median(dec_fold$mean_ub_fold), nrow(dec_fold)))
cat(sprintf("  Wilcoxon p = %.2e\n\n", wilcox_fold$p.value))

cat("--- Spearman Correlations ---\n")
cat(sprintf("  Difference vs peak density: rho = %.3f, p = %.2e\n",
            cor_density$estimate, cor_density$p.value))
cat(sprintf("  Difference vs mean K119ub fold: rho = %.3f, p = %.2e\n\n",
            cor_fold$estimate, cor_fold$p.value))

cat(sprintf("--- Significant TADs (FDR < %.2f) ---\n", TAD_FDR_THRESHOLD))
cat(sprintf("  Total: %d (Increased: %d, Decreased: %d)\n", nrow(sig_tads), nrow(sig_inc), nrow(sig_dec)))
cat(sprintf("  Peak density — Increased median: %.2f, Decreased median: %.2f\n",
            median(sig_inc$peak_density), median(sig_dec$peak_density)))
cat(sprintf("  Wilcoxon p = %.2e\n", wilcox_sig_tad_density$p.value))
sink()
cat(sprintf("  Wrote: %s\n", stats_file))

cat(sprintf("\nDone. %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
