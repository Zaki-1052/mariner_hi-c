# biomodal/downstream/scripts/viz_sections/section_26_tet_ko_comparison.R
# Section 26: TET Triple-KO Comparison (TODO 14c)
# Standalone script - sources shared config for all dependencies and data
#
# Formal quantitative comparison of BAP1-KO demethylation ratio changes
# to published TET triple-KO (TKO) data from GSE166423 (Lopez-Moyado et al.).
# TET-KO should show near-complete loss of 5hmC; BAP1-KO shows partial reduction
# consistent with indirect TET impairment. Quantifies the attenuation factor.
#
# Data sources:
#   1. demethylation_ratio_all_genes.tsv (Section 22) — BAP1-KO delta_ratio
#   2. tet_ko_gene_signal.tsv (from preprocess_tet_ko_bigwig.R) — TET-KO ratios
#
# Figures:
#   26a: Overlaid WT ratio densities — BAP1 WT vs TET-KO WT
#   26b: Delta-ratio density comparison — BAP1-KO vs TET-KO
#   26c: QQ plot — BAP1 delta vs TET-KO delta quantiles
#   26d: Effect size comparison bars — Cliff's delta, % decreased, Cohen's d
#   26e: Per-gene scatter — BAP1 delta (x) vs TET-KO delta (y)
#   26f: Chromatin-stratified boxplots — BAP1 vs TET-KO delta per state
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_26_tet_ko_comparison.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(effsize)
})

# =============================================================================
# SECTION 26 CONFIGURATION
# =============================================================================

# Input tables
BAP1_RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
TET_KO_TABLE <- file.path(BASE_DIR, "data/tet_ko_gene_signal.tsv")

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 26: TET TRIPLE-KO COMPARISON (TODO 14c)\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("demethylation_ratio_all_genes.tsv not found" = file.exists(BAP1_RATIO_TABLE))
stopifnot("tet_ko_gene_signal.tsv not found — run preprocess_tet_ko_bigwig.R on HPC first" =
            file.exists(TET_KO_TABLE))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("--- Step 1: Loading data ---\n")

bap1_df <- read.table(BAP1_RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  BAP1-KO ratios: %d genes\n", nrow(bap1_df)))

tet_df <- read.table(TET_KO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  TET-KO signal: %d genes\n", nrow(tet_df)))

# =============================================================================
# STEP 2: MERGE BY GENE NAME
# =============================================================================

cat("\n--- Step 2: Merging by gene name ---\n")

merged <- dplyr::inner_join(
  bap1_df %>% dplyr::select(gene, chr, start, end,
                             ratio_ctrl_bap1 = ratio_ctrl,
                             ratio_mut_bap1 = ratio_mut,
                             delta_ratio_bap1 = delta_ratio,
                             chromatin_state),
  tet_df %>% dplyr::select(gene,
                            ratio_wt_tet = ratio_wt,
                            ratio_ko_tet = ratio_ko,
                            delta_ratio_tet),
  by = "gene"
) %>%
  dplyr::filter(is.finite(delta_ratio_bap1) & is.finite(delta_ratio_tet))

cat(sprintf("  Matched genes: %d\n", nrow(merged)))
stopifnot("Too few matched genes — check gene name format" = nrow(merged) > 5000)

cat(sprintf("  BAP1-KO delta median: %.4f, TET-KO delta median: %.4f\n",
            median(merged$delta_ratio_bap1), median(merged$delta_ratio_tet)))

# =============================================================================
# STEP 3: COMPUTE COMPARISON STATISTICS
# =============================================================================

cat("\n--- Step 3: Computing comparison statistics ---\n")

# Cliff's delta for BAP1
bap1_cliff <- cliff.delta(merged$delta_ratio_bap1, rep(0, nrow(merged)))
cat(sprintf("  BAP1-KO Cliff's delta (vs 0): %.3f (%s)\n",
            bap1_cliff$estimate, bap1_cliff$magnitude))

# Cliff's delta for TET-KO
tet_cliff <- cliff.delta(merged$delta_ratio_tet, rep(0, nrow(merged)))
cat(sprintf("  TET-KO Cliff's delta (vs 0): %.3f (%s)\n",
            tet_cliff$estimate, tet_cliff$magnitude))

# % decreased
bap1_pct_decreased <- 100 * mean(merged$delta_ratio_bap1 < 0)
tet_pct_decreased <- 100 * mean(merged$delta_ratio_tet < 0)
cat(sprintf("  BAP1-KO %% decreased: %.1f%%\n", bap1_pct_decreased))
cat(sprintf("  TET-KO %% decreased: %.1f%%\n", tet_pct_decreased))

# Cohen's d
bap1_cohen <- cohen.d(merged$delta_ratio_bap1, rep(0, nrow(merged)))
tet_cohen <- cohen.d(merged$delta_ratio_tet, rep(0, nrow(merged)))
cat(sprintf("  BAP1-KO Cohen's d: %.3f (%s)\n",
            bap1_cohen$estimate, bap1_cohen$magnitude))
cat(sprintf("  TET-KO Cohen's d: %.3f (%s)\n",
            tet_cohen$estimate, tet_cohen$magnitude))

# KS test between delta distributions
ks_result <- ks.test(merged$delta_ratio_bap1, merged$delta_ratio_tet)
cat(sprintf("  KS test (BAP1 vs TET deltas): D=%.4f, %s\n",
            ks_result$statistic, fmt_p(ks_result$p.value)))

# Attenuation factor
bap1_median_delta <- median(merged$delta_ratio_bap1)
tet_median_delta <- median(merged$delta_ratio_tet)
attenuation <- bap1_median_delta / tet_median_delta
cat(sprintf("  Attenuation factor (BAP1/TET median deltas): %.3f\n", attenuation))
cat(sprintf("  Interpretation: BAP1-KO reproduces %.1f%% of TET-KO delta effect\n",
            100 * attenuation))

# Per-gene Spearman correlation
gene_rho <- cor.test(merged$delta_ratio_bap1, merged$delta_ratio_tet, method = "spearman")
cat(sprintf("  Per-gene Spearman rho: %.3f (%s)\n",
            gene_rho$estimate, fmt_p(gene_rho$p.value)))

# =============================================================================
# FIGURE 26a: OVERLAID WT RATIO DENSITIES
# =============================================================================

cat("\n--- Figure 26a: WT ratio density comparison ---\n")

wt_density_df <- rbind(
  data.frame(ratio = merged$ratio_ctrl_bap1, study = "BAP1 WT", stringsAsFactors = FALSE),
  data.frame(ratio = merged$ratio_wt_tet, study = "TET-KO WT", stringsAsFactors = FALSE)
)

wt_medians <- wt_density_df %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(med = median(ratio, na.rm = TRUE), .groups = "drop")

study_colors <- c("BAP1 WT" = "#2166AC", "TET-KO WT" = "#4DAF4A")

p_26a <- ggplot(wt_density_df, aes(x = ratio, fill = study, color = study)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(data = wt_medians, aes(xintercept = med, color = study),
             linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = study_colors, name = "Study") +
  scale_color_manual(values = study_colors, name = "Study") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("BAP1 WT median: %.4f\nTET-KO WT median: %.4f\nN=%s matched genes",
                           wt_medians$med[wt_medians$study == "BAP1 WT"],
                           wt_medians$med[wt_medians$study == "TET-KO WT"],
                           format(nrow(merged), big.mark = ","))) +
  labs(
    title = "WT Demethylation Ratio Distributions",
    subtitle = "Baseline 5hmC/(5mC+5hmC) in wildtype: validates comparable reference states",
    x = "Demethylation efficiency ratio (5hmC/(5mC+5hmC))",
    y = "Density"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_26a, file.path(OUTPUT_DIR, "26a_wt_ratio_comparison"), 10, 7)

# =============================================================================
# FIGURE 26b: DELTA-RATIO DENSITY COMPARISON
# =============================================================================

cat("--- Figure 26b: Delta-ratio density comparison ---\n")

delta_density_df <- rbind(
  data.frame(delta = merged$delta_ratio_bap1, study = "BAP1-KO",
             stringsAsFactors = FALSE),
  data.frame(delta = merged$delta_ratio_tet, study = "TET Triple-KO",
             stringsAsFactors = FALSE)
)

delta_medians <- delta_density_df %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(med = median(delta, na.rm = TRUE), .groups = "drop")

delta_colors <- c("BAP1-KO" = "#B2182B", "TET Triple-KO" = "#4DAF4A")

p_26b <- ggplot(delta_density_df, aes(x = delta, fill = study, color = study)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_vline(data = delta_medians, aes(xintercept = med, color = study),
             linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = delta_colors, name = "Study") +
  scale_color_manual(values = delta_colors, name = "Study") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.5,
           label = sprintf(paste0(
             "BAP1-KO median: %+.4f\n",
             "TET-KO median: %+.4f\n",
             "KS D=%.4f, %s\n",
             "Attenuation: %.1f%%"
           ),
           delta_medians$med[delta_medians$study == "BAP1-KO"],
           delta_medians$med[delta_medians$study == "TET Triple-KO"],
           ks_result$statistic, fmt_p(ks_result$p.value),
           100 * attenuation)) +
  labs(
    title = "Delta-Ratio Distributions: BAP1-KO vs TET Triple-KO",
    subtitle = "TET-KO should show much stronger negative shift (near-complete TET loss)",
    x = "Delta-ratio (KO \u2212 WT)",
    y = "Density"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_26b, file.path(OUTPUT_DIR, "26b_delta_ratio_density"), 10, 7)

# =============================================================================
# FIGURE 26c: QQ PLOT (BAP1 VS TET-KO DELTA QUANTILES)
# =============================================================================

cat("--- Figure 26c: QQ plot ---\n")

n_quantiles <- min(nrow(merged), 5000)
probs <- seq(0, 1, length.out = n_quantiles)
qq_data <- data.frame(
  bap1_q = quantile(merged$delta_ratio_bap1, probs, na.rm = TRUE),
  tet_q = quantile(merged$delta_ratio_tet, probs, na.rm = TRUE)
)

# Fit slope through origin region (attenuation factor from quantiles)
qq_fit <- lm(bap1_q ~ 0 + tet_q, data = qq_data)
qq_slope <- coef(qq_fit)["tet_q"]

p_26c <- ggplot(qq_data, aes(x = tet_q, y = bap1_q)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#333333") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = qq_slope, intercept = 0, color = "#E41A1C", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("QQ slope = %.3f\n(BAP1 = %.1f%% of TET-KO shift)\nReference: slope=1 (identity)",
                           qq_slope, 100 * qq_slope)) +
  labs(
    title = "QQ Plot: BAP1-KO vs TET Triple-KO Delta-Ratios",
    subtitle = "Slope = attenuation factor; <1 means BAP1 effect is attenuated relative to TET-KO",
    x = "TET-KO delta-ratio quantiles",
    y = "BAP1-KO delta-ratio quantiles"
  ) +
  theme_biomodal() +
  coord_equal()

save_multiformat_ggplot(p_26c, file.path(OUTPUT_DIR, "26c_qq_plot"), 9, 9)

# =============================================================================
# FIGURE 26d: EFFECT SIZE COMPARISON BARS
# =============================================================================

cat("--- Figure 26d: Effect size comparison ---\n")

effect_df <- data.frame(
  metric = rep(c("Cliff's delta", "% Decreased", "Cohen's d"), each = 2),
  study = rep(c("BAP1-KO", "TET Triple-KO"), 3),
  value = c(
    abs(bap1_cliff$estimate), abs(tet_cliff$estimate),
    bap1_pct_decreased, tet_pct_decreased,
    abs(bap1_cohen$estimate), abs(tet_cohen$estimate)
  ),
  stringsAsFactors = FALSE
)
effect_df$metric <- factor(effect_df$metric,
                           levels = c("Cliff's delta", "% Decreased", "Cohen's d"))

p_26d <- ggplot(effect_df, aes(x = study, y = value, fill = study)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.5, size = 3.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = delta_colors) +
  labs(
    title = "Effect Size Comparison: BAP1-KO vs TET Triple-KO",
    subtitle = sprintf("BAP1-KO shows attenuated effect (~%.0f%% of TET-KO magnitude)",
                       100 * attenuation),
    x = NULL, y = "Effect Size"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))

save_multiformat_ggplot(p_26d, file.path(OUTPUT_DIR, "26d_effect_size_comparison"), 12, 6)

# =============================================================================
# FIGURE 26e: PER-GENE SCATTER (BAP1 DELTA vs TET-KO DELTA)
# =============================================================================

cat("--- Figure 26e: Per-gene scatter ---\n")

# Top genes to label
top_scatter <- merged %>%
  dplyr::mutate(combined_delta = abs(delta_ratio_bap1) + abs(delta_ratio_tet)) %>%
  dplyr::arrange(dplyr::desc(combined_delta)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::bind_rows(
    merged %>% dplyr::filter(gene %in% KEY_GENES)
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

p_26e <- ggplot(merged, aes(x = delta_ratio_bap1, y = delta_ratio_tet)) +
  geom_point(alpha = 0.08, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "#E41A1C", linewidth = 0.8, se = TRUE,
              fill = "#E41A1C22") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(data = top_scatter, aes(label = gene), size = 2.8,
                  max.overlaps = 20, segment.alpha = 0.5,
                  color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("Spearman rho = %.3f\n%s\nN = %s matched genes",
                           gene_rho$estimate, fmt_p(gene_rho$p.value),
                           format(nrow(merged), big.mark = ","))) +
  labs(
    title = "Per-Gene Delta-Ratio: BAP1-KO vs TET Triple-KO",
    subtitle = "Tests whether the same genes are affected in both conditions",
    x = "BAP1-KO delta-ratio",
    y = "TET Triple-KO delta-ratio"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_26e, file.path(OUTPUT_DIR, "26e_per_gene_scatter"), 10, 9)

# =============================================================================
# FIGURE 26f: CHROMATIN-STRATIFIED BOXPLOTS
# =============================================================================

cat("--- Figure 26f: Chromatin-stratified comparison ---\n")

chrom_long <- rbind(
  merged %>%
    dplyr::select(gene, chromatin_state, delta = delta_ratio_bap1) %>%
    dplyr::mutate(study = "BAP1-KO"),
  merged %>%
    dplyr::select(gene, chromatin_state, delta = delta_ratio_tet) %>%
    dplyr::mutate(study = "TET Triple-KO")
)
chrom_long$study <- factor(chrom_long$study, levels = c("BAP1-KO", "TET Triple-KO"))

# Per-state Wilcoxon (BAP1 vs TET within each state)
state_stats <- chrom_long %>%
  dplyr::group_by(chromatin_state) %>%
  dplyr::summarise(
    n = dplyr::n() / 2,
    wilcox_p = tryCatch(
      wilcox.test(delta[study == "BAP1-KO"],
                  delta[study == "TET Triple-KO"],
                  paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    bap1_median = median(delta[study == "BAP1-KO"], na.rm = TRUE),
    tet_median = median(delta[study == "TET Triple-KO"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(wilcox_p, method = "BH"),
    label = sprintf("n=%s\np_adj=%s", format(n, big.mark = ","),
                    ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj)))
  )

cat("  Per-state comparison (BAP1 vs TET-KO delta):\n")
for (i in seq_len(nrow(state_stats))) {
  cat(sprintf("    %-20s BAP1=%.4f TET=%.4f p_adj=%s (n=%d)\n",
              state_stats$chromatin_state[i],
              state_stats$bap1_median[i], state_stats$tet_median[i],
              fmt_p(state_stats$p_adj[i]), state_stats$n[i]))
}

p_26f <- ggplot(chrom_long, aes(x = study, y = delta, fill = study)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~chromatin_state, ncol = 4, scales = "free_y") +
  geom_text(data = state_stats,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.3, size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = delta_colors) +
  scale_x_discrete(labels = c("BAP1-KO" = "BAP1", "TET Triple-KO" = "TET")) +
  labs(
    title = "Delta-Ratio by Chromatin State: BAP1-KO vs TET Triple-KO",
    subtitle = "Paired Wilcoxon test (BH-corrected) per chromatin state",
    x = NULL, y = "Delta-ratio (KO \u2212 WT)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9))

save_multiformat_ggplot(p_26f, file.path(OUTPUT_DIR, "26f_chromatin_stratified"), 14, 10)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Matched genes with both deltas
export_matched <- merged %>%
  dplyr::select(gene, chr, start, end, chromatin_state,
                ratio_ctrl_bap1, ratio_mut_bap1, delta_ratio_bap1,
                ratio_wt_tet, ratio_ko_tet, delta_ratio_tet) %>%
  dplyr::arrange(delta_ratio_bap1)

write.table(export_matched,
            file.path(TABLES_DIR, "tet_ko_comparison_matched_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved tet_ko_comparison_matched_genes.tsv (%d rows)\n", nrow(export_matched)))

# Table 2: Summary statistics
summary_df <- data.frame(
  metric = c("N matched genes",
             "BAP1-KO median delta", "TET-KO median delta",
             "Attenuation factor (BAP1/TET)",
             "BAP1-KO Cliff's delta", "TET-KO Cliff's delta",
             "BAP1-KO % decreased", "TET-KO % decreased",
             "BAP1-KO Cohen's d", "TET-KO Cohen's d",
             "KS statistic", "KS p-value",
             "Per-gene Spearman rho", "Per-gene Spearman p",
             "QQ slope (attenuation)"),
  value = c(nrow(merged),
            sprintf("%.6f", bap1_median_delta), sprintf("%.6f", tet_median_delta),
            sprintf("%.4f", attenuation),
            sprintf("%.4f", bap1_cliff$estimate), sprintf("%.4f", tet_cliff$estimate),
            sprintf("%.1f", bap1_pct_decreased), sprintf("%.1f", tet_pct_decreased),
            sprintf("%.4f", bap1_cohen$estimate), sprintf("%.4f", tet_cohen$estimate),
            sprintf("%.4f", ks_result$statistic), fmt_p(ks_result$p.value),
            sprintf("%.4f", gene_rho$estimate), fmt_p(gene_rho$p.value),
            sprintf("%.4f", qq_slope)),
  stringsAsFactors = FALSE
)

write.table(summary_df,
            file.path(TABLES_DIR, "tet_ko_comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved tet_ko_comparison_summary.tsv (%d rows)\n", nrow(summary_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 26 SUMMARY: TET TRIPLE-KO COMPARISON\n")
cat("========================================================================\n\n")

cat(sprintf("Matched genes: %d\n", nrow(merged)))
cat(sprintf("\nWT baseline ratios:\n"))
cat(sprintf("  BAP1 WT median: %.4f\n", median(merged$ratio_ctrl_bap1)))
cat(sprintf("  TET-KO WT median: %.4f\n", median(merged$ratio_wt_tet)))

cat(sprintf("\nDelta-ratio comparison:\n"))
cat(sprintf("  BAP1-KO median delta: %+.4f (%.1f%% decreased)\n",
            bap1_median_delta, bap1_pct_decreased))
cat(sprintf("  TET-KO median delta:  %+.4f (%.1f%% decreased)\n",
            tet_median_delta, tet_pct_decreased))
cat(sprintf("  Attenuation: BAP1 = %.1f%% of TET-KO effect\n", 100 * attenuation))

cat(sprintf("\nEffect sizes:\n"))
cat(sprintf("  BAP1-KO: Cliff's delta=%.3f (%s), Cohen's d=%.3f (%s)\n",
            bap1_cliff$estimate, bap1_cliff$magnitude,
            bap1_cohen$estimate, bap1_cohen$magnitude))
cat(sprintf("  TET-KO:  Cliff's delta=%.3f (%s), Cohen's d=%.3f (%s)\n",
            tet_cliff$estimate, tet_cliff$magnitude,
            tet_cohen$estimate, tet_cohen$magnitude))

cat(sprintf("\nStatistical tests:\n"))
cat(sprintf("  KS test: D=%.4f, %s\n", ks_result$statistic, fmt_p(ks_result$p.value)))
cat(sprintf("  Per-gene Spearman: rho=%.3f, %s\n",
            gene_rho$estimate, fmt_p(gene_rho$p.value)))
cat(sprintf("  QQ slope: %.3f (%.1f%% attenuation)\n", qq_slope, 100 * qq_slope))

cat(sprintf("\nInterpretation:\n"))
if (attenuation > 0 && attenuation < 1) {
  cat(sprintf("  BAP1-KO shows ATTENUATED TET impairment (%.0f%% of TET-KO effect).\n",
              100 * attenuation))
  cat("  Consistent with INDIRECT TET blockade (PRC1-mediated), not direct TET loss.\n")
  cat("  The partial nature confirms BAP1 impairs TET access rather than TET expression.\n")
}
if (gene_rho$estimate > 0 && gene_rho$p.value < 0.05) {
  cat(sprintf("  Positive per-gene correlation (rho=%.3f) confirms CONVERGENT targets:\n",
              gene_rho$estimate))
  cat("  The same genes are preferentially affected in both conditions.\n")
}

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/26{a-f}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/tet_ko_comparison_*.tsv (2 files)\n", TABLES_DIR))

cat("\n=== Section 26 complete ===\n")
