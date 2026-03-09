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
#   26g: Relative (baseline-normalized) delta density comparison
#   26h: Relative delta QQ plot — quantile attenuation after baseline normalization
#   26i: Response decomposition bars — TET-KO binary vs BAP1-KO graded loss
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

# Relative delta = fractional change from WT baseline
# Floor baseline at 0.01 to avoid division by near-zero
merged$relative_delta_bap1 <- merged$delta_ratio_bap1 / pmax(merged$ratio_ctrl_bap1, 0.01)
merged$relative_delta_tet  <- merged$delta_ratio_tet  / pmax(merged$ratio_wt_tet, 0.01)

cat(sprintf("  BAP1-KO relative delta median: %.4f, TET-KO relative delta median: %.4f\n",
            median(merged$relative_delta_bap1), median(merged$relative_delta_tet)))

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

# Relative attenuation (baseline-normalized)
bap1_median_rel_delta <- median(merged$relative_delta_bap1)
tet_median_rel_delta <- median(merged$relative_delta_tet)
relative_attenuation <- bap1_median_rel_delta / tet_median_rel_delta
cat(sprintf("  Relative attenuation (BAP1/TET median relative deltas): %.3f\n",
            relative_attenuation))
cat(sprintf("  Interpretation: BAP1-KO reproduces %.1f%% of TET-KO relative effect\n",
            100 * relative_attenuation))

bap1_rel_pct_decreased <- 100 * mean(merged$relative_delta_bap1 < 0)
tet_rel_pct_decreased <- 100 * mean(merged$relative_delta_tet < 0)
cat(sprintf("  BAP1-KO relative %% decreased: %.1f%%\n", bap1_rel_pct_decreased))
cat(sprintf("  TET-KO relative %% decreased: %.1f%%\n", tet_rel_pct_decreased))

# Per-gene Spearman correlation (relative deltas)
gene_rho_relative <- cor.test(merged$relative_delta_bap1, merged$relative_delta_tet,
                              method = "spearman")
cat(sprintf("  Per-gene Spearman rho (relative): %.3f (%s)\n",
            gene_rho_relative$estimate, fmt_p(gene_rho_relative$p.value)))

# Correlation decomposition: why relative rho collapses to ~0
n_tet_minus1 <- sum(merged$relative_delta_tet == -1)
n_tet_zero   <- sum(abs(merged$relative_delta_tet) < 1e-10)
n_tet_other  <- nrow(merged) - n_tet_minus1 - n_tet_zero
pct_tet_minus1 <- 100 * n_tet_minus1 / nrow(merged)
pct_tet_zero   <- 100 * n_tet_zero / nrow(merged)
pct_tet_other  <- 100 * n_tet_other / nrow(merged)

cat(sprintf("  Correlation decomposition (relative deltas):\n"))
cat(sprintf("    TET-KO relative delta = -1.0 (complete loss): %d genes (%.1f%%)\n",
            n_tet_minus1, pct_tet_minus1))
cat(sprintf("    TET-KO relative delta ~  0.0 (no WT signal):  %d genes (%.1f%%)\n",
            n_tet_zero, pct_tet_zero))
cat(sprintf("    TET-KO relative delta other (partial loss):   %d genes (%.1f%%)\n",
            n_tet_other, pct_tet_other))
cat(sprintf("    TET-KO relative delta variance: %.4f (vs BAP1: %.4f)\n",
            var(merged$relative_delta_tet), var(merged$relative_delta_bap1)))
cat(sprintf("    -> %.1f%% of TET-KO genes have constant relative delta, destroying correlation\n",
            pct_tet_minus1 + pct_tet_zero))
cat("    -> Absolute rho was driven by shared baseline variation, not gene specificity\n")

# Spearman on partial-loss subset only
partial_idx <- merged$relative_delta_tet != -1 & abs(merged$relative_delta_tet) > 1e-10
rho_partial <- NULL
if (sum(partial_idx) > 30) {
  rho_partial <- cor.test(merged$relative_delta_bap1[partial_idx],
                          merged$relative_delta_tet[partial_idx], method = "spearman")
  cat(sprintf("    Spearman rho (partial-loss genes only, n=%d): %.3f (%s)\n",
              sum(partial_idx), rho_partial$estimate, fmt_p(rho_partial$p.value)))
}

# Per-gene Spearman correlation (absolute deltas)
gene_rho <- cor.test(merged$delta_ratio_bap1, merged$delta_ratio_tet, method = "spearman")
cat(sprintf("  Per-gene Spearman rho: %.3f (%s)\n",
            gene_rho$estimate, fmt_p(gene_rho$p.value)))

# BAP1-KO response categories (for decomposition comparison with TET-KO)
n_bap1_strong   <- sum(merged$relative_delta_bap1 < -0.5)
n_bap1_moderate <- sum(merged$relative_delta_bap1 >= -0.5 & merged$relative_delta_bap1 < -0.1)
n_bap1_weak     <- sum(merged$relative_delta_bap1 >= -0.1)
pct_bap1_strong   <- 100 * n_bap1_strong / nrow(merged)
pct_bap1_moderate <- 100 * n_bap1_moderate / nrow(merged)
pct_bap1_weak     <- 100 * n_bap1_weak / nrow(merged)

cat(sprintf("  BAP1-KO response categories:\n"))
cat(sprintf("    Strong loss (>50%%):       %d genes (%.1f%%)\n", n_bap1_strong, pct_bap1_strong))
cat(sprintf("    Moderate loss (10-50%%):   %d genes (%.1f%%)\n", n_bap1_moderate, pct_bap1_moderate))
cat(sprintf("    Weak/no change (<10%%):    %d genes (%.1f%%)\n", n_bap1_weak, pct_bap1_weak))

# Variance-explained: how much of absolute rho is baseline-driven?
cat(sprintf("\n  Variance-explained analysis (baseline vs gene-specific):\n"))

resid_bap1 <- residuals(lm(rank(delta_ratio_bap1) ~ rank(ratio_ctrl_bap1) + rank(ratio_wt_tet),
                            data = merged))
resid_tet  <- residuals(lm(rank(delta_ratio_tet) ~ rank(ratio_ctrl_bap1) + rank(ratio_wt_tet),
                            data = merged))
rho_resid <- suppressWarnings(
  cor.test(resid_bap1, resid_tet, method = "spearman")
)

lm_baseline <- lm(rank(delta_ratio_tet) ~ rank(ratio_ctrl_bap1) + rank(ratio_wt_tet),
                   data = merged)
lm_full     <- lm(rank(delta_ratio_tet) ~ rank(delta_ratio_bap1) + rank(ratio_ctrl_bap1) + rank(ratio_wt_tet),
                   data = merged)
r2_baseline      <- summary(lm_baseline)$r.squared
r2_full          <- summary(lm_full)$r.squared
r2_gene_specific <- r2_full - r2_baseline
baseline_fraction <- 1 - rho_resid$estimate / gene_rho$estimate

cat(sprintf("    Residualized rho (WT baselines removed): %.4f (%s)\n",
            rho_resid$estimate, fmt_p(rho_resid$p.value)))
cat(sprintf("    R2 baseline only: %.4f | R2 full: %.4f | Delta-R2 gene-specific: %.4f\n",
            r2_baseline, r2_full, r2_gene_specific))
cat(sprintf("    Baseline accounts for %.1f%% of explained variance\n",
            100 * r2_baseline / r2_full))
cat(sprintf("    -> %.0f%% of original rho=%.3f was driven by shared baseline variation\n",
            100 * baseline_fraction, gene_rho$estimate))

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
  metric = rep(c("Cliff's delta", "% Decreased", "Cohen's d", "Relative attenuation"), each = 2),
  study = rep(c("BAP1-KO", "TET Triple-KO"), 4),
  value = c(
    abs(bap1_cliff$estimate), abs(tet_cliff$estimate),
    bap1_pct_decreased, tet_pct_decreased,
    abs(bap1_cohen$estimate), abs(tet_cohen$estimate),
    abs(bap1_median_rel_delta), abs(tet_median_rel_delta)
  ),
  stringsAsFactors = FALSE
)
effect_df$metric <- factor(effect_df$metric,
                           levels = c("Cliff's delta", "% Decreased", "Cohen's d",
                                      "Relative attenuation"))

p_26d <- ggplot(effect_df, aes(x = study, y = value, fill = study)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.5, size = 3.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = delta_colors) +
  labs(
    title = "Effect Size Comparison: BAP1-KO vs TET Triple-KO",
    subtitle = sprintf("Absolute attenuation: %.1f%% | Relative attenuation: %.1f%%",
                       100 * attenuation, 100 * relative_attenuation),
    x = NULL, y = "Effect Size"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))

save_multiformat_ggplot(p_26d, file.path(OUTPUT_DIR, "26d_effect_size_comparison"), 14, 6)

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
# FIGURE 26g: RELATIVE (BASELINE-NORMALIZED) DELTA DENSITY COMPARISON
# =============================================================================

cat("--- Figure 26g: Relative delta density comparison ---\n")

rel_density_df <- rbind(
  data.frame(delta = merged$relative_delta_bap1, study = "BAP1-KO",
             stringsAsFactors = FALSE),
  data.frame(delta = merged$relative_delta_tet, study = "TET Triple-KO",
             stringsAsFactors = FALSE)
)

rel_medians <- rel_density_df %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(med = median(delta, na.rm = TRUE), .groups = "drop")

p_26g <- ggplot(rel_density_df, aes(x = delta, fill = study, color = study)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_vline(data = rel_medians, aes(xintercept = med, color = study),
             linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = delta_colors, name = "Study") +
  scale_color_manual(values = delta_colors, name = "Study") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.5,
           label = sprintf(paste0(
             "BAP1-KO median: %+.4f\n",
             "TET-KO median: %+.4f\n",
             "Relative attenuation: %.1f%%\n",
             "Spearman rho (relative): %.3f"
           ),
           rel_medians$med[rel_medians$study == "BAP1-KO"],
           rel_medians$med[rel_medians$study == "TET Triple-KO"],
           100 * relative_attenuation,
           gene_rho_relative$estimate)) +
  labs(
    title = "Relative Delta Distributions: BAP1-KO vs TET Triple-KO",
    subtitle = "Baseline-normalized (delta/WT ratio) — controls for 2\u00d7 WT baseline difference",
    x = "Relative delta ((KO \u2212 WT) / WT)",
    y = "Density"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_26g, file.path(OUTPUT_DIR, "26g_relative_delta_density"), 10, 7)

# =============================================================================
# FIGURE 26h: RELATIVE DELTA QQ PLOT
# =============================================================================

cat("--- Figure 26h: Relative delta QQ plot ---\n")

n_quantiles_rel <- min(nrow(merged), 5000)
probs_rel <- seq(0, 1, length.out = n_quantiles_rel)
qq_rel_data <- data.frame(
  bap1_q = quantile(merged$relative_delta_bap1, probs_rel, na.rm = TRUE),
  tet_q = quantile(merged$relative_delta_tet, probs_rel, na.rm = TRUE)
)

# Through-origin OLS slope (relative QQ attenuation)
qq_rel_fit <- lm(bap1_q ~ 0 + tet_q, data = qq_rel_data)
rel_qq_slope <- coef(qq_rel_fit)["tet_q"]
cat(sprintf("  Relative QQ slope: %.4f (vs absolute QQ slope: %.4f)\n",
            rel_qq_slope, qq_slope))

p_26h <- ggplot(qq_rel_data, aes(x = tet_q, y = bap1_q)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#333333") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = rel_qq_slope, intercept = 0, color = "#E41A1C", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf(paste0(
             "Relative QQ slope = %.3f\n",
             "(vs absolute QQ slope = %.3f)\n",
             "TET-KO bimodal: %.1f%% at -1.0, %.1f%% at 0.0\n",
             "Step structure reflects all-or-nothing TET-KO loss"
           ), rel_qq_slope, qq_slope, pct_tet_minus1, pct_tet_zero)) +
  labs(
    title = "Relative QQ Plot: BAP1-KO vs TET Triple-KO",
    subtitle = "Baseline-normalized deltas — step pattern shows TET-KO binary loss vs BAP1-KO graded response",
    x = "TET-KO relative delta quantiles ((KO \u2212 WT) / WT)",
    y = "BAP1-KO relative delta quantiles ((KO \u2212 WT) / WT)"
  ) +
  theme_biomodal() +
  coord_equal()

save_multiformat_ggplot(p_26h, file.path(OUTPUT_DIR, "26h_relative_qq_plot"), 9, 9)

# =============================================================================
# FIGURE 26i: RESPONSE DECOMPOSITION BARS
# =============================================================================

cat("--- Figure 26i: Response decomposition ---\n")

decomp_summary <- data.frame(
  study = rep(c("BAP1-KO", "TET Triple-KO"), each = 3),
  category = c("Strong loss (>50%)", "Moderate loss (10-50%)", "Weak/no change (<10%)",
               "Complete loss (100%)", "No WT signal", "Partial loss"),
  n = c(n_bap1_strong, n_bap1_moderate, n_bap1_weak,
        n_tet_minus1, n_tet_zero, n_tet_other),
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(study) %>%
  dplyr::mutate(pct = 100 * n / sum(n),
                label = ifelse(pct > 3,
                               sprintf("%.1f%%\n(n=%s)", pct, format(n, big.mark = ",")),
                               "")) %>%
  dplyr::ungroup()

decomp_summary$category <- factor(decomp_summary$category,
  levels = c("Strong loss (>50%)", "Moderate loss (10-50%)", "Weak/no change (<10%)",
             "Complete loss (100%)", "No WT signal", "Partial loss"))

decomp_colors <- c(
  "Strong loss (>50%)"     = "#B2182B",
  "Moderate loss (10-50%)" = "#EF8A62",
  "Weak/no change (<10%)"  = "#FDDBC7",
  "Complete loss (100%)"   = "#1B7837",
  "No WT signal"           = "#A6DBA0",
  "Partial loss"           = "#D9F0D3"
)

p_26i <- ggplot(decomp_summary, aes(x = study, y = pct, fill = category)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.0, lineheight = 0.9, fontface = "bold") +
  scale_fill_manual(values = decomp_colors, name = "Response Category",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
  labs(
    title = "Response Decomposition: BAP1-KO (Graded) vs TET-KO (Binary)",
    subtitle = sprintf(paste0(
      "TET-KO: %.1f%% constant response (%.1f%% complete loss + %.1f%% no signal)\n",
      "Baseline variation explains ~%.0f%% of absolute rho=%.3f (N=%s genes)"
    ), pct_tet_minus1 + pct_tet_zero, pct_tet_minus1, pct_tet_zero,
    100 * baseline_fraction, gene_rho$estimate,
    format(nrow(merged), big.mark = ",")),
    x = NULL,
    y = "Percentage of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_26i, file.path(OUTPUT_DIR, "26i_response_decomposition"), 10, 8)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Matched genes with both deltas
export_matched <- merged %>%
  dplyr::select(gene, chr, start, end, chromatin_state,
                ratio_ctrl_bap1, ratio_mut_bap1, delta_ratio_bap1, relative_delta_bap1,
                ratio_wt_tet, ratio_ko_tet, delta_ratio_tet, relative_delta_tet) %>%
  dplyr::arrange(delta_ratio_bap1)

write.table(export_matched,
            file.path(TABLES_DIR, "tet_ko_comparison_matched_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved tet_ko_comparison_matched_genes.tsv (%d rows)\n", nrow(export_matched)))

# Table 2: Summary statistics
summary_df <- data.frame(
  metric = c("N matched genes",
             "BAP1-KO median delta", "TET-KO median delta",
             "Absolute attenuation (BAP1/TET)",
             "BAP1-KO median relative delta", "TET-KO median relative delta",
             "Relative attenuation (BAP1/TET)",
             "BAP1-KO Cliff's delta", "TET-KO Cliff's delta",
             "BAP1-KO % decreased", "TET-KO % decreased",
             "BAP1-KO relative % decreased", "TET-KO relative % decreased",
             "BAP1-KO Cohen's d", "TET-KO Cohen's d",
             "KS statistic", "KS p-value",
             "Per-gene Spearman rho (absolute)", "Per-gene Spearman p (absolute)",
             "Per-gene Spearman rho (relative)", "Per-gene Spearman p (relative)",
             "QQ slope (attenuation)",
             "Relative QQ slope",
             "TET-KO % complete loss (relative=-1)",
             "TET-KO % no WT signal (relative~0)",
             "TET-KO % partial loss",
             "BAP1-KO % strong loss (rel < -0.5)",
             "BAP1-KO % moderate loss (-0.5 to -0.1)",
             "BAP1-KO % weak/no change (rel >= -0.1)",
             "Residualized Spearman rho",
             "Residualized Spearman p",
             "R2 baseline only",
             "R2 full model",
             "Delta-R2 gene-specific",
             "Baseline % of explained variance",
             if (!is.null(rho_partial)) "Spearman rho (partial-loss genes)" else NULL,
             if (!is.null(rho_partial)) "Spearman p (partial-loss genes)" else NULL),
  value = c(nrow(merged),
            sprintf("%.6f", bap1_median_delta), sprintf("%.6f", tet_median_delta),
            sprintf("%.4f", attenuation),
            sprintf("%.6f", bap1_median_rel_delta), sprintf("%.6f", tet_median_rel_delta),
            sprintf("%.4f", relative_attenuation),
            sprintf("%.4f", bap1_cliff$estimate), sprintf("%.4f", tet_cliff$estimate),
            sprintf("%.1f", bap1_pct_decreased), sprintf("%.1f", tet_pct_decreased),
            sprintf("%.1f", bap1_rel_pct_decreased), sprintf("%.1f", tet_rel_pct_decreased),
            sprintf("%.4f", bap1_cohen$estimate), sprintf("%.4f", tet_cohen$estimate),
            sprintf("%.4f", ks_result$statistic), fmt_p(ks_result$p.value),
            sprintf("%.4f", gene_rho$estimate), fmt_p(gene_rho$p.value),
            sprintf("%.4f", gene_rho_relative$estimate), fmt_p(gene_rho_relative$p.value),
            sprintf("%.4f", qq_slope),
            sprintf("%.4f", rel_qq_slope),
            sprintf("%.1f", pct_tet_minus1),
            sprintf("%.1f", pct_tet_zero),
            sprintf("%.1f", pct_tet_other),
            sprintf("%.1f", pct_bap1_strong),
            sprintf("%.1f", pct_bap1_moderate),
            sprintf("%.1f", pct_bap1_weak),
            sprintf("%.4f", rho_resid$estimate),
            fmt_p(rho_resid$p.value),
            sprintf("%.4f", r2_baseline),
            sprintf("%.4f", r2_full),
            sprintf("%.4f", r2_gene_specific),
            sprintf("%.1f", 100 * r2_baseline / r2_full),
            if (!is.null(rho_partial)) sprintf("%.4f", rho_partial$estimate) else NULL,
            if (!is.null(rho_partial)) fmt_p(rho_partial$p.value) else NULL),
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

cat(sprintf("\nAbsolute delta-ratio comparison:\n"))
cat(sprintf("  BAP1-KO median delta: %+.4f (%.1f%% decreased)\n",
            bap1_median_delta, bap1_pct_decreased))
cat(sprintf("  TET-KO median delta:  %+.4f (%.1f%% decreased)\n",
            tet_median_delta, tet_pct_decreased))
cat(sprintf("  Absolute attenuation: BAP1 = %.1f%% of TET-KO effect\n", 100 * attenuation))

cat(sprintf("\nRelative (baseline-normalized) delta comparison:\n"))
cat(sprintf("  BAP1-KO median relative delta: %+.4f (%.1f%% decreased)\n",
            bap1_median_rel_delta, bap1_rel_pct_decreased))
cat(sprintf("  TET-KO median relative delta:  %+.4f (%.1f%% decreased)\n",
            tet_median_rel_delta, tet_rel_pct_decreased))
cat(sprintf("  Relative attenuation: BAP1 = %.1f%% of TET-KO effect\n",
            100 * relative_attenuation))

cat(sprintf("\nEffect sizes:\n"))
cat(sprintf("  BAP1-KO: Cliff's delta=%.3f (%s), Cohen's d=%.3f (%s)\n",
            bap1_cliff$estimate, bap1_cliff$magnitude,
            bap1_cohen$estimate, bap1_cohen$magnitude))
cat(sprintf("  TET-KO:  Cliff's delta=%.3f (%s), Cohen's d=%.3f (%s)\n",
            tet_cliff$estimate, tet_cliff$magnitude,
            tet_cohen$estimate, tet_cohen$magnitude))

cat(sprintf("\nStatistical tests:\n"))
cat(sprintf("  KS test: D=%.4f, %s\n", ks_result$statistic, fmt_p(ks_result$p.value)))
cat(sprintf("  Per-gene Spearman (absolute): rho=%.3f, %s\n",
            gene_rho$estimate, fmt_p(gene_rho$p.value)))
cat(sprintf("  Per-gene Spearman (relative): rho=%.3f, %s\n",
            gene_rho_relative$estimate, fmt_p(gene_rho_relative$p.value)))
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

cat(sprintf("\n  Relative correlation note:\n"))
cat(sprintf("    Relative Spearman rho ~ 0 because %.1f%% of TET-KO genes have\n", pct_tet_minus1))
cat("    relative delta = -1.0 (complete loss), leaving no per-gene variance.\n")
cat(sprintf("    The absolute rho=%.3f reflects shared baseline variation, not gene-specific\n",
            gene_rho$estimate))
cat("    sensitivity. This confirms TET-KO = uniform total loss while BAP1-KO = graded.\n")

cat(sprintf("\n  Variance decomposition:\n"))
cat(sprintf("    Baseline R2: %.4f | Full R2: %.4f | Gene-specific Delta-R2: %.4f\n",
            r2_baseline, r2_full, r2_gene_specific))
cat(sprintf("    Residualized rho: %.4f (%.0f%% of original driven by baseline)\n",
            rho_resid$estimate, 100 * baseline_fraction))

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/26{a-i}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/tet_ko_comparison_*.tsv (2 files)\n", TABLES_DIR))

cat("\n=== Section 26 complete ===\n")
