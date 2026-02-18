# abc/scripts/step10_k119ub_abc_correlation.R
# Correlate ΔABC scores with H2AK119ub ChIP-seq signal at enhancers
#
# BAP1 is an H2AK119ub deubiquitinase. In BAP1-KO, enhancers that lose ABC
# activity should show *increased* K119ub (negative correlation between
# delta_ABC and K119ub log2FC).
#
# Inputs:
#   abc/results/delta_abc_all_pairs.tsv      — 180K E-G pairs with ABC deltas
#   abc/results/k119ub_enhancer_signal.tsv   — per-replicate K119ub at enhancers (from HPC)
#   abc/results/delta_abc_annotated.tsv      — for H3K27ac stratification
#
# Outputs:
#   abc/results/figures/k119ub_correlation/  — 8 visualization panels (PDF + PNG)
#   abc/results/k119ub_abc_correlation_summary.tsv — correlation results
#   abc/results/k119ub_abc_enhancer_merged.tsv     — joined enhancer-level data
#
# Usage (local):
#   Rscript abc/scripts/step10_k119ub_abc_correlation.R

# =============================================================================
# CONFIGURATION
# =============================================================================

ABC_PAIRS_FILE   <- "results/delta_abc_all_pairs.tsv"
K119UB_FILE      <- "results/k119ub_enhancer_signal.tsv"
ANNOTATED_FILE   <- "results/delta_abc_annotated.tsv"
FIGURE_DIR       <- "results/figures/k119ub_correlation"
RESULTS_DIR      <- "results"

cat("================================================================================\n")
cat("STEP 10: K119UB-ABC ENHANCER CORRELATION ANALYSIS\n")
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
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

stopifnot(file.exists(ABC_PAIRS_FILE))
stopifnot(file.exists(K119UB_FILE))
stopifnot(file.exists(ANNOTATED_FILE))

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

k119 <- read.table(K119UB_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  K119ub signal: %d enhancers\n", nrow(k119)))

annot <- read.table(ANNOTATED_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Annotated pairs: %d rows\n", nrow(annot)))

# =============================================================================
# PART A: ENHANCER-LEVEL AGGREGATION
# =============================================================================

cat("\n--- Part A: Enhancer-level aggregation ---\n")

# delta_activity is identical for all genes from same enhancer, so take first
# delta_ABC and delta_unnorm vary per gene, so aggregate
abc$enh_key <- paste0(abc$chr, ":", abc$start, "-", abc$end)

enh_agg <- do.call(rbind, lapply(split(abc, abc$enh_key), function(df) {
  data.frame(
    enh_key          = df$enh_key[1],
    chr              = df$chr[1],
    start            = df$start[1],
    end              = df$end[1],
    delta_activity   = df$activity_base_KO[1] - df$activity_base_WT[1],
    activity_base_WT = df$activity_base_WT[1],
    activity_base_KO = df$activity_base_KO[1],
    mean_delta_abc   = mean(df$delta_ABC),
    mean_delta_unnorm = mean(df$delta_unnorm),
    max_abs_delta_unnorm = max(abs(df$delta_unnorm)),
    n_genes          = nrow(df),
    stringsAsFactors = FALSE
  )
}))
rownames(enh_agg) <- NULL

cat(sprintf("  Aggregated to %d unique enhancers\n", nrow(enh_agg)))

# H3K27ac annotation: take per-enhancer (all genes from same enhancer share it)
annot$enh_key <- paste0(annot$chr, ":", annot$start, "-", annot$end)
h3k27ac_map <- unique(annot[, c("enh_key", "has_H3K27ac")])
# Deduplicate: if any True for an enhancer, call it True
h3k27ac_agg <- tapply(h3k27ac_map$has_H3K27ac, h3k27ac_map$enh_key, function(x) any(x == "True" | x == TRUE))
h3k27ac_df <- data.frame(
  enh_key = names(h3k27ac_agg),
  has_H3K27ac = as.logical(h3k27ac_agg),
  stringsAsFactors = FALSE
)

enh_agg <- merge(enh_agg, h3k27ac_df, by = "enh_key", all.x = TRUE)

# Join with K119ub signal
k119$enh_key <- paste0(k119$chr, ":", k119$start, "-", k119$end)
enh_merged <- merge(enh_agg, k119, by = "enh_key", suffixes = c("", ".k119"))

cat(sprintf("  After joining with K119ub signal: %d enhancers\n", nrow(enh_merged)))

# Filter to quantifiable signal
enh_q <- enh_merged[enh_merged$signal_class == "quantifiable", ]
cat(sprintf("  Quantifiable K119ub signal: %d enhancers\n", nrow(enh_q)))

# Also add distance info per enhancer (median distance across gene targets)
dist_agg <- tapply(abc$distance, abc$enh_key, median)
enh_merged$median_distance <- dist_agg[enh_merged$enh_key]
enh_q$median_distance <- dist_agg[enh_q$enh_key]

# Save merged enhancer data
merged_out <- file.path(RESULTS_DIR, "k119ub_abc_enhancer_merged.tsv")
write.table(enh_merged, merged_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved merged data: %s (%d rows)\n", merged_out, nrow(enh_merged)))

# =============================================================================
# PART B: CORE CORRELATIONS
# =============================================================================

cat("\n--- Part B: Core Spearman correlations ---\n")

run_spearman <- function(x, y, label) {
  valid <- !is.na(x) & !is.na(y)
  n <- sum(valid)
  if (n < 10) {
    cat(sprintf("  %-45s: n=%d (too few)\n", label, n))
    return(data.frame(label = label, rho = NA, p_value = NA, n = n,
                      ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE))
  }
  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  # Fisher z-transform CI for Spearman
  z <- atanh(test$estimate)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1.96, 1.96) * se)
  cat(sprintf("  %-45s: rho = %+.4f [%+.4f, %+.4f], p = %.2e, n = %d\n",
              label, test$estimate, ci[1], ci[2], test$p.value, n))
  data.frame(
    label = label, rho = as.numeric(test$estimate),
    p_value = test$p.value, n = n,
    ci_lo = ci[1], ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

corr_results <- list()

corr_results[[1]] <- run_spearman(
  enh_q$delta_activity, enh_q$log2fc,
  "delta_activity vs K119ub_log2fc"
)

corr_results[[2]] <- run_spearman(
  enh_q$mean_delta_unnorm, enh_q$log2fc,
  "mean_delta_unnorm vs K119ub_log2fc"
)

corr_results[[3]] <- run_spearman(
  enh_q$mean_delta_abc, enh_q$log2fc,
  "mean_delta_ABC vs K119ub_log2fc"
)

corr_results[[4]] <- run_spearman(
  enh_q$delta_activity, enh_q$mut_mean,
  "delta_activity vs K119ub_mut_mean"
)

corr_results[[5]] <- run_spearman(
  enh_q$ctrl_mean, enh_q$mut_mean,
  "K119ub_ctrl_mean vs K119ub_mut_mean (QC)"
)

corr_table <- do.call(rbind, corr_results)
cat("\n")

# =============================================================================
# PART C: CATEGORY COMPARISON
# =============================================================================

cat("--- Part C: Category comparison ---\n")

enh_q$abc_category <- ifelse(
  enh_q$mean_delta_abc < -0.01, "Lost",
  ifelse(enh_q$mean_delta_abc > 0.01, "Gained", "Unchanged")
)
enh_q$abc_category <- factor(enh_q$abc_category, levels = c("Lost", "Unchanged", "Gained"))

cat("  Category counts:\n")
cat_counts <- table(enh_q$abc_category)
for (cat_name in names(cat_counts)) {
  cat(sprintf("    %-12s: %d\n", cat_name, cat_counts[cat_name]))
}

# Kruskal-Wallis
kw <- kruskal.test(log2fc ~ abc_category, data = enh_q)
cat(sprintf("\n  Kruskal-Wallis: chi-sq = %.2f, df = %d, p = %.2e\n",
            kw$statistic, kw$parameter, kw$p.value))

# Pairwise Wilcoxon
pairs <- list(
  c("Lost", "Unchanged"),
  c("Gained", "Unchanged"),
  c("Lost", "Gained")
)

wilcox_results <- list()
for (pr in pairs) {
  g1 <- enh_q$log2fc[enh_q$abc_category == pr[1]]
  g2 <- enh_q$log2fc[enh_q$abc_category == pr[2]]
  wt <- wilcox.test(g1, g2)

  # Rank-biserial correlation as effect size
  n1 <- length(g1)
  n2 <- length(g2)
  r_rb <- 1 - (2 * wt$statistic) / (n1 * n2)

  label <- paste(pr, collapse = " vs ")
  cat(sprintf("  Wilcoxon %s: p = %.2e, r_rb = %+.4f (n = %d vs %d)\n",
              label, wt$p.value, r_rb, n1, n2))

  wilcox_results[[label]] <- data.frame(
    comparison = label, p_value = wt$p.value,
    r_rb = as.numeric(r_rb), n1 = n1, n2 = n2,
    stringsAsFactors = FALSE
  )
}

# Median K119ub log2fc per category
cat("\n  Median K119ub log2FC by category:\n")
for (cat_name in levels(enh_q$abc_category)) {
  vals <- enh_q$log2fc[enh_q$abc_category == cat_name]
  cat(sprintf("    %-12s: median = %+.4f, mean = %+.4f, n = %d\n",
              cat_name, median(vals), mean(vals), length(vals)))
}

# =============================================================================
# PART D: STRATIFIED ANALYSES
# =============================================================================

cat("\n--- Part D: Stratified analyses ---\n")

# D1: By distance bin
cat("\n  D1: By distance bin\n")
dist_breaks <- c(0, 50e3, 200e3, 500e3, 1e6, 5e6, Inf)
dist_labels <- c("<50kb", "50-200kb", "200-500kb", "500kb-1Mb", "1-5Mb", ">5Mb")
enh_q$dist_bin <- cut(enh_q$median_distance, breaks = dist_breaks, labels = dist_labels,
                      include.lowest = TRUE, right = FALSE)

dist_corrs <- list()
for (bin in levels(enh_q$dist_bin)) {
  sub <- enh_q[enh_q$dist_bin == bin, ]
  dist_corrs[[bin]] <- run_spearman(
    sub$delta_activity, sub$log2fc,
    paste0("delta_activity vs K119ub [", bin, "]")
  )
}
dist_corr_table <- do.call(rbind, dist_corrs)

# D2: By H3K27ac status
cat("\n  D2: By H3K27ac status\n")
for (has_mark in c(TRUE, FALSE)) {
  sub <- enh_q[!is.na(enh_q$has_H3K27ac) & enh_q$has_H3K27ac == has_mark, ]
  label <- if (has_mark) "H3K27ac+" else "H3K27ac-"

  run_spearman(sub$delta_activity, sub$log2fc,
               paste0("delta_activity vs K119ub [", label, "]"))
  run_spearman(sub$mean_delta_unnorm, sub$log2fc,
               paste0("mean_delta_unnorm vs K119ub [", label, "]"))
}

# =============================================================================
# SAVE CORRELATION SUMMARY
# =============================================================================

summary_out <- file.path(RESULTS_DIR, "k119ub_abc_correlation_summary.tsv")
write.table(corr_table, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved correlation summary: %s\n", summary_out))

# =============================================================================
# PART E: VISUALIZATIONS
# =============================================================================

cat("\n--- Part E: Visualizations ---\n")

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    plot.title = element_text(size = 12, face = "bold")
  )

save_plot <- function(p, name, w = 7, h = 6) {
  pdf_path <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIGURE_DIR, paste0(name, ".png"))
  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf)
  ggsave(png_path, p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

format_rho <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  sprintf("rho = %+.3f\np = %.1e\nn = %s",
          test$estimate, test$p.value, format(sum(valid), big.mark = ","))
}

# --- Panel 1: delta_activity vs K119ub log2FC ---
p1 <- ggplot(enh_q, aes(x = delta_activity, y = log2fc)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_q$delta_activity, enh_q$log2fc),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "Enhancer activity change vs K119ub change"
  ) +
  theme_pub
save_plot(p1, "01_delta_activity_vs_k119ub")

# --- Panel 2: mean delta_unnorm vs K119ub log2FC ---
p2 <- ggplot(enh_q, aes(x = mean_delta_unnorm, y = log2fc)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_q$mean_delta_unnorm, enh_q$log2fc),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("Mean " * Delta * "(A" %*% "C) per enhancer"),
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "Unnormalized ABC change vs K119ub change"
  ) +
  theme_pub
save_plot(p2, "02_delta_unnorm_vs_k119ub")

# --- Panel 3: mean delta_ABC vs K119ub log2FC ---
p3 <- ggplot(enh_q, aes(x = mean_delta_abc, y = log2fc)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_q$mean_delta_abc, enh_q$log2fc),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("Mean " * Delta * "ABC per enhancer"),
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "Normalized ABC change vs K119ub change"
  ) +
  theme_pub
save_plot(p3, "03_delta_abc_vs_k119ub")

# --- Panel 4: Boxplot K119ub log2FC by ABC category ---
cat_colors <- c(Lost = "#2166AC", Unchanged = "grey60", Gained = "#B2182B")

# Compute p-values for bracket annotations
wilcox_p <- sapply(pairs, function(pr) {
  g1 <- enh_q$log2fc[enh_q$abc_category == pr[1]]
  g2 <- enh_q$log2fc[enh_q$abc_category == pr[2]]
  wilcox.test(g1, g2)$p.value
})

p4 <- ggplot(enh_q, aes(x = abc_category, y = log2fc, fill = abc_category)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "ABC category (based on mean ΔABC)",
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "K119ub change by ABC enhancer category",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kw$p.value)
  ) +
  theme_pub
save_plot(p4, "04_boxplot_k119ub_by_abc_category")

# --- Panel 5: Boxplot K119ub mut signal by ABC category ---
p5 <- ggplot(enh_q, aes(x = abc_category, y = mut_mean, fill = abc_category)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "ABC category (based on mean ΔABC)",
    y = "K119ub signal in KO (normalized, log1p scale)",
    title = "Absolute K119ub level in KO by ABC category"
  ) +
  theme_pub
save_plot(p5, "05_boxplot_k119ub_mut_by_category")

# --- Panel 6: Violin + box overlay ---
p6 <- ggplot(enh_q, aes(x = abc_category, y = log2fc, fill = abc_category)) +
  geom_violin(alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "ABC category (based on mean ΔABC)",
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "K119ub change distribution by ABC category"
  ) +
  theme_pub
save_plot(p6, "06_violin_k119ub_by_category")

# --- Panel 7: Faceted scatter by distance bin ---
enh_q_dist <- enh_q[!is.na(enh_q$dist_bin), ]

# Pre-compute rho per bin for annotation
rho_by_bin <- do.call(rbind, lapply(levels(enh_q_dist$dist_bin), function(bin) {
  sub <- enh_q_dist[enh_q_dist$dist_bin == bin, ]
  valid <- !is.na(sub$delta_activity) & !is.na(sub$log2fc)
  n <- sum(valid)
  if (n < 10) return(data.frame(dist_bin = bin, label = paste0("n=", n), stringsAsFactors = FALSE))
  test <- cor.test(sub$delta_activity[valid], sub$log2fc[valid], method = "spearman", exact = FALSE)
  data.frame(
    dist_bin = bin,
    label = sprintf("rho=%+.3f, n=%s", test$estimate, format(n, big.mark = ",")),
    stringsAsFactors = FALSE
  )
}))
rho_by_bin$dist_bin <- factor(rho_by_bin$dist_bin, levels = dist_labels)

p7 <- ggplot(enh_q_dist, aes(x = delta_activity, y = log2fc)) +
  geom_pointdensity(size = 0.2, alpha = 0.5) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.4, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text(data = rho_by_bin, aes(label = label),
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
            size = 2.8, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~dist_bin, nrow = 1, scales = "free_x") +
  labs(
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression("K119ub log"[2] * "FC"),
    title = "Activity-K119ub correlation by enhancer-gene distance"
  ) +
  theme_pub +
  theme(strip.text = element_text(size = 9))
save_plot(p7, "07_scatter_by_distance", w = 14, h = 4)

# --- Panel 8: Faceted scatter by H3K27ac status ---
enh_q_h3k <- enh_q[!is.na(enh_q$has_H3K27ac), ]
enh_q_h3k$h3k_label <- ifelse(enh_q_h3k$has_H3K27ac, "H3K27ac+", "H3K27ac-")

rho_by_h3k <- do.call(rbind, lapply(c("H3K27ac+", "H3K27ac-"), function(lbl) {
  sub <- enh_q_h3k[enh_q_h3k$h3k_label == lbl, ]
  valid <- !is.na(sub$delta_activity) & !is.na(sub$log2fc)
  n <- sum(valid)
  if (n < 10) return(data.frame(h3k_label = lbl, label = paste0("n=", n), stringsAsFactors = FALSE))
  test <- cor.test(sub$delta_activity[valid], sub$log2fc[valid], method = "spearman", exact = FALSE)
  data.frame(
    h3k_label = lbl,
    label = sprintf("rho=%+.3f, p=%.1e, n=%s", test$estimate, test$p.value, format(n, big.mark = ",")),
    stringsAsFactors = FALSE
  )
}))

p8 <- ggplot(enh_q_h3k, aes(x = delta_activity, y = log2fc)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = rho_by_h3k, aes(label = label),
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
            size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~h3k_label, nrow = 1) +
  labs(
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression("K119ub log"[2] * "FC (Mut / Ctrl)"),
    title = "Activity-K119ub correlation by H3K27ac status"
  ) +
  theme_pub
save_plot(p8, "08_scatter_by_h3k27ac", w = 10, h = 5)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("  Enhancers analyzed:    %d (quantifiable K119ub)\n", nrow(enh_q)))
cat(sprintf("  Figures saved to:      %s/\n", FIGURE_DIR))
cat(sprintf("  Correlation summary:   %s\n", summary_out))
cat(sprintf("  Merged enhancer data:  %s\n", merged_out))

cat("\n  Core correlations:\n")
for (i in seq_len(nrow(corr_table))) {
  cat(sprintf("    %-45s: rho = %+.4f (p = %.2e)\n",
              corr_table$label[i], corr_table$rho[i], corr_table$p_value[i]))
}

cat("\n  Biological prediction: Negative rho for ΔABC vs ΔK119ub\n")
cat("  (Lost enhancer activity → gained K119ub in BAP1-KO)\n")
cat("\nDone.\n")
