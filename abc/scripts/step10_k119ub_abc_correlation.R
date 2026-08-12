# abc/scripts/step10_k119ub_abc_correlation.R
# Correlate ΔABC scores with H2AK119ub DiffBind results at enhancers
#
# BAP1 is an H2AK119ub deubiquitinase. In BAP1-KO, enhancers that lose ABC
# activity should show *increased* K119ub (negative correlation between
# delta_ABC and K119ub Fold change).
#
# Now uses DiffBind output directly — replicate-aware normalization, proper
# statistical testing (FDR per peak), no HPC preprocessing needed.
#
# Inputs:
#   abc/results/delta_abc_all_pairs.tsv                          — 180K E-G pairs with ABC deltas
#   abc/K119ub_allATAC_diffbind_results_summit_appended_ap.txt   — DiffBind K119ub results (75K peaks)
#   abc/results/delta_abc_annotated.tsv                          — for H3K27ac stratification
#
# Outputs:
#   abc/results/figures/k119ub_correlation/  — 12 visualization panels (PDF + SVG + JPG)
#   abc/results/k119ub_abc_correlation_summary.tsv — correlation results
#   abc/results/k119ub_abc_enhancer_merged.tsv     — joined enhancer-level data
#
# Usage (local, from abc/ directory):
#   Rscript scripts/step10_k119ub_abc_correlation.R

# =============================================================================
# CONFIGURATION
# =============================================================================

ABC_PAIRS_FILE   <- "results/delta_abc_all_pairs.tsv"
DIFFBIND_FILE    <- "K119ub_allATAC_diffbind_results_summit_appended_ap.txt"
ANNOTATED_FILE   <- "results/delta_abc_annotated.tsv"
FIGURE_DIR       <- "results/figures/k119ub_correlation"
RESULTS_DIR      <- "results"
MULTIFORMAT_UTIL <- "../scripts/utils/multi_format_output.R"

cat("================================================================================\n")
cat("STEP 10: K119UB-ABC ENHANCER CORRELATION ANALYSIS (DiffBind)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggpointdensity)
})

# Load multi-format output utility
stopifnot(file.exists(MULTIFORMAT_UTIL))
source(MULTIFORMAT_UTIL)
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

stopifnot(file.exists(ABC_PAIRS_FILE))
stopifnot(file.exists(DIFFBIND_FILE))
stopifnot(file.exists(ANNOTATED_FILE))

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

diffbind <- read.table(DIFFBIND_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  DiffBind K119ub peaks: %d rows\n", nrow(diffbind)))

annot <- read.table(ANNOTATED_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Annotated pairs: %d rows\n", nrow(annot)))

# Validate DiffBind columns
required_cols <- c("Peak_Chr", "Peak_Start", "Peak_End", "Fold", "FDR", "Conc", "Conc_mut", "Conc_ctrl")
missing_cols <- setdiff(required_cols, colnames(diffbind))
if (length(missing_cols) > 0) {
  stop(sprintf("DiffBind file missing columns: %s", paste(missing_cols, collapse = ", ")))
}

# =============================================================================
# PART A: COORDINATE MATCHING VIA GENOMICRANGES
# =============================================================================

cat("\n--- Part A: Coordinate matching (GenomicRanges overlap) ---\n")

# Build unique enhancer set from ABC pairs
abc$enh_key <- paste0(abc$chr, ":", abc$start, "-", abc$end)
enh_coords <- unique(abc[, c("chr", "start", "end", "enh_key")])
cat(sprintf("  Unique ABC enhancers: %d\n", nrow(enh_coords)))

# Build GRanges for ABC enhancers
gr_abc <- GRanges(
  seqnames = enh_coords$chr,
  ranges   = IRanges(start = enh_coords$start, end = enh_coords$end)
)
gr_abc$enh_key <- enh_coords$enh_key

# Build GRanges for DiffBind peaks (using Peak coordinates for best coverage)
gr_db <- GRanges(
  seqnames = diffbind$Peak_Chr,
  ranges   = IRanges(start = diffbind$Peak_Start, end = diffbind$Peak_End)
)
gr_db$db_idx <- seq_len(nrow(diffbind))

# Find overlaps
hits <- findOverlaps(gr_abc, gr_db)
cat(sprintf("  Raw overlaps found: %d\n", length(hits)))

# Resolve multi-matches: select peak with maximum basepair overlap
hit_df <- data.frame(
  abc_idx = queryHits(hits),
  db_idx  = subjectHits(hits),
  stringsAsFactors = FALSE
)

# Compute overlap width for each hit
overlap_widths <- width(pintersect(gr_abc[hit_df$abc_idx], gr_db[hit_df$db_idx]))
hit_df$overlap_width <- overlap_widths

# For each ABC enhancer, keep the DiffBind peak with maximum overlap
hit_df <- hit_df[order(hit_df$abc_idx, -hit_df$overlap_width), ]
hit_best <- hit_df[!duplicated(hit_df$abc_idx), ]

# Count statistics
n_matched <- nrow(hit_best)
n_unmatched <- nrow(enh_coords) - n_matched
n_multi <- sum(table(hit_df$abc_idx) > 1)
match_rate <- n_matched / nrow(enh_coords) * 100

cat(sprintf("  Matched enhancers:   %d (%.1f%%)\n", n_matched, match_rate))
cat(sprintf("  Unmatched enhancers: %d\n", n_unmatched))
cat(sprintf("  Multi-match (resolved by max overlap): %d\n", n_multi))

stopifnot(n_matched > 30000)

# Build matched lookup: enh_key → DiffBind row
match_map <- data.frame(
  enh_key   = gr_abc$enh_key[hit_best$abc_idx],
  Fold      = diffbind$Fold[hit_best$db_idx],
  FDR       = diffbind$FDR[hit_best$db_idx],
  Conc      = diffbind$Conc[hit_best$db_idx],
  Conc_mut  = diffbind$Conc_mut[hit_best$db_idx],
  Conc_ctrl = diffbind$Conc_ctrl[hit_best$db_idx],
  p.value   = diffbind$p.value[hit_best$db_idx],
  stringsAsFactors = FALSE
)

# =============================================================================
# PART B: ENHANCER-LEVEL AGGREGATION
# =============================================================================

cat("\n--- Part B: Enhancer-level aggregation ---\n")

# delta_activity is identical for all genes from same enhancer, so take first
# delta_ABC and delta_unnorm vary per gene, so aggregate
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

# Join with DiffBind matched data (inner join — only matched enhancers)
enh_merged <- merge(enh_agg, match_map, by = "enh_key")

cat(sprintf("  After joining with DiffBind signal: %d enhancers\n", nrow(enh_merged)))

# Add distance info per enhancer (median distance across gene targets)
dist_agg <- tapply(abc$distance, abc$enh_key, median)
enh_merged$median_distance <- dist_agg[enh_merged$enh_key]

# Save merged enhancer data
merged_out <- file.path(RESULTS_DIR, "k119ub_abc_enhancer_merged.tsv")
write.table(enh_merged, merged_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved merged data: %s (%d rows)\n", merged_out, nrow(enh_merged)))

# =============================================================================
# PART C: CORE CORRELATIONS
# =============================================================================

cat("\n--- Part C: Core Spearman correlations ---\n")

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
  enh_merged$delta_activity, enh_merged$Fold,
  "delta_activity vs K119ub_Fold"
)

corr_results[[2]] <- run_spearman(
  enh_merged$mean_delta_unnorm, enh_merged$Fold,
  "mean_delta_unnorm vs K119ub_Fold"
)

corr_results[[3]] <- run_spearman(
  enh_merged$mean_delta_abc, enh_merged$Fold,
  "mean_delta_ABC vs K119ub_Fold"
)

corr_results[[4]] <- run_spearman(
  enh_merged$delta_activity, enh_merged$Conc_mut,
  "delta_activity vs K119ub_Conc_mut"
)

corr_results[[5]] <- run_spearman(
  enh_merged$Conc_ctrl, enh_merged$Conc_mut,
  "K119ub_Conc_ctrl vs K119ub_Conc_mut (QC)"
)

corr_table <- do.call(rbind, corr_results)
cat("\n")

# =============================================================================
# PART D: CATEGORY COMPARISON
# =============================================================================

cat("--- Part D: Category comparison ---\n")

enh_merged$abc_category <- ifelse(
  enh_merged$mean_delta_abc < -0.01, "Lost",
  ifelse(enh_merged$mean_delta_abc > 0.01, "Gained", "Unchanged")
)
enh_merged$abc_category <- factor(enh_merged$abc_category, levels = c("Lost", "Unchanged", "Gained"))

cat("  Category counts:\n")
cat_counts <- table(enh_merged$abc_category)
for (cat_name in names(cat_counts)) {
  cat(sprintf("    %-12s: %d\n", cat_name, cat_counts[cat_name]))
}

# Kruskal-Wallis
kw <- kruskal.test(Fold ~ abc_category, data = enh_merged)
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
  g1 <- enh_merged$Fold[enh_merged$abc_category == pr[1]]
  g2 <- enh_merged$Fold[enh_merged$abc_category == pr[2]]
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

# Median K119ub Fold per category
cat("\n  Median K119ub Fold by category:\n")
for (cat_name in levels(enh_merged$abc_category)) {
  vals <- enh_merged$Fold[enh_merged$abc_category == cat_name]
  cat(sprintf("    %-12s: median = %+.4f, mean = %+.4f, n = %d\n",
              cat_name, median(vals), mean(vals), length(vals)))
}

# =============================================================================
# PART E: STRATIFIED ANALYSES
# =============================================================================

cat("\n--- Part E: Stratified analyses ---\n")

# E1: By distance bin
cat("\n  E1: By distance bin\n")
dist_breaks <- c(0, 50e3, 200e3, 500e3, 1e6, 5e6, Inf)
dist_labels <- c("<50kb", "50-200kb", "200-500kb", "500kb-1Mb", "1-5Mb", ">5Mb")
enh_merged$dist_bin <- cut(enh_merged$median_distance, breaks = dist_breaks, labels = dist_labels,
                           include.lowest = TRUE, right = FALSE)

dist_corrs <- list()
for (bin in levels(enh_merged$dist_bin)) {
  sub <- enh_merged[enh_merged$dist_bin == bin, ]
  dist_corrs[[bin]] <- run_spearman(
    sub$delta_activity, sub$Fold,
    paste0("delta_activity vs K119ub [", bin, "]")
  )
}
dist_corr_table <- do.call(rbind, dist_corrs)

# E2: By H3K27ac status
cat("\n  E2: By H3K27ac status\n")
for (has_mark in c(TRUE, FALSE)) {
  sub <- enh_merged[!is.na(enh_merged$has_H3K27ac) & enh_merged$has_H3K27ac == has_mark, ]
  label <- if (has_mark) "H3K27ac+" else "H3K27ac-"

  run_spearman(sub$delta_activity, sub$Fold,
               paste0("delta_activity vs K119ub [", label, "]"))
  run_spearman(sub$mean_delta_unnorm, sub$Fold,
               paste0("mean_delta_unnorm vs K119ub [", label, "]"))
}

# =============================================================================
# PART F: DIFFBIND FDR-BASED ANALYSES
# =============================================================================

cat("\n--- Part F: DiffBind FDR-based analyses ---\n")

# Classify enhancers by K119ub significance
enh_merged$k119ub_sig <- ifelse(
  enh_merged$FDR < 0.05 & enh_merged$Fold > 0, "Sig_Up",
  ifelse(enh_merged$FDR < 0.05 & enh_merged$Fold < 0, "Sig_Down", "NS")
)
enh_merged$k119ub_sig <- factor(enh_merged$k119ub_sig, levels = c("Sig_Down", "NS", "Sig_Up"))

cat("  K119ub significance classes:\n")
sig_counts <- table(enh_merged$k119ub_sig)
for (sig_name in names(sig_counts)) {
  cat(sprintf("    %-12s: %d (%.1f%%)\n", sig_name, sig_counts[sig_name],
              sig_counts[sig_name] / nrow(enh_merged) * 100))
}

# Contingency table: K119ub significance × ABC category
ct <- table(enh_merged$abc_category, enh_merged$k119ub_sig)
cat("\n  Contingency table (ABC category × K119ub significance):\n")
print(ct)

# Fisher's exact test
fisher_res <- fisher.test(ct, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n  Fisher's exact test: p = %.4e\n", fisher_res$p.value))

# Log2 observed/expected ratios for heatmap
row_totals <- rowSums(ct)
col_totals <- colSums(ct)
grand_total <- sum(ct)
expected <- outer(row_totals, col_totals) / grand_total
log2_oe <- log2((ct + 0.5) / (expected + 0.5))
cat("\n  Log2(O/E) matrix:\n")
print(round(log2_oe, 3))

# Stratified correlation by K119ub significance
cat("\n  Correlations stratified by K119ub significance:\n")
for (sig_class in levels(enh_merged$k119ub_sig)) {
  sub <- enh_merged[enh_merged$k119ub_sig == sig_class, ]
  run_spearman(sub$delta_activity, sub$Fold,
               paste0("delta_activity vs Fold [", sig_class, "]"))
}

# =============================================================================
# SAVE CORRELATION SUMMARY
# =============================================================================

summary_out <- file.path(RESULTS_DIR, "k119ub_abc_correlation_summary.tsv")
write.table(corr_table, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved correlation summary: %s\n", summary_out))

# =============================================================================
# PART G: VISUALIZATIONS
# =============================================================================

cat("\n--- Part G: Visualizations (12 panels) ---\n")

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    plot.title = element_text(size = 12, face = "bold")
  )

format_rho <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  sprintf("rho = %+.3f\np = %.1e\nn = %s",
          test$estimate, test$p.value, format(sum(valid), big.mark = ","))
}

cat_colors <- c(Lost = "#2166AC", Unchanged = "grey60", Gained = "#B2182B")
sig_colors <- c(Sig_Down = "#2166AC", NS = "grey70", Sig_Up = "#B2182B")

# --- Panel 01: delta_activity vs K119ub Fold ---
p01 <- ggplot(enh_merged, aes(x = delta_activity, y = Fold)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_merged$delta_activity, enh_merged$Fold),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Enhancer activity change vs K119ub change"
  ) +
  theme_pub
save_multiformat_ggplot(p01, file.path(FIGURE_DIR, "01_delta_activity_vs_k119ub"), width = 7, height = 6)

# --- Panel 02: mean delta_unnorm vs K119ub Fold ---
p02 <- ggplot(enh_merged, aes(x = mean_delta_unnorm, y = Fold)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_merged$mean_delta_unnorm, enh_merged$Fold),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("Mean " * Delta * "(A" %*% "C) per enhancer"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Unnormalized ABC change vs K119ub change"
  ) +
  theme_pub
save_multiformat_ggplot(p02, file.path(FIGURE_DIR, "02_delta_unnorm_vs_k119ub"), width = 7, height = 6)

# --- Panel 03: mean delta_ABC vs K119ub Fold ---
p03 <- ggplot(enh_merged, aes(x = mean_delta_abc, y = Fold)) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_smooth(method = "lm", color = "red", linewidth = 0.5, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_merged$mean_delta_abc, enh_merged$Fold),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression("Mean " * Delta * "ABC per enhancer"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Normalized ABC change vs K119ub change"
  ) +
  theme_pub
save_multiformat_ggplot(p03, file.path(FIGURE_DIR, "03_delta_abc_vs_k119ub"), width = 7, height = 6)

# --- Panel 04: Boxplot K119ub Fold by ABC category ---
p04 <- ggplot(enh_merged, aes(x = abc_category, y = Fold, fill = abc_category)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = expression("ABC category (based on mean " * Delta * "ABC)"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "K119ub change by ABC enhancer category",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kw$p.value)
  ) +
  theme_pub
save_multiformat_ggplot(p04, file.path(FIGURE_DIR, "04_boxplot_k119ub_by_abc_category"), width = 7, height = 6)

# --- Panel 05: Boxplot K119ub Conc_mut by ABC category (log2 scale, no log1p) ---
p05 <- ggplot(enh_merged, aes(x = abc_category, y = Conc_mut, fill = abc_category)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  labs(
    x = expression("ABC category (based on mean " * Delta * "ABC)"),
    y = expression("K119ub KO concentration (DiffBind log"[2] * " scale)"),
    title = "Absolute K119ub level in KO by ABC category"
  ) +
  theme_pub
save_multiformat_ggplot(p05, file.path(FIGURE_DIR, "05_boxplot_k119ub_mut_by_category"), width = 7, height = 6)

# --- Panel 06: Violin + box overlay ---
p06 <- ggplot(enh_merged, aes(x = abc_category, y = Fold, fill = abc_category)) +
  geom_violin(alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = expression("ABC category (based on mean " * Delta * "ABC)"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "K119ub change distribution by ABC category"
  ) +
  theme_pub
save_multiformat_ggplot(p06, file.path(FIGURE_DIR, "06_violin_k119ub_by_category"), width = 7, height = 6)

# --- Panel 07: Faceted scatter by distance bin ---
enh_dist <- enh_merged[!is.na(enh_merged$dist_bin), ]

# Pre-compute rho per bin for annotation
rho_by_bin <- do.call(rbind, lapply(levels(enh_dist$dist_bin), function(bin) {
  sub <- enh_dist[enh_dist$dist_bin == bin, ]
  valid <- !is.na(sub$delta_activity) & !is.na(sub$Fold)
  n <- sum(valid)
  if (n < 10) return(data.frame(dist_bin = bin, label = paste0("n=", n), stringsAsFactors = FALSE))
  test <- cor.test(sub$delta_activity[valid], sub$Fold[valid], method = "spearman", exact = FALSE)
  data.frame(
    dist_bin = bin,
    label = sprintf("rho=%+.3f, n=%s", test$estimate, format(n, big.mark = ",")),
    stringsAsFactors = FALSE
  )
}))
rho_by_bin$dist_bin <- factor(rho_by_bin$dist_bin, levels = dist_labels)

p07 <- ggplot(enh_dist, aes(x = delta_activity, y = Fold)) +
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
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Activity-K119ub correlation by enhancer-gene distance"
  ) +
  theme_pub +
  theme(strip.text = element_text(size = 9))
save_multiformat_ggplot(p07, file.path(FIGURE_DIR, "07_scatter_by_distance"), width = 14, height = 4)

# --- Panel 08: Faceted scatter by H3K27ac status ---
enh_h3k <- enh_merged[!is.na(enh_merged$has_H3K27ac), ]
enh_h3k$h3k_label <- ifelse(enh_h3k$has_H3K27ac, "H3K27ac+", "H3K27ac-")

rho_by_h3k <- do.call(rbind, lapply(c("H3K27ac+", "H3K27ac-"), function(lbl) {
  sub <- enh_h3k[enh_h3k$h3k_label == lbl, ]
  valid <- !is.na(sub$delta_activity) & !is.na(sub$Fold)
  n <- sum(valid)
  if (n < 10) return(data.frame(h3k_label = lbl, label = paste0("n=", n), stringsAsFactors = FALSE))
  test <- cor.test(sub$delta_activity[valid], sub$Fold[valid], method = "spearman", exact = FALSE)
  data.frame(
    h3k_label = lbl,
    label = sprintf("rho=%+.3f, p=%.1e, n=%s", test$estimate, test$p.value, format(n, big.mark = ",")),
    stringsAsFactors = FALSE
  )
}))

p08 <- ggplot(enh_h3k, aes(x = delta_activity, y = Fold)) +
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
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Activity-K119ub correlation by H3K27ac status"
  ) +
  theme_pub
save_multiformat_ggplot(p08, file.path(FIGURE_DIR, "08_scatter_by_h3k27ac"), width = 10, height = 5)

# --- Panel 09: K119ub volcano at ABC enhancers ---
p09 <- ggplot(enh_merged, aes(x = Fold, y = -log10(FDR))) +
  geom_pointdensity(size = 0.3, alpha = 0.7) +
  scale_color_viridis_c(name = "Density") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("FDR<0.05: %d up, %d down\n%d NS (%.0f%%)",
                           sig_counts["Sig_Up"], sig_counts["Sig_Down"],
                           sig_counts["NS"], sig_counts["NS"] / nrow(enh_merged) * 100),
           hjust = 1.1, vjust = 1.2, size = 3.2, fontface = "bold") +
  labs(
    x = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    y = expression("-log"[10] * "(FDR)"),
    title = "K119ub differential signal at ABC enhancers"
  ) +
  theme_pub
save_multiformat_ggplot(p09, file.path(FIGURE_DIR, "09_k119ub_volcano_at_enhancers"), width = 7, height = 6)

# --- Panel 10: Contingency heatmap — ABC category × K119ub significance ---
# Reshape log2(O/E) matrix for ggplot
heatmap_df <- expand.grid(
  abc_category = rownames(log2_oe),
  k119ub_sig   = colnames(log2_oe),
  stringsAsFactors = FALSE
)
heatmap_df$log2_oe <- as.vector(log2_oe)
heatmap_df$count <- as.vector(ct)
heatmap_df$abc_category <- factor(heatmap_df$abc_category, levels = c("Lost", "Unchanged", "Gained"))
heatmap_df$k119ub_sig <- factor(heatmap_df$k119ub_sig, levels = c("Sig_Down", "NS", "Sig_Up"))

# Symmetric color limits
max_abs <- max(abs(heatmap_df$log2_oe))

p10 <- ggplot(heatmap_df, aes(x = k119ub_sig, y = abc_category, fill = log2_oe)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f\n(%d)", log2_oe, count)),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-max_abs, max_abs),
                       name = expression("log"[2] * "(O/E)")) +
  labs(
    x = "K119ub significance (DiffBind FDR < 0.05)",
    y = "ABC category",
    title = "Enrichment: ABC category vs K119ub significance",
    subtitle = sprintf("Fisher's exact test p = %.2e", fisher_res$p.value)
  ) +
  theme_pub +
  theme(panel.grid = element_blank())
save_multiformat_ggplot(p10, file.path(FIGURE_DIR, "10_contingency_heatmap"), width = 7, height = 5)

# --- Panel 11: delta_activity by K119ub significance class ---
# Kruskal-Wallis for this grouping
kw_sig <- kruskal.test(delta_activity ~ k119ub_sig, data = enh_merged)

p11 <- ggplot(enh_merged, aes(x = k119ub_sig, y = delta_activity, fill = k119ub_sig)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = sig_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "K119ub significance (DiffBind FDR < 0.05)",
    y = expression(Delta * "Activity (KO - WT)"),
    title = "Enhancer activity change by K119ub significance",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_sig$p.value)
  ) +
  theme_pub
save_multiformat_ggplot(p11, file.path(FIGURE_DIR, "11_delta_activity_by_k119ub_sig"), width = 7, height = 6)

# --- Panel 12: delta_activity vs K119ub Fold colored by significance ---
p12 <- ggplot(enh_merged, aes(x = delta_activity, y = Fold, color = k119ub_sig)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = sig_colors, name = "K119ub\nSignificance") +
  geom_smooth(data = enh_merged, aes(x = delta_activity, y = Fold),
              method = "lm", color = "black", linewidth = 0.5, se = FALSE, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf,
           label = format_rho(enh_merged$delta_activity, enh_merged$Fold),
           hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold") +
  labs(
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression("K119ub Fold (DiffBind log"[2] * "FC)"),
    title = "Activity change vs K119ub change, colored by significance"
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_pub
save_multiformat_ggplot(p12, file.path(FIGURE_DIR, "12_scatter_colored_by_significance"), width = 8, height = 6)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("  Data source:           DiffBind K119ub (%d peaks)\n", nrow(diffbind)))
cat(sprintf("  Enhancers matched:     %d / %d (%.1f%%)\n", n_matched, nrow(enh_coords), match_rate))
cat(sprintf("  Enhancers analyzed:    %d (after merge with ABC aggregates)\n", nrow(enh_merged)))
cat(sprintf("  Figures saved to:      %s/ (12 panels, PDF+SVG+JPG each)\n", FIGURE_DIR))
cat(sprintf("  Correlation summary:   %s\n", summary_out))
cat(sprintf("  Merged enhancer data:  %s\n", merged_out))

cat("\n  Core correlations:\n")
for (i in seq_len(nrow(corr_table))) {
  cat(sprintf("    %-45s: rho = %+.4f (p = %.2e)\n",
              corr_table$label[i], corr_table$rho[i], corr_table$p_value[i]))
}

cat("\n  K119ub significance at enhancers:\n")
for (sig_name in levels(enh_merged$k119ub_sig)) {
  cat(sprintf("    %-12s: %d (%.1f%%)\n", sig_name, sig_counts[sig_name],
              sig_counts[sig_name] / nrow(enh_merged) * 100))
}

cat("\n  Biological prediction: Negative rho for ΔABC vs K119ub Fold\n")
cat("  (Lost enhancer activity -> gained K119ub in BAP1-KO)\n")
cat("\nDone.\n")
