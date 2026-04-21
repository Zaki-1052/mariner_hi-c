#!/usr/bin/env Rscript
# stripes/stripenn/scripts/stripenn_visualizations.R
# Stage 7: Comprehensive Visualization for Stripenn Differential Stripe Analysis
# Adapted from stripes/quagga/scripts/stripe_visualizations.R
#
# Usage: Rscript scripts/stripenn_visualizations.R
#
# Processes both timepoints (250831/early, 250402/late) and both resolutions
# (5kb primary, 10kb validation) in a single run.

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Stripenn Visualization Analysis\n")
cat("========================================\n\n")

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(GenomicRanges)
  library(IRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(dplyr)
  library(tidyr)
  library(yaml)
  library(pheatmap)
  library(svglite)
})
cat("Packages loaded\n\n")

CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
REPO_ROOT <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c"

source(file.path(REPO_ROOT, "scripts/utils/multi_format_output.R"))

config <- yaml::read_yaml(file.path(CODE_DIR, "config", "stripenn_config.yaml"))

TIMEPOINTS <- list(
  early = list(tp = "250831", label = "Early (P12)"),
  late  = list(tp = "250402", label = "Late (Adult)")
)
RESOLUTIONS <- c(5000, 10000)

COLOR_PALETTE <- c(
  "Lost (High)"    = "#00008B",
  "Lost (Medium)"  = "#4169E1",
  "Lost (Low)"     = "#87CEEB",
  "Gained (High)"  = "#8B0000",
  "Gained (Medium)" = "#DC143C",
  "Gained (Low)"   = "#FFA07A",
  "Unchanged"      = "#D3D3D3",
  "Strengthened"   = "#006400",
  "Weakened"       = "#8B008B"
)

DIRECTION_COLORS <- c(
  "lost" = "#4575B4", "gained" = "#D73027",
  "unchanged" = "#999999", "strengthened" = "#006400", "weakened" = "#8B008B"
)

SOURCE_COLORS <- c(
  "shared" = "#999999", "control_only" = "#4575B4", "mutant_only" = "#D73027"
)

RGB_MAP <- c(
  "Lost (High)"    = "0,0,139",
  "Lost (Medium)"  = "65,105,225",
  "Lost (Low)"     = "135,206,235",
  "Gained (High)"  = "139,0,0",
  "Gained (Medium)" = "220,20,60",
  "Gained (Low)"   = "255,160,122",
  "Unchanged"      = "128,128,128"
)

assign_plot_category <- function(direction, confidence) {
  dplyr::case_when(
    direction == "lost"    & confidence == "high"   ~ "Lost (High)",
    direction == "lost"    & confidence == "medium" ~ "Lost (Medium)",
    direction == "lost"    & confidence == "low"    ~ "Lost (Low)",
    direction == "gained"  & confidence == "high"   ~ "Gained (High)",
    direction == "gained"  & confidence == "medium" ~ "Gained (Medium)",
    direction == "gained"  & confidence == "low"    ~ "Gained (Low)",
    direction == "strengthened" ~ "Strengthened",
    direction == "weakened"     ~ "Weakened",
    TRUE ~ "Unchanged"
  )
}

theme_stripe <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      axis.title    = element_text(size = base_size + 1),
      axis.text     = element_text(size = base_size - 1),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

cat("\n========================================\n")
cat("SECTION 1: Loading Data\n")
cat("========================================\n\n")

load_stripenn_data <- function(tp_id, res_kb) {
  base <- file.path(DATA_DIR, "outputs", tp_id, paste0("res_", res_kb, "kb"))
  f <- file.path(base, "05_final_differential.tsv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  df <- read.delim(f, stringsAsFactors = FALSE)
  df$stripe_length_kb <- df$stripe_length / 1000
  df$anchor_width_kb  <- df$stripe_width / 1000
  df$plot_cat <- assign_plot_category(df$direction, df$direction_confidence)
  df$FDR_safe <- ifelse(is.na(df$FDR) | df$FDR == 0, 1e-300, df$FDR)
  cat(sprintf("  Loaded %s/%dkb: %d stripes\n", tp_id, res_kb, nrow(df)))
  df
}

load_cross_res <- function(tp_id) {
  f <- file.path(DATA_DIR, "outputs", tp_id, "cross_res_merged.tsv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  df <- read.delim(f, stringsAsFactors = FALSE)
  df$stripe_length_kb <- df$stripe_length / 1000
  df$anchor_width_kb  <- df$stripe_width / 1000
  df$plot_cat <- assign_plot_category(df$direction, df$direction_confidence)
  df$FDR_safe <- ifelse(is.na(df$FDR) | df$FDR == 0, 1e-300, df$FDR)
  cat(sprintf("  Loaded %s cross-res: %d stripes\n", tp_id, nrow(df)))
  df
}

load_count_matrix <- function(tp_id, res_kb) {
  f <- file.path(DATA_DIR, "outputs", tp_id, paste0("res_", res_kb, "kb"),
                 "04_edgeR", "count_matrix.tsv")
  if (!file.exists(f)) return(NULL)
  read.delim(f, stringsAsFactors = FALSE, row.names = 1)
}

data_all <- list()
for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  data_all[[tp_name]] <- list()
  data_all[[tp_name]][["5kb"]]  <- load_stripenn_data(tp_id, 5)
  data_all[[tp_name]][["10kb"]] <- load_stripenn_data(tp_id, 10)
  data_all[[tp_name]][["cross_res"]] <- load_cross_res(tp_id)
  data_all[[tp_name]][["counts_5kb"]]  <- load_count_matrix(tp_id, 5)
  data_all[[tp_name]][["counts_10kb"]] <- load_count_matrix(tp_id, 10)
}

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  for (subdir in c("visualizations", file.path("visualizations", "enrichment"))) {
    dir.create(file.path(DATA_DIR, "outputs", tp_id, subdir),
               recursive = TRUE, showWarnings = FALSE)
  }
}
combined_dir <- file.path(DATA_DIR, "outputs", "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(combined_dir, "enrichment"), recursive = TRUE, showWarnings = FALSE)

cat("\nSection 1 complete\n")

# =============================================================================
# SECTION 2: VOLCANO PLOTS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 2: Volcano Plots\n")
cat("========================================\n\n")

create_volcano <- function(df, title_label, n_total_override = NULL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  n_total  <- nrow(df)
  n_sig    <- sum(df$significant_FDR05 == TRUE, na.rm = TRUE)
  n_lost   <- sum(df$direction == "lost", na.rm = TRUE)
  n_gained <- sum(df$direction == "gained", na.rm = TRUE)
  hi_lost  <- sum(df$direction == "lost" & df$direction_confidence == "high", na.rm = TRUE)
  hi_gain  <- sum(df$direction == "gained" & df$direction_confidence == "high", na.rm = TRUE)

  logfc_rng <- range(df$logFC, na.rm = TRUE)
  ymax <- -log10(min(df$FDR_safe, na.rm = TRUE))

  pt_alpha <- ifelse(n_total > 3000, 0.35, 0.6)
  pt_size  <- ifelse(n_total > 3000, 1.2, 2.0)

  p <- ggplot(df, aes(x = logFC, y = -log10(FDR_safe), color = plot_cat)) +
    geom_point(alpha = pt_alpha, size = pt_size) +
    scale_color_manual(values = COLOR_PALETTE, name = "Direction") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.10), linetype = "dotted", color = "gray", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    labs(
      title = sprintf("Differential Stripes: %s", title_label),
      subtitle = sprintf("Total: %d | FDR<0.05: %d (%.1f%%) | Lost: %d | Gained: %d",
                         n_total, n_sig, 100 * n_sig / n_total, n_lost, n_gained),
      x = expression(log[2]~"Fold Change (Mutant/Control)"),
      y = expression(-log[10]~"FDR")
    ) +
    annotate("text", x = logfc_rng[1] * 0.9, y = ymax * 0.95,
             label = sprintf("Lost: %d\n(High: %d)", n_lost, hi_lost),
             color = "#4169E1", size = 3.5, fontface = "bold", hjust = 0) +
    annotate("text", x = logfc_rng[2] * 0.9, y = ymax * 0.95,
             label = sprintf("Gained: %d\n(High: %d)", n_gained, hi_gain),
             color = "#DC143C", size = 3.5, fontface = "bold", hjust = 1) +
    theme_stripe()

  p
}

volcano_plots <- list()
for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  tp_label <- TIMEPOINTS[[tp_name]]$label
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")

  for (res_label in c("5kb", "10kb")) {
    df <- data_all[[tp_name]][[res_label]]
    if (is.null(df)) next
    tag <- paste0(tp_name, "_", res_label)
    p <- create_volcano(df, sprintf("%s %s", tp_label, res_label))
    if (!is.null(p)) {
      save_multiformat_ggplot(p, file.path(viz_dir, sprintf("volcano_%s_%s", tp_id, res_label)),
                              width = 10, height = 8, use_subfolders = FALSE)
      volcano_plots[[tag]] <- p
    }
  }
}

if (length(volcano_plots) == 4) {
  combined_v <- (volcano_plots$early_5kb + volcano_plots$early_10kb) /
    (volcano_plots$late_5kb + volcano_plots$late_10kb) +
    plot_annotation(title = "Differential Stripes: Stripenn Pipeline",
                    subtitle = "BAP1-KO vs Wildtype, Mouse Cerebellum (mm10)")
  save_multiformat_ggplot(combined_v, file.path(combined_dir, "volcano_combined"),
                          width = 20, height = 16, use_subfolders = FALSE)
}

cat("Section 2 complete\n")

# =============================================================================
# SECTION 3: STRIPINESS SCORE ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 3: Stripiness Score Analysis\n")
cat("========================================\n\n")

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  tp_label <- TIMEPOINTS[[tp_name]]$label
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")
  df <- data_all[[tp_name]][["5kb"]]
  if (is.null(df)) next

  cat(sprintf("Stripiness analysis: %s\n", tp_label))

  shared <- df %>% filter(!is.na(stripiness_ctrl) & !is.na(stripiness_mut))
  cat(sprintf("  Shared stripes with both scores: %d\n", nrow(shared)))

  if (nrow(shared) > 0) {
    pt_alpha <- ifelse(nrow(shared) > 2000, 0.25, 0.5)
    p_scatter <- ggplot(shared, aes(x = stripiness_ctrl, y = stripiness_mut, color = plot_cat)) +
      geom_point(alpha = pt_alpha, size = 1) +
      scale_color_manual(values = COLOR_PALETTE, name = "Direction") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "darkgray", linewidth = 0.5, linetype = "dotted") +
      labs(title = sprintf("Stripiness: Control vs Mutant (%s)", tp_label),
           subtitle = sprintf("Shared stripes only (n=%d)", nrow(shared)),
           x = "Stripiness (Control)", y = "Stripiness (Mutant)") +
      theme_stripe()

    diff_df <- df %>%
      filter(direction %in% c("lost", "gained", "unchanged")) %>%
      mutate(stripiness_val = coalesce(stripiness_ctrl, stripiness_mut))

    p_box <- ggplot(diff_df, aes(x = direction, y = stripiness_val, fill = direction)) +
      geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      labs(title = sprintf("Stripiness by Direction (%s)", tp_label),
           x = "Direction", y = "Stripiness Score") +
      theme_stripe() + theme(legend.position = "none")

    combined_s <- p_scatter + p_box + plot_layout(widths = c(2, 1))
    save_multiformat_ggplot(combined_s, file.path(viz_dir, sprintf("stripiness_%s", tp_id)),
                            width = 15, height = 7, use_subfolders = FALSE)
  }
}

cat("Section 3 complete\n")

# =============================================================================
# SECTION 4: LENGTH & WIDTH DISTRIBUTIONS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 4: Length & Width Distributions\n")
cat("========================================\n\n")

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  tp_label <- TIMEPOINTS[[tp_name]]$label
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")
  df <- data_all[[tp_name]][["5kb"]]
  if (is.null(df)) next

  cat(sprintf("Length analysis: %s\n", tp_label))

  diff_only <- df %>% filter(direction %in% c("lost", "gained"))

  if (nrow(diff_only) > 1) {
    lost_len  <- diff_only$stripe_length_kb[diff_only$direction == "lost"]
    gain_len  <- diff_only$stripe_length_kb[diff_only$direction == "gained"]
    wt <- NULL
    if (length(lost_len) > 5 && length(gain_len) > 5) {
      wt <- wilcox.test(lost_len, gain_len)
    }

    p_violin <- ggplot(diff_only, aes(x = direction, y = stripe_length_kb, fill = direction)) +
      geom_violin(alpha = 0.7, trim = TRUE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      labs(title = sprintf("Stripe Length: %s", tp_label),
           subtitle = ifelse(!is.null(wt),
                             sprintf("Lost vs Gained: Wilcoxon p = %.2e", wt$p.value), ""),
           x = "Direction", y = "Stripe Length (kb)") +
      theme_stripe() + theme(legend.position = "none")

    p_hist <- ggplot(diff_only, aes(x = stripe_length_kb, fill = direction)) +
      geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = DIRECTION_COLORS) +
      facet_wrap(~direction, ncol = 1, scales = "free_y") +
      labs(title = sprintf("Length Histogram: %s", tp_label),
           x = "Stripe Length (kb)", y = "Count") +
      theme_stripe() + theme(legend.position = "none",
                             strip.text = element_text(face = "bold"))

    width_df <- df %>%
      filter(direction %in% c("lost", "gained", "unchanged")) %>%
      mutate(anchor_width_kb = stripe_width / 1000)

    p_width <- ggplot(width_df, aes(x = direction, y = anchor_width_kb, fill = direction)) +
      geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      labs(title = sprintf("Anchor Width: %s", tp_label),
           x = "Direction", y = "Anchor Width (kb)") +
      theme_stripe() + theme(legend.position = "none")

    combined_l <- (p_violin | p_hist) / p_width + plot_layout(heights = c(2, 1))
    save_multiformat_ggplot(combined_l, file.path(viz_dir, sprintf("length_distribution_%s", tp_id)),
                            width = 14, height = 12, use_subfolders = FALSE)

    stats_df <- diff_only %>%
      group_by(direction) %>%
      summarise(n = n(),
                median_length_kb = median(stripe_length_kb, na.rm = TRUE),
                mean_length_kb   = mean(stripe_length_kb, na.rm = TRUE),
                median_width_kb  = median(stripe_width / 1000, na.rm = TRUE),
                .groups = "drop")
    write.table(stats_df,
                file.path(viz_dir, sprintf("length_statistics_%s.tsv", tp_id)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

early_df <- data_all$early[["5kb"]]
late_df  <- data_all$late[["5kb"]]
if (!is.null(early_df) && !is.null(late_df)) {
  comb_df <- bind_rows(
    early_df %>% filter(direction %in% c("lost", "gained")) %>% mutate(timepoint = "Early"),
    late_df  %>% filter(direction %in% c("lost", "gained")) %>% mutate(timepoint = "Late")
  )
  if (nrow(comb_df) > 0) {
    p_comp <- ggplot(comb_df, aes(x = direction, y = stripe_length_kb, fill = direction)) +
      geom_violin(alpha = 0.7, trim = TRUE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      facet_wrap(~timepoint) +
      labs(title = "Stripe Length: Early vs Late",
           x = "Direction", y = "Stripe Length (kb)") +
      theme_stripe() + theme(legend.position = "none",
                             strip.text = element_text(size = 12, face = "bold"))
    save_multiformat_ggplot(p_comp, file.path(combined_dir, "length_comparison"),
                            width = 12, height = 6, use_subfolders = FALSE)
  }
}

cat("Section 4 complete\n")

# =============================================================================
# SECTION 5: SOURCE & DIRECTION DISTRIBUTIONS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 5: Source & Direction Distributions\n")
cat("========================================\n\n")

summary_rows <- list()
for (tp_name in names(TIMEPOINTS)) {
  for (res_label in c("5kb", "10kb")) {
    df <- data_all[[tp_name]][[res_label]]
    if (is.null(df)) next
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      timepoint = TIMEPOINTS[[tp_name]]$label,
      resolution = res_label,
      direction = df$direction,
      source = df$source,
      direction_confidence = df$direction_confidence,
      direction_type = df$direction_type,
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_rows) > 0) {
  sum_df <- bind_rows(summary_rows)
  sum_df$dataset <- paste(sum_df$timepoint, sum_df$resolution)

  source_ct <- sum_df %>% count(dataset, source)
  p_src <- ggplot(source_ct, aes(x = dataset, y = n, fill = source)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = SOURCE_COLORS, name = "Source") +
    labs(title = "Stripe Source Distribution", x = "", y = "Count") +
    theme_stripe() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  dir_ct <- sum_df %>%
    filter(direction %in% c("lost", "gained", "unchanged")) %>%
    count(dataset, direction)
  p_dir <- ggplot(dir_ct, aes(x = dataset, y = n, fill = direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = DIRECTION_COLORS, name = "Direction") +
    labs(title = "Direction Breakdown", x = "", y = "Count") +
    theme_stripe() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  conf_df <- sum_df %>%
    filter(direction %in% c("lost", "gained")) %>%
    count(dataset, direction_confidence, direction)
  p_conf <- ggplot(conf_df, aes(x = direction_confidence, y = n, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    facet_wrap(~dataset, scales = "free_y") +
    labs(title = "Confidence Tiers (Lost/Gained Only)",
         x = "Confidence", y = "Count") +
    theme_stripe() + theme(strip.text = element_text(face = "bold"))

  dtype_df <- sum_df %>%
    filter(direction %in% c("lost", "gained")) %>%
    count(dataset, direction_type, direction)
  p_dtype <- ggplot(dtype_df, aes(x = direction_type, y = n, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    facet_wrap(~dataset, scales = "free_y") +
    labs(title = "Stripe Orientation (3'/5')", x = "Direction Type", y = "Count") +
    theme_stripe() + theme(strip.text = element_text(face = "bold"))

  combined_sd <- (p_src | p_dir) / (p_conf | p_dtype) +
    plot_annotation(title = "Stripenn: Source & Direction Summary")
  save_multiformat_ggplot(combined_sd, file.path(combined_dir, "source_direction_summary"),
                          width = 18, height = 14, use_subfolders = FALSE)

  for (tp_name in names(TIMEPOINTS)) {
    tp_id <- TIMEPOINTS[[tp_name]]$tp
    viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")
    tp_sub <- sum_df %>% filter(timepoint == TIMEPOINTS[[tp_name]]$label)
    dir_sub <- tp_sub %>%
      filter(direction %in% c("lost", "gained", "unchanged")) %>%
      count(resolution, direction)
    p_d <- ggplot(dir_sub, aes(x = resolution, y = n, fill = direction)) +
      geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
      scale_fill_manual(values = DIRECTION_COLORS) +
      labs(title = sprintf("Direction Breakdown: %s", TIMEPOINTS[[tp_name]]$label),
           x = "Resolution", y = "Count") +
      theme_stripe()
    save_multiformat_ggplot(p_d, file.path(viz_dir, sprintf("direction_breakdown_%s", tp_id)),
                            width = 8, height = 6, use_subfolders = FALSE)
  }
}

cat("Section 5 complete\n")

# =============================================================================
# SECTION 6: CROSS-RESOLUTION CONCORDANCE
# =============================================================================

cat("\n========================================\n")
cat("SECTION 6: Cross-Resolution Concordance\n")
cat("========================================\n\n")

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  tp_label <- TIMEPOINTS[[tp_name]]$label
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")
  cr <- data_all[[tp_name]][["cross_res"]]
  if (is.null(cr)) next

  cat(sprintf("Cross-res analysis: %s\n", tp_label))

  matched <- cr %>% filter(!is.na(logFC_10kb))
  if (nrow(matched) > 10) {
    r_val <- cor(matched$logFC, matched$logFC_10kb, use = "complete.obs")
    pt_alpha <- ifelse(nrow(matched) > 2000, 0.25, 0.5)

    res_colors <- c("both_concordant" = "#2ca02c", "both_discordant" = "#d62728",
                    "5kb_only" = "#1f77b4", "10kb_only" = "#ff7f0e")

    p_corr <- ggplot(matched, aes(x = logFC, y = logFC_10kb, color = resolution_support)) +
      geom_point(alpha = pt_alpha, size = 1) +
      scale_color_manual(values = res_colors, name = "Resolution\nSupport") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
      labs(title = sprintf("Cross-Resolution logFC Correlation: %s", tp_label),
           subtitle = sprintf("n=%d matched stripes, Pearson r=%.3f", nrow(matched), r_val),
           x = "logFC (5kb)", y = "logFC (10kb)") +
      theme_stripe()

    sup_ct <- cr %>% count(resolution_support)
    p_sup <- ggplot(sup_ct, aes(x = reorder(resolution_support, -n), y = n,
                                fill = resolution_support)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_manual(values = res_colors) +
      geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
      labs(title = sprintf("Resolution Support: %s", tp_label),
           x = "", y = "Count") +
      theme_stripe() + theme(legend.position = "none",
                             axis.text.x = element_text(angle = 30, hjust = 1))

    combined_cr <- p_corr + p_sup + plot_layout(widths = c(2, 1))
    save_multiformat_ggplot(combined_cr, file.path(viz_dir, sprintf("cross_res_%s", tp_id)),
                            width = 15, height = 7, use_subfolders = FALSE)
  }
}

early_cr <- data_all$early[["cross_res"]]
late_cr  <- data_all$late[["cross_res"]]
if (!is.null(early_cr) && !is.null(late_cr)) {
  sup_combined <- bind_rows(
    early_cr %>% count(resolution_support) %>% mutate(timepoint = "Early"),
    late_cr  %>% count(resolution_support) %>% mutate(timepoint = "Late")
  )
  res_colors <- c("both_concordant" = "#2ca02c", "both_discordant" = "#d62728",
                  "5kb_only" = "#1f77b4", "10kb_only" = "#ff7f0e")
  p_crc <- ggplot(sup_combined, aes(x = resolution_support, y = n, fill = resolution_support)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = res_colors) +
    facet_wrap(~timepoint, scales = "free_y") +
    labs(title = "Cross-Resolution Support: Early vs Late", x = "", y = "Count") +
    theme_stripe() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 30, hjust = 1),
                           strip.text = element_text(face = "bold"))
  save_multiformat_ggplot(p_crc, file.path(combined_dir, "cross_res_comparison"),
                          width = 12, height = 6, use_subfolders = FALSE)
}

cat("Section 6 complete\n")

# =============================================================================
# SECTION 7: REPLICATE CORRELATION HEATMAPS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 7: Replicate Correlation\n")
cat("========================================\n\n")

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  tp_label <- TIMEPOINTS[[tp_name]]$label
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")

  for (res_label in c("5kb", "10kb")) {
    counts <- data_all[[tp_name]][[paste0("counts_", res_label)]]
    if (is.null(counts)) next

    cat(sprintf("  Replicate heatmap: %s %s (%d stripes)\n", tp_label, res_label, nrow(counts)))

    cor_mat <- cor(counts, method = "spearman")
    ann_col <- data.frame(
      Condition = ifelse(grepl("^ctrl", colnames(counts)), "Control", "Mutant"),
      row.names = colnames(counts)
    )
    ann_colors <- list(Condition = c(Control = "#4575B4", Mutant = "#D73027"))

    save_multiformat_pheatmap(
      quote(pheatmap(cor_mat,
                     display_numbers = TRUE, number_format = "%.4f",
                     color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
                     breaks = seq(min(cor_mat) - 0.001, 1, length.out = 101),
                     annotation_col = ann_col, annotation_colors = ann_colors,
                     main = sprintf("Replicate Correlation: %s %s", tp_label, res_label))),
      file.path(viz_dir, sprintf("replicate_correlation_%s_%s", tp_id, res_label)),
      width = 8, height = 7, use_subfolders = FALSE
    )
  }
}

cat("Section 7 complete\n")

# =============================================================================
# SECTION 8: CHIP-SEQ ANCHOR ANNOTATION
# =============================================================================

cat("\n========================================\n")
cat("SECTION 8: ChIP-seq Anchor Annotation\n")
cat("========================================\n\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss <- resize(genes(txdb), width = 1, fix = "start")
cat(sprintf("Loaded %d TSS positions from mm10\n", length(tss)))

load_peaks <- function(peak_file) {
  if (is.null(peak_file) || !file.exists(peak_file)) return(NULL)
  pk <- tryCatch(read.table(peak_file, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(pk) || nrow(pk) == 0) return(NULL)
  GRanges(seqnames = pk$V1, ranges = IRanges(start = pk$V2 + 1, end = pk$V3))
}

annotate_stripes <- function(df, tp_name) {
  if (is.null(df) || nrow(df) == 0) return(df)

  cat(sprintf("  Annotating anchors: %s (%d stripes)\n", tp_name, nrow(df)))

  peaks_cfg <- config$chipseq_peaks[[tp_name]]
  peak_files <- lapply(peaks_cfg, function(rel) file.path(REPO_ROOT, rel))

  h3k27ac  <- load_peaks(peak_files$h3k27ac)
  h3k27me3 <- load_peaks(peak_files$h3k27me3)
  h3k4me1  <- load_peaks(peak_files$h3k4me1)
  h3k4me3  <- load_peaks(peak_files$h3k4me3)
  bivalent <- load_peaks(peak_files$bivalent)

  anchor_gr <- GRanges(seqnames = df$chr,
                       ranges = IRanges(start = df$pos1, end = df$pos2))

  nearest_tss <- distanceToNearest(anchor_gr, tss)
  df$distance_to_tss <- NA_real_
  df$nearest_gene_id <- NA_character_
  if (length(nearest_tss) > 0) {
    df$distance_to_tss[queryHits(nearest_tss)] <- mcols(nearest_tss)$distance
    df$nearest_gene_id[queryHits(nearest_tss)] <- names(tss)[subjectHits(nearest_tss)]
  }

  if (any(!is.na(df$nearest_gene_id))) {
    gm <- tryCatch(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = unique(na.omit(df$nearest_gene_id)),
                            columns = "SYMBOL", keytype = "ENTREZID"),
      error = function(e) data.frame(ENTREZID = character(), SYMBOL = character()))
    df$nearest_gene <- gm$SYMBOL[match(df$nearest_gene_id, gm$ENTREZID)]
  } else {
    df$nearest_gene <- NA_character_
  }

  df$h3k27ac  <- if (!is.null(h3k27ac))  countOverlaps(anchor_gr, h3k27ac)  > 0 else FALSE
  df$h3k27me3 <- if (!is.null(h3k27me3)) countOverlaps(anchor_gr, h3k27me3) > 0 else FALSE
  df$h3k4me1  <- if (!is.null(h3k4me1))  countOverlaps(anchor_gr, h3k4me1)  > 0 else FALSE
  df$h3k4me3  <- if (!is.null(h3k4me3))  countOverlaps(anchor_gr, h3k4me3)  > 0 else FALSE
  df$bivalent <- if (!is.null(bivalent))  countOverlaps(anchor_gr, bivalent) > 0 else FALSE

  tss_thresh <- 2000
  df$anchor_type <- dplyr::case_when(
    df$bivalent & df$distance_to_tss <= tss_thresh ~ "Bivalent_Promoter",
    df$h3k4me3 & df$h3k27ac & df$distance_to_tss <= tss_thresh ~ "Active_Promoter",
    df$h3k27me3 & !df$h3k27ac & df$distance_to_tss <= tss_thresh ~ "Repressed_Promoter",
    df$h3k27ac & df$distance_to_tss > tss_thresh ~ "Active_Enhancer",
    df$h3k4me1 & !df$h3k27ac & df$distance_to_tss > tss_thresh ~ "Poised_Enhancer",
    df$h3k27me3 & df$distance_to_tss > tss_thresh ~ "Polycomb",
    TRUE ~ "Other"
  )

  at <- table(df$anchor_type)
  cat("  Anchor types:\n")
  for (a in names(at)) cat(sprintf("    %s: %d (%.1f%%)\n", a, at[a], 100 * at[a] / nrow(df)))

  df
}

annotated <- list()
for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")

  df <- data_all[[tp_name]][["5kb"]]
  if (is.null(df)) next

  ann_df <- annotate_stripes(df, tp_name)
  annotated[[tp_name]] <- ann_df

  write.table(ann_df,
              file.path(viz_dir, sprintf("%s_annotated_stripes.tsv", tp_id)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s_annotated_stripes.tsv\n", tp_id))

  anchor_colors <- c(
    "Active_Promoter" = "#e41a1c", "Repressed_Promoter" = "#984ea3",
    "Bivalent_Promoter" = "#ff7f00", "Active_Enhancer" = "#377eb8",
    "Poised_Enhancer" = "#4daf4a", "Polycomb" = "#a65628", "Other" = "#999999"
  )

  diff_ann <- ann_df %>% filter(direction %in% c("lost", "gained", "unchanged"))
  at_pct <- diff_ann %>%
    count(direction, anchor_type) %>%
    group_by(direction) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  p_at <- ggplot(at_pct, aes(x = direction, y = pct, fill = anchor_type)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = anchor_colors, name = "Anchor Type") +
    labs(title = sprintf("Anchor Classification: %s", TIMEPOINTS[[tp_name]]$label),
         x = "Direction", y = "Percentage (%)") +
    theme_stripe()

  diff_tss <- diff_ann %>%
    filter(direction %in% c("lost", "gained") & !is.na(distance_to_tss))
  p_tss <- ggplot(diff_tss, aes(x = direction, y = distance_to_tss / 1000, fill = direction)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    scale_y_log10() +
    labs(title = sprintf("Distance to TSS: %s", TIMEPOINTS[[tp_name]]$label),
         x = "Direction", y = "Distance to TSS (kb, log10)") +
    theme_stripe() + theme(legend.position = "none")

  chip_long <- diff_ann %>%
    filter(direction %in% c("lost", "gained")) %>%
    select(direction, h3k27ac, h3k27me3, h3k4me1, h3k4me3) %>%
    pivot_longer(cols = -direction, names_to = "mark", values_to = "overlap") %>%
    group_by(direction, mark) %>%
    summarise(n_overlap = sum(overlap), n_total = n(),
              pct = 100 * n_overlap / n_total, .groups = "drop")

  p_chip <- ggplot(chip_long, aes(x = mark, y = pct, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    labs(title = sprintf("ChIP-seq Overlap: %s", TIMEPOINTS[[tp_name]]$label),
         x = "Histone Mark", y = "% of Anchors") +
    theme_stripe() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  combined_ann <- (p_at | p_tss) / p_chip + plot_layout(heights = c(2, 1)) +
    plot_annotation(title = sprintf("Anchor Annotation: %s", TIMEPOINTS[[tp_name]]$label))
  save_multiformat_ggplot(combined_ann,
                          file.path(viz_dir, sprintf("anchor_annotation_%s", tp_id)),
                          width = 14, height = 12, use_subfolders = FALSE)
}

cat("Section 8 complete\n")

# =============================================================================
# SECTION 9: BEDPE EXPORT (ANNOTATED)
# =============================================================================

cat("\n========================================\n")
cat("SECTION 9: Annotated BEDPE Export\n")
cat("========================================\n\n")

make_bedpe <- function(df, cr_df = NULL) {
  if (!is.null(cr_df)) {
    cr_cols <- cr_df %>% select(stripe_id, in_10kb, resolution_support)
    df <- df %>% left_join(cr_cols, by = "stripe_id")
  }
  if (!"in_10kb" %in% names(df)) df$in_10kb <- NA
  if (!"nearest_gene" %in% names(df)) df$nearest_gene <- NA
  if (!"distance_to_tss" %in% names(df)) df$distance_to_tss <- NA
  if (!"anchor_type" %in% names(df)) df$anchor_type <- NA
  if (!"h3k27ac" %in% names(df))  df$h3k27ac <- NA
  if (!"h3k27me3" %in% names(df)) df$h3k27me3 <- NA
  if (!"h3k4me1" %in% names(df))  df$h3k4me1 <- NA

  rgb_col <- RGB_MAP[assign_plot_category(df$direction, df$direction_confidence)]
  rgb_col[is.na(rgb_col)] <- "128,128,128"

  data.frame(
    chr1 = df$chr, x1 = df$pos1, x2 = df$pos2,
    chr2 = df$chr, y1 = df$pos3, y2 = df$pos4,
    name = df$stripe_id,
    score = pmin(300, -log10(ifelse(is.na(df$FDR) | df$FDR == 0, 1e-300, df$FDR) + 1e-300)),
    strand1 = ".", strand2 = ".",
    color = rgb_col,
    direction = df$direction,
    direction_confidence = df$direction_confidence,
    logFC = round(df$logFC, 4),
    FDR = signif(df$FDR, 4),
    source = df$source,
    pval_ctrl = signif(df$pval_ctrl, 4),
    pval_mut = signif(df$pval_mut, 4),
    detection_confidence = df$direction_confidence,
    in_10kb = df$in_10kb,
    nearest_gene = df$nearest_gene,
    distance_to_tss = df$distance_to_tss,
    anchor_type = df$anchor_type,
    h3k27ac = df$h3k27ac,
    h3k27me3 = df$h3k27me3,
    h3k4me1 = df$h3k4me1,
    stripe_length_kb = round(df$stripe_length / 1000, 1),
    anchor_width_kb = round(df$stripe_width / 1000, 1),
    stringsAsFactors = FALSE
  )
}

for (tp_name in names(TIMEPOINTS)) {
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  viz_dir <- file.path(DATA_DIR, "outputs", tp_id, "visualizations")
  ann_df <- annotated[[tp_name]]
  cr_df  <- data_all[[tp_name]][["cross_res"]]
  if (is.null(ann_df)) next

  cat(sprintf("BEDPE export: %s\n", TIMEPOINTS[[tp_name]]$label))

  highconf <- ann_df %>%
    filter(direction %in% c("lost", "gained") &
             significant_FDR05 == TRUE & direction_confidence == "high")
  allsig <- ann_df %>%
    filter(direction %in% c("lost", "gained") & significant_FDR05 == TRUE)

  if (nrow(highconf) > 0) {
    bp <- make_bedpe(highconf, cr_df)
    write.table(bp, file.path(viz_dir, sprintf("%s_stripes_highconf.bedpe", tp_id)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  highconf: %d stripes\n", nrow(bp)))

    full_start <- pmin(highconf$pos1, highconf$pos3)
    full_end   <- pmax(highconf$pos2, highconf$pos4)
    diag_bp <- bp
    diag_bp$x1 <- full_start; diag_bp$x2 <- full_end
    diag_bp$y1 <- full_start; diag_bp$y2 <- full_end
    write.table(diag_bp, file.path(viz_dir, sprintf("%s_stripes_diagonal.bedpe", tp_id)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    rect_bp <- bp
    rect_bp$y1 <- full_start; rect_bp$y2 <- full_end
    write.table(rect_bp, file.path(viz_dir, sprintf("%s_stripes_rectangle.bedpe", tp_id)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (nrow(allsig) > 0) {
    bp <- make_bedpe(allsig, cr_df)
    write.table(bp, file.path(viz_dir, sprintf("%s_stripes_allsig.bedpe", tp_id)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  allsig: %d stripes\n", nrow(bp)))
  }

  if (!is.null(cr_df)) {
    concordant <- cr_df %>%
      filter(direction %in% c("lost", "gained") & resolution_support == "both_concordant")
    if (nrow(concordant) > 0) {
      cr_ann <- annotated[[tp_name]]
      if (!is.null(cr_ann)) {
        ann_cols <- cr_ann %>%
          select(stripe_id, any_of(c("nearest_gene", "distance_to_tss", "anchor_type",
                                     "h3k27ac", "h3k27me3", "h3k4me1")))
        concordant <- concordant %>% left_join(ann_cols, by = "stripe_id")
      }
      bp <- make_bedpe(concordant)
      write.table(bp, file.path(viz_dir, sprintf("%s_stripes_concordant.bedpe", tp_id)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("  concordant: %d stripes\n", nrow(bp)))
    }
  }
}

cat("Section 9 complete\n")

# =============================================================================
# SECTION 10: GO/KEGG ENRICHMENT
# =============================================================================

cat("\n========================================\n")
cat("SECTION 10: GO/KEGG Enrichment\n")
cat("========================================\n\n")

enrichment_ok <- tryCatch({
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
  })
  TRUE
}, error = function(e) {
  cat(sprintf("  clusterProfiler not available: %s\n", e$message))
  FALSE
})

if (enrichment_ok) {
  genes_txdb <- genes(txdb)

  get_anchor_genes <- function(df, max_dist = 10000) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    anchor_gr <- GRanges(seqnames = df$chr,
                         ranges = IRanges(start = df$pos1, end = df$pos2))
    ov <- findOverlaps(anchor_gr, genes_txdb, maxgap = max_dist)
    unique(names(genes_txdb)[subjectHits(ov)])
  }

  for (tp_name in names(TIMEPOINTS)) {
    tp_id <- TIMEPOINTS[[tp_name]]$tp
    tp_label <- TIMEPOINTS[[tp_name]]$label
    enrich_dir <- file.path(DATA_DIR, "outputs", "combined", "enrichment")

    ann_df <- annotated[[tp_name]]
    if (is.null(ann_df)) next

    sig_df <- ann_df %>% filter(significant_FDR05 == TRUE & direction %in% c("lost", "gained"))
    lost_df  <- sig_df %>% filter(direction == "lost")
    gained_df <- sig_df %>% filter(direction == "gained")

    lost_genes  <- get_anchor_genes(lost_df)
    gained_genes <- get_anchor_genes(gained_df)
    bg_genes    <- get_anchor_genes(ann_df)

    cat(sprintf("  %s: lost=%d genes, gained=%d genes, background=%d\n",
                tp_label, length(lost_genes), length(gained_genes), length(bg_genes)))

    gene_list <- list()
    if (length(lost_genes) >= 15)  gene_list$Lost  <- lost_genes
    if (length(gained_genes) >= 15) gene_list$Gained <- gained_genes

    if (length(gene_list) == 0) {
      cat(sprintf("  Skipping enrichment for %s: too few genes\n", tp_label))
      next
    }

    gene_df <- data.frame(
      gene_id = c(lost_genes, gained_genes),
      direction = c(rep("lost", length(lost_genes)), rep("gained", length(gained_genes))),
      stringsAsFactors = FALSE
    )
    sym_map <- tryCatch(
      AnnotationDbi::select(org.Mm.eg.db, keys = unique(gene_df$gene_id),
                            columns = "SYMBOL", keytype = "ENTREZID"),
      error = function(e) data.frame(ENTREZID = character(), SYMBOL = character()))
    gene_df$symbol <- sym_map$SYMBOL[match(gene_df$gene_id, sym_map$ENTREZID)]
    write.table(gene_df,
                file.path(enrich_dir, sprintf("stripe_anchor_genes_%s.tsv", tp_name)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    for (ont in c("BP", "CC", "MF")) {
      cat(sprintf("  GO %s...\n", ont))
      go_res <- tryCatch(
        compareCluster(geneCluster = gene_list, fun = "enrichGO",
                       OrgDb = org.Mm.eg.db, ont = ont,
                       universe = bg_genes,
                       pAdjustMethod = "BH", pvalueCutoff = 0.05),
        error = function(e) { cat(sprintf("    Error: %s\n", e$message)); NULL })

      if (!is.null(go_res) && nrow(go_res@compareClusterResult) > 0) {
        p_go <- dotplot(go_res, showCategory = 15) +
          labs(title = sprintf("GO %s: %s Stripes", ont, tp_label)) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        save_multiformat_ggplot(p_go,
                                file.path(enrich_dir, sprintf("go_%s_%s", tolower(ont), tp_name)),
                                width = 10, height = 8, use_subfolders = FALSE)
      } else {
        cat(sprintf("    No significant GO %s terms\n", ont))
      }
    }

    cat("  KEGG...\n")
    kegg_res <- tryCatch(
      compareCluster(geneCluster = gene_list, fun = "enrichKEGG",
                     organism = "mmu", universe = bg_genes,
                     pAdjustMethod = "BH", pvalueCutoff = 0.05),
      error = function(e) { cat(sprintf("    Error: %s\n", e$message)); NULL })

    if (!is.null(kegg_res) && nrow(kegg_res@compareClusterResult) > 0) {
      p_kegg <- dotplot(kegg_res, showCategory = 15) +
        labs(title = sprintf("KEGG: %s Stripes", tp_label)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      save_multiformat_ggplot(p_kegg,
                              file.path(enrich_dir, sprintf("kegg_%s", tp_name)),
                              width = 10, height = 8, use_subfolders = FALSE)
    } else {
      cat("    No significant KEGG pathways\n")
    }
  }
}

cat("Section 10 complete\n")

# =============================================================================
# SECTION 11: COMBINED SUMMARY & COMPARISON
# =============================================================================

cat("\n========================================\n")
cat("SECTION 11: Summary & Comparison\n")
cat("========================================\n\n")

count_rows <- list()
for (tp_name in names(TIMEPOINTS)) {
  tp_label <- TIMEPOINTS[[tp_name]]$label
  for (res_label in c("5kb", "10kb")) {
    df <- data_all[[tp_name]][[res_label]]
    if (is.null(df)) next
    count_rows[[length(count_rows) + 1]] <- data.frame(
      Dataset = paste(tp_label, res_label),
      Total = nrow(df),
      Significant = sum(df$significant_FDR05 == TRUE, na.rm = TRUE),
      Pct_sig = round(100 * sum(df$significant_FDR05 == TRUE, na.rm = TRUE) / nrow(df), 1),
      Lost = sum(df$direction == "lost", na.rm = TRUE),
      Gained = sum(df$direction == "gained", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
if (length(count_rows) > 0) {
  count_df <- bind_rows(count_rows)

  p_sig <- ggplot(count_df, aes(x = Dataset, y = Pct_sig, fill = Dataset)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Significant, Pct_sig)),
              vjust = -0.2, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Significance Rate by Dataset",
         x = "", y = "% Significant (FDR<0.05)") +
    theme_stripe() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 30, hjust = 1)) +
    ylim(0, max(count_df$Pct_sig) * 1.3)

  count_long <- count_df %>%
    select(Dataset, Lost, Gained) %>%
    pivot_longer(cols = c(Lost, Gained), names_to = "Direction", values_to = "Count")
  p_lg <- ggplot(count_long, aes(x = Dataset, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Lost" = "#4575B4", "Gained" = "#D73027")) +
    labs(title = "Lost vs Gained by Dataset", x = "", y = "Count") +
    theme_stripe() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  combined_comp <- p_sig + p_lg +
    plot_annotation(title = "Stripenn Differential Stripe Summary")
  save_multiformat_ggplot(combined_comp,
                          file.path(combined_dir, "comparison_summary"),
                          width = 14, height = 6, use_subfolders = FALSE)
}

sink_file <- file.path(combined_dir, "combined_summary.txt")
sink(sink_file)
cat("=========================================\n")
cat("STRIPENN DIFFERENTIAL STRIPE ANALYSIS SUMMARY\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat("=========================================\n\n")

for (tp_name in names(TIMEPOINTS)) {
  tp_label <- TIMEPOINTS[[tp_name]]$label
  tp_id <- TIMEPOINTS[[tp_name]]$tp
  cat(sprintf("--- %s (%s) ---\n", toupper(tp_label), tp_id))

  for (res_label in c("5kb", "10kb")) {
    df <- data_all[[tp_name]][[res_label]]
    if (is.null(df)) next
    n_sig <- sum(df$significant_FDR05 == TRUE, na.rm = TRUE)
    cat(sprintf("\n  Resolution: %s\n", res_label))
    cat(sprintf("  Total: %d | Significant: %d (%.1f%%)\n", nrow(df), n_sig, 100*n_sig/nrow(df)))
    cat(sprintf("  Lost: %d | Gained: %d | Unchanged: %d\n",
                sum(df$direction == "lost"), sum(df$direction == "gained"),
                sum(df$direction == "unchanged")))
    cat(sprintf("  High confidence (lost/gained): %d / %d\n",
                sum(df$direction == "lost" & df$direction_confidence == "high"),
                sum(df$direction == "gained" & df$direction_confidence == "high")))
    cat(sprintf("  Directional consistency: %.1f%%\n",
                100 * mean(df$direction_consistent == TRUE, na.rm = TRUE)))
    cat(sprintf("  Median |logFC| (sig): %.3f\n",
                median(abs(df$logFC[df$significant_FDR05 == TRUE]), na.rm = TRUE)))
  }

  cr <- data_all[[tp_name]][["cross_res"]]
  if (!is.null(cr)) {
    matched <- cr %>% filter(!is.na(logFC_10kb))
    cat(sprintf("\n  Cross-resolution:\n"))
    cat(sprintf("    Matched: %d\n", nrow(matched)))
    if (nrow(matched) > 0) {
      cat(sprintf("    logFC correlation: %.3f\n",
                  cor(matched$logFC, matched$logFC_10kb, use = "complete.obs")))
    }
    cat(sprintf("    both_concordant: %d\n", sum(cr$resolution_support == "both_concordant", na.rm = TRUE)))
  }

  ann <- annotated[[tp_name]]
  if (!is.null(ann) && "anchor_type" %in% names(ann)) {
    cat("\n  Anchor types (5kb):\n")
    at <- table(ann$anchor_type)
    for (a in names(sort(at, decreasing = TRUE))) {
      cat(sprintf("    %s: %d (%.1f%%)\n", a, at[a], 100*at[a]/sum(at)))
    }
  }

  cat("\n")
}

cat("=========================================\n")
cat("OUTPUT FILES:\n")
cat("  Per timepoint: visualizations/{volcano,stripiness,length,direction,cross_res,replicate,anchor}.*\n")
cat("  Per timepoint: visualizations/{tp}_stripes_{highconf,allsig,concordant}.bedpe\n")
cat("  Per timepoint: visualizations/{tp}_stripes_{diagonal,rectangle}.bedpe\n")
cat("  Per timepoint: visualizations/{tp}_annotated_stripes.tsv\n")
cat("  Combined: combined/{volcano,length,source_direction,cross_res,comparison}_*\n")
cat("  Enrichment: combined/enrichment/go_{bp,cc,mf}_{tp}.* kegg_{tp}.*\n")
cat("=========================================\n")
sink()
cat(sprintf("Summary written: %s\n", sink_file))

cat("\n========================================\n")
cat("STRIPENN VISUALIZATION COMPLETE\n")
cat("========================================\n\n")
