# abc/scripts/step12_activity_contact_scatter.R
# Scatter plots decomposing ΔABC into activity vs contact components.
#
# Each enhancer-gene pair's ABC change can be driven by changes in enhancer
# activity (ATAC × H3K27ac), Hi-C contact frequency, or both. Plotting
# δactivity vs δcontact, colored by enhancer class, reveals whether specific
# enhancer classes are driven more by activity vs contact changes.
#
# Inputs (relative to abc/ working directory):
#   results/delta_abc_all_pairs.tsv                           — 180K E-G pairs
#   results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv — 21.5K classified enhancers
#
# Outputs:
#   results/figures/activity_contact_scatter/   — 4 scatter plots (PDF+SVG+JPG)
#   results/figures/activity_contact_scatter/activity_contact_summary.tsv
#
# Usage:
#   cd abc && Rscript scripts/step12_activity_contact_scatter.R

# =============================================================================
# CONFIGURATION
# =============================================================================

ABC_PAIRS_FILE  <- "results/delta_abc_all_pairs.tsv"
ENH_CLASS_FILE  <- "results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv"
OUTPUT_DIR      <- "results/figures/activity_contact_scatter"

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

cat("================================================================================\n")
cat("STEP 12: ACTIVITY vs CONTACT SCATTER PLOTS\n")
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
  library(dplyr)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
stopifnot(file.exists(ABC_PAIRS_FILE))
cat(sprintf("  [OK] ABC pairs: %s\n", ABC_PAIRS_FILE))
stopifnot(file.exists(ENH_CLASS_FILE))
cat(sprintf("  [OK] Enhancer classes: %s\n", ENH_CLASS_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

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

UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

enh_class <- read.table(ENH_CLASS_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
cat(sprintf("  Classified enhancers: %d\n", nrow(enh_class)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_class$enhancer_class == cls)))
}

# =============================================================================
# COMPUTE DELTAS
# =============================================================================

cat("\nComputing delta columns...\n")

abc$delta_activity <- abc$activity_base_KO - abc$activity_base_WT
abc$delta_contact  <- abc$hic_contact_pl_scaled_adj_KO - abc$hic_contact_pl_scaled_adj_WT

# Log2FC with pseudocount (half of minimum nonzero value per column)
pseudo_act <- min(c(
  abc$activity_base_WT[abc$activity_base_WT > 0],
  abc$activity_base_KO[abc$activity_base_KO > 0]
)) / 2
pseudo_con <- min(c(
  abc$hic_contact_pl_scaled_adj_WT[abc$hic_contact_pl_scaled_adj_WT > 0],
  abc$hic_contact_pl_scaled_adj_KO[abc$hic_contact_pl_scaled_adj_KO > 0]
)) / 2

cat(sprintf("  Pseudocounts: activity = %.6f, contact = %.6f\n",
            pseudo_act, pseudo_con))

abc$log2fc_activity <- log2((abc$activity_base_KO + pseudo_act) /
                            (abc$activity_base_WT + pseudo_act))
abc$log2fc_contact  <- log2((abc$hic_contact_pl_scaled_adj_KO + pseudo_con) /
                            (abc$hic_contact_pl_scaled_adj_WT + pseudo_con))

cat(sprintf("  delta_activity range: [%.4f, %.4f]\n",
            min(abc$delta_activity), max(abc$delta_activity)))
cat(sprintf("  delta_contact  range: [%.4f, %.4f]\n",
            min(abc$delta_contact), max(abc$delta_contact)))
cat(sprintf("  log2fc_activity range: [%.4f, %.4f]\n",
            min(abc$log2fc_activity, na.rm = TRUE),
            max(abc$log2fc_activity, na.rm = TRUE)))
cat(sprintf("  log2fc_contact  range: [%.4f, %.4f]\n",
            min(abc$log2fc_contact, na.rm = TRUE),
            max(abc$log2fc_contact, na.rm = TRUE)))

# =============================================================================
# JOIN ENHANCER CLASSES
# =============================================================================

cat("\nJoining enhancer classes via overlap...\n")

# Coordinates differ between files (ABC uses 500bp padded peaks,
# class file uses raw narrowPeak coords) — use overlap-based join
abc_gr <- GRanges(seqnames = abc$chr,
                  ranges = IRanges(start = abc$start, end = abc$end))
class_gr <- GRanges(seqnames = enh_class$chr,
                    ranges = IRanges(start = enh_class$start,
                                    end = enh_class$end))

# Find overlaps (each ABC enhancer → best-matching class enhancer)
hits <- findOverlaps(abc_gr, class_gr)

# For ABC enhancers with multiple class hits, take first (classes are exclusive)
first_hits <- hits[!duplicated(queryHits(hits))]

abc$enhancer_class <- NA_character_
abc$enhancer_class[queryHits(first_hits)] <-
  enh_class$enhancer_class[subjectHits(first_hits)]

n_classified <- sum(!is.na(abc$enhancer_class))
n_unclassified <- sum(is.na(abc$enhancer_class))
cat(sprintf("  Classified pairs: %d (%.1f%%)\n", n_classified,
            100 * n_classified / nrow(abc)))
cat(sprintf("  Unclassified pairs: %d (%.1f%%)\n", n_unclassified,
            100 * n_unclassified / nrow(abc)))

# Label unclassified
abc$enhancer_class[is.na(abc$enhancer_class)] <- "Unclassified"

# Factor with Unclassified last
abc$enhancer_class <- factor(abc$enhancer_class,
                             levels = c(CLASS_ORDER, "Unclassified"))

# Split for plotting
abc_classified <- abc[abc$enhancer_class != "Unclassified", ]
abc_classified$enhancer_class <- droplevels(abc_classified$enhancer_class)
cat(sprintf("  Classified subset: %d pairs\n", nrow(abc_classified)))

# =============================================================================
# HELPER: QUADRANT LABEL POSITIONS
# =============================================================================

# Assign quadrant labels to each point for annotation counts
assign_quadrant <- function(dx, dy) {
  ifelse(dx > 0 & dy > 0, "Q1_up_up",
  ifelse(dx < 0 & dy > 0, "Q2_down_up",
  ifelse(dx < 0 & dy < 0, "Q3_down_down",
  ifelse(dx > 0 & dy < 0, "Q4_up_down", "origin"))))
}

# Build quadrant annotation for a scatter plot
make_quadrant_labels <- function(x, y, x_range, y_range) {
  quads <- assign_quadrant(x, y)
  counts <- table(quads)
  total <- length(quads)

  # Position labels at corners of plot range
  pad_x <- diff(x_range) * 0.04
  pad_y <- diff(y_range) * 0.04

  data.frame(
    label = c(
      sprintf("n=%d\n(%.0f%%)", counts["Q2_down_up"], 100 * counts["Q2_down_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q1_up_up"], 100 * counts["Q1_up_up"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q3_down_down"], 100 * counts["Q3_down_down"] / total),
      sprintf("n=%d\n(%.0f%%)", counts["Q4_up_down"], 100 * counts["Q4_up_down"] / total)
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

# =============================================================================
# PLOT 1: RAW DELTA — CLASSIFIED ONLY
# =============================================================================

cat("\n--- Plot 1: Raw delta, classified only ---\n")

# Compute per-class Spearman correlations
class_cors <- abc_classified %>%
  group_by(enhancer_class) %>%
  summarise(
    rho = cor(delta_activity, delta_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  )
cat("  Per-class Spearman(δact, δcontact):\n")
for (i in seq_len(nrow(class_cors))) {
  cat(sprintf("    %s: rho=%.3f (n=%d)\n",
              class_cors$enhancer_class[i], class_cors$rho[i], class_cors$n[i]))
}

# Clip outliers for axis limits (99.5th percentile symmetric)
clip_raw <- quantile(c(abs(abc_classified$delta_activity),
                       abs(abc_classified$delta_contact)),
                     0.995, na.rm = TRUE)
x_lim_raw <- c(-clip_raw, clip_raw)
y_lim_raw <- c(-clip_raw, clip_raw)

# Quadrant annotations
q_labels_raw <- make_quadrant_labels(
  abc_classified$delta_activity, abc_classified$delta_contact,
  x_lim_raw, y_lim_raw
)

# Correlation label
rho_all_raw <- cor(abc_classified$delta_activity,
                   abc_classified$delta_contact, method = "spearman")

p1 <- ggplot(abc_classified,
             aes(x = delta_activity, y = delta_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_text(data = q_labels_raw, aes(x = x, y = y, label = label,
            hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 2.8, color = "grey30") +
  scale_color_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_raw, ylim = y_lim_raw) +
  labs(
    title = "Activity vs Contact Changes (Classified Enhancers)",
    subtitle = sprintf("Spearman rho = %.3f | n = %s E-G pairs",
                       rho_all_raw, format(nrow(abc_classified), big.mark = ",")),
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression(Delta * "Contact (KO - WT)")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p1, "01_raw_delta_classified", w = 7, h = 6.5)

# =============================================================================
# PLOT 2: RAW DELTA — ALL PAIRS (UNCLASSIFIED AS GRAY REFERENCE)
# =============================================================================

cat("\n--- Plot 2: Raw delta, all pairs ---\n")

# Plot unclassified first, then classified on top
abc_ordered <- rbind(
  abc[abc$enhancer_class == "Unclassified", ],
  abc[abc$enhancer_class != "Unclassified", ]
)

all_colors <- c(CLASS_COLORS, Unclassified = "grey80")

clip_all <- quantile(c(abs(abc$delta_activity), abs(abc$delta_contact)),
                     0.995, na.rm = TRUE)
x_lim_all <- c(-clip_all, clip_all)
y_lim_all <- c(-clip_all, clip_all)

q_labels_all <- make_quadrant_labels(
  abc$delta_activity, abc$delta_contact, x_lim_all, y_lim_all
)

rho_all <- cor(abc$delta_activity, abc$delta_contact, method = "spearman")

p2 <- ggplot(abc_ordered,
             aes(x = delta_activity, y = delta_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_text(data = q_labels_all, aes(x = x, y = y, label = label,
            hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 2.8, color = "grey30") +
  scale_color_manual(values = all_colors, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_all, ylim = y_lim_all) +
  labs(
    title = "Activity vs Contact Changes (All E-G Pairs)",
    subtitle = sprintf("Spearman rho = %.3f | n = %s total (%s classified)",
                       rho_all,
                       format(nrow(abc), big.mark = ","),
                       format(nrow(abc_classified), big.mark = ",")),
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression(Delta * "Contact (KO - WT)")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p2, "02_raw_delta_all_pairs", w = 7, h = 6.5)

# =============================================================================
# PLOT 3: LOG2FC — CLASSIFIED ONLY
# =============================================================================

cat("\n--- Plot 3: Log2FC, classified only ---\n")

# Remove infinite/NA log2fc values
valid_lfc <- abc_classified[is.finite(abc_classified$log2fc_activity) &
                            is.finite(abc_classified$log2fc_contact), ]
cat(sprintf("  Valid log2FC pairs (classified): %d / %d\n",
            nrow(valid_lfc), nrow(abc_classified)))

class_cors_lfc <- valid_lfc %>%
  group_by(enhancer_class) %>%
  summarise(
    rho = cor(log2fc_activity, log2fc_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  )
cat("  Per-class Spearman(log2FC_act, log2FC_contact):\n")
for (i in seq_len(nrow(class_cors_lfc))) {
  cat(sprintf("    %s: rho=%.3f (n=%d)\n",
              class_cors_lfc$enhancer_class[i],
              class_cors_lfc$rho[i], class_cors_lfc$n[i]))
}

clip_lfc <- quantile(c(abs(valid_lfc$log2fc_activity),
                       abs(valid_lfc$log2fc_contact)),
                     0.995, na.rm = TRUE)
x_lim_lfc <- c(-clip_lfc, clip_lfc)
y_lim_lfc <- c(-clip_lfc, clip_lfc)

q_labels_lfc <- make_quadrant_labels(
  valid_lfc$log2fc_activity, valid_lfc$log2fc_contact,
  x_lim_lfc, y_lim_lfc
)

rho_lfc_class <- cor(valid_lfc$log2fc_activity, valid_lfc$log2fc_contact,
                     method = "spearman")

p3 <- ggplot(valid_lfc,
             aes(x = log2fc_activity, y = log2fc_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_text(data = q_labels_lfc, aes(x = x, y = y, label = label,
            hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 2.8, color = "grey30") +
  scale_color_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_lfc, ylim = y_lim_lfc) +
  labs(
    title = "Log2FC Activity vs Contact (Classified Enhancers)",
    subtitle = sprintf("Spearman rho = %.3f | n = %s E-G pairs",
                       rho_lfc_class, format(nrow(valid_lfc), big.mark = ",")),
    x = expression(log[2] * "FC Activity (KO/WT)"),
    y = expression(log[2] * "FC Contact (KO/WT)")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p3, "03_log2fc_classified", w = 7, h = 6.5)

# =============================================================================
# PLOT 4: LOG2FC — ALL PAIRS
# =============================================================================

cat("\n--- Plot 4: Log2FC, all pairs ---\n")

valid_all <- abc[is.finite(abc$log2fc_activity) &
                 is.finite(abc$log2fc_contact), ]
cat(sprintf("  Valid log2FC pairs (all): %d / %d\n",
            nrow(valid_all), nrow(abc)))

valid_all_ordered <- rbind(
  valid_all[valid_all$enhancer_class == "Unclassified", ],
  valid_all[valid_all$enhancer_class != "Unclassified", ]
)

clip_lfc_all <- quantile(c(abs(valid_all$log2fc_activity),
                           abs(valid_all$log2fc_contact)),
                         0.995, na.rm = TRUE)
x_lim_lfc_all <- c(-clip_lfc_all, clip_lfc_all)
y_lim_lfc_all <- c(-clip_lfc_all, clip_lfc_all)

q_labels_lfc_all <- make_quadrant_labels(
  valid_all$log2fc_activity, valid_all$log2fc_contact,
  x_lim_lfc_all, y_lim_lfc_all
)

rho_lfc_all <- cor(valid_all$log2fc_activity, valid_all$log2fc_contact,
                   method = "spearman")

p4 <- ggplot(valid_all_ordered,
             aes(x = log2fc_activity, y = log2fc_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_text(data = q_labels_lfc_all, aes(x = x, y = y, label = label,
            hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 2.8, color = "grey30") +
  scale_color_manual(values = all_colors, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_lfc_all, ylim = y_lim_lfc_all) +
  labs(
    title = "Log2FC Activity vs Contact (All E-G Pairs)",
    subtitle = sprintf("Spearman rho = %.3f | n = %s total (%s classified)",
                       rho_lfc_all,
                       format(nrow(valid_all), big.mark = ","),
                       format(sum(valid_all$enhancer_class != "Unclassified"),
                              big.mark = ",")),
    x = expression(log[2] * "FC Activity (KO/WT)"),
    y = expression(log[2] * "FC Contact (KO/WT)")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p4, "04_log2fc_all_pairs", w = 7, h = 6.5)

# =============================================================================
# FACETED PANELS BY CLASS
# =============================================================================

cat("\n--- Faceted plots by enhancer class ---\n")

# Faceted raw delta
p5 <- ggplot(abc_classified,
             aes(x = delta_activity, y = delta_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.2, size = 0.4) +
  facet_wrap(~enhancer_class, scales = "free") +
  scale_color_manual(values = CLASS_COLORS, guide = "none") +
  labs(
    title = "Activity vs Contact Changes by Enhancer Class",
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression(Delta * "Contact (KO - WT)")
  ) +
  theme_pub

# Add per-facet correlation labels
cor_labels_raw <- class_cors %>%
  mutate(label = sprintf("rho = %.3f\nn = %s", rho, format(n, big.mark = ",")))

p5 <- p5 +
  geom_text(data = cor_labels_raw,
            aes(label = label), x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.2, size = 3, color = "black",
            inherit.aes = FALSE)

save_plot(p5, "05_raw_delta_faceted", w = 10, h = 8)

# Faceted log2FC
p6 <- ggplot(valid_lfc,
             aes(x = log2fc_activity, y = log2fc_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.2, size = 0.4) +
  facet_wrap(~enhancer_class, scales = "free") +
  scale_color_manual(values = CLASS_COLORS, guide = "none") +
  labs(
    title = "Log2FC Activity vs Contact by Enhancer Class",
    x = expression(log[2] * "FC Activity (KO/WT)"),
    y = expression(log[2] * "FC Contact (KO/WT)")
  ) +
  theme_pub

cor_labels_lfc <- class_cors_lfc %>%
  mutate(label = sprintf("rho = %.3f\nn = %s", rho, format(n, big.mark = ",")))

p6 <- p6 +
  geom_text(data = cor_labels_lfc,
            aes(label = label), x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.2, size = 3, color = "black",
            inherit.aes = FALSE)

save_plot(p6, "06_log2fc_faceted", w = 10, h = 8)

# =============================================================================
# MARGINAL DENSITY RIDGES BY CLASS
# =============================================================================

cat("\n--- Marginal density plots ---\n")

# Activity marginal
p_act_dens <- ggplot(abc_classified,
                     aes(x = delta_activity, fill = enhancer_class)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_raw) +
  labs(title = "Distribution of Activity Changes by Class",
       x = expression(Delta * "Activity (KO - WT)"), y = "Density") +
  theme_pub

# Contact marginal
p_con_dens <- ggplot(abc_classified,
                     aes(x = delta_contact, fill = enhancer_class)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = y_lim_raw) +
  labs(title = "Distribution of Contact Changes by Class",
       x = expression(Delta * "Contact (KO - WT)"), y = "Density") +
  theme_pub

p_marginals <- p_act_dens / p_con_dens +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_marginals, "07_marginal_densities", w = 8, h = 8)

# =============================================================================
# STATISTICAL SUMMARY TABLE
# =============================================================================

cat("\n--- Computing summary statistics ---\n")

compute_class_stats <- function(df, x_col, y_col, class_col) {
  classes <- unique(df[[class_col]])
  results <- lapply(classes, function(cls) {
    sub <- df[df[[class_col]] == cls, ]
    x <- sub[[x_col]]
    y <- sub[[y_col]]
    quads <- assign_quadrant(x, y)
    qtab <- table(factor(quads, levels = c("Q1_up_up", "Q2_down_up",
                                           "Q3_down_down", "Q4_up_down",
                                           "origin")))
    n <- nrow(sub)
    data.frame(
      enhancer_class = cls,
      n_pairs = n,
      median_delta_x = median(x, na.rm = TRUE),
      median_delta_y = median(y, na.rm = TRUE),
      spearman_rho = cor(x, y, method = "spearman", use = "complete.obs"),
      pct_Q1_up_up     = 100 * qtab["Q1_up_up"] / n,
      pct_Q2_down_up   = 100 * qtab["Q2_down_up"] / n,
      pct_Q3_down_down = 100 * qtab["Q3_down_down"] / n,
      pct_Q4_up_down   = 100 * qtab["Q4_up_down"] / n,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

# Raw delta stats
stats_raw <- compute_class_stats(abc, "delta_activity", "delta_contact",
                                 "enhancer_class")
stats_raw$metric <- "raw_delta"

# Log2FC stats
stats_lfc <- compute_class_stats(
  valid_all, "log2fc_activity", "log2fc_contact", "enhancer_class"
)
stats_lfc$metric <- "log2fc"

# Combine and rename columns for clarity
stats_all <- rbind(stats_raw, stats_lfc)
names(stats_all) <- gsub("delta_x", "activity", gsub("delta_y", "contact",
                          names(stats_all)))

out_tsv <- file.path(OUTPUT_DIR, "activity_contact_summary.tsv")
write.table(stats_all, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Summary table saved: %s\n", out_tsv))

# Print summary
cat("\n--- Summary: Raw delta ---\n")
print(stats_raw[, c("enhancer_class", "n_pairs", "median_delta_x",
                     "median_delta_y", "spearman_rho",
                     "pct_Q1_up_up", "pct_Q3_down_down")],
      row.names = FALSE)

cat("\n--- Summary: Log2FC ---\n")
print(stats_lfc[, c("enhancer_class", "n_pairs", "median_delta_x",
                     "median_delta_y", "spearman_rho",
                     "pct_Q1_up_up", "pct_Q3_down_down")],
      row.names = FALSE)

cat("\n================================================================================\n")
cat("STEP 12 COMPLETE\n")
cat(sprintf("Outputs in: %s\n", OUTPUT_DIR))
cat("================================================================================\n")
