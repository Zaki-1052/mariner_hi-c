# abc/scripts/step12b_promoter_distal_scatter.R
#
# Stratify activity-vs-contact scatter by ABC `class` column (promoter vs distal).
# The ABC class column (promoter/genic/intergenic) reflects genomic location of
# the enhancer element, orthogonal to the phenotypic enhancer_class from step10.
# Self-promoters already removed in step6; remaining "promoter" = ATAC peak at
# one gene's promoter acting as enhancer for a DIFFERENT gene.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
})

cat("=== Step 12b: Promoter vs Distal Activity-Contact Scatter ===\n")
cat(sprintf("  Started: %s\n\n", Sys.time()))

# =============================================================================
# CONFIGURATION
# =============================================================================

ABC_PAIRS_FILE <- "results/delta_abc_all_pairs.tsv"
ENH_CLASS_FILE <- "results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv"
OUTPUT_DIR     <- "results/figures/activity_contact_scatter"
UTIL_PATH      <- "../scripts/utils/multi_format_output.R"

stopifnot(file.exists(ABC_PAIRS_FILE))
stopifnot(file.exists(ENH_CLASS_FILE))
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared aesthetics (identical to step12)
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# Quadrant assignment (reused from step12)
assign_quadrant <- function(dx, dy) {
  ifelse(dx > 0 & dy > 0, "Q1_up_up",
  ifelse(dx < 0 & dy > 0, "Q2_down_up",
  ifelse(dx < 0 & dy < 0, "Q3_down_down",
  ifelse(dx > 0 & dy < 0, "Q4_up_down", "origin"))))
}

make_quadrant_labels <- function(x, y, x_range, y_range) {
  quads <- assign_quadrant(x, y)
  counts <- table(quads)
  total <- length(quads)

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
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

enh_class <- read.table(ENH_CLASS_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
cat(sprintf("  Classified enhancers: %d\n", nrow(enh_class)))

# ABC class distribution
cat("\n  ABC class column distribution:\n")
for (cls in sort(unique(abc$class))) {
  cat(sprintf("    %s: %s\n", cls, format(sum(abc$class == cls), big.mark = ",")))
}

# =============================================================================
# COMPUTE DELTAS (same as step12)
# =============================================================================

cat("\nComputing delta columns...\n")

abc$delta_activity <- abc$activity_base_KO - abc$activity_base_WT
abc$delta_contact  <- abc$hic_contact_pl_scaled_adj_KO - abc$hic_contact_pl_scaled_adj_WT

# Log2FC with pseudocount (half of minimum nonzero value)
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

# =============================================================================
# JOIN ENHANCER CLASSES (same overlap approach as step12)
# =============================================================================

cat("\nJoining enhancer classes via overlap...\n")

abc_gr <- GRanges(seqnames = abc$chr,
                  ranges = IRanges(start = abc$start, end = abc$end))
class_gr <- GRanges(seqnames = enh_class$chr,
                    ranges = IRanges(start = enh_class$start,
                                    end = enh_class$end))

hits <- findOverlaps(abc_gr, class_gr)
first_hits <- hits[!duplicated(queryHits(hits))]

abc$enhancer_class <- NA_character_
abc$enhancer_class[queryHits(first_hits)] <-
  enh_class$enhancer_class[subjectHits(first_hits)]

n_classified <- sum(!is.na(abc$enhancer_class))
cat(sprintf("  Classified pairs: %d (%.1f%%)\n", n_classified,
            100 * n_classified / nrow(abc)))

abc$enhancer_class[is.na(abc$enhancer_class)] <- "Unclassified"
abc$enhancer_class <- factor(abc$enhancer_class,
                             levels = c(CLASS_ORDER, "Unclassified"))

# =============================================================================
# BINARY LOCATION: PROMOTER vs DISTAL
# =============================================================================

cat("\nCreating binary location variable...\n")

abc$location <- ifelse(abc$class == "promoter", "Promoter", "Distal")
abc$location <- factor(abc$location, levels = c("Promoter", "Distal"))

cat(sprintf("  Promoter: %s | Distal: %s\n",
            format(sum(abc$location == "Promoter"), big.mark = ","),
            format(sum(abc$location == "Distal"), big.mark = ",")))

# Classified subset only
abc_classified <- abc[abc$enhancer_class != "Unclassified", ]
abc_classified$enhancer_class <- droplevels(abc_classified$enhancer_class)
cat(sprintf("  Classified subset: %d pairs\n", nrow(abc_classified)))
cat(sprintf("    Promoter classified: %d | Distal classified: %d\n",
            sum(abc_classified$location == "Promoter"),
            sum(abc_classified$location == "Distal")))

# =============================================================================
# PLOT 08: RAW DELTA — FACETED BY PROMOTER/DISTAL
# =============================================================================

cat("\n--- Plot 08: Raw delta, promoter vs distal ---\n")

# Per-location + per-class correlations
loc_class_cors_raw <- abc_classified %>%
  group_by(location, enhancer_class) %>%
  summarise(
    rho = cor(delta_activity, delta_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  )
cat("  Per-location × class Spearman(δact, δcontact):\n")
for (i in seq_len(nrow(loc_class_cors_raw))) {
  cat(sprintf("    %s | %s: rho=%.3f (n=%d)\n",
              loc_class_cors_raw$location[i],
              loc_class_cors_raw$enhancer_class[i],
              loc_class_cors_raw$rho[i],
              loc_class_cors_raw$n[i]))
}

# Per-facet overall rho for subtitles
loc_cors_raw <- abc_classified %>%
  group_by(location) %>%
  summarise(
    rho = cor(delta_activity, delta_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%s: rho = %.3f (n = %s)",
                         location, rho, format(n, big.mark = ",")))

# Clip outliers (99.5th percentile symmetric)
clip_raw <- quantile(c(abs(abc_classified$delta_activity),
                       abs(abc_classified$delta_contact)),
                     0.995, na.rm = TRUE)
x_lim_raw <- c(-clip_raw, clip_raw)
y_lim_raw <- c(-clip_raw, clip_raw)

p8 <- ggplot(abc_classified,
             aes(x = delta_activity, y = delta_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.6) +
  facet_wrap(~ location, ncol = 2) +
  scale_color_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_raw, ylim = y_lim_raw) +
  labs(
    title = "Activity vs Contact Changes by Genomic Location",
    subtitle = paste(loc_cors_raw$label, collapse = " | "),
    x = expression(Delta * "Activity (KO - WT)"),
    y = expression(Delta * "Contact (KO - WT)")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p8, "08_raw_delta_promoter_distal", w = 12, h = 6.5)

# =============================================================================
# PLOT 09: LOG2FC — FACETED BY PROMOTER/DISTAL
# =============================================================================

cat("\n--- Plot 09: Log2FC, promoter vs distal ---\n")

# Filter infinite/NA log2FC values
abc_log2 <- abc_classified[is.finite(abc_classified$log2fc_activity) &
                           is.finite(abc_classified$log2fc_contact), ]

loc_class_cors_log <- abc_log2 %>%
  group_by(location, enhancer_class) %>%
  summarise(
    rho = cor(log2fc_activity, log2fc_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  )
cat("  Per-location × class Spearman(log2FC_act, log2FC_contact):\n")
for (i in seq_len(nrow(loc_class_cors_log))) {
  cat(sprintf("    %s | %s: rho=%.3f (n=%d)\n",
              loc_class_cors_log$location[i],
              loc_class_cors_log$enhancer_class[i],
              loc_class_cors_log$rho[i],
              loc_class_cors_log$n[i]))
}

loc_cors_log <- abc_log2 %>%
  group_by(location) %>%
  summarise(
    rho = cor(log2fc_activity, log2fc_contact, method = "spearman"),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%s: rho = %.3f (n = %s)",
                         location, rho, format(n, big.mark = ",")))

clip_log <- quantile(c(abs(abc_log2$log2fc_activity),
                       abs(abc_log2$log2fc_contact)),
                     0.995, na.rm = TRUE)
x_lim_log <- c(-clip_log, clip_log)
y_lim_log <- c(-clip_log, clip_log)

p9 <- ggplot(abc_log2,
             aes(x = log2fc_activity, y = log2fc_contact,
                 color = enhancer_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.6) +
  facet_wrap(~ location, ncol = 2) +
  scale_color_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  coord_cartesian(xlim = x_lim_log, ylim = y_lim_log) +
  labs(
    title = "Activity vs Contact Changes (log2FC) by Genomic Location",
    subtitle = paste(loc_cors_log$label, collapse = " | "),
    x = expression(log[2] * "FC Activity"),
    y = expression(log[2] * "FC Contact")
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p9, "09_log2fc_promoter_distal", w = 12, h = 6.5)

# =============================================================================
# PLOT 10: CLASS COMPOSITION BAR — BY ABC CLASS
# =============================================================================

cat("\n--- Plot 10: Enhancer class composition by ABC class ---\n")

# Use all 3 ABC classes (promoter/genic/intergenic), not binary
abc_for_comp <- abc[abc$enhancer_class != "Unclassified", ]
abc_for_comp$enhancer_class <- droplevels(abc_for_comp$enhancer_class)

comp_df <- abc_for_comp %>%
  group_by(class, enhancer_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(class) %>%
  mutate(
    total = sum(n),
    pct = 100 * n / total,
    label = sprintf("%.0f%%", pct)
  ) %>%
  ungroup()

# Order ABC classes for display
comp_df$class <- factor(comp_df$class,
                        levels = c("promoter", "genic", "intergenic"))

# Fisher's test: enhancer_class distribution differs by location?
cont_tbl <- table(abc_for_comp$location, abc_for_comp$enhancer_class)
fisher_loc <- fisher.test(cont_tbl, simulate.p.value = TRUE, B = 10000)
cat(sprintf("  Fisher's exact (location × enhancer_class): %s\n",
            fmt_p(fisher_loc$p.value)))

# Print composition table
cat("  Composition:\n")
for (cls in levels(comp_df$class)) {
  sub <- comp_df[comp_df$class == cls, ]
  cat(sprintf("    %s (n=%s):", cls,
              format(sub$total[1], big.mark = ",")))
  for (j in seq_len(nrow(sub))) {
    cat(sprintf(" %s=%.1f%%", sub$enhancer_class[j], sub$pct[j]))
  }
  cat("\n")
}

p10 <- ggplot(comp_df, aes(x = class, y = pct, fill = enhancer_class)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = CLASS_COLORS, name = "Enhancer Class") +
  labs(
    title = "Enhancer Class Composition by ABC Genomic Class",
    subtitle = sprintf("Fisher's exact (Promoter vs Distal): %s",
                       fmt_p(fisher_loc$p.value)),
    x = "ABC Class", y = "Percentage"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 11))

save_plot(p10, "10_class_composition_bar", w = 7, h = 6)

# =============================================================================
# SAVE CORRELATION TABLE
# =============================================================================

cat("\nSaving correlation summary...\n")

cors_combined <- bind_rows(
  loc_class_cors_raw %>% mutate(metric = "raw_delta"),
  loc_class_cors_log %>% mutate(metric = "log2fc")
)

write.table(cors_combined,
            file.path(OUTPUT_DIR, "promoter_distal_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Saved: %d rows to promoter_distal_correlations.tsv\n",
            nrow(cors_combined)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("  Total classified pairs: %d\n", nrow(abc_classified)))
cat(sprintf("    Promoter: %d | Distal: %d\n",
            sum(abc_classified$location == "Promoter"),
            sum(abc_classified$location == "Distal")))
cat(sprintf("  Fisher (location × class): %s\n", fmt_p(fisher_loc$p.value)))
cat(sprintf("  Plots saved: 08, 09, 10 in %s\n", OUTPUT_DIR))

cat(sprintf("\nCompleted: %s\n", Sys.time()))
cat("=== Step 12b complete ===\n")
