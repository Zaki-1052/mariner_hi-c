#!/usr/bin/env Rscript
# scripts/compartment_genome_percentage.R
# % of Genome with Differential PC1 Compartment Shifts
# Author: Zakir Alibhai
# Date: 2026-03-14
#
# Purpose:
#   Compute and visualize the percentage of the mm10 genome that undergoes
#   A/B compartment shifts in BAP1-KO vs wildtype. Reads pre-computed
#   significant compartment region TSVs from compartment_volcano_plot.R.
#
#   Paper panel: Figure 1D — "% of genome with differential PC1"
#
# Output:
#   - Stacked bar chart: % genome A->B vs B->A vs unchanged (both thresholds)
#   - Pie chart (binary): compartment shift proportions per threshold
#   - Pie chart (7-category): weakened/strengthened/flipped breakdown per threshold
#   - Per-chromosome bar chart: % differential by chromosome (standard threshold)
#   - Summary statistics text file
#   - Summary table TSV
#
# Usage:
#   Rscript scripts/compartment_genome_percentage.R
#
# =============================================================================
# CONFIGURATION
# =============================================================================

cat("\n")
cat("================================================\n")
cat("% of Genome with Differential PC1\n")
cat("================================================\n\n")

# mm10 genome size (assembled autosomes + chrX)
MM10_GENOME_SIZE <- 2730871774

# Input files (from compartment_volcano_plot.R)
# Original: BASE_DIR <- "tads/tad-pc-analysis/output/compartment_analysis"
TSV_DIR <- "data/tsvs/figure_1_tads_boundaries_compartments"
PLOT_DIR <- "data/plots/figure_1_tads_boundaries_compartments"
INPUT_FILES <- list(
  # Original: file.path(BASE_DIR, "compartment_significant_standard.tsv")
  standard = file.path(TSV_DIR, "1D_compartment_significant_standard.tsv"),
  # Original: file.path(BASE_DIR, "compartment_significant_relaxed.tsv")
  relaxed  = file.path(TSV_DIR, "1D_compartment_significant_relaxed.tsv"),
  # Original: file.path(BASE_DIR, "compartment_all_annotated.tsv")
  all      = file.path(TSV_DIR, "1D_compartment_all_annotated.tsv")
)
# Original: OUTPUT_DIR <- BASE_DIR
dir.create(TSV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Validate all inputs exist
for (name in names(INPUT_FILES)) {
  if (!file.exists(INPUT_FILES[[name]])) {
    stop(sprintf("Input file not found: %s", INPUT_FILES[[name]]))
  }
}

# =============================================================================
# LOAD PACKAGES AND UTILITIES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Original: source("scripts/utils/multi_format_output.R")
source("data/scripts/_shared/multi_format_output.R")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading compartment data...\n")

sig_standard <- read.delim(INPUT_FILES$standard, stringsAsFactors = FALSE)
sig_relaxed  <- read.delim(INPUT_FILES$relaxed, stringsAsFactors = FALSE)
all_regions  <- read.delim(INPUT_FILES$all, stringsAsFactors = FALSE)

# Validate expected columns
required_cols <- c("Region_size", "direction", "Chr", "ctrl_avg_PC1", "mut_avg_PC1")
stopifnot("Missing columns in standard TSV" = all(required_cols %in% names(sig_standard)))
stopifnot("Missing columns in relaxed TSV"  = all(required_cols %in% names(sig_relaxed)))
stopifnot("Missing columns in all TSV"      = all(c("Region_size", "Chr", "ctrl_avg_PC1", "mut_avg_PC1") %in% names(all_regions)))

cat(sprintf("  All regions analyzed: %d\n", nrow(all_regions)))
cat(sprintf("  Significant (standard): %d\n", nrow(sig_standard)))
cat(sprintf("  Significant (relaxed):  %d\n", nrow(sig_relaxed)))

# =============================================================================
# COMPUTE GENOME PERCENTAGES
# =============================================================================

cat("\nComputing genome percentages...\n")

compute_pct <- function(sig_df, label) {
  a_to_b_bp <- sum(sig_df$Region_size[sig_df$direction == "A_to_B_in_Mutant"])
  b_to_a_bp <- sum(sig_df$Region_size[sig_df$direction == "B_to_A_in_Mutant"])
  total_sig_bp <- a_to_b_bp + b_to_a_bp
  unchanged_bp <- MM10_GENOME_SIZE - total_sig_bp

  data.frame(
    threshold     = label,
    category      = c("A->B (More Inactive)", "B->A (More Active)", "Unchanged"),
    bp            = c(a_to_b_bp, b_to_a_bp, unchanged_bp),
    pct_genome    = c(a_to_b_bp, b_to_a_bp, unchanged_bp) / MM10_GENOME_SIZE * 100,
    n_regions     = c(
      sum(sig_df$direction == "A_to_B_in_Mutant"),
      sum(sig_df$direction == "B_to_A_in_Mutant"),
      nrow(all_regions) - nrow(sig_df)
    ),
    stringsAsFactors = FALSE
  )
}

pct_standard <- compute_pct(sig_standard, "Standard\n(FDR<0.05, |Diff|>0.30)")
pct_relaxed  <- compute_pct(sig_relaxed, "Relaxed\n(FDR<0.15, |Diff|>0.15)")
pct_combined <- rbind(pct_standard, pct_relaxed)

# Set factor order for plotting
pct_combined$category <- factor(pct_combined$category,
                                levels = c("Unchanged", "A->B (More Inactive)", "B->A (More Active)"))
pct_combined$threshold <- factor(pct_combined$threshold,
                                 levels = c("Standard\n(FDR<0.05, |Diff|>0.30)",
                                            "Relaxed\n(FDR<0.15, |Diff|>0.15)"))

# Color scheme
shift_colors <- c(
  "A->B (More Inactive)" = "#4393C3",
  "B->A (More Active)"   = "#D6604D",
  "Unchanged"            = "#E0E0E0"
)

# Print results
for (thresh in c("Standard\n(FDR<0.05, |Diff|>0.30)", "Relaxed\n(FDR<0.15, |Diff|>0.15)")) {
  sub <- pct_combined[pct_combined$threshold == thresh, ]
  cat(sprintf("\n  %s:\n", gsub("\n", " ", thresh)))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %s: %.1f Mb (%.2f%% of genome, %d regions)\n",
                sub$category[i],
                sub$bp[i] / 1e6,
                sub$pct_genome[i],
                sub$n_regions[i]))
  }
}

# =============================================================================
# FIGURE 1: STACKED BAR CHART (primary paper figure)
# =============================================================================

cat("\nGenerating stacked bar chart...\n")

# Subset to significant categories only (exclude "Unchanged" from bar, show as annotation)
pct_sig <- pct_combined[pct_combined$category != "Unchanged", ]

p_bar <- ggplot(pct_sig, aes(x = threshold, y = pct_genome, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d regions)", pct_genome, n_regions)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  # Add total % annotation above bars
  stat_summary(aes(label = sprintf("Total: %.1f%%", after_stat(y))),
               fun = sum, geom = "text", vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = shift_colors, name = "Compartment Shift") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "% of Genome with Differential Compartment (PC1)",
    subtitle = "BAP1-KO vs Wildtype | mm10 genome (2.73 Gb)",
    x = NULL,
    y = "% of Genome"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0, color = "grey40"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

save_multiformat_ggplot(p_bar,
                        file.path(PLOT_DIR, "1D_compartment_genome_pct_bar"),  # Original: file.path(OUTPUT_DIR, "compartment_genome_pct_bar")
                        width = 7, height = 7)

# =============================================================================
# FIGURE 2: PIE CHARTS (one per threshold)
# =============================================================================

cat("Generating pie charts...\n")

make_pie <- function(pct_df, title_suffix) {
  # Only show significant slices with labels; unchanged gets a simple label
  pct_df$label <- ifelse(
    pct_df$category == "Unchanged",
    sprintf("%.1f%%", pct_df$pct_genome),
    sprintf("%s\n%.1f%%\n(%d regions)", pct_df$category, pct_df$pct_genome, pct_df$n_regions)
  )

  ggplot(pct_df, aes(x = "", y = pct_genome, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.8) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5), size = 3.2) +
    scale_fill_manual(values = shift_colors, name = "Category") +
    labs(title = sprintf("Compartment Shifts — %s", title_suffix)) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
}

p_pie_std <- make_pie(pct_standard, "Standard (FDR<0.05, |Diff|>0.30)")
p_pie_rlx <- make_pie(pct_relaxed, "Relaxed (FDR<0.15, |Diff|>0.15)")

# Combine into one figure with patchwork
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_pie_combined <- p_pie_std + p_pie_rlx +
    plot_annotation(
      title = "% of Genome with Differential PC1",
      subtitle = "BAP1-KO vs Wildtype | mm10 genome (2.73 Gb)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40")
      )
    )
  save_multiformat_ggplot(p_pie_combined,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie"),  # Original: file.path(OUTPUT_DIR, "compartment_genome_pct_pie")
                          width = 12, height = 6)
} else {
  save_multiformat_ggplot(p_pie_std,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_standard"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 7, height = 6)
  save_multiformat_ggplot(p_pie_rlx,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_relaxed"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 7, height = 6)
}

# =============================================================================
# FIGURE 2B: 7-CATEGORY PIE CHARTS (weakened/strengthened/flipped)
# =============================================================================

cat("Generating 7-category pie charts...\n")

# Classify each significant region by PC1 sign and direction of change
# Classify by whether compartment identity flipped (sign change), weakened, or strengthened
# Convention: PC1 >= 0 = A compartment, PC1 < 0 = B compartment
classify_compartment <- function(ctrl_pc1, mut_pc1) {
  ctrl_is_A <- ctrl_pc1 >= 0
  mut_is_A  <- mut_pc1 >= 0
  sign_flip <- ctrl_is_A != mut_is_A
  dplyr::case_when(
    sign_flip & ctrl_is_A              ~ "Flipped A->B",
    sign_flip & !ctrl_is_A             ~ "Flipped B->A",
    !sign_flip & ctrl_is_A & mut_pc1 < ctrl_pc1  ~ "Weakened A",
    !sign_flip & ctrl_is_A & mut_pc1 >= ctrl_pc1 ~ "Strengthened A",
    !sign_flip & !ctrl_is_A & mut_pc1 > ctrl_pc1 ~ "Weakened B",
    !sign_flip & !ctrl_is_A & mut_pc1 <= ctrl_pc1 ~ "Strengthened B",
    TRUE ~ "Unclassified"
  )
}

compute_pct_7cat <- function(sig_df, all_df, label) {
  sig_df$compartment_category <- classify_compartment(sig_df$ctrl_avg_PC1,
                                                       sig_df$mut_avg_PC1)

  cat_levels <- c("Weakened A", "Strengthened A", "Flipped A->B",
                   "Weakened B", "Strengthened B", "Flipped B->A",
                   "No Change")

  # Tally bp and region counts per category
  cat_summary <- sig_df %>%
    group_by(compartment_category) %>%
    summarise(bp = sum(Region_size), n_regions = n(), .groups = "drop")

  # Build full table including "No Change"
  sig_bp <- sum(cat_summary$bp)
  unchanged_bp <- MM10_GENOME_SIZE - sig_bp
  unchanged_n  <- nrow(all_df) - nrow(sig_df)

  result <- data.frame(
    threshold = label,
    category = cat_summary$compartment_category,
    bp = cat_summary$bp,
    n_regions = cat_summary$n_regions,
    stringsAsFactors = FALSE
  )
  result <- rbind(result, data.frame(
    threshold = label,
    category = "No Change",
    bp = unchanged_bp,
    n_regions = unchanged_n,
    stringsAsFactors = FALSE
  ))

  result$pct_genome <- result$bp / MM10_GENOME_SIZE * 100
  result$category <- factor(result$category, levels = cat_levels)
  result
}

pct_7cat_std <- compute_pct_7cat(sig_standard, all_regions,
                                  "Standard (FDR<0.05, |Diff|>0.30)")
pct_7cat_rlx <- compute_pct_7cat(sig_relaxed, all_regions,
                                  "Relaxed (FDR<0.15, |Diff|>0.15)")

# Print 7-category results
for (pct_df in list(pct_7cat_std, pct_7cat_rlx)) {
  cat(sprintf("\n  7-Category — %s:\n", pct_df$threshold[1]))
  for (i in seq_len(nrow(pct_df))) {
    cat(sprintf("    %-18s: %8.1f Mb (%6.2f%%, %5d regions)\n",
                as.character(pct_df$category[i]),
                pct_df$bp[i] / 1e6,
                pct_df$pct_genome[i],
                pct_df$n_regions[i]))
  }
}

# Color scheme for 7 categories
cat_colors <- c(
  "Weakened A"      = "#FCBBA1",
  "Strengthened A"  = "#A50F15",
  "Flipped A->B"    = "#2166AC",
  "Weakened B"      = "#9ECAE1",
  "Strengthened B"  = "#08306B",
  "Flipped B->A"    = "#B2182B",
  "No Change"       = "#E0E0E0"
)

make_pie_7cat <- function(pct_df, title_suffix) {
  # Label significant slices with category + stats; "No Change" gets simple label
  pct_df$label <- ifelse(
    pct_df$category == "No Change",
    sprintf("%.1f%%", pct_df$pct_genome),
    sprintf("%s\n%.2f%%\n(%d)", pct_df$category, pct_df$pct_genome, pct_df$n_regions)
  )
  # Suppress labels for tiny slices
  pct_df$label[pct_df$pct_genome < 0.1 & pct_df$category != "No Change"] <-
    sprintf("%s\n%.2f%%", pct_df$category[pct_df$pct_genome < 0.1 & pct_df$category != "No Change"],
            pct_df$pct_genome[pct_df$pct_genome < 0.1 & pct_df$category != "No Change"])

  ggplot(pct_df, aes(x = "", y = pct_genome, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5), size = 2.8) +
    scale_fill_manual(values = cat_colors, name = "Category", drop = FALSE) +
    labs(title = sprintf("Compartment Shifts — %s", title_suffix)) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(fill = guide_legend(nrow = 2))
}

p_pie7_std <- make_pie_7cat(pct_7cat_std, "Standard (FDR<0.05, |Diff|>0.30)")
p_pie7_rlx <- make_pie_7cat(pct_7cat_rlx, "Relaxed (FDR<0.15, |Diff|>0.15)")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_pie7_combined <- p_pie7_std + p_pie7_rlx +
    plot_annotation(
      title = "PC1 Compartment Changes: Weakened vs Strengthened vs Flipped",
      subtitle = "BAP1-KO vs Wildtype | mm10 genome (2.73 Gb)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
      )
    )
  save_multiformat_ggplot(p_pie7_combined,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 14, height = 7)
} else {
  save_multiformat_ggplot(p_pie7_std,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat_standard"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 8, height = 7)
  save_multiformat_ggplot(p_pie7_rlx,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat_relaxed"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 8, height = 7)
}

# Sig-only version (exclude No Change for readability)
make_pie_7cat_sigonly <- function(pct_df, title_suffix) {
  pct_sig <- pct_df[pct_df$category != "No Change", ]
  pct_sig$pct_of_sig <- pct_sig$bp / sum(pct_sig$bp) * 100
  pct_sig$label <- sprintf("%s\n%.1f%%\n(%d regions)",
                           pct_sig$category, pct_sig$pct_of_sig, pct_sig$n_regions)

  total_pct <- sum(pct_sig$pct_genome)

  ggplot(pct_sig, aes(x = "", y = bp, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = cat_colors, name = "Category", drop = FALSE) +
    labs(title = sprintf("Significant Compartment Shifts -- %s", title_suffix),
         subtitle = sprintf("%.1f%% of genome affected", total_pct)) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(fill = guide_legend(nrow = 1))
}

p_pie7_sig_std <- make_pie_7cat_sigonly(pct_7cat_std, "Standard (FDR<0.05, |Diff|>0.30)")
p_pie7_sig_rlx <- make_pie_7cat_sigonly(pct_7cat_rlx, "Relaxed (FDR<0.15, |Diff|>0.15)")

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_pie7_sig_combined <- p_pie7_sig_std + p_pie7_sig_rlx +
    plot_annotation(
      title = "Differential PC1 Categories (Significant Regions Only)",
      subtitle = "BAP1-KO vs Wildtype | Proportion of significant compartment changes",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
      )
    )
  save_multiformat_ggplot(p_pie7_sig_combined,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat_sigonly"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 14, height = 7)
} else {
  save_multiformat_ggplot(p_pie7_sig_std,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat_sigonly_standard"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 8, height = 7)
  save_multiformat_ggplot(p_pie7_sig_rlx,
                          file.path(PLOT_DIR, "1D_compartment_genome_pct_pie_7cat_sigonly_relaxed"),  # Original: file.path(OUTPUT_DIR, ...)
                          width = 8, height = 7)
}

# Save 7-category table
pct_7cat_export <- rbind(
  transform(pct_7cat_std, threshold = as.character(threshold)),
  transform(pct_7cat_rlx, threshold = as.character(threshold))
)
pct_7cat_export$category <- as.character(pct_7cat_export$category)
# Original: file.path(OUTPUT_DIR, "compartment_genome_pct_7cat_table.tsv")
write.table(pct_7cat_export,
            file.path(TSV_DIR, "1D_compartment_genome_pct_7cat_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(TSV_DIR, "1D_compartment_genome_pct_7cat_table.tsv")))

# =============================================================================
# FIGURE 3: PER-CHROMOSOME BREAKDOWN (standard threshold)
# =============================================================================

cat("Generating per-chromosome breakdown...\n")

# mm10 chromosome sizes (assembled)
mm10_chr_sizes <- c(
  chr1 = 195471971, chr2 = 182113224, chr3 = 160039680, chr4 = 156508116,
  chr5 = 151834684, chr6 = 149736546, chr7 = 145441459, chr8 = 129401213,
  chr9 = 124595110, chr10 = 130694993, chr11 = 122082543, chr12 = 120129022,
  chr13 = 120421639, chr14 = 124902244, chr15 = 104043685, chr16 = 98207768,
  chr17 = 94987271, chr18 = 90702639, chr19 = 61431566, chrX = 171031299
)

# Tally significant bp per chromosome per direction
chr_summary <- sig_standard %>%
  group_by(Chr, direction) %>%
  summarise(sig_bp = sum(Region_size), n_regions = n(), .groups = "drop") %>%
  mutate(
    chr_size = mm10_chr_sizes[Chr],
    pct_chr = sig_bp / chr_size * 100,
    direction_label = ifelse(direction == "A_to_B_in_Mutant",
                             "A->B (More Inactive)", "B->A (More Active)")
  )

# Order chromosomes numerically
chr_order <- c(paste0("chr", 1:19), "chrX")
chr_summary$Chr <- factor(chr_summary$Chr, levels = chr_order)

p_chr <- ggplot(chr_summary, aes(x = Chr, y = pct_chr, fill = direction_label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("A->B (More Inactive)" = "#4393C3",
                                "B->A (More Active)" = "#D6604D"),
                    name = "Compartment Shift") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Differential Compartment Shifts by Chromosome",
    subtitle = "Standard thresholds (FDR<0.05, |Diff|>0.30) | % of each chromosome affected",
    x = NULL,
    y = "% of Chromosome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

save_multiformat_ggplot(p_chr,
                        file.path(PLOT_DIR, "1D_compartment_genome_pct_by_chr"),  # Original: file.path(OUTPUT_DIR, "compartment_genome_pct_by_chr")
                        width = 10, height = 6)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("Writing summary statistics...\n")

total_analyzed_bp <- sum(all_regions$Region_size)

summary_lines <- c(
  "% of Genome with Differential PC1 — Summary Statistics",
  "=======================================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("mm10 genome size: %s bp (%.2f Gb)", format(MM10_GENOME_SIZE, big.mark = ","), MM10_GENOME_SIZE / 1e9),
  sprintf("Total regions analyzed: %d (%s bp = %.2f Gb)",
          nrow(all_regions), format(total_analyzed_bp, big.mark = ","), total_analyzed_bp / 1e9),
  ""
)

for (thresh_label in c("Standard", "Relaxed")) {
  if (thresh_label == "Standard") {
    sig_df <- sig_standard
    thresh_desc <- "FDR<0.05, |Diff|>0.30"
  } else {
    sig_df <- sig_relaxed
    thresh_desc <- "FDR<0.15, |Diff|>0.15"
  }

  a2b_bp <- sum(sig_df$Region_size[sig_df$direction == "A_to_B_in_Mutant"])
  b2a_bp <- sum(sig_df$Region_size[sig_df$direction == "B_to_A_in_Mutant"])
  total_bp <- a2b_bp + b2a_bp

  summary_lines <- c(summary_lines,
    "------------------------------------------------------------",
    sprintf("%s Thresholds (%s)", toupper(thresh_label), thresh_desc),
    "------------------------------------------------------------",
    sprintf("  Total significant regions: %d", nrow(sig_df)),
    sprintf("  Total significant bp: %s (%.1f Mb)", format(total_bp, big.mark = ","), total_bp / 1e6),
    sprintf("  %% of genome: %.2f%%", total_bp / MM10_GENOME_SIZE * 100),
    sprintf("  %% of analyzed regions: %.2f%%", nrow(sig_df) / nrow(all_regions) * 100),
    "",
    sprintf("  A->B (More Inactive in Mutant):"),
    sprintf("    Regions: %d", sum(sig_df$direction == "A_to_B_in_Mutant")),
    sprintf("    bp: %s (%.1f Mb)", format(a2b_bp, big.mark = ","), a2b_bp / 1e6),
    sprintf("    %% of genome: %.2f%%", a2b_bp / MM10_GENOME_SIZE * 100),
    "",
    sprintf("  B->A (More Active in Mutant):"),
    sprintf("    Regions: %d", sum(sig_df$direction == "B_to_A_in_Mutant")),
    sprintf("    bp: %s (%.1f Mb)", format(b2a_bp, big.mark = ","), b2a_bp / 1e6),
    sprintf("    %% of genome: %.2f%%", b2a_bp / MM10_GENOME_SIZE * 100),
    ""
  )
}

# Per-chromosome summary
summary_lines <- c(summary_lines,
  "------------------------------------------------------------",
  "Per-Chromosome Breakdown (Standard Thresholds)",
  "------------------------------------------------------------",
  sprintf("%-6s  %12s  %10s  %8s  %8s",
          "Chr", "Chr Size", "Sig bp", "% Chr", "Regions")
)

for (chr in chr_order) {
  chr_sig <- sig_standard[sig_standard$Chr == chr, ]
  if (nrow(chr_sig) == 0) next
  chr_bp <- sum(chr_sig$Region_size)
  summary_lines <- c(summary_lines,
    sprintf("%-6s  %12s  %10s  %7.2f%%  %7d",
            chr,
            format(mm10_chr_sizes[chr], big.mark = ","),
            format(chr_bp, big.mark = ","),
            chr_bp / mm10_chr_sizes[chr] * 100,
            nrow(chr_sig)))
}

# 7-category breakdown
for (pct_7 in list(pct_7cat_std, pct_7cat_rlx)) {
  thresh_name <- as.character(pct_7$threshold[1])
  summary_lines <- c(summary_lines,
    "",
    "------------------------------------------------------------",
    sprintf("7-Category Breakdown — %s", thresh_name),
    "------------------------------------------------------------",
    sprintf("  %-18s  %10s  %8s  %8s", "Category", "bp (Mb)", "% Genome", "Regions")
  )
  for (i in seq_len(nrow(pct_7))) {
    summary_lines <- c(summary_lines,
      sprintf("  %-18s  %10.1f  %7.2f%%  %7d",
              as.character(pct_7$category[i]),
              pct_7$bp[i] / 1e6,
              pct_7$pct_genome[i],
              pct_7$n_regions[i]))
  }
}

# Original: summary_file <- file.path(OUTPUT_DIR, "compartment_genome_percentage_summary.txt")
summary_file <- file.path(TSV_DIR, "1D_compartment_genome_percentage_summary.txt")
writeLines(summary_lines, summary_file)
cat(sprintf("  Saved: %s\n", summary_file))

# Save summary table as TSV
# Original: table_file <- file.path(OUTPUT_DIR, "compartment_genome_percentage_table.tsv")
table_file <- file.path(TSV_DIR, "1D_compartment_genome_percentage_table.tsv")
pct_export <- pct_combined
pct_export$threshold <- gsub("\n", " ", pct_export$threshold)
write.table(pct_export, table_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", table_file))

# =============================================================================
# DONE
# =============================================================================

cat("\n================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================\n\n")

cat(sprintf("Output directory (TSVs): %s\n", TSV_DIR))  # Original: OUTPUT_DIR
cat(sprintf("Output directory (plots): %s\n\n", PLOT_DIR))
cat("Generated files:\n")
cat("  - compartment_genome_pct_bar/{pdf,svg,jpg}          (stacked bar chart)\n")
cat("  - compartment_genome_pct_pie/{pdf,svg,jpg}          (binary pie charts)\n")
cat("  - compartment_genome_pct_pie_7cat/{pdf,svg,jpg}     (7-category pie charts)\n")
cat("  - compartment_genome_pct_by_chr/{pdf,svg,jpg}       (per-chromosome breakdown)\n")
cat("  - compartment_genome_percentage_summary.txt          (summary statistics)\n")
cat("  - compartment_genome_percentage_table.tsv            (binary summary table)\n")
cat("  - compartment_genome_pct_7cat_table.tsv              (7-category summary table)\n\n")

# Print key numbers for quick reference
sig_std_bp <- sum(sig_standard$Region_size)
sig_rlx_bp <- sum(sig_relaxed$Region_size)
cat("Key results:\n")
cat(sprintf("  Standard: %.1f Mb (%.2f%% of genome) — %d regions\n",
            sig_std_bp / 1e6, sig_std_bp / MM10_GENOME_SIZE * 100, nrow(sig_standard)))
cat(sprintf("  Relaxed:  %.1f Mb (%.2f%% of genome) — %d regions\n",
            sig_rlx_bp / 1e6, sig_rlx_bp / MM10_GENOME_SIZE * 100, nrow(sig_relaxed)))
cat("\n")
