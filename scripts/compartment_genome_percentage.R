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
#   - Pie chart: compartment shift proportions per threshold
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
BASE_DIR <- "tads/tad-pc-analysis/output/compartment_analysis"
INPUT_FILES <- list(
  standard = file.path(BASE_DIR, "compartment_significant_standard.tsv"),
  relaxed  = file.path(BASE_DIR, "compartment_significant_relaxed.tsv"),
  all      = file.path(BASE_DIR, "compartment_all_annotated.tsv")
)

OUTPUT_DIR <- BASE_DIR

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

source("scripts/utils/multi_format_output.R")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading compartment data...\n")

sig_standard <- read.delim(INPUT_FILES$standard, stringsAsFactors = FALSE)
sig_relaxed  <- read.delim(INPUT_FILES$relaxed, stringsAsFactors = FALSE)
all_regions  <- read.delim(INPUT_FILES$all, stringsAsFactors = FALSE)

# Validate expected columns
required_cols <- c("Region_size", "direction", "Chr")
stopifnot("Missing columns in standard TSV" = all(required_cols %in% names(sig_standard)))
stopifnot("Missing columns in relaxed TSV"  = all(required_cols %in% names(sig_relaxed)))
stopifnot("Missing columns in all TSV"      = all(c("Region_size", "Chr") %in% names(all_regions)))

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
                        file.path(OUTPUT_DIR, "compartment_genome_pct_bar"),
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
                          file.path(OUTPUT_DIR, "compartment_genome_pct_pie"),
                          width = 12, height = 6)
} else {
  save_multiformat_ggplot(p_pie_std,
                          file.path(OUTPUT_DIR, "compartment_genome_pct_pie_standard"),
                          width = 7, height = 6)
  save_multiformat_ggplot(p_pie_rlx,
                          file.path(OUTPUT_DIR, "compartment_genome_pct_pie_relaxed"),
                          width = 7, height = 6)
}

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
                        file.path(OUTPUT_DIR, "compartment_genome_pct_by_chr"),
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

summary_file <- file.path(OUTPUT_DIR, "compartment_genome_percentage_summary.txt")
writeLines(summary_lines, summary_file)
cat(sprintf("  Saved: %s\n", summary_file))

# Save summary table as TSV
table_file <- file.path(OUTPUT_DIR, "compartment_genome_percentage_table.tsv")
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

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))
cat("Generated files:\n")
cat("  - compartment_genome_pct_bar/{pdf,svg,jpg}     (stacked bar chart)\n")
cat("  - compartment_genome_pct_pie/{pdf,svg,jpg}     (pie charts)\n")
cat("  - compartment_genome_pct_by_chr/{pdf,svg,jpg}  (per-chromosome breakdown)\n")
cat("  - compartment_genome_percentage_summary.txt     (summary statistics)\n")
cat("  - compartment_genome_percentage_table.tsv       (summary table)\n\n")

# Print key numbers for quick reference
sig_std_bp <- sum(sig_standard$Region_size)
sig_rlx_bp <- sum(sig_relaxed$Region_size)
cat("Key results:\n")
cat(sprintf("  Standard: %.1f Mb (%.2f%% of genome) — %d regions\n",
            sig_std_bp / 1e6, sig_std_bp / MM10_GENOME_SIZE * 100, nrow(sig_standard)))
cat(sprintf("  Relaxed:  %.1f Mb (%.2f%% of genome) — %d regions\n",
            sig_rlx_bp / 1e6, sig_rlx_bp / MM10_GENOME_SIZE * 100, nrow(sig_relaxed)))
cat("\n")
