# scripts/compartment_genome_pct_pie_fdr_only.R
# 7-Category Sig-Only Pie Chart — FDR < 0.05 (no |Diff| threshold)
# Author: Zakir Alibhai
# Date: 2026-03-20
#
# Purpose:
#   Generate a sig-only 7-category pie chart for compartment shifts using
#   FDR < 0.05 only (no minimum difference cutoff). Single panel, no relaxed.
#
# Output:
#   - compartment_genome_pct_pie_7cat_sigonly_fdr05/{pdf,svg,jpg}
#
# Usage:
#   Rscript scripts/compartment_genome_pct_pie_fdr_only.R

cat("\n")
cat("================================================\n")
cat("7-Cat Sig-Only Pie — FDR<0.05 (no Diff cutoff)\n")
cat("================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

MM10_GENOME_SIZE <- 2730871774

BASE_DIR <- "tads/tad-pc-analysis/output/compartment_analysis"
INPUT_FILE <- file.path(BASE_DIR, "compartment_all_annotated.tsv")
OUTPUT_DIR <- BASE_DIR

FDR_THRESHOLD <- 0.05

if (!file.exists(INPUT_FILE)) {
  stop(sprintf("Input file not found: %s", INPUT_FILE))
}

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("scripts/utils/multi_format_output.R")

# =============================================================================
# LOAD AND FILTER DATA
# =============================================================================

cat("Loading compartment data...\n")
all_regions <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)

required_cols <- c("Region_size", "Chr", "ctrl_avg_PC1", "mut_avg_PC1", "adj_pvalue")
stopifnot("Missing columns in input TSV" = all(required_cols %in% names(all_regions)))

sig_regions <- all_regions[all_regions$adj_pvalue < FDR_THRESHOLD, ]

cat(sprintf("  All regions: %d\n", nrow(all_regions)))
cat(sprintf("  Significant (FDR < %.2f): %d\n", FDR_THRESHOLD, nrow(sig_regions)))

# =============================================================================
# CLASSIFY INTO 7 CATEGORIES
# =============================================================================

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

sig_regions$compartment_category <- classify_compartment(
  sig_regions$ctrl_avg_PC1, sig_regions$mut_avg_PC1
)

cat_levels <- c("Weakened A", "Strengthened A", "Flipped A->B",
                "Weakened B", "Strengthened B", "Flipped B->A")

cat_summary <- sig_regions %>%
  group_by(compartment_category) %>%
  summarise(bp = sum(Region_size), n_regions = n(), .groups = "drop")

pct_7cat <- data.frame(
  category = cat_summary$compartment_category,
  bp = cat_summary$bp,
  n_regions = cat_summary$n_regions,
  stringsAsFactors = FALSE
)
pct_7cat$pct_genome <- pct_7cat$bp / MM10_GENOME_SIZE * 100
pct_7cat$category <- factor(pct_7cat$category, levels = cat_levels)

# Print results
cat(sprintf("\n  7-Category — FDR < %.2f (no Diff cutoff):\n", FDR_THRESHOLD))
for (i in seq_len(nrow(pct_7cat))) {
  cat(sprintf("    %-18s: %8.1f Mb (%6.2f%%, %5d regions)\n",
              as.character(pct_7cat$category[i]),
              pct_7cat$bp[i] / 1e6,
              pct_7cat$pct_genome[i],
              pct_7cat$n_regions[i]))
}
total_pct <- sum(pct_7cat$pct_genome)
cat(sprintf("    TOTAL:              %8.1f Mb (%6.2f%%, %5d regions)\n",
            sum(pct_7cat$bp) / 1e6, total_pct, sum(pct_7cat$n_regions)))

# =============================================================================
# GENERATE SIG-ONLY PIE CHART
# =============================================================================

cat("\nGenerating sig-only 7-category pie chart...\n")

cat_colors <- c(
  "Weakened A"      = "#FCBBA1",
  "Strengthened A"  = "#A50F15",
  "Flipped A->B"    = "#2166AC",
  "Weakened B"      = "#9ECAE1",
  "Strengthened B"  = "#08306B",
  "Flipped B->A"    = "#B2182B"
)

pct_sig <- pct_7cat[!is.na(pct_7cat$category), ]
pct_sig$pct_of_sig <- pct_sig$bp / sum(pct_sig$bp) * 100
pct_sig$label <- sprintf("%s\n%.1f%%\n(%d regions)",
                         pct_sig$category, pct_sig$pct_of_sig, pct_sig$n_regions)

p_pie <- ggplot(pct_sig, aes(x = "", y = bp, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = cat_colors, name = "Category", drop = FALSE) +
  labs(title = "Significant Compartment Shifts (FDR < 0.05)",
       subtitle = sprintf("BAP1-KO vs Wildtype | %.1f%% of genome affected (%d regions)",
                          total_pct, sum(pct_sig$n_regions))) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1))

save_multiformat_ggplot(p_pie,
                        file.path(OUTPUT_DIR, "compartment_genome_pct_pie_7cat_sigonly_fdr05"),
                        width = 8, height = 7)

cat("\n================================================\n")
cat("DONE\n")
cat("================================================\n\n")
