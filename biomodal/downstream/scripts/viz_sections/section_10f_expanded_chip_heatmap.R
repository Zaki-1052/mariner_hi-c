# biomodal/downstream/scripts/viz_sections/section_10f_expanded_chip_heatmap.R
# Section 10f: Expanded ChIP-seq / Chromatin Mark Overlap Heatmap (10 marks)
# Standalone script - sources shared config for all dependencies and data
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_10f_expanded_chip_heatmap.R

source("scripts/viz_sections/_shared_config.R")

cat("================================================================================\n")
cat("SECTION 10f: EXPANDED CHROMATIN MARK OVERLAP HEATMAP\n")
cat("================================================================================\n\n")

# =============================================================================
# Load all 10 marks as GRanges
# =============================================================================

cat("Loading ChIP-seq / chromatin mark peak files...\n\n")

cat("  Original 6 marks:\n")
chip_peaks <- list(
  CTCF     = load_chip_peaks(CHIP_PEAK_FILES$ctcf,     "CTCF"),
  H3K27ac  = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac,  "H3K27ac"),
  H3K27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  H3K4me1  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1,  "H3K4me1"),
  H3K4me3  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3,  "H3K4me3"),
  Bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent,  "Bivalent")
)

cat("\n  New marks (union of condition-specific peaks):\n")

k119ub_ctrl <- load_chip_peaks(K119UB_FILES$ctrl, "H2AK119ub ctrl")
k119ub_mut  <- load_chip_peaks(K119UB_FILES$mut,  "H2AK119ub mut")
chip_peaks$H2AK119ub <- c(k119ub_ctrl, k119ub_mut)
cat(sprintf("    H2AK119ub union: %d peaks\n", length(chip_peaks$H2AK119ub)))

atac_ctrl <- load_chip_peaks(ATAC_FILES$consensus_ctrl, "ATAC consensus ctrl")
atac_mut  <- load_chip_peaks(ATAC_FILES$consensus_mut,  "ATAC consensus mut")
chip_peaks$ATAC <- c(atac_ctrl, atac_mut)
cat(sprintf("    ATAC union: %d peaks\n", length(chip_peaks$ATAC)))

me2_up   <- load_chip_peaks(H3K36ME2_FILES$up,   "H3K36me2 up")
me2_down <- load_chip_peaks(H3K36ME2_FILES$down,  "H3K36me2 down")
chip_peaks$H3K36me2 <- c(me2_up, me2_down)
cat(sprintf("    H3K36me2 union: %d peaks\n", length(chip_peaks$H3K36me2)))

me3_up   <- load_chip_peaks(H3K36ME3_FILES$up,   "H3K36me3 up")
me3_down <- load_chip_peaks(H3K36ME3_FILES$down,  "H3K36me3 down")
chip_peaks$H3K36me3 <- c(me3_up, me3_down)
cat(sprintf("    H3K36me3 union: %d peaks\n", length(chip_peaks$H3K36me3)))

failed <- names(chip_peaks)[sapply(chip_peaks, is.null)]
if (length(failed) > 0) {
  stop(sprintf("Peak files missing for: %s", paste(failed, collapse = ", ")))
}

# =============================================================================
# Compute overlaps for all 10 marks against significant mC DMRs
# =============================================================================

mc_sig <- mc_dmr %>% dplyr::filter(significant)
mc_gr  <- dmr_to_granges(mc_sig)
cat(sprintf("\n%d significant mC DMRs\n", length(mc_gr)))

MARK_ORDER <- names(chip_peaks)

cat("\nComputing overlaps:\n")
for (mark in MARK_ORDER) {
  col <- paste0(mark, "_overlap")
  mc_sig[[col]] <- countOverlaps(mc_gr, chip_peaks[[mark]]) > 0
  cat(sprintf("  %-12s: %5d DMRs (%.1f%%)\n",
              mark, sum(mc_sig[[col]]), 100 * mean(mc_sig[[col]])))
}

# =============================================================================
# Build heatmap data (rows = direction, columns = marks)
# =============================================================================

overlap_cols <- paste0(MARK_ORDER, "_overlap")

mark_by_direction <- mc_sig %>%
  group_by(direction) %>%
  summarise(across(all_of(overlap_cols), ~ 100 * mean(.x)), .groups = "drop") %>%
  pivot_longer(cols = all_of(overlap_cols),
               names_to = "Mark", values_to = "Percentage") %>%
  mutate(Mark = sub("_overlap$", "", Mark))

mark_overall <- mc_sig %>%
  summarise(across(all_of(overlap_cols), ~ 100 * mean(.x))) %>%
  mutate(direction = "All Significant") %>%
  pivot_longer(cols = all_of(overlap_cols),
               names_to = "Mark", values_to = "Percentage") %>%
  mutate(Mark = sub("_overlap$", "", Mark))

heatmap_df <- bind_rows(mark_overall, mark_by_direction)
heatmap_df$direction <- factor(heatmap_df$direction,
                               levels = c("All Significant", "Hypermethylated", "Hypomethylated"))
heatmap_df$Mark <- factor(heatmap_df$Mark, levels = MARK_ORDER)

# =============================================================================
# Draw heatmap
# =============================================================================

cat("\nCreating expanded chromatin mark overlap heatmap...\n")

n_hyper <- sum(mc_sig$direction == "Hypermethylated")
n_hypo  <- sum(mc_sig$direction == "Hypomethylated")

p_10f <- ggplot(heatmap_df, aes(x = Mark, y = direction, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 3.5) +
  scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                       midpoint = 50, name = "% Overlap", limits = c(0, 100)) +
  labs(
    title = "Chromatin Mark Overlap at Significant mC DMRs",
    subtitle = sprintf("10 marks | %d significant DMRs (hyper=%d, hypo=%d)",
                       nrow(mc_sig), n_hyper, n_hypo),
    x = "Chromatin Mark", y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x   = element_text(face = "bold", size = 10, angle = 30, hjust = 1),
    axis.text.y   = element_text(size = 11),
    panel.grid    = element_blank()
  )

save_multiformat_ggplot(p_10f, file.path(OUTPUT_DIR, "10f_expanded_chip_heatmap"),
                        width = 14, height = 5)

# =============================================================================
# Save overlap table
# =============================================================================

write.table(heatmap_df, file.path(TABLES_DIR, "10f_expanded_chip_overlap_percentages.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 10f_expanded_chip_overlap_percentages.tsv\n")

cat("\nSection 10f complete.\n\n")
