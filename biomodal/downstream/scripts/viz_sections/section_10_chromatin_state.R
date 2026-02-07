# biomodal/downstream/scripts/viz_sections/section_10_chromatin_state.R
# Section 10: ChIP-seq Chromatin State Analysis
# Standalone script - sources shared config for all dependencies and data
# NOTE: This section depends on the 'coordinated' data frame from Section 5.
# The coordinated analysis is re-computed here so this script is self-contained.

# Source shared config (handles both Rscript and source())
local({
  script_dir <- NULL
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", args, value = TRUE)
  if (length(f) > 0) script_dir <- dirname(normalizePath(sub("--file=", "", f)))
  if (is.null(script_dir)) for (i in seq_len(sys.nframe())) {
    fi <- sys.frame(i)$ofile
    if (!is.null(fi)) { script_dir <- dirname(normalizePath(fi)); break }
  }
  if (is.null(script_dir)) script_dir <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections"
  source(file.path(script_dir, "_shared_config.R"))
})

# =============================================================================
# Re-compute coordinated changes (from Section 5) for chromatin analysis
# =============================================================================

mc_sig_coord <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_ctrl = mean_mod_group1, mc_mut = mean_mod_group2)

hmc_sig_coord <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2)

coordinated <- inner_join(mc_sig_coord, hmc_sig_coord, by = "gene")

coordinated <- coordinated %>%
  mutate(
    coordinated_pattern = (mc_diff > 0 & hmc_diff < 0),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  ) %>%
  arrange(desc(combined_effect))

# =============================================================================
# SECTION 10: ChIP-seq CHROMATIN STATE ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 10: ChIP-seq CHROMATIN STATE ANALYSIS\n")
cat("================================================================================\n\n")

cat("Loading ChIP-seq peak files...\n")

# Load all ChIP-seq peak files
chip_peaks <- list(
  ctcf = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  h3k27ac = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  h3k4me1 = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  h3k4me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)

# Check if all peaks loaded successfully
if (any(sapply(chip_peaks, is.null))) {
  cat("  WARNING: Some ChIP-seq peak files not found. Skipping Section 10.\n")
} else {
  cat("\nCreating GRanges from mC DMR data...\n")

  # Create GRanges from mC DMR data (significant genes only for analysis)
  mc_sig <- mc_dmr %>% dplyr::filter(significant)
  mc_gr <- dmr_to_granges(mc_sig)
  cat(sprintf("  %d significant mC DMRs\n", length(mc_gr)))

  # Compute TSS distances
  cat("\nComputing TSS distances...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss_gr <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d genes\n", length(genes)))

  nearest <- distanceToNearest(mc_gr, tss_gr)
  distance_to_tss <- rep(NA_real_, length(mc_gr))
  if (length(nearest) > 0) {
    distance_to_tss[queryHits(nearest)] <- mcols(nearest)$distance
  }

  cat(sprintf("  DMRs <= 2kb from TSS: %d (%.1f%%)\n",
              sum(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE),
              100 * mean(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE)))

  # Compute ChIP-seq overlaps
  cat("\nComputing ChIP-seq overlaps...\n")
  mc_overlaps <- compute_chip_overlaps(mc_gr, chip_peaks)

  cat("  Overlap summary:\n")
  cat(sprintf("    CTCF:     %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$CTCF_overlap), 100 * mean(mc_overlaps$CTCF_overlap)))
  cat(sprintf("    H3K27ac:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K27ac_overlap), 100 * mean(mc_overlaps$H3K27ac_overlap)))
  cat(sprintf("    H3K27me3: %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K27me3_overlap), 100 * mean(mc_overlaps$H3K27me3_overlap)))
  cat(sprintf("    H3K4me1:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K4me1_overlap), 100 * mean(mc_overlaps$H3K4me1_overlap)))
  cat(sprintf("    H3K4me3:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K4me3_overlap), 100 * mean(mc_overlaps$H3K4me3_overlap)))
  cat(sprintf("    Bivalent: %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$Bivalent_overlap), 100 * mean(mc_overlaps$Bivalent_overlap)))

  # Classify chromatin states
  cat("\nClassifying chromatin states (7 categories)...\n")
  mc_sig$chromatin_state <- classify_chromatin_state(mc_overlaps, distance_to_tss, TSS_THRESHOLD)

  # Add overlap columns to dataframe
  mc_sig$CTCF_overlap <- mc_overlaps$CTCF_overlap
  mc_sig$H3K27ac_overlap <- mc_overlaps$H3K27ac_overlap
  mc_sig$H3K27me3_overlap <- mc_overlaps$H3K27me3_overlap
  mc_sig$H3K4me1_overlap <- mc_overlaps$H3K4me1_overlap
  mc_sig$H3K4me3_overlap <- mc_overlaps$H3K4me3_overlap
  mc_sig$Bivalent_overlap <- mc_overlaps$Bivalent_overlap
  mc_sig$distance_to_tss <- distance_to_tss

  cat("\n  Chromatin state distribution (significant mC DMRs):\n")
  for (state in CHROMATIN_STATE_ORDER) {
    count <- sum(mc_sig$chromatin_state == state)
    pct <- 100 * count / nrow(mc_sig)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", state, count, pct))
  }

  # -----------------------------------------------------------------------
  # FIGURE 10a: Overall chromatin state distribution
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 10a: Chromatin state distribution...\n")

  state_summary <- mc_sig %>%
    group_by(chromatin_state, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction) %>%
    mutate(percentage = 100 * count / sum(count))

  # Bar chart
  p_10a_bar <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.3, size = 3) +
    facet_wrap(~direction, ncol = 2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.15), expand = c(0, 0)) +
    labs(
      title = "Chromatin State Distribution of Significant mC DMRs",
      subtitle = sprintf("Hypermethylated: n=%d | Hypomethylated: n=%d",
                         sum(mc_sig$direction == "Hypermethylated"),
                         sum(mc_sig$direction == "Hypomethylated")),
      x = "Chromatin State", y = "Percentage (%)"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Pie chart - overall
  state_overall <- mc_sig %>%
    group_by(chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  p_10a_pie <- ggplot(state_overall, aes(x = "", y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    labs(title = "Overall Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    ) +
    geom_text(aes(label = ifelse(percentage > 4, sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.5)

  p_10a <- p_10a_bar + p_10a_pie +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "Chromatin State Analysis of Differentially Methylated Genes",
      subtitle = "Based on ChIP-seq peak overlaps (Late timepoint)",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_10a, file.path(OUTPUT_DIR, "10a_chromatin_state_distribution"),
                          width = 14, height = 7)

  # -----------------------------------------------------------------------
  # FIGURE 10b: Chromatin state by methylation direction comparison
  # -----------------------------------------------------------------------
  cat("Creating Figure 10b: Chromatin state by direction...\n")

  # Stacked bar comparison
  p_10b_stacked <- ggplot(state_summary, aes(x = direction, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Chromatin State by Methylation Direction",
      x = "Direction", y = "Percentage (%)"
    ) +
    theme_biomodal()

  # Grouped bar comparison
  p_10b_grouped <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = COLORS$direction, name = "Direction") +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(
      title = "Comparison by Chromatin State",
      subtitle = "Do hypermethylated genes have specific chromatin signatures?",
      x = "Chromatin State", y = "Percentage (%)"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  p_10b <- p_10b_stacked / p_10b_grouped +
    plot_annotation(
      title = "Methylation Direction vs Chromatin State",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

  save_multiformat_ggplot(p_10b, file.path(OUTPUT_DIR, "10b_chromatin_by_methylation_direction"),
                          width = 12, height = 10)

  # -----------------------------------------------------------------------
  # FIGURE 10c: ChIP mark overlap heatmap
  # -----------------------------------------------------------------------
  cat("Creating Figure 10c: ChIP mark overlap heatmap...\n")

  # Calculate overlap percentages by direction
  mark_overlap <- mc_sig %>%
    group_by(direction) %>%
    summarise(
      CTCF = 100 * mean(CTCF_overlap),
      H3K27ac = 100 * mean(H3K27ac_overlap),
      H3K27me3 = 100 * mean(H3K27me3_overlap),
      H3K4me1 = 100 * mean(H3K4me1_overlap),
      H3K4me3 = 100 * mean(H3K4me3_overlap),
      Bivalent = 100 * mean(Bivalent_overlap),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -direction, names_to = "Mark", values_to = "Percentage")

  # Add overall row
  mark_overall <- mc_sig %>%
    summarise(
      CTCF = 100 * mean(CTCF_overlap),
      H3K27ac = 100 * mean(H3K27ac_overlap),
      H3K27me3 = 100 * mean(H3K27me3_overlap),
      H3K4me1 = 100 * mean(H3K4me1_overlap),
      H3K4me3 = 100 * mean(H3K4me3_overlap),
      Bivalent = 100 * mean(Bivalent_overlap)
    ) %>%
    mutate(direction = "All Significant") %>%
    pivot_longer(cols = -direction, names_to = "Mark", values_to = "Percentage")

  mark_overlap <- bind_rows(mark_overall, mark_overlap)
  mark_overlap$direction <- factor(mark_overlap$direction,
                                   levels = c("All Significant", "Hypermethylated", "Hypomethylated"))

  p_10c <- ggplot(mark_overlap, aes(x = Mark, y = direction, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 4) +
    scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                         midpoint = 50, name = "% Overlap", limits = c(0, 100)) +
    labs(
      title = "ChIP-seq Mark Overlap at Significant mC DMRs",
      subtitle = "% of DMRs overlapping each ChIP-seq mark",
      x = "ChIP-seq Mark", y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank()
    )

  save_multiformat_ggplot(p_10c, file.path(OUTPUT_DIR, "10c_chip_mark_overlap_heatmap"),
                          width = 10, height = 5)

  # -----------------------------------------------------------------------
  # FIGURE 10d: Chromatin state of coordinated genes
  # -----------------------------------------------------------------------
  cat("Creating Figure 10d: Coordinated genes chromatin state...\n")

  # Get chromatin states for coordinated genes (mC up / hmC down)
  coordinated_mc <- mc_sig %>%
    dplyr::filter(gene %in% coordinated$gene[coordinated$coordinated_pattern])

  if (nrow(coordinated_mc) > 0) {
    coord_state_summary <- coordinated_mc %>%
      group_by(chromatin_state) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(percentage = 100 * count / sum(count))

    p_10d_bar <- ggplot(coord_state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, count)), vjust = -0.1, size = 3) +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
      scale_y_continuous(limits = c(0, max(coord_state_summary$percentage) * 1.2), expand = c(0, 0)) +
      labs(
        title = "Chromatin State of Coordinated mC\u2191/hmC\u2193 Genes",
        subtitle = sprintf("n = %d genes with coordinated methylation changes", nrow(coordinated_mc)),
        x = "Chromatin State", y = "Percentage (%)"
      ) +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    p_10d_pie <- ggplot(coord_state_summary, aes(x = "", y = percentage, fill = chromatin_state)) +
      geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
      labs(title = "Distribution") +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "right"
      ) +
      geom_text(aes(label = ifelse(percentage > 5, sprintf("%.0f%%", percentage), "")),
                position = position_stack(vjust = 0.5), size = 3)

    p_10d <- p_10d_bar + p_10d_pie +
      plot_layout(widths = c(2, 1)) +
      plot_annotation(
        title = "Chromatin State of Genes with Coordinated Methylation Changes",
        subtitle = "Genes showing mC increase + hmC decrease pattern",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
        )
      )

    save_multiformat_ggplot(p_10d, file.path(OUTPUT_DIR, "10d_coordinated_genes_chromatin"),
                            width = 13, height = 6)
  }

  # -----------------------------------------------------------------------
  # FIGURE 10e: Top genes with chromatin annotations
  # -----------------------------------------------------------------------
  cat("Creating Figure 10e: Top genes with chromatin annotations...\n")

  # Get top 20 coordinated genes with chromatin state
  top_coord_chrom <- coordinated %>%
    dplyr::filter(coordinated_pattern) %>%
    head(20) %>%
    left_join(mc_sig %>% dplyr::select(gene, chromatin_state), by = "gene")

  if (nrow(top_coord_chrom) > 0) {
    top_coord_chrom <- top_coord_chrom %>%
      mutate(gene = factor(gene, levels = rev(unique(gene))))

    p_10e <- ggplot(top_coord_chrom, aes(x = gene, y = combined_effect * 100, fill = chromatin_state)) +
      geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
      geom_text(aes(label = as.character(chromatin_state)),
                hjust = -0.1, size = 3, color = "grey30") +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
      coord_flip() +
      labs(
        title = "Top 20 Coordinated Genes with Chromatin State Annotation",
        subtitle = "Sorted by combined effect (|mC change| + |hmC change|)",
        x = "", y = "Combined Effect (%)"
      ) +
      theme_biomodal() +
      theme(legend.position = "bottom")

    save_multiformat_ggplot(p_10e, file.path(OUTPUT_DIR, "10e_top_genes_chromatin_annotation"),
                            width = 12, height = 10)
  }

  # -----------------------------------------------------------------------
  # Save chromatin state annotations
  # -----------------------------------------------------------------------
  cat("\nSaving chromatin state annotation tables...\n")

  # Full annotation table
  write.table(mc_sig, file.path(TABLES_DIR, "dmr_chromatin_state_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: dmr_chromatin_state_annotation.tsv\n")

  # Summary table
  chromatin_summary <- mc_sig %>%
    group_by(chromatin_state, direction) %>%
    summarise(
      count = n(),
      mean_mc_change = mean(mod_difference) * 100,
      median_mc_change = median(mod_difference) * 100,
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = direction,
                values_from = c(count, mean_mc_change, median_mc_change))

  write.table(chromatin_summary, file.path(TABLES_DIR, "chromatin_state_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: chromatin_state_summary.tsv\n")
}

cat("\nSection 10 complete.\n\n")
