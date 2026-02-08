# biomodal/downstream/scripts/viz_sections/section_01_qc_overview.R
# Section 1: QC & Data Overview
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 1: QC & DATA OVERVIEW
# =============================================================================

cat("================================================================================\n")
cat("SECTION 1: QC & DATA OVERVIEW\n")
cat("================================================================================\n\n")

# Extract QC metrics from bioqc JSON
if (!is.null(bioqc)) {
  # Parse QC metrics from the JSON structure
  qc_df <- data.frame(
    sample = SAMPLES$short_name,
    condition = SAMPLES$condition,
    stringsAsFactors = FALSE
  )

  # Create QC overview plot
  cat("Creating QC overview plot...\n")

  # If upstream data available, create sequencing metrics plot
  if (!is.null(upstream)) {
    # Extract key metrics (using actual column names from CSV)
    seq_metrics <- upstream %>%
      dplyr::select(sample_id, bamlet_mapped_reads_genome, bamlet_mapped_bases_cigar_genome,
                    bamlet_prop_duplicated_reads_genome, bamlet_mean_phred_genome) %>%
      dplyr::mutate(
        Mapped_Reads_M = bamlet_mapped_reads_genome / 1e6,
        Mapped_Bases_B = bamlet_mapped_bases_cigar_genome / 1e9,
        Duplication_Rate = bamlet_prop_duplicated_reads_genome,
        Mean_Phred_Score = bamlet_mean_phred_genome,
        sample_short = gsub("evoC-Bap1-", "", sample_id),
        condition = ifelse(grepl("ctrl", sample_id), "Control", "Mutant")
      )

    # Mapped reads bar chart - use points+bars to show actual values clearly
    p1 <- ggplot(seq_metrics, aes(x = sample_short, y = Mapped_Reads_M, fill = condition)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%.0f", Mapped_Reads_M)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(0, max(seq_metrics$Mapped_Reads_M) * 1.15)) +
      labs(title = "Mapped Reads", x = "", y = "Reads (millions)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Mapped bases bar chart
    p2 <- ggplot(seq_metrics, aes(x = sample_short, y = Mapped_Bases_B, fill = condition)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%.1f", Mapped_Bases_B)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(0, max(seq_metrics$Mapped_Bases_B) * 1.15)) +
      labs(title = "Mapped Bases", x = "", y = "Bases (billions)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Duplication rate - use lollipop plot with zoomed y-axis to show differences
    dup_rate_pct <- seq_metrics$Duplication_Rate * 100
    p3 <- ggplot(seq_metrics, aes(x = sample_short, y = Duplication_Rate * 100,
                                   color = condition)) +
      geom_segment(aes(xend = sample_short, y = min(dup_rate_pct) - 1,
                       yend = Duplication_Rate * 100), linewidth = 2) +
      geom_point(size = 5) +
      geom_text(aes(label = sprintf("%.1f%%", Duplication_Rate * 100)),
                vjust = -1, size = 3, color = "black") +
      scale_color_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(min(dup_rate_pct) - 2, max(dup_rate_pct) + 2)) +
      labs(title = "Duplication Rate",
           subtitle = sprintf("Range: %.1f%% - %.1f%%", min(dup_rate_pct), max(dup_rate_pct)),
           x = "", y = "Rate (%)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Mean Phred score - use lollipop plot with zoomed y-axis
    p4 <- ggplot(seq_metrics, aes(x = sample_short, y = Mean_Phred_Score,
                                   color = condition)) +
      geom_segment(aes(xend = sample_short,
                       y = min(seq_metrics$Mean_Phred_Score) - 0.5,
                       yend = Mean_Phred_Score), linewidth = 2) +
      geom_point(size = 5) +
      geom_text(aes(label = sprintf("%.2f", Mean_Phred_Score)),
                vjust = -1, size = 3, color = "black") +
      scale_color_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(min(seq_metrics$Mean_Phred_Score) - 1,
                                    max(seq_metrics$Mean_Phred_Score) + 0.5)) +
      labs(title = "Mean Phred Score",
           subtitle = sprintf("All samples: ~%.1f (excellent)",
                             mean(seq_metrics$Mean_Phred_Score)),
           x = "", y = "Phred Score") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    qc_combined <- (p1 | p2) / (p3 | p4) +
      plot_annotation(
        title = "Sequencing Quality Metrics - Biomodal DUET evoC",
        subtitle = "BAP1-KO vs Wildtype Control",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                      plot.subtitle = element_text(hjust = 0.5, size = 12))
      )

    save_multiformat_ggplot(qc_combined, file.path(OUTPUT_DIR, "01_qc_overview"),
                            width = 12, height = 10)
  } else {
    cat("  Skipping QC plot (upstream data not available)\n")
  }
}

cat("Section 1 complete.\n\n")
