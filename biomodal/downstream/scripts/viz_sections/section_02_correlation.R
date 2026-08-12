# biomodal/downstream/scripts/viz_sections/section_02_correlation.R
# Section 2: Sample Correlation Analysis
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 2: SAMPLE CORRELATION ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 2: SAMPLE CORRELATION ANALYSIS\n")
cat("================================================================================\n\n")

# Extract correlation matrices from BioQC JSON
if (!is.null(bioqc)) {
  cat("Extracting correlation matrices...\n")

  # Parse correlation data from JSON (it's a data frame with nested data)
  corr_data <- bioqc$correlation_data

  # Find indices for 5mC and 5hmC
  mc_idx <- which(corr_data$methylation_type == "5mC")
  hmc_idx <- which(corr_data$methylation_type == "5hmC")

  if (length(mc_idx) > 0 && length(hmc_idx) > 0) {
    # Extract correlation matrices (stored as matrices with NA for upper triangle)
    mc_corr_mat <- corr_data$correlation_matrix$data[[mc_idx]]
    hmc_corr_mat <- corr_data$correlation_matrix$data[[hmc_idx]]

    # Get sample names from JSON
    sample_names <- corr_data$correlation_matrix$index[[mc_idx]]
    sample_names_short <- gsub("evoC-Bap1-", "", sample_names)

    # Fill upper triangle by symmetry (matrices have NA in upper triangle)
    n_samples <- length(sample_names)
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        if (is.na(mc_corr_mat[i, j]) && !is.na(mc_corr_mat[j, i])) {
          mc_corr_mat[i, j] <- mc_corr_mat[j, i]
        }
        if (is.na(hmc_corr_mat[i, j]) && !is.na(hmc_corr_mat[j, i])) {
          hmc_corr_mat[i, j] <- hmc_corr_mat[j, i]
        }
      }
    }

    rownames(mc_corr_mat) <- colnames(mc_corr_mat) <- sample_names_short
    rownames(hmc_corr_mat) <- colnames(hmc_corr_mat) <- sample_names_short

    # Create annotation for heatmap (match order from JSON)
    annotation_df <- data.frame(
      Condition = ifelse(grepl("ctrl", sample_names_short), "Control", "Mutant"),
      Sex = ifelse(grepl("-F", sample_names_short), "Female", "Male"),
      row.names = sample_names_short
    )

    annotation_colors <- list(
      Condition = COLORS$condition,
      Sex = COLORS$sex
    )

    cat("Creating 5mC correlation heatmap...\n")
    save_multiformat_pheatmap(
      quote(pheatmap(
        mc_corr_mat,
        main = "5mC Sample Correlation (Gene Bodies)",
        display_numbers = TRUE,
        number_format = "%.2f",
        color = colorRampPalette(brewer.pal(9, "Blues"))(100),
        breaks = seq(0.5, 1, length.out = 101),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        annotation_row = annotation_df,
        annotation_colors = annotation_colors,
        fontsize = 12,
        fontsize_number = 10
      )),
      file.path(OUTPUT_DIR, "02a_mc_correlation_heatmap"),
      width = 8, height = 7
    )

    cat("Creating 5hmC correlation heatmap...\n")
    save_multiformat_pheatmap(
      quote(pheatmap(
        hmc_corr_mat,
        main = "5hmC Sample Correlation (Gene Bodies)",
        display_numbers = TRUE,
        number_format = "%.2f",
        color = colorRampPalette(brewer.pal(9, "Greens"))(100),
        breaks = seq(0.3, 1, length.out = 101),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        annotation_row = annotation_df,
        annotation_colors = annotation_colors,
        fontsize = 12,
        fontsize_number = 10
      )),
      file.path(OUTPUT_DIR, "02b_hmc_correlation_heatmap"),
      width = 8, height = 7
    )

    # Combined correlation comparison
    cat("Creating combined correlation comparison...\n")

    # Build all pairwise comparisons dynamically
    ctrl_samples <- sample_names_short[grepl("ctrl", sample_names_short)]
    mut_samples <- sample_names_short[grepl("mut", sample_names_short)]

    corr_rows <- list()
    for (mod_type in c("5mC", "5hmC")) {
      mat <- if (mod_type == "5mC") mc_corr_mat else hmc_corr_mat
      # Within-Control
      if (length(ctrl_samples) >= 2) {
        for (pair in combn(ctrl_samples, 2, simplify = FALSE)) {
          corr_rows[[length(corr_rows) + 1]] <- data.frame(
            Type = mod_type, Comparison = paste(pair[1], "vs", pair[2]),
            Correlation = mat[pair[1], pair[2]], Group = "Within-Control",
            stringsAsFactors = FALSE)
        }
      }
      # Within-Mutant
      if (length(mut_samples) >= 2) {
        for (pair in combn(mut_samples, 2, simplify = FALSE)) {
          corr_rows[[length(corr_rows) + 1]] <- data.frame(
            Type = mod_type, Comparison = paste(pair[1], "vs", pair[2]),
            Correlation = mat[pair[1], pair[2]], Group = "Within-Mutant",
            stringsAsFactors = FALSE)
        }
      }
      # Between
      for (c_s in ctrl_samples) {
        for (m_s in mut_samples) {
          corr_rows[[length(corr_rows) + 1]] <- data.frame(
            Type = mod_type, Comparison = paste(c_s, "vs", m_s),
            Correlation = mat[c_s, m_s], Group = "Between",
            stringsAsFactors = FALSE)
        }
      }
    }
    corr_summary <- do.call(rbind, corr_rows)

    p_corr <- ggplot(corr_summary, aes(x = Comparison, y = Correlation, fill = Type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = c(0.5, 0.75), linetype = "dashed", color = "grey50", alpha = 0.5) +
      scale_fill_manual(values = COLORS$methylation) +
      facet_wrap(~Group, scales = "free_x") +
      labs(
        title = "Sample Correlations: 5mC vs 5hmC",
        subtitle = sprintf("%d samples — within-group vs between-group correlations", n_samples),
        x = "", y = "Pearson Correlation"
      ) +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    save_multiformat_ggplot(p_corr, file.path(OUTPUT_DIR, "02c_correlation_comparison"),
                            width = 16, height = 7)
  }
}

cat("Section 2 complete.\n\n")
