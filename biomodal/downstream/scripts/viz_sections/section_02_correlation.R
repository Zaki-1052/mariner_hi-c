# biomodal/downstream/scripts/viz_sections/section_02_correlation.R
# Section 2: Sample Correlation Analysis
# Standalone script - sources shared config for all dependencies and data

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
    n_samples <- 4
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
      Sex = ifelse(grepl("F$", sample_names_short), "Female", "Male"),
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

    # Extract correlations using named indexing for robustness
    get_corr <- function(mat, s1, s2) {
      mat[s1, s2]
    }

    corr_summary <- data.frame(
      Type = c(rep("5mC", 6), rep("5hmC", 6)),
      Comparison = rep(c("ctrl-F vs ctrl-M", "ctrl-F vs mut-F", "ctrl-F vs mut-M",
                         "ctrl-M vs mut-F", "ctrl-M vs mut-M", "mut-F vs mut-M"), 2),
      Correlation = c(
        get_corr(mc_corr_mat, "ctrl-F", "ctrl-M"),
        get_corr(mc_corr_mat, "ctrl-F", "mut-F"),
        get_corr(mc_corr_mat, "ctrl-F", "mut-M"),
        get_corr(mc_corr_mat, "ctrl-M", "mut-F"),
        get_corr(mc_corr_mat, "ctrl-M", "mut-M"),
        get_corr(mc_corr_mat, "mut-F", "mut-M"),
        get_corr(hmc_corr_mat, "ctrl-F", "ctrl-M"),
        get_corr(hmc_corr_mat, "ctrl-F", "mut-F"),
        get_corr(hmc_corr_mat, "ctrl-F", "mut-M"),
        get_corr(hmc_corr_mat, "ctrl-M", "mut-F"),
        get_corr(hmc_corr_mat, "ctrl-M", "mut-M"),
        get_corr(hmc_corr_mat, "mut-F", "mut-M")
      ),
      Group = rep(c("Within-Control", "Between", "Between",
                    "Between", "Between", "Within-Mutant"), 2)
    )

    p_corr <- ggplot(corr_summary, aes(x = Comparison, y = Correlation, fill = Type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = c(0.5, 0.75), linetype = "dashed", color = "grey50", alpha = 0.5) +
      scale_fill_manual(values = COLORS$methylation) +
      facet_wrap(~Group, scales = "free_x") +
      labs(
        title = "Sample Correlations: 5mC vs 5hmC",
        subtitle = "5mC shows higher within-group correlation than 5hmC",
        x = "", y = "Pearson Correlation"
      ) +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    save_multiformat_ggplot(p_corr, file.path(OUTPUT_DIR, "02c_correlation_comparison"),
                            width = 14, height = 7)
  }
}

cat("Section 2 complete.\n\n")
