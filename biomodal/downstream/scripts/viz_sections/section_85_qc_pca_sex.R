# scripts/viz_sections/section_85_qc_pca_sex.R
# Section 85: QC PCA with Sex Stratification
#
# PCA of per-sample gene-body methylation fractions for 5mC and 5hmC,
# showing PC1=condition and PC2=sex. Demonstrates that the sex covariate
# is properly handled and samples cluster by condition, not sex.
#
# Panels:
#   85a: 5mC PCA (colored by condition, shaped by sex)
#   85b: 5hmC PCA (same layout)
#   85_composite: side-by-side 5mC | 5hmC
#
# Data: per-sample regional fractions from EXTRACT_PATHS
# Output: plots/visualizations/85_qc_pca_sex/
# =============================================================================

source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

SEC85_DIR <- file.path(OUTPUT_DIR, "85_qc_pca_sex")
dir.create(SEC85_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC85_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

cat("\n=== Section 85: QC PCA with Sex Stratification ===\n\n")

mc_frac_file  <- EXTRACT_PATHS$mc_regional_frac
hmc_frac_file <- EXTRACT_PATHS$hmc_regional_frac
if (!file.exists(mc_frac_file) || !file.exists(hmc_frac_file)) {
  stop("Section 85: per-sample regional-fraction file(s) missing:\n  ",
       mc_frac_file, "\n  ", hmc_frac_file)
}

mc_frac  <- utils::read.table(mc_frac_file,  header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)
hmc_frac <- utils::read.table(hmc_frac_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)

build_sample_meta <- function(frac_cols, suffix) {
  ids <- sub(suffix, "", frac_cols)
  meta <- SAMPLES[match(ids, SAMPLES$sample_id), ]
  if (anyNA(meta$sample_id)) {
    stop("Section 85: regional-fraction columns did not map to SAMPLES$sample_id: ",
         paste(frac_cols[is.na(meta$sample_id)], collapse = ", "))
  }
  data.frame(
    col_name  = frac_cols,
    sample_id = meta$sample_id,
    condition = factor(meta$condition, levels = c("Control", "Mutant")),
    sex       = factor(meta$sex, levels = c("Female", "Male")),
    batch     = factor(paste("Batch", meta$batch), levels = c("Batch 1", "Batch 2")),
    short     = meta$short_name,
    stringsAsFactors = FALSE
  )
}

mc_cols  <- colnames(mc_frac)[grepl("_num_mc_region_frac$",  colnames(mc_frac))]
hmc_cols <- colnames(hmc_frac)[grepl("_num_hmc_region_frac$", colnames(hmc_frac))]
if (length(mc_cols) != 8 || length(hmc_cols) != 8) {
  stop("Section 85: expected 8 per-sample columns each; got ",
       length(mc_cols), " (mC) / ", length(hmc_cols), " (hmC).")
}

mc_meta  <- build_sample_meta(mc_cols,  "_num_mc_region_frac")
hmc_meta <- build_sample_meta(hmc_cols, "_num_hmc_region_frac")

run_methylation_pca <- function(frac, sample_cols, meta, meth_label) {
  mat <- as.matrix(frac[, sample_cols, drop = FALSE])
  for (i in seq_len(nrow(mat))) {
    na_idx <- is.na(mat[i, ])
    if (any(na_idx) && !all(na_idx)) mat[i, na_idx] <- mean(mat[i, !na_idx])
  }
  gene_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat <- mat[!is.na(gene_vars) & gene_vars > 1e-10, , drop = FALSE]

  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_explained <- summary(pca)$importance[2, ] * 100

  pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    condition = meta$condition, sex = meta$sex,
    batch = meta$batch, label = meta$short
  )

  cat(sprintf("  %s PCA: PC1 = %.1f%% variance, PC2 = %.1f%% variance\n",
              meth_label, var_explained[1], var_explained[2]))

  ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = sex)) +
    geom_point(size = 2.5, stroke = 0.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.2,
                             show.legend = FALSE, max.overlaps = Inf,
                             segment.color = NA) +
    scale_colour_manual(values = COLORS$condition, name = "Condition") +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
    labs(
      title = sprintf("%s PCA", meth_label),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_pub()
}

pca_mc  <- run_methylation_pca(mc_frac,  mc_cols,  mc_meta,  "5mC")
pca_hmc <- run_methylation_pca(hmc_frac, hmc_cols, hmc_meta, "5hmC")

save_plot(pca_mc,  "85a_5mc_pca",  w = 7, h = 6)
save_plot(pca_hmc, "85b_5hmc_pca", w = 7, h = 6)

sec85_composite <- (pca_mc | pca_hmc) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

sec85_composite <- sec85_composite +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Section 85: Per-sample methylation PCA (condition + sex)",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

save_plot(sec85_composite, "85_composite", w = 14, h = 6)

cat("\nSection 85 complete.\n")
