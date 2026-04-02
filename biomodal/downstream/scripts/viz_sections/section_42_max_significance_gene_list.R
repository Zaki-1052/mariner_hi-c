# biomodal/downstream/scripts/viz_sections/section_42_max_significance_gene_list.R
# Section 42: Extract genes at maximum significance (q-value floor)
# Identifies genes whose q-values hit the numerical floor (~0), appearing at
# the -log10(q) = 300 ceiling in volcano plots. Exports merged mC/hmC table.
# Also produces PCA of per-sample methylation across these genes.

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 42: MAXIMUM SIGNIFICANCE GENE LIST
# =============================================================================

cat("================================================================================\n")
cat("SECTION 42: MAXIMUM SIGNIFICANCE GENE LIST\n")
cat("================================================================================\n\n")

# Threshold: q-values that reach the numerical floor used by load_dmr_bed()
# load_dmr_bed() caps at pmax(q, 1e-300) => neg_log10_q = 300 for q < 1e-300
Q_FLOOR <- 1e-300

# ---- Extract floor-significance genes from each modality --------------------

mc_floor <- mc_dmr %>%
  dplyr::filter(dmr_qvalue < Q_FLOOR) %>%
  dplyr::select(chr, start, end, gene, num_contexts, mean_coverage,
                mean_mod_group1, mean_mod_group2,
                mod_difference, dmr_qvalue, direction) %>%
  dplyr::rename(
    mc_difference  = mod_difference,
    mc_qvalue      = dmr_qvalue,
    mc_direction   = direction,
    mc_ctrl_mean   = mean_mod_group1,
    mc_mut_mean    = mean_mod_group2
  )

hmc_floor <- hmc_dmr %>%
  dplyr::filter(dmr_qvalue < Q_FLOOR) %>%
  dplyr::select(chr, start, end, gene, mod_difference, dmr_qvalue, direction,
                mean_mod_group1, mean_mod_group2) %>%
  dplyr::rename(
    hmc_difference = mod_difference,
    hmc_qvalue     = dmr_qvalue,
    hmc_direction  = direction,
    hmc_ctrl_mean  = mean_mod_group1,
    hmc_mut_mean   = mean_mod_group2
  )

cat(sprintf("  5mC genes at q-value floor: %d\n", nrow(mc_floor)))
cat(sprintf("  5hmC genes at q-value floor: %d\n", nrow(hmc_floor)))

# ---- Merge on gene name ----------------------------------------------------

# Join on gene + coordinates to handle genes with multiple isoforms (e.g. Arhgap26)
merged <- dplyr::full_join(mc_floor, hmc_floor, by = c("gene", "chr", "start", "end"))

# Flag which modalities hit the floor
merged$mc_at_floor  <- !is.na(merged$mc_difference)
merged$hmc_at_floor <- !is.na(merged$hmc_difference)

# For genes at floor in both: classify coordinated pattern
merged$coordinated <- dplyr::case_when(
  merged$mc_at_floor & merged$hmc_at_floor &
    merged$mc_direction == "Hypermethylated" &
    merged$hmc_direction == "Hypomethylated" ~ "mC_up_hmC_down",
  merged$mc_at_floor & merged$hmc_at_floor &
    merged$mc_direction == "Hypomethylated" &
    merged$hmc_direction == "Hypermethylated" ~ "mC_down_hmC_up",
  merged$mc_at_floor & merged$hmc_at_floor ~ "same_direction",
  TRUE ~ "single_modality"
)

# Convert differences to percentage for readability
merged$mc_diff_pct  <- round(merged$mc_difference * 100, 2)
merged$hmc_diff_pct <- round(merged$hmc_difference * 100, 2)

# Sort by absolute mC effect size (largest first), then hmC
merged <- merged %>%
  dplyr::arrange(desc(abs(mc_diff_pct)), desc(abs(hmc_diff_pct)))

# ---- Summary statistics -----------------------------------------------------

n_both  <- sum(merged$mc_at_floor & merged$hmc_at_floor)
n_mc_only  <- sum(merged$mc_at_floor & !merged$hmc_at_floor)
n_hmc_only <- sum(!merged$mc_at_floor & merged$hmc_at_floor)
n_coordinated <- sum(merged$coordinated == "mC_up_hmC_down", na.rm = TRUE)

cat(sprintf("\n  At floor in both mC and hmC: %d\n", n_both))
cat(sprintf("  At floor in mC only: %d\n", n_mc_only))
cat(sprintf("  At floor in hmC only: %d\n", n_hmc_only))
cat(sprintf("  Coordinated mC-up/hmC-down (TET block): %d / %d (%.1f%%)\n",
            n_coordinated, n_both,
            ifelse(n_both > 0, n_coordinated / n_both * 100, 0)))

# ---- Export tables -----------------------------------------------------------

# Full merged table
output_merged <- merged %>%
  dplyr::select(gene, chr, start, end,
                mc_diff_pct, mc_direction, mc_ctrl_mean, mc_mut_mean, mc_at_floor,
                hmc_diff_pct, hmc_direction, hmc_ctrl_mean, hmc_mut_mean, hmc_at_floor,
                coordinated, num_contexts, mean_coverage)

write.table(output_merged,
            file.path(TABLES_DIR, "max_significance_genes_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_merged.tsv"),
            nrow(output_merged)))

# mC-only list (sorted by effect size)
mc_only_table <- mc_floor %>%
  dplyr::mutate(mc_diff_pct = round(mc_difference * 100, 2)) %>%
  dplyr::arrange(desc(mc_diff_pct)) %>%
  dplyr::select(gene, chr, start, end, mc_diff_pct, mc_direction,
                mc_ctrl_mean, mc_mut_mean, num_contexts, mean_coverage)

write.table(mc_only_table,
            file.path(TABLES_DIR, "max_significance_genes_mc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_mc.tsv"),
            nrow(mc_only_table)))

# hmC-only list (sorted by effect size)
hmc_only_table <- hmc_floor %>%
  dplyr::mutate(hmc_diff_pct = round(hmc_difference * 100, 2)) %>%
  dplyr::arrange(hmc_diff_pct) %>%
  dplyr::select(gene, hmc_diff_pct, hmc_direction,
                hmc_ctrl_mean, hmc_mut_mean)

write.table(hmc_only_table,
            file.path(TABLES_DIR, "max_significance_genes_hmc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_hmc.tsv"),
            nrow(hmc_only_table)))

# Gene-only list (just names, one per line) for pathway tools / GREAT / etc.
gene_list <- sort(unique(merged$gene))
writeLines(gene_list,
           file.path(TABLES_DIR, "max_significance_gene_names.txt"))
cat(sprintf("  Saved: %s (%d unique gene names)\n",
            file.path(TABLES_DIR, "max_significance_gene_names.txt"),
            length(gene_list)))

# =============================================================================
# PCA OF MAX-SIGNIFICANCE GENES (per-sample methylation)
# =============================================================================

cat("\n--- PCA of max-significance genes ---\n\n")

# Per-sample regional fraction files from modality feature extraction
EXTRACT_DIR <- file.path(BASE_DIR,
  "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260329_201021")

mc_frac_file  <- file.path(EXTRACT_DIR, "Extract_mc_regional-frac_20260329_201021.tsv.gz")
hmc_frac_file <- file.path(EXTRACT_DIR, "Extract_hmc_regional-frac_20260329_201021.tsv.gz")

stopifnot(file.exists(mc_frac_file), file.exists(hmc_frac_file))

# Load per-sample matrices
mc_frac  <- read.table(mc_frac_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hmc_frac <- read.table(hmc_frac_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Sample columns are columns 4:11 (between coordinates and annotation)
sample_cols <- colnames(mc_frac)[4:11]

# Build short sample labels and metadata from column names
sample_meta <- data.frame(
  col_name  = sample_cols,
  condition = ifelse(grepl("ctrl", sample_cols), "Control", "Mutant"),
  sex       = ifelse(grepl("_F", sample_cols), "Female", "Male"),
  batch     = ifelse(grepl("B2", sample_cols), "Batch 2", "Batch 1"),
  short     = gsub("evoC.Bap1.", "", sample_cols) %>%
    gsub("_num_mc_region_frac|_num_hmc_region_frac", "", .),
  stringsAsFactors = FALSE
)

# Subset to max-significance genes
mc_sub  <- mc_frac[mc_frac$Name %in% gene_list, ]
hmc_sub <- hmc_frac[hmc_frac$Name %in% gene_list, ]

cat(sprintf("  Genes matched in mC matrix: %d / %d\n", nrow(mc_sub), length(gene_list)))
cat(sprintf("  Genes matched in hmC matrix: %d / %d\n", nrow(hmc_sub), length(gene_list)))

# Helper: run PCA on a gene x sample matrix and return a ggplot
run_pca_plot <- function(df, sample_cols, meta, meth_label) {
  # Genes as rows, samples as columns => transpose for prcomp (samples as rows)
  mat <- as.matrix(df[, sample_cols])
  rownames(mat) <- df$Name

  # Remove genes with zero variance (constant across samples)
  var_filter <- apply(mat, 1, var, na.rm = TRUE) > 0
  mat <- mat[var_filter, , drop = FALSE]

  # Replace any NAs with row means
  for (i in seq_len(nrow(mat))) {
    na_idx <- is.na(mat[i, ])
    if (any(na_idx)) mat[i, na_idx] <- mean(mat[i, !na_idx])
  }

  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

  var_explained <- summary(pca)$importance[2, ] * 100  # % variance

  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = meta$condition,
    sex       = meta$sex,
    batch     = meta$batch,
    label     = meta$short
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
    geom_point(size = 4, stroke = 1.2) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, show.legend = FALSE) +
    scale_color_manual(values = COLORS$condition, name = "Condition") +
    scale_shape_manual(values = c("Batch 1" = 16, "Batch 2" = 17), name = "Batch") +
    labs(
      title = sprintf("PCA: %s at Max-Significance Genes", meth_label),
      subtitle = sprintf("%d genes | q < 1e-300", nrow(mat)),
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
      y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
    ) +
    theme_biomodal() +
    theme(legend.position = "right")
}

# Build PCA plots
p_pca_mc  <- run_pca_plot(mc_sub,  sample_cols, sample_meta, "5mC")

# hmC sample columns have different suffix
hmc_sample_cols <- colnames(hmc_frac)[4:11]
p_pca_hmc <- run_pca_plot(hmc_sub, hmc_sample_cols, sample_meta, "5hmC")

# Combined panel
p_pca_combined <- p_pca_mc | p_pca_hmc
p_pca_combined <- p_pca_combined +
  patchwork::plot_annotation(
    title = "PCA of Max-Significance Genes: 5mC vs 5hmC",
    subtitle = "Per-sample regional methylation fraction | BAP1-KO vs Control",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

# Save
SECTION_DIR <- file.path(OUTPUT_DIR, "42_max_significance")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

save_multiformat_ggplot(p_pca_combined,
                        file.path(SECTION_DIR, "42_pca_max_significance_genes"),
                        width = 16, height = 7)

cat(sprintf("  Saved PCA plot to: %s\n", SECTION_DIR))

cat("\n✅ Section 42 complete.\n")
