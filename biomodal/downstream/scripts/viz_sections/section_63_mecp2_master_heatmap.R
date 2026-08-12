# biomodal/downstream/scripts/viz_sections/section_63_mecp2_master_heatmap.R
# Section 63: Master Multi-Mark Heatmap at MeCP2 Binding Sites
#
# Clustered heatmap showing all chromatin marks at MeCP2-significant peaks,
# annotated by MeCP2 direction, chromatin state, and CG concordance class.
# Reveals whether MeCP2-Up peaks cluster into distinct Polycomb-enriched profiles.
#
# Figures:
#   63a: Full multi-mark heatmap (ctrl signal, Z-scored, pheatmap)
#   63b: Per-cluster mean mark profile bar chart
#   63c: Cluster composition stacked bar (chromatin state per cluster)
#   63d: Log2FC heatmap (mut/ctrl change per mark)
#
# Input:
#   62_mecp2_peak_chromatin_signal.tsv (from section 62)
#   56_mecp2_peak_annotated.tsv (for cg_class)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_63_mecp2_master_heatmap.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(pheatmap)
  library(patchwork)
})

# =============================================================================
# SECTION 63 CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 63: MASTER MULTI-MARK HEATMAP AT MeCP2 BINDING SITES\n")
cat("================================================================================\n\n")

SEC63_DIR <- file.path(OUTPUT_DIR, "63_mecp2_master_heatmap")
dir.create(SEC63_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNAL_PATH <- file.path(TABLES_DIR, "62_mecp2_peak_chromatin_signal.tsv")
SEC56_PATH  <- file.path(TABLES_DIR, "56_mecp2_peak_annotated.tsv")

N_CLUSTERS <- 4
MAX_PEAKS_HEATMAP <- 5000

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC63_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n")

if (!file.exists(SIGNAL_PATH)) stop("Run section 62 first: ", SIGNAL_PATH)
peaks <- read.table(SIGNAL_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Signal table: %d peaks × %d columns\n", nrow(peaks), ncol(peaks)))

peaks$mecp2_class <- dplyr::case_when(
  peaks$FDR < 0.05 & peaks$Fold > 0 ~ "MeCP2 Up",
  peaks$FDR < 0.05 & peaks$Fold < 0 ~ "MeCP2 Down",
  TRUE ~ "Not Significant"
)

# Merge CG class from section 56
if (file.exists(SEC56_PATH)) {
  sec56 <- read.table(SEC56_PATH, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  peaks <- merge(peaks, sec56[, c("Chr", "Start", "End", "cg_class")],
                 by = c("Chr", "Start", "End"), all.x = TRUE)
  peaks$cg_class[is.na(peaks$cg_class)] <- "Unknown"
} else {
  peaks$cg_class <- "Unknown"
}

sig_peaks <- peaks[peaks$mecp2_class != "Not Significant", ]
cat(sprintf("  Significant MeCP2 peaks: %d (Up=%d, Down=%d)\n",
            nrow(sig_peaks),
            sum(sig_peaks$mecp2_class == "MeCP2 Up"),
            sum(sig_peaks$mecp2_class == "MeCP2 Down")))

# =============================================================================
# CLASSIFY CHROMATIN STATE
# =============================================================================

cat("\n--- Classifying chromatin state ---\n")

sig_gr <- GRanges(seqnames = sig_peaks$Chr,
                  ranges = IRanges(start = sig_peaks$Start, end = sig_peaks$End))

chip_overlaps <- data.frame(
  H3K27ac_overlap  = countOverlaps(sig_gr, load_chip_peaks(CHIP_PEAK_FILES$h3k27ac)) > 0,
  H3K27me3_overlap = countOverlaps(sig_gr, load_chip_peaks(CHIP_PEAK_FILES$h3k27me3)) > 0,
  H3K4me1_overlap  = countOverlaps(sig_gr, load_chip_peaks(CHIP_PEAK_FILES$h3k4me1)) > 0,
  H3K4me3_overlap  = countOverlaps(sig_gr, load_chip_peaks(CHIP_PEAK_FILES$h3k4me3)) > 0,
  Bivalent_overlap = countOverlaps(sig_gr, load_chip_peaks(CHIP_PEAK_FILES$bivalent)) > 0
)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
tss_gr <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
dist_to_tss <- GenomicRanges::distanceToNearest(sig_gr, tss_gr)
sig_peaks$dist_to_tss <- NA_real_
sig_peaks$dist_to_tss[queryHits(dist_to_tss)] <- mcols(dist_to_tss)$distance

sig_peaks$chromatin_state <- classify_chromatin_state(chip_overlaps, sig_peaks$dist_to_tss)

state_ct <- table(sig_peaks$chromatin_state)
cat("  Chromatin states:\n")
for (s in names(sort(state_ct, decreasing = TRUE))) {
  cat(sprintf("    %-20s %d\n", s, state_ct[s]))
}

# =============================================================================
# 63a: CTRL SIGNAL HEATMAP
# =============================================================================

cat("\n--- 63a: Multi-mark heatmap (ctrl signal) ---\n")

ctrl_cols <- c("mecp2_ctrl", "k27me3_ctrl", "k119ub_ctrl", "k27me1_ctrl",
               "k27ac_ctrl", "atac_ctrl", "k4me3_ctrl", "k36me3_ctrl",
               "cg_mc_ctrl", "cg_hmc_ctrl")

display_names <- c("MeCP2", "H3K27me3", "H2AK119ub", "H3K27me1",
                   "H3K27ac", "ATAC-seq", "H3K4me3", "H3K36me3",
                   "CG 5mC", "CG 5hmC")

available_cols <- intersect(ctrl_cols, colnames(sig_peaks))
if (length(available_cols) < 5) {
  cat("  WARNING: Too few mark columns available, skipping heatmap\n")
} else {
  mat <- as.matrix(sig_peaks[, available_cols])

  # Log2-transform
  for (j in seq_len(ncol(mat))) {
    vals <- mat[, j]
    nz <- vals[!is.na(vals) & vals > 0]
    pseudo <- if (length(nz) > 0) max(quantile(nz, 0.05), 1e-6) else 1e-6
    mat[, j] <- log2(vals + pseudo)
  }

  # Z-score per column
  for (j in seq_len(ncol(mat))) {
    vals <- mat[, j]
    mu <- mean(vals, na.rm = TRUE)
    sd_val <- sd(vals, na.rm = TRUE)
    if (!is.na(sd_val) && sd_val > 0) {
      mat[, j] <- (vals - mu) / sd_val
    }
  }

  # Replace NAs with 0 for clustering
  mat[is.na(mat)] <- 0

  # Subsample if too many peaks
  if (nrow(mat) > MAX_PEAKS_HEATMAP) {
    set.seed(42)
    sample_idx <- sample(nrow(mat), MAX_PEAKS_HEATMAP)
    mat_plot <- mat[sample_idx, ]
    ann_idx <- sample_idx
    cat(sprintf("  Subsampled %d → %d peaks for heatmap\n", nrow(mat), MAX_PEAKS_HEATMAP))
  } else {
    mat_plot <- mat
    ann_idx <- seq_len(nrow(mat))
  }

  # Rename columns
  col_idx <- match(available_cols, ctrl_cols)
  colnames(mat_plot) <- display_names[col_idx]

  # Row annotation
  ann_row <- data.frame(
    MeCP2 = sig_peaks$mecp2_class[ann_idx],
    State = sig_peaks$chromatin_state[ann_idx],
    row.names = paste0("p", seq_len(nrow(mat_plot)))
  )
  rownames(mat_plot) <- rownames(ann_row)

  ann_colors <- list(
    MeCP2 = c("MeCP2 Up" = "#D73027", "MeCP2 Down" = "#4575B4"),
    State = CHROMATIN_STATE_COLORS[names(CHROMATIN_STATE_COLORS) %in% unique(ann_row$State)]
  )

  # Clip extreme Z-scores
  mat_plot[mat_plot > 3] <- 3
  mat_plot[mat_plot < -3] <- -3

  heatmap_colors <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)

  cat("  Rendering heatmap...\n")

  # Run once silently to get the tree for cluster extraction
  ph_obj <- pheatmap(
    mat_plot,
    color = heatmap_colors,
    cluster_rows = TRUE, cluster_cols = FALSE,
    clustering_method = "ward.D2",
    cutree_rows = N_CLUSTERS,
    show_rownames = FALSE,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    fontsize = 10, fontsize_col = 11,
    main = sprintf("Multi-mark signal at MeCP2 binding sites (n=%d, Z-scored)", nrow(mat_plot)),
    silent = TRUE
  )

  # save_multiformat_pheatmap expects a quote()'d expression
  pheatmap_call <- quote(pheatmap(
    mat_plot,
    color = heatmap_colors,
    cluster_rows = TRUE, cluster_cols = FALSE,
    clustering_method = "ward.D2",
    cutree_rows = N_CLUSTERS,
    show_rownames = FALSE,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    fontsize = 10, fontsize_col = 11,
    main = sprintf("Multi-mark signal at MeCP2 binding sites (n=%d, Z-scored)", nrow(mat_plot))
  ))

  save_multiformat_pheatmap(
    pheatmap_call,
    base_path = file.path(SEC63_DIR, "63a_multimark_heatmap"),
    width = 14, height = 12, dpi = 300, verbose = TRUE
  )

  # Extract cluster assignments from the silent run
  row_order <- ph_obj$tree_row$order
  clusters <- cutree(ph_obj$tree_row, k = N_CLUSTERS)

  # =============================================================================
  # 63b: PER-CLUSTER MEAN MARK PROFILE
  # =============================================================================

  cat("\n--- 63b: Per-cluster mark profiles ---\n")

  cluster_df <- data.frame(
    peak_id = names(clusters),
    cluster = paste0("Cluster ", clusters),
    stringsAsFactors = FALSE
  )
  cluster_df <- cbind(cluster_df, as.data.frame(mat_plot[cluster_df$peak_id, ]))

  cluster_long <- tidyr::pivot_longer(
    cluster_df, cols = all_of(colnames(mat_plot)),
    names_to = "mark", values_to = "zscore"
  )
  cluster_long$mark <- factor(cluster_long$mark, levels = colnames(mat_plot))

  cluster_means <- cluster_long %>%
    dplyr::group_by(cluster, mark) %>%
    dplyr::summarize(mean_z = mean(zscore, na.rm = TRUE),
                     se = sd(zscore, na.rm = TRUE) / sqrt(dplyr::n()),
                     .groups = "drop")

  p_63b <- ggplot(cluster_means, aes(x = mark, y = mean_z, fill = mark)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_z - se, ymax = mean_z + se), width = 0.3) +
    facet_wrap(~ cluster, ncol = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c(
      "MeCP2" = "#333333", "H3K27me3" = "#984EA3", "H2AK119ub" = "#FF7F00",
      "H3K27me1" = "#A65628", "H3K27ac" = "#4DAF4A", "ATAC-seq" = "#377EB8",
      "H3K4me3" = "#E41A1C", "H3K36me3" = "#F781BF",
      "CG 5mC" = "#D6604D", "CG 5hmC" = "#FDB863"
    )) +
    labs(title = "Mean mark profile per cluster",
         x = NULL, y = "Mean Z-score") +
    theme_biomodal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p_63b, "63b_cluster_profiles", w = 12, h = 9)

  # =============================================================================
  # 63c: CLUSTER COMPOSITION (chromatin state per cluster)
  # =============================================================================

  cat("\n--- 63c: Cluster composition ---\n")

  comp_df <- data.frame(
    cluster = cluster_df$cluster,
    state = ann_row$State,
    mecp2 = ann_row$MeCP2,
    stringsAsFactors = FALSE
  )

  # Chromatin state by cluster
  state_by_cluster <- comp_df %>%
    dplyr::count(cluster, state) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()

  p_63c_state <- ggplot(state_by_cluster, aes(x = cluster, y = frac, fill = state)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Chromatin state composition per cluster",
         x = NULL, y = "Fraction") +
    theme_biomodal()

  # MeCP2 direction by cluster
  dir_by_cluster <- comp_df %>%
    dplyr::count(cluster, mecp2) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()

  p_63c_dir <- ggplot(dir_by_cluster, aes(x = cluster, y = frac, fill = mecp2)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("MeCP2 Up" = "#D73027", "MeCP2 Down" = "#4575B4"),
                      name = "MeCP2 Direction") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "MeCP2 direction per cluster",
         x = NULL, y = "Fraction") +
    theme_biomodal()

  p_63c <- p_63c_state + p_63c_dir +
    plot_annotation(title = "Cluster composition analysis")

  save_plot(p_63c, "63c_cluster_composition", w = 14, h = 6)

  # Print cluster summary
  cat("  Cluster summary:\n")
  for (cl in sort(unique(cluster_df$cluster))) {
    cl_mask <- cluster_df$cluster == cl
    n_cl <- sum(cl_mask)
    n_up <- sum(comp_df$mecp2[cl_mask] == "MeCP2 Up")
    n_down <- sum(comp_df$mecp2[cl_mask] == "MeCP2 Down")
    top_state <- names(sort(table(comp_df$state[cl_mask]), decreasing = TRUE))[1]
    cat(sprintf("    %s: n=%d  Up=%d  Down=%d  top_state=%s\n",
                cl, n_cl, n_up, n_down, top_state))
  }

  # =============================================================================
  # 63d: LOG2FC HEATMAP
  # =============================================================================

  cat("\n--- 63d: Log2FC heatmap ---\n")

  fc_cols <- c("mecp2_log2fc", "k27me3_log2fc", "k119ub_log2fc", "k27me1_log2fc",
               "k27ac_log2fc", "atac_log2fc", "k4me3_log2fc", "k36me3_log2fc",
               "cg_mc_log2fc", "cg_hmc_log2fc")
  fc_display <- c("MeCP2 Δ", "K27me3 Δ", "K119ub Δ", "K27me1 Δ",
                  "K27ac Δ", "ATAC Δ", "K4me3 Δ", "K36me3 Δ",
                  "CG mC Δ", "CG hmC Δ")

  available_fc <- intersect(fc_cols, colnames(sig_peaks))

  if (length(available_fc) >= 5) {
    mat_fc <- as.matrix(sig_peaks[ann_idx, available_fc])
    mat_fc[is.na(mat_fc)] <- 0

    # Clip extremes
    mat_fc[mat_fc > 3] <- 3
    mat_fc[mat_fc < -3] <- -3

    fc_idx <- match(available_fc, fc_cols)
    colnames(mat_fc) <- fc_display[fc_idx]
    rownames(mat_fc) <- rownames(mat_plot)

    pheatmap_fc_call <- quote(pheatmap(
      mat_fc,
      color = heatmap_colors,
      cluster_rows = TRUE, cluster_cols = FALSE,
      clustering_method = "ward.D2",
      cutree_rows = N_CLUSTERS,
      show_rownames = FALSE,
      annotation_row = ann_row,
      annotation_colors = ann_colors,
      fontsize = 10, fontsize_col = 11,
      main = sprintf("Log2FC (mut/ctrl) at MeCP2 binding sites (n=%d)", nrow(mat_fc))
    ))

    save_multiformat_pheatmap(
      pheatmap_fc_call,
      base_path = file.path(SEC63_DIR, "63d_log2fc_heatmap"),
      width = 14, height = 12, dpi = 300, verbose = TRUE
    )
  }

  # =============================================================================
  # SAVE TABLES
  # =============================================================================

  cat("\n--- Saving tables ---\n")

  signal_out <- data.frame(sig_peaks[ann_idx, c("Chr", "Start", "End", "mecp2_class",
                                                 "chromatin_state", "cg_class")],
                           mat_plot, check.names = FALSE)
  write.table(signal_out, file.path(TABLES_DIR, "63_peak_signal_matrix.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved: 63_peak_signal_matrix.tsv\n")

  cluster_out <- data.frame(
    Chr = sig_peaks$Chr[ann_idx],
    Start = sig_peaks$Start[ann_idx],
    End = sig_peaks$End[ann_idx],
    cluster = cluster_df$cluster,
    mecp2_class = ann_row$MeCP2,
    chromatin_state = ann_row$State,
    stringsAsFactors = FALSE
  )
  write.table(cluster_out, file.path(TABLES_DIR, "63_cluster_assignments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved: 63_cluster_assignments.tsv\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 63 SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("  Significant MeCP2 peaks: %d\n", nrow(sig_peaks)))
cat(sprintf("  Heatmap peaks: %d (subsampled if >%d)\n",
            min(nrow(sig_peaks), MAX_PEAKS_HEATMAP), MAX_PEAKS_HEATMAP))
cat(sprintf("  Clusters: %d\n", N_CLUSTERS))

cat("\n  Plots saved to:", SEC63_DIR, "\n")
cat("    63a: Multi-mark ctrl signal heatmap (pheatmap)\n")
cat("    63b: Per-cluster mark profile bars\n")
cat("    63c: Cluster composition (chromatin state + MeCP2 direction)\n")
cat("    63d: Log2FC heatmap\n")

cat("\nSection 63 complete.\n")
