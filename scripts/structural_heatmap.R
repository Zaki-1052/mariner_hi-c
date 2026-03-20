# scripts/structural_heatmap.R
#
# Top 50 Genes Heatmap by Combined Structural Score (Figure 5D)
#
# Produces two heatmaps:
#   A) Top 50 genes ranked by combined structural score (boundary + loop + ABC)
#   B) Top 50 genes ranked by |ΔABC| only (enhancer-gene disruption)
#
# Usage:
#   Rscript scripts/structural_heatmap.R
#
# Input:
#   output/network_analysis/late/tables/gene_structural_profile_all.tsv
#   abc/results/gene_level_summary.tsv
#
# Output:
#   output/structural_heatmap/
#     combined_score_heatmap/   (PDF/SVG/JPG)
#     abc_only_heatmap/         (PDF/SVG/JPG)
#     tables/top50_combined_score.tsv
#     tables/top50_abc_only.tsv
#     tables/heatmap_summary.txt

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

PROFILE_FILE   <- "output/network_analysis/late/tables/gene_structural_profile_all.tsv"
ABC_FILE       <- "abc/results/gene_level_summary.tsv"
OUTPUT_DIR     <- "output/structural_heatmap"
PLOT_DIR       <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR      <- file.path(OUTPUT_DIR, "tables")

TOP_N          <- 50
BASEMEAN_MIN   <- 10
PSEUDOCOUNT    <- 1e-6

# Color palette (blue-white-red, matching reference style)
HEATMAP_COLORS <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Row annotation palettes
LAYER_COLORS <- c("1" = "#FED976", "2" = "#FD8D3C", "3" = "#E31A1C")
DYSREG_COLORS <- c("TRUE" = "#31A354", "FALSE" = "#D9D9D9")

# ==============================================================================
# 3. LOAD DATA
# ==============================================================================

load_data <- function() {
  stopifnot(file.exists(PROFILE_FILE))
  stopifnot(file.exists(ABC_FILE))

  cat("Loading gene structural profiles...\n")
  profile <- read_tsv(PROFILE_FILE, show_col_types = FALSE)
  cat(sprintf("  %d genes loaded\n", nrow(profile)))

  cat("Loading ABC gene-level summary...\n")
  abc <- read_tsv(ABC_FILE, show_col_types = FALSE) %>%
    select(TargetGene, ABC.Score_WT, ABC.Score_KO)
  cat(sprintf("  %d genes loaded\n", nrow(abc)))

  merged <- profile %>%
    left_join(abc, by = c("gene" = "TargetGene"))

  stopifnot(nrow(merged) == nrow(profile))
  cat(sprintf("  Joined: %d genes with ABC scores\n",
              sum(!is.na(merged$ABC.Score_WT))))

  merged
}

# ==============================================================================
# 4. COMPUTE logFC AxC
# ==============================================================================

add_logfc_axc <- function(df) {
  df %>%
    mutate(
      logFC_AxC = case_when(
        is.na(ABC.Score_WT) | is.na(ABC.Score_KO) ~ NA_real_,
        TRUE ~ log2((ABC.Score_KO + PSEUDOCOUNT) / (ABC.Score_WT + PSEUDOCOUNT))
      )
    )
}

# ==============================================================================
# 5. BUILD HEATMAP MATRIX (z-scored columns)
# ==============================================================================

build_zscore_matrix <- function(df, columns, col_labels = NULL) {
  mat <- df %>% select(all_of(columns)) %>% as.matrix()
  rownames(mat) <- df$gene

  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    col_sd <- sd(col_vals, na.rm = TRUE)
    col_mean <- mean(col_vals, na.rm = TRUE)
    if (!is.na(col_sd) && col_sd > 0) {
      mat[, j] <- (col_vals - col_mean) / col_sd
    }
  }

  # Cap at +/- 3 SD for visual clarity
  mat[mat > 3] <- 3
  mat[mat < -3] <- -3

  if (!is.null(col_labels)) {
    colnames(mat) <- col_labels
  }

  mat
}

# ==============================================================================
# 6. RENDER HEATMAP
# ==============================================================================

render_heatmap <- function(mat, annotation_row, annotation_colors,
                           title, base_path, width, height) {

  max_abs <- max(abs(mat), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)

  hm_call <- bquote(
    pheatmap(
      .(mat),
      main              = .(title),
      color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks            = .(breaks),
      cluster_rows      = TRUE,
      cluster_cols      = FALSE,
      clustering_method = "ward.D2",
      show_rownames     = TRUE,
      show_colnames     = TRUE,
      fontsize_row      = 8,
      fontsize_col      = 11,
      fontsize          = 10,
      annotation_row    = .(annotation_row),
      annotation_colors = .(annotation_colors),
      border_color      = NA,
      cellwidth         = 40,
      cellheight        = 10,
      angle_col         = 45
    )
  )

  save_multiformat_pheatmap(hm_call, base_path, width = width, height = height)
  cat(sprintf("  Saved: %s\n", base_path))
}

# ==============================================================================
# 7. HEATMAP A — TOP 50 BY COMBINED STRUCTURAL SCORE
# ==============================================================================

make_combined_score_heatmap <- function(df) {
  cat("\n=== Heatmap A: Top 50 by Combined Structural Score ===\n")

  top50 <- df %>%
    filter(baseMean >= BASEMEAN_MIN) %>%
    filter(!is.na(combined_score)) %>%
    arrange(desc(combined_score)) %>%
    slice_head(n = TOP_N) %>%
    mutate(
      loop_logFC_signed = replace_na(mean_loop_logFC, 0),
      boundary_flag     = as.numeric(has_boundary)
    )

  cat(sprintf("  Top gene: %s (combined_score=%.2f, n_layers=%d)\n",
              top50$gene[1], top50$combined_score[1], top50$n_layers[1]))
  cat(sprintf("  Score range: %.2f to %.2f\n",
              min(top50$combined_score), max(top50$combined_score)))

  columns <- c("log2FC", "logFC_AxC", "loop_logFC_signed", "boundary_flag")
  col_labels <- c("mRNA logFC", "AxC logFC", "Loop logFC", "Boundary")

  mat <- build_zscore_matrix(top50, columns, col_labels)

  annotation_row <- data.frame(
    Layers       = factor(top50$n_layers),
    Dysregulated = factor(as.character(top50$dysregulated)),
    row.names    = top50$gene
  )

  annotation_colors <- list(
    Layers       = LAYER_COLORS,
    Dysregulated = DYSREG_COLORS
  )

  render_heatmap(
    mat, annotation_row, annotation_colors,
    title     = "Top 50 Genes by Combined Structural Score",
    base_path = file.path(PLOT_DIR, "combined_score_heatmap"),
    width     = 10, height = 14
  )

  write_tsv(top50, file.path(TABLE_DIR, "top50_combined_score.tsv"))
  cat(sprintf("  Saved table: %d genes\n", nrow(top50)))

  top50
}

# ==============================================================================
# 8. HEATMAP B — TOP 50 BY ABC ONLY
# ==============================================================================

make_abc_only_heatmap <- function(df) {
  cat("\n=== Heatmap B: Top 50 by |delta ABC| ===\n")

  top50 <- df %>%
    filter(baseMean >= BASEMEAN_MIN) %>%
    filter(!is.na(max_delta_unnorm)) %>%
    arrange(desc(abs(max_delta_unnorm))) %>%
    slice_head(n = TOP_N)

  cat(sprintf("  Top gene: %s (|ΔABC|=%.4f, logFC_AxC=%.2f)\n",
              top50$gene[1], abs(top50$max_delta_unnorm[1]),
              top50$logFC_AxC[1]))

  columns <- c("log2FC", "logFC_AxC")
  col_labels <- c("mRNA logFC", "AxC logFC")

  mat <- build_zscore_matrix(top50, columns, col_labels)

  n_enh_breaks <- cut(top50$n_enhancers,
                      breaks = c(0, 5, 10, 15, Inf),
                      labels = c("1-5", "6-10", "11-15", ">15"))

  annotation_row <- data.frame(
    Dysregulated = factor(as.character(top50$dysregulated)),
    Enhancers    = n_enh_breaks,
    row.names    = top50$gene
  )

  enh_colors <- c("1-5" = "#DEEBF7", "6-10" = "#9ECAE1",
                   "11-15" = "#4292C6", ">15" = "#08519C")

  annotation_colors <- list(
    Dysregulated = DYSREG_COLORS,
    Enhancers    = enh_colors
  )

  render_heatmap(
    mat, annotation_row, annotation_colors,
    title     = "Top 50 Genes by |AxC Change|",
    base_path = file.path(PLOT_DIR, "abc_only_heatmap"),
    width     = 8, height = 14
  )

  write_tsv(top50, file.path(TABLE_DIR, "top50_abc_only.tsv"))
  cat(sprintf("  Saved table: %d genes\n", nrow(top50)))

  top50
}

# ==============================================================================
# 9. SUMMARY REPORT
# ==============================================================================

write_summary <- function(df, top_combined, top_abc) {
  sink(file.path(TABLE_DIR, "heatmap_summary.txt"))

  cat("=== Structural Heatmap Summary (Figure 5D) ===\n\n")
  cat(sprintf("Total genes in profile: %d\n", nrow(df)))
  cat(sprintf("After baseMean >= %d filter: %d\n",
              BASEMEAN_MIN, sum(df$baseMean >= BASEMEAN_MIN, na.rm = TRUE)))

  cat("\n--- Heatmap A: Combined Structural Score ---\n")
  cat(sprintf("Ranking: combined_score (z_boundary + z_loop_count + z_loop_strength + z_abc)\n"))
  cat(sprintf("Top 50 score range: %.2f to %.2f\n",
              min(top_combined$combined_score), max(top_combined$combined_score)))
  cat(sprintf("Layer distribution in top 50:\n"))
  print(table(top_combined$n_layers))
  cat(sprintf("Dysregulated in top 50: %d / %d\n",
              sum(top_combined$dysregulated), nrow(top_combined)))

  overlap_genes <- intersect(top_combined$gene, top_abc$gene)
  cat(sprintf("\nOverlap with ABC-only top 50: %d genes\n", length(overlap_genes)))
  if (length(overlap_genes) > 0) {
    cat(sprintf("  Shared: %s\n", paste(overlap_genes, collapse = ", ")))
  }

  cat("\n--- Heatmap B: ABC-Only ---\n")
  cat(sprintf("Ranking: |max_delta_unnorm|\n"))
  cat(sprintf("Top 50 |ΔABC| range: %.4f to %.4f\n",
              min(abs(top_abc$max_delta_unnorm)),
              max(abs(top_abc$max_delta_unnorm))))
  cat(sprintf("Dysregulated in top 50: %d / %d\n",
              sum(top_abc$dysregulated), nrow(top_abc)))

  cat("\n--- Columns ---\n")
  cat("Heatmap A: mRNA logFC, AxC logFC, Loop logFC, Boundary (all z-scored)\n")
  cat("Heatmap B: mRNA logFC, AxC logFC (z-scored)\n")
  cat("logFC_AxC = log2((ABC.Score_KO + 1e-6) / (ABC.Score_WT + 1e-6))\n")

  sink()
  cat("  Saved summary report\n")
}

# ==============================================================================
# 10. MAIN
# ==============================================================================

main <- function() {
  cat("=== Structural Heatmap (Figure 5D) ===\n\n")

  dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  df <- load_data() %>%
    add_logfc_axc()

  cat(sprintf("  logFC_AxC computed for %d genes\n",
              sum(!is.na(df$logFC_AxC))))

  top_combined <- make_combined_score_heatmap(df)
  top_abc      <- make_abc_only_heatmap(df)

  write_summary(df, top_combined, top_abc)

  cat("\n=== Done ===\n")
  cat(sprintf("Output: %s\n", OUTPUT_DIR))
}

main()
