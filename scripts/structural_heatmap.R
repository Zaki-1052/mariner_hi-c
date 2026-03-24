# scripts/structural_heatmap.R
#
# Top 100 Genes Heatmap by Combined Structural Score (Figure 5D)
#
# Produces two heatmaps:
#   A) Top 100 genes ranked by combined structural score (boundary + loop + ABC)
#      with modification fold-change columns (H2AK119ub, H3K27me3, H3K27ac, ATAC)
#   B) Top 100 genes ranked by |ΔABC| only (enhancer-gene disruption)
#
# Usage:
#   Rscript scripts/structural_heatmap.R
#
# Input:
#   output/network_analysis/late/tables/gene_structural_profile_all.tsv
#   abc/results/gene_level_summary.tsv
#   peaks/diffbind/*.txt (DiffBind results for 4 marks)
#   abc/results/delta_abc_all_pairs.tsv (enhancer-gene linkage)
#
# Output:
#   output/structural_heatmap/
#     plots/combined_score_heatmap/  (PDF/SVG/JPG)
#     plots/abc_only_heatmap/        (PDF/SVG/JPG)
#     tables/top100_combined_score.tsv
#     tables/top100_abc_only.tsv
#     tables/heatmap_summary.txt

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
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

TOP_N          <- 100
BASEMEAN_MIN   <- 10
PSEUDOCOUNT    <- 1e-6

# Color palette (blue-white-red, matching reference style)
HEATMAP_COLORS <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Row annotation palettes
LAYER_COLORS <- c("1" = "#FED976", "2" = "#FD8D3C", "3" = "#E31A1C")
DYSREG_COLORS <- c("TRUE" = "#31A354", "FALSE" = "#D9D9D9")

# --- DiffBind files (per-peak fold-changes, all same column format) ---
DIFFBIND_FILES <- list(
  H2AK119ub = "peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt",
  H3K27me3  = "peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt",
  H3K27ac   = "peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt",
  ATAC      = "peaks/diffbind/ATAC_allATAC_diffbind_results_summit_appended_ap.txt"
)

# Peak-to-gene linkage
ABC_PAIRS_FILE     <- "abc/results/delta_abc_all_pairs.tsv"
LOOP_ASSIGN_FILE   <- "abc/results/loops_with_gene_assignments.tsv"
USE_ABC            <- TRUE
USE_LOOP_FALLBACK  <- FALSE

TSS_REGION <- c(-2000, 2000)
TXDB       <- TxDb.Mmusculus.UCSC.mm10.knownGene

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
# 5. MODIFICATION-GENE LINKAGE
# ==============================================================================

load_diffbind_file <- function(mark_name, file_path) {
  stopifnot(file.exists(file_path))

  df <- read_tsv(file_path, show_col_types = FALSE)

  required_cols <- c("Summit_Chr", "Summit_Start", "Summit_End",
                     "Peak_Chr", "Peak_Start", "Peak_End", "Fold", "FDR")
  stopifnot(all(required_cols %in% names(df)))

  df <- df %>%
    select(all_of(required_cols)) %>%
    filter(!is.na(Summit_Start) & !is.na(Summit_End)) %>%
    filter(!is.na(Peak_Start) & !is.na(Peak_End)) %>%
    mutate(mark = mark_name)

  cat(sprintf("  %s: %d peaks loaded\n", mark_name, nrow(df)))
  df
}

annotate_peaks_promoter_distal <- function(peaks_df) {
  cat("  Annotating peaks with ChIPseeker...\n")

  peaks_df <- peaks_df %>% mutate(.peak_idx = row_number())

  gr <- GRanges(
    seqnames = peaks_df$Summit_Chr,
    ranges   = IRanges(start = peaks_df$Summit_Start, end = peaks_df$Summit_End)
  )
  names(gr) <- peaks_df$.peak_idx

  anno <- annotatePeak(gr, tssRegion = TSS_REGION, TxDb = TXDB,
                       annoDb = "org.Mm.eg.db", verbose = FALSE)
  anno_df <- as.data.frame(anno) %>%
    mutate(.peak_idx = as.integer(names(gr)[seq_len(nrow(.))])) %>%
    select(.peak_idx, annotation, SYMBOL) %>%
    # Keep first mapping per peak if 1:many
    distinct(.peak_idx, .keep_all = TRUE)

  peaks_df <- peaks_df %>%
    left_join(anno_df, by = ".peak_idx") %>%
    mutate(
      is_promoter    = grepl("Promoter", annotation),
      chipseeker_gene = SYMBOL,
      is_promoter     = replace_na(is_promoter, FALSE)
    ) %>%
    select(-.peak_idx, -annotation, -SYMBOL)

  n_prom <- sum(peaks_df$is_promoter, na.rm = TRUE)
  cat(sprintf("  %d promoter / %d distal peaks\n", n_prom, nrow(peaks_df) - n_prom))

  peaks_df
}

link_promoter_peaks_to_genes <- function(peaks_df) {
  peaks_df %>%
    filter(is_promoter) %>%
    filter(!is.na(chipseeker_gene)) %>%
    mutate(assigned_gene = chipseeker_gene)
}

link_distal_peaks_to_genes <- function(peaks_df) {
  distal <- peaks_df %>%
    filter(!is_promoter)

  if (nrow(distal) == 0) return(distal %>% mutate(assigned_gene = character(0)))

  assigned_rows <- list()

  if (USE_ABC) {
    stopifnot(file.exists(ABC_PAIRS_FILE))
    cat("  Linking distal peaks via ABC enhancer-gene pairs...\n")

    abc_pairs <- read_tsv(ABC_PAIRS_FILE, show_col_types = FALSE) %>%
      select(chr, start, end, TargetGene) %>%
      distinct()

    abc_gr <- GRanges(
      seqnames = abc_pairs$chr,
      ranges   = IRanges(start = abc_pairs$start, end = abc_pairs$end)
    )
    abc_gr$TargetGene <- abc_pairs$TargetGene

    peak_gr <- GRanges(
      seqnames = distal$Peak_Chr,
      ranges   = IRanges(start = distal$Peak_Start, end = distal$Peak_End)
    )

    hits <- findOverlaps(peak_gr, abc_gr)

    if (length(hits) > 0) {
      abc_assigned <- tibble(
        peak_idx   = queryHits(hits),
        assigned_gene = abc_gr$TargetGene[subjectHits(hits)]
      ) %>%
        left_join(
          distal %>% mutate(peak_idx = row_number()),
          by = "peak_idx"
        ) %>%
        select(-peak_idx)

      assigned_rows <- append(assigned_rows, list(abc_assigned))

      abc_peak_idxs <- unique(queryHits(hits))
      cat(sprintf("    ABC linked %d / %d distal peaks\n",
                  length(abc_peak_idxs), nrow(distal)))

      remaining_idx <- setdiff(seq_len(nrow(distal)), abc_peak_idxs)
    } else {
      remaining_idx <- seq_len(nrow(distal))
    }
  } else {
    remaining_idx <- seq_len(nrow(distal))
  }

  if (USE_LOOP_FALLBACK && length(remaining_idx) > 0) {
    stopifnot(file.exists(LOOP_ASSIGN_FILE))
    cat("  Linking remaining distal peaks via loop anchors...\n")

    loops <- read_tsv(LOOP_ASSIGN_FILE, show_col_types = FALSE)

    remaining <- distal[remaining_idx, ]
    remain_gr <- GRanges(
      seqnames = remaining$Peak_Chr,
      ranges   = IRanges(start = remaining$Peak_Start, end = remaining$Peak_End)
    )

    a1_gr <- GRanges(
      seqnames = loops$anchor1_chr,
      ranges   = IRanges(start = loops$anchor1_start, end = loops$anchor1_end)
    )
    a1_gr$gene <- loops$anchor2_gene

    a2_gr <- GRanges(
      seqnames = loops$anchor2_chr,
      ranges   = IRanges(start = loops$anchor2_start, end = loops$anchor2_end)
    )
    a2_gr$gene <- loops$anchor1_gene

    hits1 <- findOverlaps(remain_gr, a1_gr)
    hits2 <- findOverlaps(remain_gr, a2_gr)

    loop_genes <- bind_rows(
      tibble(peak_idx = queryHits(hits1), assigned_gene = a1_gr$gene[subjectHits(hits1)]),
      tibble(peak_idx = queryHits(hits2), assigned_gene = a2_gr$gene[subjectHits(hits2)])
    ) %>%
      filter(!is.na(assigned_gene) & assigned_gene != "0") %>%
      distinct()

    if (nrow(loop_genes) > 0) {
      loop_assigned <- loop_genes %>%
        left_join(
          remaining %>% mutate(peak_idx = row_number()),
          by = "peak_idx"
        ) %>%
        select(-peak_idx)

      assigned_rows <- append(assigned_rows, list(loop_assigned))
      cat(sprintf("    Loop fallback linked %d additional peaks\n",
                  n_distinct(loop_genes$peak_idx)))
    }
  }

  # Nearest-gene fallback: assign remaining distal peaks to ChIPseeker gene.
  # Biologically correct for repressive marks (H3K27me3) that act locally
  # and won't overlap ABC enhancers (built from ATAC + H3K27ac).
  if (length(remaining_idx) > 0) {
    nearest_fallback <- distal[remaining_idx, ] %>%
      filter(!is.na(chipseeker_gene)) %>%
      mutate(assigned_gene = chipseeker_gene)

    if (nrow(nearest_fallback) > 0) {
      assigned_rows <- append(assigned_rows, list(nearest_fallback))
      cat(sprintf("    Nearest-gene fallback linked %d additional peaks\n",
                  nrow(nearest_fallback)))
    }
  }

  if (length(assigned_rows) == 0) {
    return(distal %>% slice(0) %>% mutate(assigned_gene = character(0)))
  }

  bind_rows(assigned_rows)
}

aggregate_peak_fc_per_gene <- function(peaks_with_genes) {
  peaks_with_genes %>%
    filter(!is.na(assigned_gene)) %>%
    group_by(mark, assigned_gene) %>%
    slice_max(abs(Fold), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(mark, gene = assigned_gene, Fold) %>%
    pivot_wider(names_from = mark, values_from = Fold, names_glue = "{mark}_FC")
}

build_modification_fc_table <- function() {
  cat("\n=== Modification-Gene Linkage ===\n")

  all_peaks <- list()
  for (mark_name in names(DIFFBIND_FILES)) {
    all_peaks[[mark_name]] <- load_diffbind_file(mark_name, DIFFBIND_FILES[[mark_name]])
  }

  combined <- bind_rows(all_peaks)
  cat(sprintf("  Total: %d peaks across %d marks\n", nrow(combined), length(all_peaks)))

  combined <- annotate_peaks_promoter_distal(combined)

  promoter_linked <- link_promoter_peaks_to_genes(combined)
  distal_linked   <- link_distal_peaks_to_genes(combined)

  all_linked <- bind_rows(promoter_linked, distal_linked)
  cat(sprintf("  Total linked: %d peak-gene pairs\n", nrow(all_linked)))

  mod_fc <- aggregate_peak_fc_per_gene(all_linked)
  cat(sprintf("  Gene-level FC table: %d genes\n", nrow(mod_fc)))

  for (mark_name in names(DIFFBIND_FILES)) {
    col_name <- paste0(mark_name, "_FC")
    if (col_name %in% names(mod_fc)) {
      n_genes <- sum(!is.na(mod_fc[[col_name]]))
      cat(sprintf("    %s: %d genes with FC\n", mark_name, n_genes))
    }
  }

  mod_fc
}

# ==============================================================================
# 6. BUILD HEATMAP MATRIX (z-scored columns)
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
# 7. RENDER HEATMAP
# ==============================================================================

render_heatmap <- function(mat, annotation_row, annotation_colors,
                           title, base_path, width, height,
                           gaps_col = NULL) {

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
      fontsize_row      = 6,
      fontsize_col      = 10,
      fontsize          = 10,
      annotation_row    = .(annotation_row),
      annotation_colors = .(annotation_colors),
      border_color      = NA,
      cellwidth         = 36,
      cellheight        = 8,
      angle_col         = 45,
      na_col            = "grey90",
      gaps_col          = .(gaps_col)
    )
  )

  save_multiformat_pheatmap(hm_call, base_path, width = width, height = height)
  cat(sprintf("  Saved: %s\n", base_path))
}

# ==============================================================================
# 8. HEATMAP A — TOP 100 BY COMBINED STRUCTURAL SCORE
# ==============================================================================

make_combined_score_heatmap <- function(df, mod_fc) {
  cat("\n=== Heatmap A: Top 100 by Combined Structural Score ===\n")

  top <- df %>%
    filter(baseMean >= BASEMEAN_MIN) %>%
    filter(!is.na(combined_score)) %>%
    arrange(desc(combined_score)) %>%
    slice_head(n = TOP_N) %>%
    mutate(
      loop_logFC_signed = replace_na(mean_loop_logFC, 0),
      boundary_flag     = as.numeric(has_boundary)
    )

  cat(sprintf("  Top gene: %s (combined_score=%.2f, n_layers=%d)\n",
              top$gene[1], top$combined_score[1], top$n_layers[1]))
  cat(sprintf("  Score range: %.2f to %.2f\n",
              min(top$combined_score), max(top$combined_score)))

  # Join modification FCs
  top <- top %>% left_join(mod_fc, by = "gene")

  # Build column lists: structural + available modification columns
  struct_columns <- c("log2FC", "logFC_AxC", "loop_logFC_signed", "boundary_flag")
  struct_labels  <- c("mRNA logFC", "AxC logFC", "Loop logFC", "Boundary")

  mod_col_map <- c(
    H2AK119ub_FC = "H2AK119ub FC",
    H3K27me3_FC  = "H3K27me3 FC",
    H3K27ac_FC   = "H3K27ac FC",
    ATAC_FC      = "ATAC-seq FC"
  )
  available_mod_cols <- intersect(names(mod_col_map), names(top))

  columns    <- c(struct_columns, available_mod_cols)
  col_labels <- c(struct_labels, mod_col_map[available_mod_cols])

  n_mod <- length(available_mod_cols)
  cat(sprintf("  Columns: %d structural + %d modification\n",
              length(struct_columns), n_mod))

  mat <- build_zscore_matrix(top, columns, col_labels)

  annotation_row <- data.frame(
    Layers       = factor(top$n_layers),
    Dysregulated = factor(as.character(top$dysregulated)),
    row.names    = top$gene
  )

  annotation_colors <- list(
    Layers       = LAYER_COLORS,
    Dysregulated = DYSREG_COLORS
  )

  hm_width <- 10 + 2 * n_mod

  render_heatmap(
    mat, annotation_row, annotation_colors,
    title     = "Top 100 Genes by Combined Structural Score",
    base_path = file.path(PLOT_DIR, "combined_score_heatmap"),
    width     = hm_width, height = 22,
    gaps_col  = length(struct_columns)
  )

  write_tsv(top, file.path(TABLE_DIR, "top100_combined_score.tsv"))
  cat(sprintf("  Saved table: %d genes\n", nrow(top)))

  top
}

# ==============================================================================
# 9. HEATMAP B — TOP 100 BY ABC ONLY
# ==============================================================================

make_abc_only_heatmap <- function(df) {
  cat("\n=== Heatmap B: Top 100 by |delta ABC| ===\n")

  top <- df %>%
    filter(baseMean >= BASEMEAN_MIN) %>%
    filter(!is.na(max_delta_unnorm)) %>%
    arrange(desc(abs(max_delta_unnorm))) %>%
    slice_head(n = TOP_N)

  cat(sprintf("  Top gene: %s (|ΔABC|=%.4f, logFC_AxC=%.2f)\n",
              top$gene[1], abs(top$max_delta_unnorm[1]),
              top$logFC_AxC[1]))

  columns <- c("log2FC", "logFC_AxC")
  col_labels <- c("mRNA logFC", "AxC logFC")

  mat <- build_zscore_matrix(top, columns, col_labels)

  n_enh_breaks <- cut(top$n_enhancers,
                      breaks = c(0, 5, 10, 15, Inf),
                      labels = c("1-5", "6-10", "11-15", ">15"))

  annotation_row <- data.frame(
    Dysregulated = factor(as.character(top$dysregulated)),
    Enhancers    = n_enh_breaks,
    row.names    = top$gene
  )

  # Purple palette to avoid conflict with blue-white-red heatmap body
  enh_colors <- c("1-5" = "#EFEDF5", "6-10" = "#BCBDDC",
                   "11-15" = "#807DBA", ">15" = "#4A1486")

  annotation_colors <- list(
    Dysregulated = DYSREG_COLORS,
    Enhancers    = enh_colors
  )

  render_heatmap(
    mat, annotation_row, annotation_colors,
    title     = "Top 100 Genes by |AxC Change|",
    base_path = file.path(PLOT_DIR, "abc_only_heatmap"),
    width     = 8, height = 22
  )

  write_tsv(top, file.path(TABLE_DIR, "top100_abc_only.tsv"))
  cat(sprintf("  Saved table: %d genes\n", nrow(top)))

  top
}

# ==============================================================================
# 10. SUMMARY REPORT
# ==============================================================================

write_summary <- function(df, top_combined, top_abc, mod_fc) {
  sink(file.path(TABLE_DIR, "heatmap_summary.txt"))

  cat("=== Structural Heatmap Summary (Figure 5D) ===\n\n")
  cat(sprintf("Total genes in profile: %d\n", nrow(df)))
  cat(sprintf("After baseMean >= %d filter: %d\n",
              BASEMEAN_MIN, sum(df$baseMean >= BASEMEAN_MIN, na.rm = TRUE)))

  cat("\n--- Heatmap A: Combined Structural Score ---\n")
  cat(sprintf("Ranking: combined_score (z_boundary + z_loop_count + z_loop_strength + z_abc)\n"))
  cat(sprintf("Top %d score range: %.2f to %.2f\n", TOP_N,
              min(top_combined$combined_score), max(top_combined$combined_score)))
  cat(sprintf("Layer distribution in top %d:\n", TOP_N))
  print(table(top_combined$n_layers))
  cat(sprintf("Dysregulated in top %d: %d / %d\n", TOP_N,
              sum(top_combined$dysregulated), nrow(top_combined)))

  overlap_genes <- intersect(top_combined$gene, top_abc$gene)
  cat(sprintf("\nOverlap with ABC-only top %d: %d genes\n", TOP_N, length(overlap_genes)))
  if (length(overlap_genes) > 0) {
    cat(sprintf("  Shared: %s\n", paste(overlap_genes, collapse = ", ")))
  }

  cat("\n--- Heatmap B: ABC-Only ---\n")
  cat(sprintf("Ranking: |max_delta_unnorm|\n"))
  cat(sprintf("Top %d |ΔABC| range: %.4f to %.4f\n", TOP_N,
              min(abs(top_abc$max_delta_unnorm)),
              max(abs(top_abc$max_delta_unnorm))))
  cat(sprintf("Dysregulated in top %d: %d / %d\n", TOP_N,
              sum(top_abc$dysregulated), nrow(top_abc)))

  cat("\n--- Modification-Gene Linkage ---\n")
  cat(sprintf("Linkage mode: ABC=%s, Loop fallback=%s\n", USE_ABC, USE_LOOP_FALLBACK))
  cat(sprintf("Total genes with modification FC: %d\n", nrow(mod_fc)))
  for (mark_name in names(DIFFBIND_FILES)) {
    col_name <- paste0(mark_name, "_FC")
    if (col_name %in% names(mod_fc)) {
      n_genes <- sum(!is.na(mod_fc[[col_name]]))
      n_in_top <- sum(!is.na(top_combined[[col_name]]))
      cat(sprintf("  %s: %d total genes, %d in top %d\n",
                  mark_name, n_genes, n_in_top, TOP_N))
    }
  }

  cat("\n--- Columns ---\n")
  cat("Heatmap A: mRNA logFC, AxC logFC, Loop logFC, Boundary | H2AK119ub FC, H3K27me3 FC, H3K27ac FC, ATAC-seq FC (all z-scored)\n")
  cat("Heatmap B: mRNA logFC, AxC logFC (z-scored)\n")
  cat("logFC_AxC = log2((ABC.Score_KO + 1e-6) / (ABC.Score_WT + 1e-6))\n")
  cat("Modification FC = DiffBind Fold (log2 mut/ctrl), strongest peak per gene\n")
  cat("Peak-gene linkage: promoter peaks -> ChIPseeker gene; distal peaks -> ABC TargetGene\n")

  sink()
  cat("  Saved summary report\n")
}

# ==============================================================================
# 11. MAIN
# ==============================================================================

main <- function() {
  cat("=== Structural Heatmap (Figure 5D) ===\n\n")

  dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  df <- load_data() %>%
    add_logfc_axc()

  cat(sprintf("  logFC_AxC computed for %d genes\n",
              sum(!is.na(df$logFC_AxC))))

  mod_fc <- build_modification_fc_table()

  top_combined <- make_combined_score_heatmap(df, mod_fc)
  top_abc      <- make_abc_only_heatmap(df)

  write_summary(df, top_combined, top_abc, mod_fc)

  cat("\n=== Done ===\n")
  cat(sprintf("Output: %s\n", OUTPUT_DIR))
}

main()
