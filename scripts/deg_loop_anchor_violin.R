# scripts/deg_loop_anchor_violin.R
#
# Creates SATB2-style violin plots showing log2FC distribution of DEGs
# proximal to differential chromatin loop anchors (lost vs gained)
#
# Usage:
#   Rscript scripts/deg_loop_anchor_violin.R --timepoint late
#   Rscript scripts/deg_loop_anchor_violin.R --timepoint early
#   Rscript scripts/deg_loop_anchor_violin.R --timepoint both
#
# Input:
#   - Characterized loops: {timepoint}_outputs/merged_loops/characterized_loops.tsv
#   - RNA-seq: tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx (late)
#              tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx (early)
#
# Output:
#   - output/deg_loop_violin/{timepoint}/deg_loop_anchor_violin.pdf
#   - output/deg_loop_violin/{timepoint}/deg_anchor_genes.tsv
#   - output/deg_loop_violin/{timepoint}/deg_loop_statistics.txt
#   Plus distance-stratified and chromatin-stratified analyses

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(GenomicRanges)
  library(IRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggpubr)
})

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# Base directory (run from project root)
BASE_DIR <- getwd()

# Parameters
DEG_PADJ_THRESHOLD <- 0.05      # Significance threshold for DEGs
DEG_LFC_THRESHOLD <- 0.3        # Minimum |log2FoldChange| for DEGs

# GREAT-style gene-anchor association parameters
GREAT_UPSTREAM <- 5000          # 5kb upstream of TSS
GREAT_DOWNSTREAM <- 1000        # 1kb downstream of TSS
GREAT_MAX_EXTENSION <- 100000   # 100kb maximum extension

# Distance thresholds for stratification
SHORT_RANGE_THRESHOLD <- 500000  # 500kb boundary for short vs long range

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    loops_file = file.path(BASE_DIR, "25042-late_outputs/merged_loops/characterized_loops.tsv"),
    rna_file = file.path(BASE_DIR, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"),
    output_dir = file.path(BASE_DIR, "output/deg_loop_violin/late"),
    label = "Late (Adult)"
  ),
  early = list(
    loops_file = file.path(BASE_DIR, "250831-early_outputs/merged_loops/characterized_loops.tsv"),
    rna_file = file.path(BASE_DIR, "tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),
    output_dir = file.path(BASE_DIR, "output/deg_loop_violin/early"),
    label = "Early (Young)"
  )
)

# Color scheme matching SATB2 Figure G and existing pipeline
DIRECTION_COLORS <- c("lost" = "#4575b4", "gained" = "#e0a730")

# Chromatin state colors
CHROMATIN_COLORS <- c(
  "Active_Promoter" = "#FF0000",
  "Repressed_Promoter" = "#808080",
  "Bivalent_Promoter" = "#A020F0",
  "Polycomb" = "#0000FF",
  "Active_Enhancer" = "#FFA500",
  "Poised_Enhancer" = "#FFFF00",
  "Other" = "#C0C0C0"
)

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

#' Load characterized loops from differential analysis
#' @param loops_file Path to characterized_loops.tsv
#' @return tibble of differential loops
load_characterized_loops <- function(loops_file) {
  if (!file.exists(loops_file)) {
    stop(sprintf("Loops file not found: %s", loops_file))
  }

  loops_df <- read_tsv(loops_file, show_col_types = FALSE)

  # Ensure required columns exist
  required_cols <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
                     "direction", "logFC", "FDR", "loop_distance", "distance_category",
                     "anchor1_type", "anchor2_type")
  missing_cols <- setdiff(required_cols, colnames(loops_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Filter for differential loops with valid direction
  diff_loops <- loops_df %>%
    dplyr::filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  cat(sprintf("  Loaded %d differential loops\n", nrow(diff_loops)))
  cat(sprintf("    - Lost (down_in_mutant): %d\n", sum(diff_loops$direction == "down_in_mutant")))
  cat(sprintf("    - Gained (up_in_mutant): %d\n", sum(diff_loops$direction == "up_in_mutant")))

  return(diff_loops)
}

#' Load RNA-seq results and filter for significant DEGs
#' @param rna_file Path to RNA-seq Excel file
#' @return tibble of significant DEGs with gene symbols
load_rnaseq_degs <- function(rna_file) {
  if (!file.exists(rna_file)) {
    stop(sprintf("RNA-seq file not found: %s", rna_file))
  }

  rna_df <- read_excel(rna_file)

  # The column is named ensembl_gene_id but actually contains gene symbols
  if (!"ensembl_gene_id" %in% colnames(rna_df)) {
    stop("Expected column 'ensembl_gene_id' not found in RNA-seq file")
  }

  # Filter for significant DEGs: padj < 0.05 AND |log2FC| > 0.3
  deg_df <- rna_df %>%
    dplyr::filter(!is.na(padj) & padj < DEG_PADJ_THRESHOLD) %>%
    dplyr::filter(abs(log2FoldChange) > DEG_LFC_THRESHOLD) %>%
    dplyr::select(gene_symbol = ensembl_gene_id, log2FoldChange, padj, baseMean)

  cat(sprintf("  Loaded %d significant DEGs (padj < %g, |log2FC| > %g)\n",
              nrow(deg_df), DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD))

  return(deg_df)
}

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Vector of gene symbols
#' @return tibble mapping gene symbols to Entrez IDs
convert_symbols_to_entrez <- function(gene_symbols) {
  mapping <- tryCatch({
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(gene_symbols),
      columns = c("ENTREZID"),
      keytype = "SYMBOL"
    )
  }, error = function(e) {
    warning(sprintf("Error converting gene symbols: %s", e$message))
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })

  # Remove NAs and duplicates (keep first mapping)
  mapping <- mapping %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)

  cat(sprintf("  Converted %d / %d gene symbols to Entrez\n",
              nrow(mapping), length(unique(gene_symbols))))

  return(mapping)
}

#' Calculate GREAT-style regulatory domains for all genes
#' @return tibble with gene regulatory domains
calculate_great_domains <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)

  # Get TSS and strand information
  tss_pos <- ifelse(as.character(strand(genes_gr)) == "-",
                    end(genes_gr),
                    start(genes_gr))
  gene_strand <- as.character(strand(genes_gr))
  gene_chr <- as.character(seqnames(genes_gr))
  gene_ids <- names(genes_gr)

  gene_info <- data.frame(
    entrez_id = gene_ids,
    chr = gene_chr,
    tss = tss_pos,
    strand = gene_strand,
    stringsAsFactors = FALSE
  )

  # Calculate GREAT-style regulatory domains
  gene_info <- gene_info %>%
    dplyr::arrange(chr, tss) %>%
    dplyr::mutate(
      # Basal domain boundaries (strand-aware)
      basal_start = ifelse(strand == "+", tss - GREAT_UPSTREAM, tss - GREAT_DOWNSTREAM),
      basal_end = ifelse(strand == "+", tss + GREAT_DOWNSTREAM, tss + GREAT_UPSTREAM),
      # Maximum possible extension
      max_start = tss - GREAT_MAX_EXTENSION,
      max_end = tss + GREAT_MAX_EXTENSION
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      # Get neighboring genes' basal domains
      prev_basal_end = dplyr::lag(basal_end, default = -Inf),
      next_basal_start = dplyr::lead(basal_start, default = Inf),
      # Extended domain: extend up to max, but stop at neighbor's basal domain
      reg_start = pmax(max_start, prev_basal_end, 1),
      reg_end = pmin(max_end, next_basal_start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      reg_start = ifelse(reg_end < reg_start, basal_start, reg_start),
      reg_end = ifelse(reg_end < reg_start, basal_end, reg_end),
      reg_end = pmax(reg_end, reg_start)
    ) %>%
    dplyr::select(entrez_id, chr, tss, strand, reg_start, reg_end)

  return(gene_info)
}

#' Extract loop anchors as GRanges
#' @param loops_df Differential loops tibble
#' @return GRanges with anchor coordinates and metadata
extract_loop_anchors <- function(loops_df) {
  # Extract anchor 1
  anchor1 <- loops_df %>%
    dplyr::select(
      loop_id, chr = chr1, start = start1, end = end1,
      direction, logFC, FDR, loop_distance, distance_category,
      anchor_type = anchor1_type
    ) %>%
    dplyr::mutate(anchor_num = 1)

  # Extract anchor 2
  anchor2 <- loops_df %>%
    dplyr::select(
      loop_id, chr = chr2, start = start2, end = end2,
      direction, logFC, FDR, loop_distance, distance_category,
      anchor_type = anchor2_type
    ) %>%
    dplyr::mutate(anchor_num = 2)

  # Combine anchors
  all_anchors <- bind_rows(anchor1, anchor2) %>%
    dplyr::mutate(
      anchor_id = paste0(loop_id, "_anchor", anchor_num),
      # Map direction to lost/gained terminology
      boundary_class = case_when(
        direction == "down_in_mutant" ~ "lost",
        direction == "up_in_mutant" ~ "gained",
        TRUE ~ NA_character_
      ),
      # Distance classification
      is_long_range = loop_distance > SHORT_RANGE_THRESHOLD,
      distance_class = ifelse(is_long_range, "long_range", "short_range")
    ) %>%
    dplyr::filter(!is.na(boundary_class))

  cat(sprintf("  Extracted %d anchors from %d loops\n",
              nrow(all_anchors), nrow(loops_df)))

  return(all_anchors)
}

#' Associate anchors with genes using GREAT-style regulatory domains
#' @param anchors_df Anchors tibble with genomic coordinates
#' @param gene_domains Gene regulatory domains tibble
#' @return tibble of anchor-gene associations
associate_anchors_with_genes <- function(anchors_df, gene_domains) {
  # Create GRanges for gene regulatory domains
  gene_domains_gr <- GRanges(
    seqnames = gene_domains$chr,
    ranges = IRanges(start = gene_domains$reg_start, end = gene_domains$reg_end),
    entrez_id = gene_domains$entrez_id
  )

  # Create GRanges for anchors (use midpoint for association)
  anchor_gr <- GRanges(
    seqnames = anchors_df$chr,
    ranges = IRanges(start = anchors_df$start, end = anchors_df$end),
    anchor_id = anchors_df$anchor_id,
    boundary_class = anchors_df$boundary_class,
    anchor_type = anchors_df$anchor_type,
    loop_distance = anchors_df$loop_distance,
    distance_class = anchors_df$distance_class
  )

  # Find anchors that fall within gene regulatory domains
  overlaps <- findOverlaps(anchor_gr, gene_domains_gr, ignore.strand = TRUE)

  if (length(overlaps) == 0) {
    warning("No gene-anchor overlaps found")
    return(tibble(
      entrez_id = character(),
      anchor_id = character(),
      boundary_class = character(),
      anchor_type = character(),
      loop_distance = numeric(),
      distance_class = character(),
      chr = character()
    ))
  }

  # Build anchor-gene association table
  anchor_gene_df <- tibble(
    entrez_id = gene_domains_gr$entrez_id[subjectHits(overlaps)],
    anchor_id = anchor_gr$anchor_id[queryHits(overlaps)],
    boundary_class = anchor_gr$boundary_class[queryHits(overlaps)],
    anchor_type = anchor_gr$anchor_type[queryHits(overlaps)],
    loop_distance = anchor_gr$loop_distance[queryHits(overlaps)],
    distance_class = anchor_gr$distance_class[queryHits(overlaps)],
    chr = as.character(seqnames(anchor_gr)[queryHits(overlaps)])
  )

  cat(sprintf("  Found %d gene-anchor associations (GREAT-style)\n", nrow(anchor_gene_df)))
  cat(sprintf("    - Unique genes: %d\n", length(unique(anchor_gene_df$entrez_id))))
  cat(sprintf("    - Parameters: %dkb upstream, %dkb downstream, %dkb max extension\n",
              GREAT_UPSTREAM/1000, GREAT_DOWNSTREAM/1000, GREAT_MAX_EXTENSION/1000))

  return(anchor_gene_df)
}

#' Merge anchor-gene associations with RNA-seq data
#' @param anchor_gene_df Anchor-gene associations
#' @param deg_df Significant DEGs with gene symbols
#' @param id_mapping Gene symbol to Entrez ID mapping
#' @return tibble ready for plotting
merge_deg_anchors <- function(anchor_gene_df, deg_df, id_mapping) {
  # Add Entrez IDs to DEG data
  deg_with_entrez <- deg_df %>%
    dplyr::inner_join(id_mapping, by = c("gene_symbol" = "SYMBOL"))

  # Merge with anchor-gene associations
  plot_data <- anchor_gene_df %>%
    dplyr::inner_join(deg_with_entrez, by = c("entrez_id" = "ENTREZID")) %>%
    dplyr::select(entrez_id, symbol = gene_symbol, boundary_class, log2FoldChange, padj,
                  chr, anchor_id, anchor_type, loop_distance, distance_class) %>%
    # For genes near multiple anchors, keep distinct genes per boundary class
    dplyr::distinct(entrez_id, boundary_class, .keep_all = TRUE)

  cat(sprintf("  Merged data: %d gene-anchor pairs\n", nrow(plot_data)))
  cat(sprintf("    - Near lost anchors: %d genes\n", sum(plot_data$boundary_class == "lost")))
  cat(sprintf("    - Near gained anchors: %d genes\n", sum(plot_data$boundary_class == "gained")))

  return(plot_data)
}

#' Create SATB2-style violin plot
#' @param plot_data Merged DEG-anchor data
#' @param title Plot title
#' @param group_col Column to use for grouping (default: boundary_class)
#' @param colors Named vector of colors for groups
#' @return list with plot object and test result
create_violin_plot <- function(plot_data, title = "DEGs proximal to differential loop anchors",
                               group_col = "boundary_class", colors = DIRECTION_COLORS) {

  # Ensure group column exists and has valid data
  if (!group_col %in% colnames(plot_data)) {
    stop(sprintf("Group column '%s' not found in data", group_col))
  }

  plot_data <- plot_data %>%
    dplyr::filter(!is.na(.data[[group_col]]))

  groups <- unique(plot_data[[group_col]])

  # Handle case with only one group or no data
  if (length(groups) < 2) {
    warning(sprintf("Need at least 2 groups for comparison, got: %s", paste(groups, collapse = ", ")))
    return(list(
      plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void(),
      test = NULL,
      n_per_group = table(plot_data[[group_col]])
    ))
  }

  # Calculate counts per group
  n_per_group <- table(plot_data[[group_col]])

  # Perform Mann-Whitney U test (for 2 groups)
  if (length(groups) == 2) {
    test_result <- wilcox.test(
      as.formula(paste("log2FoldChange ~", group_col)),
      data = plot_data,
      alternative = "two.sided"
    )
    p_value <- test_result$p.value
  } else {
    # Kruskal-Wallis for multiple groups
    test_result <- kruskal.test(
      as.formula(paste("log2FoldChange ~", group_col)),
      data = plot_data
    )
    p_value <- test_result$p.value
  }

  # Format p-value
  p_formatted <- if (p_value < 0.00001) {
    sprintf("p < 0.00001")
  } else if (p_value < 0.001) {
    sprintf("p = %.2e", p_value)
  } else {
    sprintf("p = %.4f", p_value)
  }

  # Determine significance stars
  sig_stars <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )

  # Calculate y positions for annotation
  y_max <- max(plot_data$log2FoldChange, na.rm = TRUE)
  y_min <- min(plot_data$log2FoldChange, na.rm = TRUE)
  y_range <- y_max - y_min

  # Create plot
  p <- ggplot(plot_data, aes(x = .data[[group_col]], y = log2FoldChange, fill = .data[[group_col]])) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_fill_manual(values = colors) +
    annotate("text", x = (length(groups) + 1) / 2, y = y_max + y_range * 0.12,
             label = p_formatted, size = 4, fontface = "plain") +
    labs(
      title = title,
      x = NULL,
      y = expression(log[2]*"FC BAP1-KO vs WT")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )

  # Add significance bracket for 2-group comparison
  if (length(groups) == 2) {
    p <- p +
      annotate("segment", x = 1, xend = 2,
               y = y_max + y_range * 0.02,
               yend = y_max + y_range * 0.02,
               color = "black", linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_max + y_range * 0.05,
               label = sig_stars, size = 5)
  }

  return(list(plot = p, test = test_result, n_per_group = n_per_group, p_value = p_value))
}

#' Generate summary statistics text
#' @param plot_data Merged DEG-anchor data
#' @param test_result Statistical test result
#' @param timepoint Timepoint label
#' @param analysis_type Type of analysis performed
#' @return Character vector of summary lines
generate_statistics <- function(plot_data, test_result, timepoint, analysis_type = "basic") {
  lost_data <- plot_data %>% dplyr::filter(boundary_class == "lost")
  gained_data <- plot_data %>% dplyr::filter(boundary_class == "gained")

  stats_lines <- c(
    "===========================================",
    sprintf("DEG-Loop Anchor Analysis: %s", timepoint),
    sprintf("Analysis Type: %s", analysis_type),
    "===========================================",
    "",
    sprintf("Date: %s", Sys.time()),
    "",
    "--- METHODOLOGY ---",
    "Gene-anchor association: GREAT-style regulatory domains",
    sprintf("  Basal domain: %dkb upstream, %dkb downstream of TSS", GREAT_UPSTREAM/1000, GREAT_DOWNSTREAM/1000),
    sprintf("  Max extension: %dkb (stops at neighboring gene's basal domain)", GREAT_MAX_EXTENSION/1000),
    sprintf("DEG thresholds: padj < %g AND |log2FC| > %g", DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD),
    "",
    "--- DIRECTION INTERPRETATION ---",
    "  'lost' = down_in_mutant = loop weaker in BAP1-KO",
    "  'gained' = up_in_mutant = loop stronger in BAP1-KO",
    "",
    "--- GENE COUNTS ---",
    sprintf("Total DEGs near differential anchors: %d", nrow(plot_data)),
    sprintf("  DEGs near LOST anchors: %d", nrow(lost_data)),
    sprintf("  DEGs near GAINED anchors: %d", nrow(gained_data)),
    ""
  )

  # Add statistics for each group if they exist
  if (nrow(lost_data) > 0) {
    stats_lines <- c(stats_lines,
      "--- LOG2FC STATISTICS: LOST ANCHORS ---",
      sprintf("  Median log2FC: %.3f", median(lost_data$log2FoldChange)),
      sprintf("  Mean log2FC: %.3f", mean(lost_data$log2FoldChange)),
      sprintf("  SD: %.3f", sd(lost_data$log2FoldChange)),
      sprintf("  Range: [%.3f, %.3f]", min(lost_data$log2FoldChange), max(lost_data$log2FoldChange)),
      ""
    )
  }

  if (nrow(gained_data) > 0) {
    stats_lines <- c(stats_lines,
      "--- LOG2FC STATISTICS: GAINED ANCHORS ---",
      sprintf("  Median log2FC: %.3f", median(gained_data$log2FoldChange)),
      sprintf("  Mean log2FC: %.3f", mean(gained_data$log2FoldChange)),
      sprintf("  SD: %.3f", sd(gained_data$log2FoldChange)),
      sprintf("  Range: [%.3f, %.3f]", min(gained_data$log2FoldChange), max(gained_data$log2FoldChange)),
      ""
    )
  }

  # Add test result if available
  if (!is.null(test_result)) {
    stats_lines <- c(stats_lines,
      "--- STATISTICAL TEST ---",
      "Mann-Whitney U test (Wilcoxon rank-sum):",
      sprintf("  W statistic: %.0f", test_result$statistic),
      sprintf("  p-value: %.2e", test_result$p.value),
      sprintf("  Significant (p < 0.05): %s", ifelse(test_result$p.value < 0.05, "YES", "NO")),
      ""
    )
  }

  # Add interpretation
  if (nrow(lost_data) > 0 && nrow(gained_data) > 0) {
    median_diff <- median(lost_data$log2FoldChange) - median(gained_data$log2FoldChange)
    stats_lines <- c(stats_lines,
      "--- INTERPRETATION ---",
      sprintf("Median difference (lost - gained): %.3f", median_diff),
      ifelse(median(lost_data$log2FoldChange) < median(gained_data$log2FoldChange),
             "Genes near LOST anchors show LOWER expression in mutant (as expected from SATB2)",
             "Genes near LOST anchors do NOT show lower expression in mutant"),
      "==========================================="
    )
  }

  return(stats_lines)
}

#' Save outputs in multiple formats
#' @param plot_result List containing plot, test, and counts
#' @param plot_data Merged DEG-anchor data
#' @param output_dir Output directory
#' @param base_name Base filename (without extension)
#' @param timepoint Timepoint name
#' @param analysis_type Analysis type for statistics file
save_outputs <- function(plot_result, plot_data, output_dir, base_name, timepoint, analysis_type = "basic") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Save plot in multiple formats
  pdf_file <- file.path(output_dir, sprintf("%s.pdf", base_name))
  svg_file <- file.path(output_dir, sprintf("%s.svg", base_name))
  jpg_file <- file.path(output_dir, sprintf("%s.jpg", base_name))

  ggsave(pdf_file, plot_result$plot, width = 5, height = 6, dpi = 300)
  ggsave(svg_file, plot_result$plot, width = 5, height = 6, dpi = 300)
  ggsave(jpg_file, plot_result$plot, width = 5, height = 6, dpi = 300)

  cat(sprintf("  Saved: %s\n", pdf_file))

  # Save gene list
  gene_file <- file.path(output_dir, sprintf("%s_genes.tsv", base_name))
  write_tsv(plot_data, gene_file)
  cat(sprintf("  Saved: %s\n", gene_file))

  # Save statistics
  stats_file <- file.path(output_dir, sprintf("%s_statistics.txt", base_name))
  stats_lines <- generate_statistics(plot_data, plot_result$test, timepoint, analysis_type)
  writeLines(stats_lines, stats_file)
  cat(sprintf("  Saved: %s\n", stats_file))
}

# ==============================================================================
# 4. MAIN ANALYSIS FUNCTIONS
# ==============================================================================

#' Task 2b: Basic violin plot - DEGs near lost vs gained anchors
#' @param plot_data Merged DEG-anchor data
#' @param output_dir Output directory
#' @param timepoint Timepoint label
run_basic_analysis <- function(plot_data, output_dir, timepoint) {
  cat("\n=== Task 2b: Basic Lost vs Gained Analysis ===\n")

  plot_result <- create_violin_plot(
    plot_data,
    title = sprintf("DEGs proximal to differential\nloop anchors (%s)", timepoint),
    group_col = "boundary_class",
    colors = DIRECTION_COLORS
  )

  save_outputs(plot_result, plot_data, output_dir, "deg_loop_anchor_violin", timepoint, "basic")

  cat(sprintf("  Mann-Whitney p-value: %.2e\n", plot_result$p_value))

  return(plot_result)
}

#' Task 2c-i: Focused comparison - Long-range lost vs Short-range gained
#' @param plot_data Merged DEG-anchor data with distance info
#' @param output_dir Output directory
#' @param timepoint Timepoint label
run_distance_focused_analysis <- function(plot_data, output_dir, timepoint) {
  cat("\n=== Task 2c-i: Long-Range Lost vs Short-Range Gained ===\n")

  # Create focused comparison categories
  focused_data <- plot_data %>%
    dplyr::mutate(
      focused_category = case_when(
        boundary_class == "lost" & distance_class == "long_range" ~ "Long-range Lost",
        boundary_class == "gained" & distance_class == "short_range" ~ "Short-range Gained",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(focused_category))

  cat(sprintf("  Long-range lost (>500kb): %d genes\n", sum(focused_data$focused_category == "Long-range Lost")))
  cat(sprintf("  Short-range gained (<=500kb): %d genes\n", sum(focused_data$focused_category == "Short-range Gained")))

  if (nrow(focused_data) < 5) {
    cat("  Insufficient data for focused analysis, skipping...\n")
    return(NULL)
  }

  focused_colors <- c("Long-range Lost" = "#4575b4", "Short-range Gained" = "#e0a730")

  plot_result <- create_violin_plot(
    focused_data,
    title = sprintf("Long-range Lost vs Short-range Gained\n(%s)", timepoint),
    group_col = "focused_category",
    colors = focused_colors
  )

  save_outputs(plot_result, focused_data, output_dir,
               "deg_loop_violin_longrange_lost_vs_shortrange_gained", timepoint, "distance_focused")

  return(plot_result)
}

#' Task 2c-ii: Full 2x2 stratification by distance and direction
#' @param plot_data Merged DEG-anchor data with distance info
#' @param output_dir Output directory
#' @param timepoint Timepoint label
run_distance_stratified_analysis <- function(plot_data, output_dir, timepoint) {
  cat("\n=== Task 2c-ii: Full Distance Stratification (2x2) ===\n")

  # Create 4-category stratification
  stratified_data <- plot_data %>%
    dplyr::mutate(
      stratified_category = case_when(
        boundary_class == "lost" & distance_class == "short_range" ~ "Lost Short-range",
        boundary_class == "lost" & distance_class == "long_range" ~ "Lost Long-range",
        boundary_class == "gained" & distance_class == "short_range" ~ "Gained Short-range",
        boundary_class == "gained" & distance_class == "long_range" ~ "Gained Long-range",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(stratified_category))

  # Count per category
  category_counts <- table(stratified_data$stratified_category)
  cat("  Category counts:\n")
  for (cat_name in names(category_counts)) {
    cat(sprintf("    - %s: %d genes\n", cat_name, category_counts[cat_name]))
  }

  if (nrow(stratified_data) < 10 || length(unique(stratified_data$stratified_category)) < 2) {
    cat("  Insufficient data for stratified analysis, skipping...\n")
    return(NULL)
  }

  # Define colors for 4 categories
  stratified_colors <- c(
    "Lost Short-range" = "#6baed6",    # Light blue
    "Lost Long-range" = "#08519c",     # Dark blue
    "Gained Short-range" = "#fdae6b",  # Light orange
    "Gained Long-range" = "#d94801"    # Dark orange
  )

  # Order categories logically
  stratified_data$stratified_category <- factor(
    stratified_data$stratified_category,
    levels = c("Lost Short-range", "Lost Long-range", "Gained Short-range", "Gained Long-range")
  )

  # Perform Kruskal-Wallis test
  kw_test <- kruskal.test(log2FoldChange ~ stratified_category, data = stratified_data)

  # Calculate pairwise Wilcoxon tests
  pairwise_tests <- pairwise.wilcox.test(
    stratified_data$log2FoldChange,
    stratified_data$stratified_category,
    p.adjust.method = "BH"
  )

  # Format p-value
  p_formatted <- if (kw_test$p.value < 0.00001) {
    "p < 0.00001"
  } else if (kw_test$p.value < 0.001) {
    sprintf("p = %.2e", kw_test$p.value)
  } else {
    sprintf("p = %.4f", kw_test$p.value)
  }

  # Calculate y positions
  y_max <- max(stratified_data$log2FoldChange, na.rm = TRUE)
  y_range <- y_max - min(stratified_data$log2FoldChange, na.rm = TRUE)

  # Create faceted violin plot
  p <- ggplot(stratified_data, aes(x = stratified_category, y = log2FoldChange, fill = stratified_category)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_fill_manual(values = stratified_colors) +
    annotate("text", x = 2.5, y = y_max + y_range * 0.1,
             label = sprintf("Kruskal-Wallis %s", p_formatted), size = 3.5) +
    labs(
      title = sprintf("DEGs by Loop Direction and Distance\n(%s)", timepoint),
      subtitle = sprintf("Distance threshold: %dkb", SHORT_RANGE_THRESHOLD/1000),
      x = NULL,
      y = expression(log[2]*"FC BAP1-KO vs WT")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )

  # Save plot
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, "deg_loop_violin_distance_2x2_stratified.pdf")
  svg_file <- file.path(output_dir, "deg_loop_violin_distance_2x2_stratified.svg")
  jpg_file <- file.path(output_dir, "deg_loop_violin_distance_2x2_stratified.jpg")

  ggsave(pdf_file, p, width = 8, height = 6, dpi = 300)
  ggsave(svg_file, p, width = 8, height = 6, dpi = 300)
  ggsave(jpg_file, p, width = 8, height = 6, dpi = 300)
  cat(sprintf("  Saved: %s\n", pdf_file))

  # Save gene list
  gene_file <- file.path(output_dir, "deg_loop_distance_stratified_genes.tsv")
  write_tsv(stratified_data, gene_file)
  cat(sprintf("  Saved: %s\n", gene_file))

  # Save statistics with pairwise comparisons
  stats_file <- file.path(output_dir, "deg_loop_violin_distance_statistics.txt")
  stats_lines <- c(
    "===========================================",
    sprintf("Distance-Stratified DEG Analysis: %s", timepoint),
    "===========================================",
    "",
    sprintf("Date: %s", Sys.time()),
    sprintf("Distance threshold: %d kb", SHORT_RANGE_THRESHOLD/1000),
    "",
    "--- CATEGORY COUNTS ---"
  )
  for (cat_name in levels(stratified_data$stratified_category)) {
    n <- sum(stratified_data$stratified_category == cat_name)
    med <- median(stratified_data$log2FoldChange[stratified_data$stratified_category == cat_name])
    stats_lines <- c(stats_lines, sprintf("  %s: n=%d, median log2FC=%.3f", cat_name, n, med))
  }
  stats_lines <- c(stats_lines,
    "",
    "--- KRUSKAL-WALLIS TEST ---",
    sprintf("  Chi-squared: %.2f", kw_test$statistic),
    sprintf("  df: %d", kw_test$parameter),
    sprintf("  p-value: %.2e", kw_test$p.value),
    "",
    "--- PAIRWISE WILCOXON TESTS (BH-adjusted) ---"
  )

  # Add pairwise comparison results
  pw_mat <- pairwise_tests$p.value
  for (i in 1:nrow(pw_mat)) {
    for (j in 1:ncol(pw_mat)) {
      if (!is.na(pw_mat[i, j])) {
        stats_lines <- c(stats_lines,
                         sprintf("  %s vs %s: p = %.4e", rownames(pw_mat)[i], colnames(pw_mat)[j], pw_mat[i, j]))
      }
    }
  }

  writeLines(stats_lines, stats_file)
  cat(sprintf("  Saved: %s\n", stats_file))

  return(list(plot = p, test = kw_test, pairwise = pairwise_tests))
}

#' Task 2d: Polycomb-focused analysis
#' @param plot_data Merged DEG-anchor data with chromatin state info
#' @param output_dir Output directory
#' @param timepoint Timepoint label
run_polycomb_analysis <- function(plot_data, output_dir, timepoint) {
  cat("\n=== Task 2d: Polycomb-Focused Analysis ===\n")

  # Filter for Polycomb-marked anchors
  polycomb_data <- plot_data %>%
    dplyr::filter(anchor_type == "Polycomb")

  cat(sprintf("  Total Polycomb-anchored DEGs: %d\n", nrow(polycomb_data)))
  cat(sprintf("    - Near lost Polycomb anchors: %d\n", sum(polycomb_data$boundary_class == "lost")))
  cat(sprintf("    - Near gained Polycomb anchors: %d\n", sum(polycomb_data$boundary_class == "gained")))

  if (nrow(polycomb_data) < 5 || length(unique(polycomb_data$boundary_class)) < 2) {
    cat("  Insufficient Polycomb data for analysis, skipping...\n")
    return(NULL)
  }

  plot_result <- create_violin_plot(
    polycomb_data,
    title = sprintf("DEGs near Polycomb-marked anchors\n(%s)", timepoint),
    group_col = "boundary_class",
    colors = DIRECTION_COLORS
  )

  save_outputs(plot_result, polycomb_data, output_dir,
               "deg_loop_anchor_violin_polycomb", timepoint, "polycomb")

  return(plot_result)
}

#' Task 2d extended: Analysis by all chromatin states
#' @param plot_data Merged DEG-anchor data with chromatin state info
#' @param output_dir Output directory
#' @param timepoint Timepoint label
run_chromatin_state_analysis <- function(plot_data, output_dir, timepoint) {
  cat("\n=== Task 2d Extended: All Chromatin States ===\n")

  # Count genes per chromatin state
  state_counts <- plot_data %>%
    dplyr::group_by(anchor_type, boundary_class) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = boundary_class, values_from = n, values_fill = 0)

  cat("  Genes per chromatin state:\n")
  print(state_counts)

  # Filter for states with sufficient data (at least 5 genes in each direction)
  valid_states <- plot_data %>%
    dplyr::group_by(anchor_type, boundary_class) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    dplyr::group_by(anchor_type) %>%
    dplyr::filter(all(n >= 3)) %>%  # At least 3 genes in each direction
    dplyr::pull(anchor_type) %>%
    unique()

  if (length(valid_states) < 1) {
    cat("  No chromatin states have sufficient data for comparison, skipping...\n")
    return(NULL)
  }

  cat(sprintf("  States with sufficient data: %s\n", paste(valid_states, collapse = ", ")))

  # Filter data
  chromatin_data <- plot_data %>%
    dplyr::filter(anchor_type %in% valid_states)

  # Calculate median log2FC per state and direction
  state_summary <- chromatin_data %>%
    dplyr::group_by(anchor_type, boundary_class) %>%
    dplyr::summarize(
      n = n(),
      median_log2FC = median(log2FoldChange),
      mean_log2FC = mean(log2FoldChange),
      .groups = "drop"
    )

  # Create faceted plot
  p <- ggplot(chromatin_data, aes(x = boundary_class, y = log2FoldChange, fill = boundary_class)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    facet_wrap(~ anchor_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = DIRECTION_COLORS) +
    labs(
      title = sprintf("DEGs by Anchor Chromatin State (%s)", timepoint),
      x = NULL,
      y = expression(log[2]*"FC BAP1-KO vs WT")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray90", color = NA),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )

  # Save plot
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, "deg_loop_anchor_violin_by_chromatin_state.pdf")
  svg_file <- file.path(output_dir, "deg_loop_anchor_violin_by_chromatin_state.svg")
  jpg_file <- file.path(output_dir, "deg_loop_anchor_violin_by_chromatin_state.jpg")

  # Adjust width based on number of panels
  n_panels <- length(valid_states)
  plot_width <- min(12, max(6, n_panels * 3))

  ggsave(pdf_file, p, width = plot_width, height = 6, dpi = 300)
  ggsave(svg_file, p, width = plot_width, height = 6, dpi = 300)
  ggsave(jpg_file, p, width = plot_width, height = 6, dpi = 300)
  cat(sprintf("  Saved: %s\n", pdf_file))

  # Save summary
  summary_file <- file.path(output_dir, "deg_loop_chromatin_state_summary.tsv")
  write_tsv(state_summary, summary_file)
  cat(sprintf("  Saved: %s\n", summary_file))

  # Perform per-state Wilcoxon tests
  state_tests <- chromatin_data %>%
    dplyr::group_by(anchor_type) %>%
    dplyr::summarize(
      n_lost = sum(boundary_class == "lost"),
      n_gained = sum(boundary_class == "gained"),
      median_lost = median(log2FoldChange[boundary_class == "lost"]),
      median_gained = median(log2FoldChange[boundary_class == "gained"]),
      p_value = tryCatch({
        wilcox.test(log2FoldChange ~ boundary_class, data = pick(everything()))$p.value
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    )

  tests_file <- file.path(output_dir, "deg_loop_chromatin_state_tests.tsv")
  write_tsv(state_tests, tests_file)
  cat(sprintf("  Saved: %s\n", tests_file))

  return(list(plot = p, summary = state_summary, tests = state_tests))
}

# ==============================================================================
# 5. PROCESS TIMEPOINT
# ==============================================================================

#' Process a single timepoint through all analyses
#' @param timepoint Timepoint name ("late" or "early")
process_timepoint <- function(timepoint) {
  cat("\n")
  cat("==================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("==================================================\n")

  config <- TIMEPOINT_CONFIG[[timepoint]]

  if (is.null(config)) {
    stop(sprintf("Unknown timepoint: %s", timepoint))
  }

  # Step 1: Load characterized loops
  cat("\n[Step 1] Loading characterized loops...\n")
  loops_df <- load_characterized_loops(config$loops_file)

  # Step 2: Load RNA-seq DEGs
  cat("\n[Step 2] Loading RNA-seq DEGs...\n")
  deg_df <- load_rnaseq_degs(config$rna_file)

  # Step 3: Convert gene symbols to Entrez IDs
  cat("\n[Step 3] Converting gene symbols to Entrez IDs...\n")
  id_mapping <- convert_symbols_to_entrez(deg_df$gene_symbol)

  # Step 4: Calculate GREAT-style regulatory domains
  cat("\n[Step 4] Calculating GREAT-style regulatory domains...\n")
  gene_domains <- calculate_great_domains()
  cat(sprintf("  Calculated domains for %d genes\n", nrow(gene_domains)))

  # Step 5: Extract loop anchors
  cat("\n[Step 5] Extracting loop anchors...\n")
  anchors_df <- extract_loop_anchors(loops_df)

  # Step 6: Associate anchors with genes
  cat("\n[Step 6] Associating anchors with genes (GREAT-style)...\n")
  anchor_gene_df <- associate_anchors_with_genes(anchors_df, gene_domains)

  # Step 7: Merge with RNA-seq data
  cat("\n[Step 7] Merging with RNA-seq expression data...\n")
  plot_data <- merge_deg_anchors(anchor_gene_df, deg_df, id_mapping)

  # Check if we have data to analyze
  if (nrow(plot_data) == 0) {
    warning(sprintf("No DEGs found near differential loop anchors for %s timepoint", timepoint))
    return(NULL)
  }

  # Save complete merged dataset
  full_output_file <- file.path(config$output_dir, "deg_anchor_genes.tsv")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(plot_data, full_output_file)
  cat(sprintf("\n[Saved] Complete gene list: %s\n", full_output_file))

  results <- list()

  # Task 2b: Basic analysis
  results$basic <- run_basic_analysis(plot_data, config$output_dir, config$label)

  # Task 2c-i: Distance-focused analysis (long-range lost vs short-range gained)
  results$distance_focused <- run_distance_focused_analysis(plot_data, config$output_dir, config$label)

  # Task 2c-ii: Full distance stratification
  results$distance_stratified <- run_distance_stratified_analysis(plot_data, config$output_dir, config$label)

  # Task 2d: Polycomb-focused analysis
  results$polycomb <- run_polycomb_analysis(plot_data, config$output_dir, config$label)

  # Task 2d extended: All chromatin states
  results$chromatin_states <- run_chromatin_state_analysis(plot_data, config$output_dir, config$label)

  cat("\n")
  cat(sprintf("Completed %s timepoint\n", timepoint))
  cat(sprintf("Output directory: %s\n", config$output_dir))

  return(results)
}

# ==============================================================================
# 6. COMMAND LINE INTERFACE
# ==============================================================================

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default to both timepoints
  timepoint_arg <- "both"

  # Parse --timepoint argument
  if (length(args) >= 2) {
    for (i in seq(1, length(args) - 1)) {
      if (args[i] == "--timepoint") {
        timepoint_arg <- args[i + 1]
        break
      }
    }
  } else if (length(args) == 1) {
    timepoint_arg <- args[1]
  }

  cat("====================================================\n")
  cat("DEG-Loop Anchor Violin Plot Analysis\n")
  cat("====================================================\n")
  cat(sprintf("Base directory: %s\n", BASE_DIR))
  cat(sprintf("Timepoint: %s\n", timepoint_arg))
  cat(sprintf("Gene association: GREAT-style (%dkb up, %dkb down, %dkb max)\n",
              GREAT_UPSTREAM/1000, GREAT_DOWNSTREAM/1000, GREAT_MAX_EXTENSION/1000))
  cat(sprintf("DEG threshold: padj < %g AND |log2FC| > %g\n", DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD))
  cat(sprintf("Distance threshold for stratification: %dkb\n", SHORT_RANGE_THRESHOLD/1000))

  # Determine which timepoints to process
  if (tolower(timepoint_arg) == "both") {
    timepoints <- c("late", "early")
  } else if (tolower(timepoint_arg) %in% c("late", "early")) {
    timepoints <- tolower(timepoint_arg)
  } else {
    stop(sprintf("Invalid timepoint: %s. Use 'late', 'early', or 'both'", timepoint_arg))
  }

  # Process each timepoint
  all_results <- list()
  for (tp in timepoints) {
    all_results[[tp]] <- tryCatch({
      process_timepoint(tp)
    }, error = function(e) {
      cat(sprintf("\nERROR processing %s: %s\n", tp, e$message))
      return(NULL)
    })
  }

  cat("\n")
  cat("====================================================\n")
  cat("Analysis Complete\n")
  cat("====================================================\n")

  # Summary
  for (tp in names(all_results)) {
    if (!is.null(all_results[[tp]])) {
      config <- TIMEPOINT_CONFIG[[tp]]
      cat(sprintf("\n%s timepoint outputs:\n", toupper(tp)))
      cat(sprintf("  %s/\n", config$output_dir))

      # List generated files
      files <- list.files(config$output_dir, pattern = "\\.(pdf|tsv|txt)$")
      for (f in files) {
        cat(sprintf("    - %s\n", f))
      }
    }
  }
}

# Run if called directly
if (!interactive()) {
  main()
}
