# scripts/deg_tad_violin.R
#
# Creates SATB2-style violin plots showing log2FC distribution of DEGs
# proximal to differential TAD boundaries (lost vs gained)
#
# Usage:
#   Rscript scripts/deg_tad_violin.R --timepoint late
#   Rscript scripts/deg_tad_violin.R --timepoint early
#   Rscript scripts/deg_tad_violin.R --timepoint both
#
# Input:
#   - TADCompare results: results/{timepoint}/final/tadcompare_final_filtered.tsv
#   - RNA-seq: adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx (late)
#              young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx (early)
#
# Output:
#   - results/visualizations/{timepoint}/deg_violin/deg_tad_violin_{timepoint}.pdf
#   - results/visualizations/{timepoint}/deg_violin/deg_boundary_genes_{timepoint}.tsv
#   - results/visualizations/{timepoint}/deg_violin/deg_tad_statistics_{timepoint}.txt

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(GenomicRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ggpubr)
})

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# Base directory (run from tads/)
BASE_DIR <- getwd()

# Parameters
PROXIMITY_WINDOW <- 10000  # 10kb window around boundary
DEG_PADJ_THRESHOLD <- 0.05  # Significance threshold for DEGs

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    tad_file = file.path(BASE_DIR, "results/late/final/tadcompare_final_filtered.tsv"),
    rna_file = file.path(BASE_DIR, "adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"),
    output_dir = file.path(BASE_DIR, "results/visualizations/late/deg_violin"),
    label = "Late (Adult)"
  ),
  early = list(
    tad_file = file.path(BASE_DIR, "results/early/final/tadcompare_final_filtered.tsv"),
    rna_file = file.path(BASE_DIR, "young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),
    output_dir = file.path(BASE_DIR, "results/visualizations/early/deg_violin"),
    label = "Early (Young)"
  )
)

# Color scheme matching SATB2 Figure G
BOUNDARY_COLORS <- c("lost" = "#4575b4", "gained" = "#e0a730")

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

#' Load TADCompare differential boundaries
#' @param tad_file Path to TADCompare filtered results
#' @return tibble of differential boundaries
load_tad_boundaries <- function(tad_file) {
  if (!file.exists(tad_file)) {
    stop(sprintf("TADCompare file not found: %s", tad_file))
  }

  tad_df <- read_tsv(tad_file, show_col_types = FALSE)

  # Filter for differential boundaries only
  diff_boundaries <- tad_df %>%
    dplyr::filter(Differential == "Differential") %>%
    dplyr::select(chr, Boundary, Gap_Score, Enriched_In, Type) %>%
    dplyr::mutate(
      # Map enrichment to lost/gained terminology
      # Matrix 1 = Control-enriched = LOST in mutant
      # Matrix 2 = Mutant-enriched = GAINED in mutant
      boundary_class = case_when(
        Enriched_In == "Matrix 1" ~ "lost",
        Enriched_In == "Matrix 2" ~ "gained",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(boundary_class))

  cat(sprintf("  Loaded %d differential boundaries\n", nrow(diff_boundaries)))
  cat(sprintf("    - Lost (Control-enriched): %d\n", sum(diff_boundaries$boundary_class == "lost")))
  cat(sprintf("    - Gained (Mutant-enriched): %d\n", sum(diff_boundaries$boundary_class == "gained")))

  return(diff_boundaries)
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

  # Filter for significant DEGs and rename column to reflect actual content
  deg_df <- rna_df %>%
    dplyr::filter(!is.na(padj) & padj < DEG_PADJ_THRESHOLD) %>%
    dplyr::select(gene_symbol = ensembl_gene_id, log2FoldChange, padj, baseMean)

  cat(sprintf("  Loaded %d significant DEGs (padj < %g)\n", nrow(deg_df), DEG_PADJ_THRESHOLD))

  return(deg_df)
}

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Vector of gene symbols
#' @return tibble mapping gene symbols to Entrez IDs
convert_symbols_to_entrez <- function(gene_symbols) {
  # Use AnnotationDbi::select explicitly to avoid dplyr conflict
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

#' Find genes within proximity window of TAD boundaries
#' @param boundaries tibble of differential boundaries
#' @param proximity_window Window size in bp (default 10kb)
#' @return tibble of gene-boundary associations
find_boundary_genes <- function(boundaries, proximity_window = PROXIMITY_WINDOW) {
  # Create GRanges for boundaries with flanking regions
  boundary_gr <- GRanges(
    seqnames = boundaries$chr,
    ranges = IRanges(
      start = pmax(1, boundaries$Boundary - proximity_window),
      end = boundaries$Boundary + proximity_window
    ),
    boundary_class = boundaries$boundary_class,
    boundary_pos = boundaries$Boundary
  )

  # Get gene TSS positions from TxDb
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)
  tss_gr <- resize(genes_gr, width = 1, fix = "start")

  # Find overlaps between boundary windows and gene TSS
  overlaps <- findOverlaps(boundary_gr, tss_gr, ignore.strand = TRUE)

  if (length(overlaps) == 0) {
    warning("No gene-boundary overlaps found")
    return(tibble(
      entrez_id = character(),
      boundary_class = character(),
      boundary_pos = integer(),
      chr = character()
    ))
  }

  # Build gene-boundary association table
  gene_boundary_df <- tibble(
    entrez_id = names(tss_gr)[subjectHits(overlaps)],
    boundary_class = boundary_gr$boundary_class[queryHits(overlaps)],
    boundary_pos = boundary_gr$boundary_pos[queryHits(overlaps)],
    chr = as.character(seqnames(boundary_gr)[queryHits(overlaps)])
  )

  cat(sprintf("  Found %d gene-boundary associations\n", nrow(gene_boundary_df)))
  cat(sprintf("    - Unique genes: %d\n", length(unique(gene_boundary_df$entrez_id))))

  return(gene_boundary_df)
}

#' Merge gene-boundary associations with RNA-seq data
#' @param gene_boundary_df Gene-boundary associations
#' @param deg_df Significant DEGs with gene symbols
#' @param id_mapping Gene symbol to Entrez ID mapping
#' @return tibble ready for plotting
merge_deg_boundaries <- function(gene_boundary_df, deg_df, id_mapping) {
  # Add Entrez IDs to DEG data
  deg_with_entrez <- deg_df %>%
    dplyr::inner_join(id_mapping, by = c("gene_symbol" = "SYMBOL"))

  # Merge with gene-boundary associations
  plot_data <- gene_boundary_df %>%
    dplyr::inner_join(deg_with_entrez, by = c("entrez_id" = "ENTREZID")) %>%
    dplyr::select(entrez_id, symbol = gene_symbol, boundary_class, log2FoldChange, padj, chr, boundary_pos) %>%
    # For genes near multiple boundaries, take the one with highest |Gap_Score| effect
    # Simplified: just keep distinct genes per boundary class
    dplyr::distinct(entrez_id, boundary_class, .keep_all = TRUE)

  cat(sprintf("  Merged data: %d gene-boundary pairs\n", nrow(plot_data)))
  cat(sprintf("    - Near lost boundaries: %d genes\n", sum(plot_data$boundary_class == "lost")))
  cat(sprintf("    - Near gained boundaries: %d genes\n", sum(plot_data$boundary_class == "gained")))

  return(plot_data)
}

#' Create SATB2-style violin plot
#' @param plot_data Merged DEG-boundary data
#' @param title Plot title
#' @return ggplot object
create_violin_plot <- function(plot_data, title = "DEGs proximal to differential TAD boundaries") {
  # Calculate statistics for annotation
  n_lost <- sum(plot_data$boundary_class == "lost")
  n_gained <- sum(plot_data$boundary_class == "gained")

  # Perform Mann-Whitney U test
  test_result <- wilcox.test(
    log2FoldChange ~ boundary_class,
    data = plot_data,
    alternative = "two.sided"
  )

  # Format p-value
  p_formatted <- if (test_result$p.value < 0.00001) {
    sprintf("p < 0.00001")
  } else if (test_result$p.value < 0.001) {
    sprintf("p = %.2e", test_result$p.value)
  } else {
    sprintf("p = %.4f", test_result$p.value)
  }

  # Create plot matching SATB2 Figure G style
  p <- ggplot(plot_data, aes(x = boundary_class, y = log2FoldChange, fill = boundary_class)) +
    # Violin with density
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    # Inner boxplot
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    # Reference line at 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Colors matching SATB2 figure
    scale_fill_manual(values = BOUNDARY_COLORS) +
    # X-axis labels
    scale_x_discrete(labels = c("lost" = "lost", "gained" = "gained")) +
    # Add p-value annotation
    annotate("text", x = 1.5, y = max(plot_data$log2FoldChange) * 0.95,
             label = p_formatted, size = 4, fontface = "plain") +
    # Add significance bracket
    annotate("segment", x = 1, xend = 2,
             y = max(plot_data$log2FoldChange) * 0.85,
             yend = max(plot_data$log2FoldChange) * 0.85,
             color = "black", linewidth = 0.5) +
    annotate("text", x = 1.5, y = max(plot_data$log2FoldChange) * 0.88,
             label = "***", size = 5) +
    # Labels
    labs(
      title = title,
      x = NULL,
      y = expression(log[2]*"FC cKO vs floxed")
    ) +
    # Theme matching SATB2 style
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )

  return(list(plot = p, test = test_result, n_lost = n_lost, n_gained = n_gained))
}

#' Generate summary statistics
#' @param plot_data Merged DEG-boundary data
#' @param test_result Wilcoxon test result
#' @param timepoint Timepoint label
#' @return Character vector of summary lines
generate_statistics <- function(plot_data, test_result, timepoint) {
  lost_data <- plot_data %>% dplyr::filter(boundary_class == "lost")
  gained_data <- plot_data %>% dplyr::filter(boundary_class == "gained")

  stats_lines <- c(
    "===========================================",
    sprintf("DEG-TAD Boundary Analysis: %s", timepoint),
    "===========================================",
    "",
    sprintf("Date: %s", Sys.time()),
    sprintf("Proximity window: %d kb", PROXIMITY_WINDOW / 1000),
    sprintf("DEG significance threshold: padj < %g", DEG_PADJ_THRESHOLD),
    "",
    "--- GENE COUNTS ---",
    sprintf("Total genes near differential boundaries: %d", nrow(plot_data)),
    sprintf("  Genes near LOST boundaries: %d", nrow(lost_data)),
    sprintf("  Genes near GAINED boundaries: %d", nrow(gained_data)),
    "",
    "--- LOG2FC STATISTICS ---",
    "",
    "LOST boundaries (Control-enriched):",
    sprintf("  Median log2FC: %.3f", median(lost_data$log2FoldChange)),
    sprintf("  Mean log2FC: %.3f", mean(lost_data$log2FoldChange)),
    sprintf("  SD: %.3f", sd(lost_data$log2FoldChange)),
    sprintf("  Range: [%.3f, %.3f]", min(lost_data$log2FoldChange), max(lost_data$log2FoldChange)),
    "",
    "GAINED boundaries (Mutant-enriched):",
    sprintf("  Median log2FC: %.3f", median(gained_data$log2FoldChange)),
    sprintf("  Mean log2FC: %.3f", mean(gained_data$log2FoldChange)),
    sprintf("  SD: %.3f", sd(gained_data$log2FoldChange)),
    sprintf("  Range: [%.3f, %.3f]", min(gained_data$log2FoldChange), max(gained_data$log2FoldChange)),
    "",
    "--- STATISTICAL TEST ---",
    "Mann-Whitney U test (Wilcoxon rank-sum):",
    sprintf("  W statistic: %.0f", test_result$statistic),
    sprintf("  p-value: %.2e", test_result$p.value),
    sprintf("  Significant (p < 0.05): %s", ifelse(test_result$p.value < 0.05, "YES", "NO")),
    "",
    "--- INTERPRETATION ---",
    ifelse(median(lost_data$log2FoldChange) < median(gained_data$log2FoldChange),
           "Genes near LOST boundaries show LOWER expression in mutant (as expected from SATB2)",
           "Genes near LOST boundaries do NOT show lower expression in mutant"),
    "==========================================="
  )

  return(stats_lines)
}

#' Save outputs for a single timepoint
#' @param plot_result List containing plot, test, and counts
#' @param plot_data Merged DEG-boundary data
#' @param output_dir Output directory
#' @param timepoint Timepoint name
save_outputs <- function(plot_result, plot_data, output_dir, timepoint) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Save plot in multiple formats
  pdf_file <- file.path(output_dir, sprintf("deg_tad_violin_%s.pdf", timepoint))
  svg_file <- file.path(output_dir, sprintf("deg_tad_violin_%s.svg", timepoint))

  ggsave(pdf_file, plot_result$plot, width = 5, height = 6, dpi = 300)
  ggsave(svg_file, plot_result$plot, width = 5, height = 6, dpi = 300)

  cat(sprintf("  Saved: %s\n", pdf_file))
  cat(sprintf("  Saved: %s\n", svg_file))

  # Save gene list
  gene_file <- file.path(output_dir, sprintf("deg_boundary_genes_%s.tsv", timepoint))
  write_tsv(plot_data, gene_file)
  cat(sprintf("  Saved: %s\n", gene_file))

  # Save statistics
  stats_file <- file.path(output_dir, sprintf("deg_tad_statistics_%s.txt", timepoint))
  stats_lines <- generate_statistics(plot_data, plot_result$test, timepoint)
  writeLines(stats_lines, stats_file)
  cat(sprintf("  Saved: %s\n", stats_file))
}

# ==============================================================================
# 4. MAIN PROCESSING FUNCTION
# ==============================================================================

#' Process a single timepoint
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

  # Step 1: Load TAD boundaries
  cat("\n[Step 1] Loading TAD boundaries...\n")
  boundaries <- load_tad_boundaries(config$tad_file)

  # Step 2: Load RNA-seq DEGs
  cat("\n[Step 2] Loading RNA-seq DEGs...\n")
  deg_df <- load_rnaseq_degs(config$rna_file)

  # Step 3: Convert gene symbols to Entrez IDs
  cat("\n[Step 3] Converting gene symbols to Entrez IDs...\n")
  id_mapping <- convert_symbols_to_entrez(deg_df$gene_symbol)

  # Step 4: Find genes near boundaries
  cat("\n[Step 4] Finding genes within %d kb of boundaries...\n", PROXIMITY_WINDOW / 1000)
  gene_boundary_df <- find_boundary_genes(boundaries)

  # Step 5: Merge with RNA-seq data
  cat("\n[Step 5] Merging with RNA-seq expression data...\n")
  plot_data <- merge_deg_boundaries(gene_boundary_df, deg_df, id_mapping)

  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    warning(sprintf("No genes found near differential boundaries for %s timepoint", timepoint))
    return(NULL)
  }

  # Step 6: Create violin plot
  cat("\n[Step 6] Creating violin plot...\n")
  plot_result <- create_violin_plot(plot_data,
                                    title = sprintf("DEGs proximal to differential\nTAD boundaries (%s)", config$label))

  # Step 7: Save outputs
  cat("\n[Step 7] Saving outputs...\n")
  save_outputs(plot_result, plot_data, config$output_dir, timepoint)

  cat("\n")
  cat(sprintf("Completed %s timepoint\n", timepoint))
  cat(sprintf("  Mann-Whitney p-value: %.2e\n", plot_result$test$p.value))
  cat(sprintf("  Genes near lost: %d, Genes near gained: %d\n",
              plot_result$n_lost, plot_result$n_gained))

  return(plot_result)
}

# ==============================================================================
# 5. COMMAND LINE INTERFACE
# ==============================================================================

main <- function() {
  # Parse command line arguments
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
  cat("DEG-TAD Boundary Violin Plot Analysis\n")
  cat("====================================================\n")
  cat(sprintf("Base directory: %s\n", BASE_DIR))
  cat(sprintf("Timepoint: %s\n", timepoint_arg))
  cat(sprintf("Proximity window: %d kb\n", PROXIMITY_WINDOW / 1000))
  cat(sprintf("DEG threshold: padj < %g\n", DEG_PADJ_THRESHOLD))

  # Determine which timepoints to process
  if (tolower(timepoint_arg) == "both") {
    timepoints <- c("late", "early")
  } else if (tolower(timepoint_arg) %in% c("late", "early")) {
    timepoints <- tolower(timepoint_arg)
  } else {
    stop(sprintf("Invalid timepoint: %s. Use 'late', 'early', or 'both'", timepoint_arg))
  }

  # Process each timepoint
  results <- list()
  for (tp in timepoints) {
    results[[tp]] <- tryCatch({
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
  for (tp in names(results)) {
    if (!is.null(results[[tp]])) {
      config <- TIMEPOINT_CONFIG[[tp]]
      cat(sprintf("\n%s timepoint outputs:\n", toupper(tp)))
      cat(sprintf("  %s/\n", config$output_dir))
    }
  }
}

# Run if called directly
if (!interactive()) {
  main()
}
