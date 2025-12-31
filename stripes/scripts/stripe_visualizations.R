#!/usr/bin/env Rscript
# stripes/scripts/stripe_visualizations.R
# Comprehensive Visualization Analysis for Differential Chromatin Stripes
# Author: Zakir Alibhai
# Date: 2024-12-30
#
# Purpose:
#   Generate publication-quality visualizations from differential stripe analysis:
#   1. Medium/high confidence differential stripe lists
#   2. Volcano plots (EnhancedVolcano)
#   3. Length distribution analysis
#   4. ChIP-seq-based anchor annotation (enhancer/promoter/Polycomb)
#   5. Summary statistics and exports
#
# Usage:
#   Rscript scripts/stripe_visualizations.R
#
# Processes both early and late timepoints in a single run

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Stripe Visualization Analysis\n")
cat("========================================\n\n")

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  # Core visualization
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(EnhancedVolcano)

  # Genomic annotation
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)

  # Data wrangling
  library(tidyverse)
  library(yaml)
})

cat("Packages loaded\n\n")

# Load configuration
config_path <- "config/stripe_config.yaml"
if (!file.exists(config_path)) {
  stop("ERROR: Configuration file not found: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Set working directory
base_dir <- config$project$base_dir
if (dir.exists(base_dir)) {
  setwd(base_dir)
  cat(sprintf("Working directory: %s\n", base_dir))
} else {
  cat(sprintf("Using current directory: %s\n", getwd()))
}

# Create output directories
output_base <- file.path(config$outputs$base_dir, "visualizations")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "early"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "late"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "combined"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output directory: %s\n\n", output_base))

# Define timepoints (directory names match these directly)
timepoints <- c("early", "late")

# =============================================================================
# SECTION 1: LOAD & FILTER DATA
# =============================================================================

cat("\n========================================\n")
cat("SECTION 1: Loading & Filtering Data\n")
cat("========================================\n\n")

# Function to load stripe data for a timepoint
load_stripe_data <- function(timepoint_name) {
  # Directory structure: outputs/{timepoint}/
  input_file <- file.path(config$outputs$base_dir, timepoint_name,
                          "04_final_differential_stripes.tsv")

  if (!file.exists(input_file)) {
    warning(sprintf("Stripe file not found for %s: %s", timepoint_name, input_file))
    return(NULL)
  }

  df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df$timepoint <- timepoint_name

  cat(sprintf("  %s: Loaded %d stripes\n", timepoint_name, nrow(df)))
  cat(sprintf("    Lost: %d | Gained: %d | Shared: %d\n",
              sum(df$direction == "lost", na.rm = TRUE),
              sum(df$direction == "gained", na.rm = TRUE),
              sum(df$direction %in% c("unchanged", "strengthened", "weakened"), na.rm = TRUE)))

  return(df)
}

# Load data for both timepoints
stripes_list <- list()
for (tp_name in timepoints) {
  stripes_list[[tp_name]] <- load_stripe_data(tp_name)
}

# Check if data loaded successfully
if (all(sapply(stripes_list, is.null))) {
  stop("ERROR: No stripe data could be loaded. Check file paths.")
}

# Filter to medium/high confidence stripes
filter_medium_high <- function(df) {
  if (is.null(df)) return(NULL)

  filtered <- df %>%
    filter(direction %in% c("lost", "gained") &
             direction_confidence %in% c("high", "medium"))

  cat(sprintf("  Filtered to medium/high confidence: %d stripes\n", nrow(filtered)))
  cat(sprintf("    High: %d | Medium: %d\n",
              sum(filtered$direction_confidence == "high", na.rm = TRUE),
              sum(filtered$direction_confidence == "medium", na.rm = TRUE)))

  return(filtered)
}

# Apply filtering
cat("\nFiltering to medium/high confidence stripes:\n")
medium_high_list <- list()
for (tp_name in names(stripes_list)) {
  cat(sprintf("\n%s timepoint:\n", tp_name))
  medium_high_list[[tp_name]] <- filter_medium_high(stripes_list[[tp_name]])
}

# Export medium/high confidence lists
cat("\nExporting filtered differential lists...\n")
for (tp_name in names(medium_high_list)) {
  if (!is.null(medium_high_list[[tp_name]]) && nrow(medium_high_list[[tp_name]]) > 0) {
    output_file <- file.path(output_base, tp_name,
                             sprintf("%s_medium_high_confidence_stripes.tsv", tp_name))
    write.table(medium_high_list[[tp_name]], output_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d stripes)\n", basename(output_file),
                nrow(medium_high_list[[tp_name]])))
  }
}

cat("\nSection 1 complete\n")

# =============================================================================
# SECTION 2: VOLCANO PLOTS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 2: Volcano Plots\n")
cat("========================================\n\n")

# Function to create publication-quality volcano plot for stripes
create_stripe_volcano <- function(stripes_df, timepoint_name, output_path) {
  if (is.null(stripes_df) || nrow(stripes_df) == 0) {
    cat(sprintf("  Skipping volcano plot for %s (no data)\n", timepoint_name))
    return(NULL)
  }

  cat(sprintf("Creating volcano plot: %s\n", timepoint_name))

  # Prepare data - ensure required columns exist
  df <- stripes_df %>%
    mutate(
      # Handle any missing FDR values
      FDR = ifelse(is.na(FDR) | FDR == 0, 1e-300, FDR),
      # Create display labels
      significance_label = case_when(
        direction_confidence == "high" ~ "High confidence",
        direction_confidence == "medium" ~ "Medium confidence",
        direction_confidence == "low" ~ "Low confidence",
        TRUE ~ "Not significant"
      )
    )

  # Calculate counts for annotations
  n_total <- nrow(df)
  n_lost <- sum(df$direction == "lost", na.rm = TRUE)
  n_gained <- sum(df$direction == "gained", na.rm = TRUE)
  n_medium_high_lost <- sum(df$direction == "lost" &
                              df$direction_confidence %in% c("high", "medium"), na.rm = TRUE)
  n_medium_high_gained <- sum(df$direction == "gained" &
                                df$direction_confidence %in% c("high", "medium"), na.rm = TRUE)

  cat(sprintf("  Total: %d | Lost: %d (M/H: %d) | Gained: %d (M/H: %d)\n",
              n_total, n_lost, n_medium_high_lost, n_gained, n_medium_high_gained))

  # Get data ranges for annotation positioning
  logfc_range <- range(df$logFC, na.rm = TRUE)
  fdr_range <- range(df$FDR[df$FDR > 0], na.rm = TRUE)
  neg_log10_fdr_max <- -log10(min(fdr_range))

  # Custom colors based on direction and confidence
  df <- df %>%
    mutate(
      plot_color = case_when(
        direction == "lost" & direction_confidence == "high" ~ "Lost (High)",
        direction == "lost" & direction_confidence == "medium" ~ "Lost (Medium)",
        direction == "lost" & direction_confidence == "low" ~ "Lost (Low)",
        direction == "gained" & direction_confidence == "high" ~ "Gained (High)",
        direction == "gained" & direction_confidence == "medium" ~ "Gained (Medium)",
        direction == "gained" & direction_confidence == "low" ~ "Gained (Low)",
        TRUE ~ "Unchanged"
      )
    )

  # Define color palette
  color_palette <- c(
    "Lost (High)" = "#00008B",       # Dark blue
    "Lost (Medium)" = "#4169E1",     # Royal blue
    "Lost (Low)" = "#87CEEB",        # Sky blue
    "Gained (High)" = "#8B0000",     # Dark red
    "Gained (Medium)" = "#DC143C",   # Crimson
    "Gained (Low)" = "#FFA07A",      # Light salmon
    "Unchanged" = "#808080"          # Gray
  )

  # Create base volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = plot_color)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = color_palette, name = "Direction") +
    # Threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.10), linetype = "dotted", color = "gray", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    # Labels
    labs(
      title = sprintf("Differential Stripes: %s Timepoint", tools::toTitleCase(timepoint_name)),
      subtitle = sprintf("Total: %d stripes | Lost: %d | Gained: %d", n_total, n_lost, n_gained),
      x = expression(log[2]~"Fold Change (Mutant/Control)"),
      y = expression(-log[10]~"FDR")
    ) +
    # Annotations for counts
    annotate("text", x = logfc_range[1] + 0.1, y = neg_log10_fdr_max * 0.95,
             label = sprintf("Lost: %d\n(M/H: %d)", n_lost, n_medium_high_lost),
             color = "#4169E1", size = 4, fontface = "bold", hjust = 0) +
    annotate("text", x = logfc_range[2] - 0.1, y = neg_log10_fdr_max * 0.95,
             label = sprintf("Gained: %d\n(M/H: %d)", n_gained, n_medium_high_gained),
             color = "#DC143C", size = 4, fontface = "bold", hjust = 1) +
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  # Save plot
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_path)))

  return(p)
}

# Generate volcano plots for each timepoint
volcano_plots <- list()
for (tp_name in names(stripes_list)) {
  if (!is.null(stripes_list[[tp_name]])) {
    output_file <- file.path(output_base, tp_name, sprintf("volcano_%s.pdf", tp_name))
    volcano_plots[[tp_name]] <- create_stripe_volcano(stripes_list[[tp_name]], tp_name, output_file)
  }
}

# Create combined volcano plot (side by side)
if (length(volcano_plots) == 2 && !any(sapply(volcano_plots, is.null))) {
  cat("\nCreating combined volcano plot...\n")

  combined_volcano <- volcano_plots$early + volcano_plots$late +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Differential Chromatin Stripes: BAP1-KO vs Wildtype",
      subtitle = "Comparing early and late timepoints",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )

  combined_output <- file.path(output_base, "combined", "volcano_combined.pdf")
  ggsave(combined_output, combined_volcano, width = 18, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(combined_output)))
}

cat("\nSection 2 complete: Volcano plots generated\n")

# =============================================================================
# SECTION 3: LENGTH DISTRIBUTION ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 3: Length Distribution Analysis\n")
cat("========================================\n\n")

# Function to analyze length distribution for a timepoint
analyze_length_distribution <- function(stripes_df, timepoint_name, output_dir) {
  if (is.null(stripes_df) || nrow(stripes_df) == 0) {
    cat(sprintf("  Skipping length analysis for %s (no data)\n", timepoint_name))
    return(NULL)
  }

  cat(sprintf("Analyzing length distribution: %s\n", timepoint_name))

  # Prepare data
  df <- stripes_df %>%
    mutate(
      stripe_length_kb = stripe_length / 1000,
      anchor_width_kb = anchor_width / 1000,
      direction_label = case_when(
        direction == "lost" ~ "Lost",
        direction == "gained" ~ "Gained",
        TRUE ~ "Shared"
      )
    )

  # Summary statistics
  length_stats <- df %>%
    group_by(direction_label) %>%
    summarise(
      n = n(),
      median_length_kb = median(stripe_length_kb, na.rm = TRUE),
      mean_length_kb = mean(stripe_length_kb, na.rm = TRUE),
      sd_length_kb = sd(stripe_length_kb, na.rm = TRUE),
      min_length_kb = min(stripe_length_kb, na.rm = TRUE),
      max_length_kb = max(stripe_length_kb, na.rm = TRUE),
      median_anchor_kb = median(anchor_width_kb, na.rm = TRUE),
      .groups = "drop"
    )

  cat("  Length statistics by direction:\n")
  print(length_stats)

  # Statistical test: Lost vs Gained
  lost_lengths <- df$stripe_length[df$direction == "lost"]
  gained_lengths <- df$stripe_length[df$direction == "gained"]

  wilcox_result <- NULL
  if (length(lost_lengths) > 0 && length(gained_lengths) > 0) {
    wilcox_result <- wilcox.test(lost_lengths, gained_lengths)
    cat(sprintf("  Wilcoxon test (Lost vs Gained): p = %.3e\n", wilcox_result$p.value))
  }

  # Color palette
  direction_colors <- c("Lost" = "#4575b4", "Gained" = "#d73027", "Shared" = "#808080")

  # Plot 1: Violin plot of stripe lengths
  p_violin <- ggplot(df, aes(x = direction_label, y = stripe_length_kb, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_jitter(alpha = 0.3, width = 0.1, size = 1) +
    scale_fill_manual(values = direction_colors) +
    labs(
      title = sprintf("Stripe Length Distribution: %s", tools::toTitleCase(timepoint_name)),
      subtitle = ifelse(!is.null(wilcox_result),
                        sprintf("Wilcoxon test (Lost vs Gained): p = %.3e", wilcox_result$p.value),
                        ""),
      x = "Direction",
      y = "Stripe Length (kb)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  # Plot 2: Histogram overlay
  p_hist <- ggplot(df, aes(x = stripe_length_kb, fill = direction_label)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = direction_colors) +
    facet_wrap(~direction_label, ncol = 1, scales = "free_y") +
    labs(
      title = sprintf("Stripe Length Histogram: %s", tools::toTitleCase(timepoint_name)),
      x = "Stripe Length (kb)",
      y = "Count",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold")
    )

  # Plot 3: Anchor width distribution
  p_anchor <- ggplot(df, aes(x = direction_label, y = anchor_width_kb, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = direction_colors) +
    labs(
      title = sprintf("Anchor Width Distribution: %s", tools::toTitleCase(timepoint_name)),
      x = "Direction",
      y = "Anchor Width (kb)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  # Combine plots
  combined <- (p_violin | p_hist) / p_anchor +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = sprintf("Length Analysis: %s Timepoint", tools::toTitleCase(timepoint_name)),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  # Save
  output_file <- file.path(output_dir, sprintf("length_distribution_%s.pdf", timepoint_name))
  ggsave(output_file, combined, width = 14, height = 12, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_file)))

  # Save statistics
  stats_file <- file.path(output_dir, sprintf("length_statistics_%s.tsv", timepoint_name))
  write.table(length_stats, stats_file, sep = "\t", quote = FALSE, row.names = FALSE)

  return(list(
    stats = length_stats,
    wilcox = wilcox_result,
    plots = combined
  ))
}

# Analyze length distribution for each timepoint
length_results <- list()
for (tp_name in names(stripes_list)) {
  output_dir <- file.path(output_base, tp_name)
  length_results[[tp_name]] <- analyze_length_distribution(stripes_list[[tp_name]], tp_name, output_dir)
}

# Create combined length comparison plot
if (length(stripes_list) == 2) {
  cat("\nCreating combined length comparison...\n")

  # Combine data from both timepoints
  combined_df <- bind_rows(
    stripes_list$early %>% mutate(timepoint = "Early"),
    stripes_list$late %>% mutate(timepoint = "Late")
  ) %>%
    filter(!is.na(stripe_length)) %>%
    mutate(
      stripe_length_kb = stripe_length / 1000,
      direction_label = case_when(
        direction == "lost" ~ "Lost",
        direction == "gained" ~ "Gained",
        TRUE ~ "Shared"
      )
    )

  direction_colors <- c("Lost" = "#4575b4", "Gained" = "#d73027", "Shared" = "#808080")

  p_comparison <- ggplot(combined_df, aes(x = direction_label, y = stripe_length_kb, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = direction_colors) +
    facet_wrap(~timepoint) +
    labs(
      title = "Stripe Length Comparison: Early vs Late",
      x = "Direction",
      y = "Stripe Length (kb)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold")
    )

  combined_output <- file.path(output_base, "combined", "length_comparison.pdf")
  ggsave(combined_output, p_comparison, width = 12, height = 6, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(combined_output)))
}

cat("\nSection 3 complete: Length distribution analyzed\n")

# =============================================================================
# SECTION 4: ChIP-seq ANCHOR ANNOTATION
# =============================================================================

cat("\n========================================\n")
cat("SECTION 4: ChIP-seq Anchor Annotation\n")
cat("========================================\n\n")

# Load TxDb for TSS annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss <- resize(genes(txdb), width = 1, fix = "start")
cat(sprintf("Loaded %d TSS positions from mm10\n", length(tss)))

# Function to load ChIP-seq peaks
load_peaks <- function(peak_file, peak_type) {
  if (is.null(peak_file) || !file.exists(peak_file)) {
    warning(sprintf("Peak file not found: %s", peak_file))
    return(NULL)
  }

  cat(sprintf("  Loading %s: %s\n", peak_type, basename(peak_file)))

  # Read as table first to determine format
  peak_df <- tryCatch({
    read.table(peak_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(sprintf("Could not read peak file: %s", e$message))
    return(NULL)
  })

  if (is.null(peak_df) || nrow(peak_df) == 0) {
    return(NULL)
  }

  # Create GRanges (assuming BED format: chr, start, end, ...)
  peaks_gr <- GRanges(
    seqnames = peak_df$V1,
    ranges = IRanges(start = peak_df$V2 + 1, end = peak_df$V3)  # BED is 0-based
  )

  cat(sprintf("    Loaded %d peaks\n", length(peaks_gr)))
  return(peaks_gr)
}

# Function to annotate stripe anchors for a timepoint
annotate_stripe_anchors <- function(stripes_df, timepoint_name, config) {
  if (is.null(stripes_df) || nrow(stripes_df) == 0) {
    cat(sprintf("  Skipping annotation for %s (no data)\n", timepoint_name))
    return(stripes_df)
  }

  cat(sprintf("\nAnnotating anchors: %s timepoint\n", timepoint_name))

  # Load timepoint-specific peaks
  if (timepoint_name == "early") {
    h3k27ac_peaks <- load_peaks(config$chipseq_peaks$h3k27ac_early, "H3K27ac")
    h3k27me3_peaks <- load_peaks(config$chipseq_peaks$h3k27me3_early, "H3K27me3")
    h3k4me1_peaks <- NULL
  } else {
    h3k27ac_peaks <- load_peaks(config$chipseq_peaks$h3k27ac_late, "H3K27ac")
    h3k4me1_peaks <- load_peaks(config$chipseq_peaks$h3k4me1_late, "H3K4me1")
    h3k27me3_peaks <- NULL
  }

  # Create GRanges for stripe anchors
  anchor_gr <- GRanges(
    seqnames = stripes_df$chr,
    ranges = IRanges(start = stripes_df$anchor_x1, end = stripes_df$anchor_x2)
  )

  # Compute distance to nearest TSS
  cat("  Computing distance to nearest TSS...\n")
  nearest_tss <- distanceToNearest(anchor_gr, tss)

  stripes_df$distance_to_tss <- NA
  stripes_df$nearest_gene_id <- NA

  if (length(nearest_tss) > 0) {
    stripes_df$distance_to_tss[queryHits(nearest_tss)] <- mcols(nearest_tss)$distance
    stripes_df$nearest_gene_id[queryHits(nearest_tss)] <- names(tss)[subjectHits(nearest_tss)]
  }

  # Convert Entrez IDs to gene symbols
  if (any(!is.na(stripes_df$nearest_gene_id))) {
    gene_map <- tryCatch({
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = unique(na.omit(stripes_df$nearest_gene_id)),
                            columns = "SYMBOL",
                            keytype = "ENTREZID")
    }, error = function(e) {
      data.frame(ENTREZID = character(), SYMBOL = character())
    })

    stripes_df$nearest_gene_symbol <- gene_map$SYMBOL[match(stripes_df$nearest_gene_id, gene_map$ENTREZID)]
  } else {
    stripes_df$nearest_gene_symbol <- NA
  }

  # Compute ChIP-seq overlaps
  cat("  Computing ChIP-seq overlaps...\n")

  stripes_df$h3k27ac_overlap <- FALSE
  stripes_df$h3k27me3_overlap <- FALSE
  stripes_df$h3k4me1_overlap <- FALSE

  if (!is.null(h3k27ac_peaks)) {
    stripes_df$h3k27ac_overlap <- countOverlaps(anchor_gr, h3k27ac_peaks) > 0
    cat(sprintf("    H3K27ac overlaps: %d (%.1f%%)\n",
                sum(stripes_df$h3k27ac_overlap),
                100 * mean(stripes_df$h3k27ac_overlap)))
  }

  if (!is.null(h3k27me3_peaks)) {
    stripes_df$h3k27me3_overlap <- countOverlaps(anchor_gr, h3k27me3_peaks) > 0
    cat(sprintf("    H3K27me3 overlaps: %d (%.1f%%)\n",
                sum(stripes_df$h3k27me3_overlap),
                100 * mean(stripes_df$h3k27me3_overlap)))
  }

  if (!is.null(h3k4me1_peaks)) {
    stripes_df$h3k4me1_overlap <- countOverlaps(anchor_gr, h3k4me1_peaks) > 0
    cat(sprintf("    H3K4me1 overlaps: %d (%.1f%%)\n",
                sum(stripes_df$h3k4me1_overlap),
                100 * mean(stripes_df$h3k4me1_overlap)))
  }

  # Classify anchor type based on ChIP-seq overlaps
  cat("  Classifying anchor types...\n")
  tss_threshold <- 2000  # 2kb

  stripes_df$anchor_type <- case_when(
    # Promoter: H3K27ac+ AND within 2kb of TSS
    stripes_df$h3k27ac_overlap & stripes_df$distance_to_tss <= tss_threshold ~ "Promoter",
    # Active Enhancer: H3K27ac+ AND > 2kb from TSS
    stripes_df$h3k27ac_overlap & stripes_df$distance_to_tss > tss_threshold ~ "Active_Enhancer",
    # Poised Enhancer: H3K4me1+ (no H3K27ac) - late only
    stripes_df$h3k4me1_overlap & !stripes_df$h3k27ac_overlap &
      stripes_df$distance_to_tss > tss_threshold ~ "Poised_Enhancer",
    # Polycomb Domain: H3K27me3+ - early only
    stripes_df$h3k27me3_overlap ~ "Polycomb_Domain",
    # Other: No marks
    TRUE ~ "Other"
  )

  # Print summary
  anchor_summary <- table(stripes_df$anchor_type)
  cat("\n  Anchor type distribution:\n")
  for (atype in names(anchor_summary)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", atype, anchor_summary[atype],
                100 * anchor_summary[atype] / nrow(stripes_df)))
  }

  return(stripes_df)
}

# Annotate stripes for each timepoint
annotated_stripes <- list()
for (tp_name in names(stripes_list)) {
  annotated_stripes[[tp_name]] <- annotate_stripe_anchors(stripes_list[[tp_name]], tp_name, config)
}

# Function to create anchor classification plots
create_anchor_plots <- function(stripes_df, timepoint_name, output_dir) {
  if (is.null(stripes_df) || nrow(stripes_df) == 0 || !"anchor_type" %in% colnames(stripes_df)) {
    return(NULL)
  }

  cat(sprintf("\nCreating anchor classification plots: %s\n", timepoint_name))

  # Prepare data
  df <- stripes_df %>%
    mutate(
      direction_label = case_when(
        direction == "lost" ~ "Lost",
        direction == "gained" ~ "Gained",
        TRUE ~ "Shared"
      )
    )

  # Color palette for anchor types
  anchor_colors <- c(
    "Promoter" = "#e41a1c",
    "Active_Enhancer" = "#ff7f00",
    "Poised_Enhancer" = "#4daf4a",
    "Polycomb_Domain" = "#984ea3",
    "Other" = "#999999"
  )

  direction_colors <- c("Lost" = "#4575b4", "Gained" = "#d73027", "Shared" = "#808080")

  # Plot 1: Stacked bar - Anchor types by direction
  anchor_by_direction <- df %>%
    group_by(direction_label, anchor_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction_label) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()

  p_stacked <- ggplot(anchor_by_direction, aes(x = direction_label, y = percentage, fill = anchor_type)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = anchor_colors, name = "Anchor Type") +
    labs(
      title = sprintf("Anchor Classification by Direction: %s", tools::toTitleCase(timepoint_name)),
      x = "Stripe Direction",
      y = "Percentage (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right",
      axis.title = element_text(size = 12)
    )

  # Plot 2: Distance to TSS by direction
  p_tss_dist <- ggplot(df %>% filter(!is.na(distance_to_tss)),
                        aes(x = direction_label, y = distance_to_tss / 1000, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = direction_colors) +
    scale_y_log10() +
    labs(
      title = sprintf("Distance to Nearest TSS: %s", tools::toTitleCase(timepoint_name)),
      x = "Stripe Direction",
      y = "Distance to TSS (kb, log10)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  # Plot 3: ChIP-seq overlap summary
  chip_summary <- df %>%
    summarise(
      `H3K27ac+` = sum(h3k27ac_overlap, na.rm = TRUE),
      `H3K27me3+` = sum(h3k27me3_overlap, na.rm = TRUE),
      `H3K4me1+` = sum(h3k4me1_overlap, na.rm = TRUE),
      Total = n()
    ) %>%
    pivot_longer(cols = -Total, names_to = "Mark", values_to = "Count") %>%
    mutate(Percentage = 100 * Count / Total)

  p_chip <- ggplot(chip_summary %>% filter(Count > 0),
                    aes(x = Mark, y = Percentage, fill = Mark)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
              vjust = -0.3, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = sprintf("ChIP-seq Mark Overlaps: %s", tools::toTitleCase(timepoint_name)),
      x = "Histone Mark",
      y = "% of Anchors"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    ylim(0, max(chip_summary$Percentage) * 1.3)

  # Combine plots
  combined <- (p_stacked | p_tss_dist) / p_chip +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = sprintf("Anchor Annotation: %s Timepoint", tools::toTitleCase(timepoint_name)),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  # Save
  output_file <- file.path(output_dir, sprintf("anchor_classification_%s.pdf", timepoint_name))
  ggsave(output_file, combined, width = 14, height = 10, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_file)))

  return(combined)
}

# Create anchor plots for each timepoint
for (tp_name in names(annotated_stripes)) {
  output_dir <- file.path(output_base, tp_name)
  create_anchor_plots(annotated_stripes[[tp_name]], tp_name, output_dir)
}

# Export annotated stripe lists
cat("\nExporting annotated stripe lists...\n")
for (tp_name in names(annotated_stripes)) {
  if (!is.null(annotated_stripes[[tp_name]])) {
    output_file <- file.path(output_base, tp_name,
                             sprintf("%s_annotated_stripes.tsv", tp_name))
    write.table(annotated_stripes[[tp_name]], output_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d stripes)\n", basename(output_file),
                nrow(annotated_stripes[[tp_name]])))
  }
}

cat("\nSection 4 complete: ChIP-seq anchor annotation done\n")

# =============================================================================
# SECTION 4B: EXPORT ANNOTATED CONFIDENT STRIPES (TSV + BEDPE)
# =============================================================================

cat("\n========================================\n")
cat("SECTION 4B: Annotated Confident Stripe Export\n")
cat("========================================\n\n")

# Function to create BEDPE from annotated stripes
create_annotated_bedpe <- function(stripes_df, output_path) {
  if (is.null(stripes_df) || nrow(stripes_df) == 0) {
    return(NULL)
  }

  # Create BEDPE format with annotations
  # Standard BEDPE: chr1, x1, x2, chr2, y1, y2, name, score, strand1, strand2
  # Plus annotation columns
  bedpe_df <- data.frame(
    chr1 = stripes_df$chr,
    x1 = stripes_df$anchor_x1,
    x2 = stripes_df$anchor_x2,
    chr2 = stripes_df$chr,
    y1 = stripes_df$span_y1,
    y2 = stripes_df$span_y2,
    name = stripes_df$stripe_id,
    score = ifelse(is.na(stripes_df$FDR) | stripes_df$FDR == 0,
                   300, -log10(stripes_df$FDR + 1e-300)),
    strand1 = ".",
    strand2 = ".",
    # Color by direction and confidence (RGB format for JuiceBox)
    color = case_when(
      stripes_df$direction == "lost" & stripes_df$direction_confidence == "high" ~ "0,0,139",
      stripes_df$direction == "lost" & stripes_df$direction_confidence == "medium" ~ "65,105,225",
      stripes_df$direction == "lost" ~ "135,206,235",
      stripes_df$direction == "gained" & stripes_df$direction_confidence == "high" ~ "139,0,0",
      stripes_df$direction == "gained" & stripes_df$direction_confidence == "medium" ~ "220,20,60",
      stripes_df$direction == "gained" ~ "255,160,122",
      TRUE ~ "128,128,128"
    ),
    # Key analysis columns
    direction = stripes_df$direction,
    direction_confidence = stripes_df$direction_confidence,
    logFC = round(stripes_df$logFC, 4),
    FDR = signif(stripes_df$FDR, 4),
    source = stripes_df$source,
    # Quagga detection p-values
    pval_ctrl = signif(stripes_df$pval_ctrl, 4),
    pval_mut = signif(stripes_df$pval_mut, 4),
    # Detection confidence
    detection_confidence = stripes_df$confidence,
    in_10kb = stripes_df$in_10kb,
    # Genomic annotations
    nearest_gene = stripes_df$nearest_gene_symbol,
    distance_to_tss = stripes_df$distance_to_tss,
    anchor_type = stripes_df$anchor_type,
    # ChIP-seq overlaps
    h3k27ac = stripes_df$h3k27ac_overlap,
    h3k27me3 = stripes_df$h3k27me3_overlap,
    h3k4me1 = stripes_df$h3k4me1_overlap,
    # Geometry
    stripe_length_kb = round(stripes_df$stripe_length / 1000, 1),
    anchor_width_kb = round(stripes_df$anchor_width / 1000, 1),
    stringsAsFactors = FALSE
  )

  write.table(bedpe_df, output_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)

  return(bedpe_df)
}

# Filter annotated stripes to medium/high confidence and export
cat("Exporting annotated confident stripes (TSV + BEDPE)...\n\n")

for (tp_name in names(annotated_stripes)) {
  df <- annotated_stripes[[tp_name]]
  if (is.null(df)) next

  # Filter to medium/high confidence lost/gained
  confident_df <- df %>%
    filter(direction %in% c("lost", "gained") &
             direction_confidence %in% c("high", "medium"))

  if (nrow(confident_df) == 0) {
    cat(sprintf("  %s: No medium/high confidence stripes\n", tp_name))
    next
  }

  cat(sprintf("%s timepoint: %d confident stripes\n", tp_name, nrow(confident_df)))
  cat(sprintf("  Lost: %d | Gained: %d\n",
              sum(confident_df$direction == "lost"),
              sum(confident_df$direction == "gained")))
  cat(sprintf("  High: %d | Medium: %d\n",
              sum(confident_df$direction_confidence == "high"),
              sum(confident_df$direction_confidence == "medium")))

  output_dir <- file.path(output_base, tp_name)

  # Export annotated TSV (overwrites previous non-annotated version)
  tsv_file <- file.path(output_dir,
                        sprintf("%s_medium_high_confidence_stripes.tsv", tp_name))
  write.table(confident_df, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved TSV: %s\n", basename(tsv_file)))

  # Export BEDPE with annotations
  bedpe_file <- file.path(output_dir,
                          sprintf("%s_medium_high_confidence_stripes.bedpe", tp_name))
  create_annotated_bedpe(confident_df, bedpe_file)
  cat(sprintf("  Saved BEDPE: %s\n", basename(bedpe_file)))

  cat("\n")
}

cat("Section 4B complete: Annotated confident stripe exports done\n")

# =============================================================================
# SECTION 5: SUMMARY STATISTICS & FINAL EXPORT
# =============================================================================

cat("\n========================================\n")
cat("SECTION 5: Summary Statistics\n")
cat("========================================\n\n")

# Generate comprehensive summary
summary_lines <- c(
  "=========================================",
  "DIFFERENTIAL STRIPE ANALYSIS SUMMARY",
  "=========================================",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Working directory: %s", getwd()),
  ""
)

for (tp_name in names(annotated_stripes)) {
  df <- annotated_stripes[[tp_name]]
  if (is.null(df)) next

  summary_lines <- c(summary_lines,
    sprintf("--- %s TIMEPOINT ---", toupper(tp_name)),
    sprintf("Total stripes: %d", nrow(df)),
    "",
    "Direction breakdown:",
    sprintf("  Lost (control_only): %d", sum(df$direction == "lost", na.rm = TRUE)),
    sprintf("  Gained (mutant_only): %d", sum(df$direction == "gained", na.rm = TRUE)),
    sprintf("  Shared: %d", sum(df$direction %in% c("unchanged", "strengthened", "weakened"), na.rm = TRUE)),
    "",
    "Confidence tiers (lost/gained only):",
    sprintf("  High: %d", sum(df$direction_confidence == "high" & df$direction %in% c("lost", "gained"), na.rm = TRUE)),
    sprintf("  Medium: %d", sum(df$direction_confidence == "medium" & df$direction %in% c("lost", "gained"), na.rm = TRUE)),
    sprintf("  Low: %d", sum(df$direction_confidence == "low" & df$direction %in% c("lost", "gained"), na.rm = TRUE)),
    "",
    "Stripe length statistics:",
    sprintf("  Median: %.1f kb", median(df$stripe_length, na.rm = TRUE) / 1000),
    sprintf("  Mean: %.1f kb", mean(df$stripe_length, na.rm = TRUE) / 1000),
    sprintf("  Range: %.1f - %.1f kb", min(df$stripe_length, na.rm = TRUE) / 1000,
            max(df$stripe_length, na.rm = TRUE) / 1000),
    ""
  )

  if ("anchor_type" %in% colnames(df)) {
    anchor_table <- table(df$anchor_type)
    summary_lines <- c(summary_lines,
      "Anchor type distribution:",
      paste0("  ", names(anchor_table), ": ", anchor_table, " (",
             round(100 * anchor_table / sum(anchor_table), 1), "%)"),
      ""
    )
  }

  summary_lines <- c(summary_lines,
    "Directional consistency:",
    sprintf("  Consistent: %d (%.1f%%)",
            sum(df$direction_consistent == TRUE, na.rm = TRUE),
            100 * mean(df$direction_consistent == TRUE, na.rm = TRUE)),
    sprintf("  Inconsistent: %d", sum(df$direction_consistent == FALSE, na.rm = TRUE)),
    "",
    "Statistical significance:",
    sprintf("  FDR < 0.05: %d", sum(df$FDR < 0.05, na.rm = TRUE)),
    sprintf("  FDR < 0.10: %d", sum(df$FDR < 0.10, na.rm = TRUE)),
    "",
    "========================================",
    ""
  )
}

# Add output files summary
summary_lines <- c(summary_lines,
  "OUTPUT FILES GENERATED:",
  "",
  "Per timepoint (early/late):",
  "  - volcano_{timepoint}.pdf",
  "  - length_distribution_{timepoint}.pdf",
  "  - length_statistics_{timepoint}.tsv",
  "  - anchor_classification_{timepoint}.pdf",
  "  - {timepoint}_annotated_stripes.tsv (all stripes)",
  "  - {timepoint}_medium_high_confidence_stripes.tsv (confident only)",
  "  - {timepoint}_medium_high_confidence_stripes.bedpe (confident, for JuiceBox)",
  "",
  "Combined:",
  "  - volcano_combined.pdf",
  "  - length_comparison.pdf",
  "  - summary_statistics.txt",
  "",
  "========================================="
)

# Write summary
summary_file <- file.path(output_base, "combined", "summary_statistics.txt")
writeLines(summary_lines, summary_file)
cat("Summary statistics saved to:", summary_file, "\n")

# Print summary to console
cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n")

# =============================================================================
# COMPLETION MESSAGE
# =============================================================================

cat("\n========================================\n")
cat("STRIPE VISUALIZATION ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Output directory: %s\n\n", output_base))

cat("Generated files:\n\n")

cat("Early timepoint (outputs/visualizations/early/):\n")
cat("  - volcano_early.pdf\n")
cat("  - length_distribution_early.pdf\n")
cat("  - anchor_classification_early.pdf\n")
cat("  - early_annotated_stripes.tsv (all stripes)\n")
cat("  - early_medium_high_confidence_stripes.tsv (confident only)\n")
cat("  - early_medium_high_confidence_stripes.bedpe (for JuiceBox)\n\n")

cat("Late timepoint (outputs/visualizations/late/):\n")
cat("  - volcano_late.pdf\n")
cat("  - length_distribution_late.pdf\n")
cat("  - anchor_classification_late.pdf\n")
cat("  - late_annotated_stripes.tsv (all stripes)\n")
cat("  - late_medium_high_confidence_stripes.tsv (confident only)\n")
cat("  - late_medium_high_confidence_stripes.bedpe (for JuiceBox)\n\n")

cat("Combined (outputs/visualizations/combined/):\n")
cat("  - volcano_combined.pdf\n")
cat("  - length_comparison.pdf\n")
cat("  - summary_statistics.txt\n\n")

cat("========================================\n\n")
