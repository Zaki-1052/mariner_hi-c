# biomodal/downstream/scripts/biomodal_visualizations.R
# Comprehensive visualization script for biomodal DUET evoC differential methylation analysis
# BAP1-KO mutant vs wildtype control - Mouse mm10/GRCm38
# Author: Zakir Alibhai
# Date: 2026-01-21

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base paths - determine script location
get_script_dir <- function() {
  # Try multiple methods to find script directory
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # Fallback to known path
  return("/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts")
}

SCRIPT_DIR <- get_script_dir()
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)

# Data file paths
DATA_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_mc_control__mutant_20260121_172049.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_hmc_control__mutant_20260121_172049.bed"),
  bioqc_json = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/BioQC_20260121_170327/biological_qc_report_4_samples_20260121_170327.json"),
  upstream_csv = file.path(BASE_DIR, "../upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv"),
  # Additional region DMR files for region comparison
  cpg_islands_mc = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260121_171455/DMR_mc_control__mutant_20260121_171455.bed"),
  cpg_shores_mc = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260121_171910/DMR_mc_control__mutant_20260121_171910.bed"),
  cpg_shelves_mc = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260121_171718/DMR_mc_control__mutant_20260121_171718.bed"),
  promoters_mc = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260121_172220/DMR_mc_control__mutant_20260121_172220.bed"),
  tss_mc = file.path(BASE_DIR, "modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260121_172333/DMR_mc_control__mutant_20260121_172333.bed")
)

# ChIP-seq peak file paths (from peaks/beds/)
# Use Late timepoint peaks to match the adult BAP1-KO analysis
CHIP_PEAK_FILES <- list(
  ctcf     = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/CTCF.bed",
  h3k27ac  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/H3K27acCerebellumLate2.bed",
  h3k27me3 = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/H3K27me3CerebellumLate1.bed",
  h3k4me1  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/H3K4me1CerebellumLate1.bed",
  h3k4me3  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/H3K4me3CerebellumLate2.bed",
  bivalent = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/beds/Bivalent_Cerebellum_Late.bed"
)

# Output directory
OUTPUT_DIR <- file.path(BASE_DIR, "plots/visualizations")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Key genes to highlight
KEY_GENES <- c("Syt1", "Zbtb20", "Trpm3", "Epha3", "Mcu", "Cntnap2", "Lpp", "Dlgap1", "Arhgap26", "Cdh8")

# Statistical thresholds
Q_THRESHOLD <- 0.05

# Sample metadata
SAMPLES <- data.frame(
  sample_id = c("evoC-Bap1-ctrl-F", "evoC-Bap1-ctrl-M", "evoC-Bap1-mut-F", "evoC-Bap1-mut-M"),
  condition = c("Control", "Control", "Mutant", "Mutant"),
  sex = c("Female", "Male", "Female", "Male"),
  short_name = c("ctrl-F", "ctrl-M", "mut-F", "mut-M"),
  stringsAsFactors = FALSE
)

# Color schemes
COLORS <- list(
  condition = c("Control" = "#2166AC", "Mutant" = "#B2182B"),
  sex = c("Female" = "#E377C2", "Male" = "#17BECF"),
  direction = c("Hypermethylated" = "#D7191C", "Hypomethylated" = "#2C7BB6"),
  methylation = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
  significant = c("Significant" = "#E41A1C", "Not Significant" = "grey70")
)

# Chromatin state classification (consistent with annotate_loops_extended.R)
CHROMATIN_STATE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                           "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")

CHROMATIN_STATE_COLORS <- c(
  "Active_Promoter" = "#e41a1c",     # Red - active transcription
  "Repressed_Promoter" = "#756bb1",  # Purple - Polycomb-silenced promoter
  "Bivalent_Promoter" = "#984ea3",   # Magenta - developmental poised
  "Polycomb" = "#4daf4a",            # Green - distal repressive
  "Active_Enhancer" = "#377eb8",     # Blue - active enhancer
  "Poised_Enhancer" = "#ff7f00",     # Orange - primed enhancer
  "Other" = "#999999"                # Gray - unmarked
)

# TSS threshold for promoter classification
TSS_THRESHOLD <- 2000  # 2kb

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\n")

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
  library(jsonlite)
  library(pheatmap)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(ggVennDiagram)
  # ChIP-seq analysis packages (Section 10)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# Source multi-format output utility
util_path <- file.path(dirname(BASE_DIR), "scripts/utils/multi_format_output.R")
if (!file.exists(util_path)) {
  util_path <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/scripts/utils/multi_format_output.R"
}
source(util_path)

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Packages loaded successfully.\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load DMR BED file with proper column names
load_dmr_bed <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(NULL)
  }

  # Read with header (file has header row)
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "")

  # Standardize column names based on number of columns
  base_cols <- c("chr", "start", "end", "num_contexts", "mean_coverage",
                 "mean_mod_group1", "mean_mod_group2", "mod_fold_change",
                 "mod_difference", "test_statistic", "dmr_pvalue", "dmr_qvalue",
                 "annotation")

  if (ncol(df) == 14) {
    colnames(df) <- c(base_cols, "gene")
  } else if (ncol(df) == 13) {
    colnames(df) <- base_cols
    df$gene <- df$annotation  # Use annotation as gene placeholder
  } else {
    warning("Unexpected number of columns: ", ncol(df))
    return(NULL)
  }

  # Add chr prefix if missing
  if (nrow(df) > 0 && !grepl("^chr", df$chr[1])) {
    df$chr <- paste0("chr", df$chr)
  }

  # Ensure numeric columns are properly typed
  numeric_cols <- c("start", "end", "num_contexts", "mean_coverage",
                    "mean_mod_group1", "mean_mod_group2", "mod_fold_change",
                    "mod_difference", "test_statistic", "dmr_pvalue", "dmr_qvalue")
  for (col in numeric_cols) {
    df[[col]] <- as.numeric(df[[col]])
  }

  # Add significance column
  df$significant <- df$dmr_qvalue < Q_THRESHOLD
  df$direction <- ifelse(df$mod_difference > 0, "Hypermethylated", "Hypomethylated")

  # Handle very small q-values for plotting (cap at 1e-300)
  df$neg_log10_q <- -log10(pmax(df$dmr_qvalue, 1e-300, na.rm = TRUE))
  df$neg_log10_q[is.infinite(df$neg_log10_q)] <- 300

  return(df)
}

#' Create common ggplot theme
theme_biomodal <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

# =============================================================================
# ChIP-seq HELPER FUNCTIONS (Section 10)
# =============================================================================

#' Load ChIP-seq peaks from BED file as GRanges
#' @param bed_path Path to BED file
#' @param mark_name Name of the mark (for logging)
#' @return GRanges object
load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    warning(sprintf("%s BED file not found: %s", mark_name, bed_path))
    return(NULL)
  }

  df <- read.table(bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("  %s: %d peaks\n", mark_name, length(gr)))
  return(gr)
}

#' Create GRanges from DMR data frame
#' @param dmr_df Data frame with chr, start, end columns
#' @return GRanges object
dmr_to_granges <- function(dmr_df) {
  GRanges(
    seqnames = dmr_df$chr,
    ranges = IRanges(start = dmr_df$start, end = dmr_df$end),
    gene = dmr_df$gene,
    mod_difference = dmr_df$mod_difference,
    dmr_qvalue = dmr_df$dmr_qvalue,
    significant = dmr_df$significant,
    direction = dmr_df$direction
  )
}

#' Compute ChIP-seq overlaps for DMR GRanges
#' @param dmr_gr GRanges with DMR coordinates
#' @param chip_peaks List of GRanges for each ChIP mark
#' @return data.frame with overlap columns
compute_chip_overlaps <- function(dmr_gr, chip_peaks) {
  overlaps <- data.frame(
    CTCF_overlap = countOverlaps(dmr_gr, chip_peaks$ctcf) > 0,
    H3K27ac_overlap = countOverlaps(dmr_gr, chip_peaks$h3k27ac) > 0,
    H3K27me3_overlap = countOverlaps(dmr_gr, chip_peaks$h3k27me3) > 0,
    H3K4me1_overlap = countOverlaps(dmr_gr, chip_peaks$h3k4me1) > 0,
    H3K4me3_overlap = countOverlaps(dmr_gr, chip_peaks$h3k4me3) > 0,
    Bivalent_overlap = countOverlaps(dmr_gr, chip_peaks$bivalent) > 0
  )
  return(overlaps)
}

#' Classify chromatin state using 7-category priority system
#'
#' Priority order (consistent with annotate_loops_extended.R):
#'   1. Active_Promoter:    H3K4me3+ AND NOT H3K27me3 AND <=2kb from TSS
#'   2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND <=2kb from TSS
#'   3. Bivalent_Promoter:  K4me3+K27me3 overlap (pre-computed)
#'   4. Polycomb:           H3K27me3+ AND >2kb from TSS
#'   5. Active_Enhancer:    H3K27ac+ AND >2kb from TSS
#'   6. Poised_Enhancer:    H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
#'   7. Other:              Default (no marks)
#'
#' @param overlaps data.frame with overlap columns
#' @param distance_to_tss Numeric vector of distances to nearest TSS
#' @param tss_threshold Distance threshold for promoter (default 2000bp)
#' @return Character vector with chromatin state classifications
classify_chromatin_state <- function(overlaps, distance_to_tss, tss_threshold = 2000) {
  n <- nrow(overlaps)
  chromatin_state <- rep("Other", n)

  # Extract overlap columns
  h3k27ac <- overlaps$H3K27ac_overlap
  h3k27me3 <- overlaps$H3K27me3_overlap
  h3k4me1 <- overlaps$H3K4me1_overlap
  h3k4me3 <- overlaps$H3K4me3_overlap
  bivalent <- overlaps$Bivalent_overlap

  # 1. Active_Promoter: H3K4me3+ AND NOT H3K27me3 AND <=2kb from TSS
  is_active_promoter <- h3k4me3 & !h3k27me3 &
    !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  chromatin_state[is_active_promoter] <- "Active_Promoter"

  # 2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND <=2kb from TSS
  is_repressed_promoter <- !is_active_promoter &
    h3k27me3 & !h3k27ac &
    !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  chromatin_state[is_repressed_promoter] <- "Repressed_Promoter"

  # 3. Bivalent_Promoter: K4me3+K27me3 overlap (pre-computed)
  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent
  chromatin_state[is_bivalent] <- "Bivalent_Promoter"

  # 4. Polycomb: H3K27me3+ AND >2kb from TSS
  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
    h3k27me3 & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_polycomb] <- "Polycomb"

  # 5. Active_Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb &
    h3k27ac & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_active_enhancer] <- "Active_Enhancer"

  # 6. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb & !is_active_enhancer &
    h3k4me1 & !h3k27ac & !h3k27me3 &
    (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_poised_enhancer] <- "Poised_Enhancer"

  # 7. Other: Default (no marks)
  return(factor(chromatin_state, levels = CHROMATIN_STATE_ORDER))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("================================================================================\n")
cat("LOADING DATA\n")
cat("================================================================================\n\n")

# Load gene body DMR data
cat("Loading gene body mC DMRs...\n")
mc_dmr <- load_dmr_bed(DATA_PATHS$mc_dmr)
cat(sprintf("  Loaded %d genes\n", nrow(mc_dmr)))

cat("Loading gene body hmC DMRs...\n")
hmc_dmr <- load_dmr_bed(DATA_PATHS$hmc_dmr)
cat(sprintf("  Loaded %d genes\n", nrow(hmc_dmr)))

# Load BioQC JSON
cat("Loading Biological QC data...\n")
bioqc <- fromJSON(DATA_PATHS$bioqc_json)
cat("  Loaded QC data successfully\n")

# Load upstream summary
cat("Loading upstream sequencing metrics...\n")
if (file.exists(DATA_PATHS$upstream_csv)) {
  upstream <- read.csv(DATA_PATHS$upstream_csv, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d samples\n", nrow(upstream)))
} else {
  cat("  WARNING: Upstream CSV not found, skipping QC plots\n")
  upstream <- NULL
}

# Load other region DMRs for comparison
cat("\nLoading regional DMR data for comparison...\n")
region_dmrs <- list(
  "Gene bodies" = mc_dmr,
  "CpG islands" = load_dmr_bed(DATA_PATHS$cpg_islands_mc),
  "CpG shores" = load_dmr_bed(DATA_PATHS$cpg_shores_mc),
  "CpG shelves" = load_dmr_bed(DATA_PATHS$cpg_shelves_mc),
  "Promoters" = load_dmr_bed(DATA_PATHS$promoters_mc),
  "TSS regions" = load_dmr_bed(DATA_PATHS$tss_mc)
)

for (region in names(region_dmrs)) {
  if (!is.null(region_dmrs[[region]])) {
    n_sig <- sum(region_dmrs[[region]]$significant)
    n_total <- nrow(region_dmrs[[region]])
    cat(sprintf("  %s: %d significant / %d total (%.1f%%)\n",
                region, n_sig, n_total, 100 * n_sig / n_total))
  }
}

cat("\nData loading complete.\n\n")

# =============================================================================
# SECTION 1: QC & DATA OVERVIEW
# =============================================================================

cat("================================================================================\n")
cat("SECTION 1: QC & DATA OVERVIEW\n")
cat("================================================================================\n\n")

# Extract QC metrics from bioqc JSON
if (!is.null(bioqc)) {
  # Parse QC metrics from the JSON structure
  qc_df <- data.frame(
    sample = SAMPLES$short_name,
    condition = SAMPLES$condition,
    stringsAsFactors = FALSE
  )

  # Create QC overview plot
  cat("Creating QC overview plot...\n")

  # If upstream data available, create sequencing metrics plot
  if (!is.null(upstream)) {
    # Extract key metrics (using actual column names from CSV)
    seq_metrics <- upstream %>%
      dplyr::select(sample_id, bamlet_mapped_reads_genome, bamlet_mapped_bases_cigar_genome,
                    bamlet_prop_duplicated_reads_genome, bamlet_mean_phred_genome) %>%
      dplyr::mutate(
        Mapped_Reads_M = bamlet_mapped_reads_genome / 1e6,
        Mapped_Bases_B = bamlet_mapped_bases_cigar_genome / 1e9,
        Duplication_Rate = bamlet_prop_duplicated_reads_genome,
        Mean_Phred_Score = bamlet_mean_phred_genome,
        sample_short = gsub("evoC-Bap1-", "", sample_id),
        condition = ifelse(grepl("ctrl", sample_id), "Control", "Mutant")
      )

    # Mapped reads bar chart - use points+bars to show actual values clearly
    p1 <- ggplot(seq_metrics, aes(x = sample_short, y = Mapped_Reads_M, fill = condition)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%.0f", Mapped_Reads_M)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(0, max(seq_metrics$Mapped_Reads_M) * 1.15)) +
      labs(title = "Mapped Reads", x = "", y = "Reads (millions)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Mapped bases bar chart
    p2 <- ggplot(seq_metrics, aes(x = sample_short, y = Mapped_Bases_B, fill = condition)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%.1f", Mapped_Bases_B)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(0, max(seq_metrics$Mapped_Bases_B) * 1.15)) +
      labs(title = "Mapped Bases", x = "", y = "Bases (billions)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Duplication rate - use lollipop plot with zoomed y-axis to show differences
    dup_rate_pct <- seq_metrics$Duplication_Rate * 100
    p3 <- ggplot(seq_metrics, aes(x = sample_short, y = Duplication_Rate * 100,
                                   color = condition)) +
      geom_segment(aes(xend = sample_short, y = min(dup_rate_pct) - 1,
                       yend = Duplication_Rate * 100), linewidth = 2) +
      geom_point(size = 5) +
      geom_text(aes(label = sprintf("%.1f%%", Duplication_Rate * 100)),
                vjust = -1, size = 3, color = "black") +
      scale_color_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(min(dup_rate_pct) - 2, max(dup_rate_pct) + 2)) +
      labs(title = "Duplication Rate",
           subtitle = sprintf("Range: %.1f%% - %.1f%%", min(dup_rate_pct), max(dup_rate_pct)),
           x = "", y = "Rate (%)") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Mean Phred score - use lollipop plot with zoomed y-axis
    p4 <- ggplot(seq_metrics, aes(x = sample_short, y = Mean_Phred_Score,
                                   color = condition)) +
      geom_segment(aes(xend = sample_short,
                       y = min(seq_metrics$Mean_Phred_Score) - 0.5,
                       yend = Mean_Phred_Score), linewidth = 2) +
      geom_point(size = 5) +
      geom_text(aes(label = sprintf("%.2f", Mean_Phred_Score)),
                vjust = -1, size = 3, color = "black") +
      scale_color_manual(values = COLORS$condition) +
      scale_y_continuous(limits = c(min(seq_metrics$Mean_Phred_Score) - 1,
                                    max(seq_metrics$Mean_Phred_Score) + 0.5)) +
      labs(title = "Mean Phred Score",
           subtitle = sprintf("All samples: ~%.1f (excellent)",
                             mean(seq_metrics$Mean_Phred_Score)),
           x = "", y = "Phred Score") +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    qc_combined <- (p1 | p2) / (p3 | p4) +
      plot_annotation(
        title = "Sequencing Quality Metrics - Biomodal DUET evoC",
        subtitle = "BAP1-KO vs Wildtype Control",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                      plot.subtitle = element_text(hjust = 0.5, size = 12))
      )

    save_multiformat_ggplot(qc_combined, file.path(OUTPUT_DIR, "01_qc_overview"),
                            width = 12, height = 10)
  } else {
    cat("  Skipping QC plot (upstream data not available)\n")
  }
}

cat("Section 1 complete.\n\n")

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

# =============================================================================
# SECTION 3: DMR STATISTICS BY GENOMIC REGION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 3: DMR STATISTICS BY GENOMIC REGION\n")
cat("================================================================================\n\n")

# Create region comparison data frame
region_stats <- data.frame(
  Region = character(),
  Total = integer(),
  Significant = integer(),
  Percentage = numeric(),
  stringsAsFactors = FALSE
)

for (region in names(region_dmrs)) {
  if (!is.null(region_dmrs[[region]])) {
    region_stats <- rbind(region_stats, data.frame(
      Region = region,
      Total = nrow(region_dmrs[[region]]),
      Significant = sum(region_dmrs[[region]]$significant),
      Percentage = 100 * sum(region_dmrs[[region]]$significant) / nrow(region_dmrs[[region]])
    ))
  }
}

# Order by significant count
region_stats <- region_stats %>%
  mutate(Region = factor(Region, levels = Region[order(Significant, decreasing = TRUE)]))

cat("Creating region comparison bar chart...\n")

p_region <- ggplot(region_stats, aes(x = Region, y = Significant, fill = Region)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Significant, Percentage)),
            vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Significant mC DMRs by Genomic Region (CG Context)",
    subtitle = sprintf("q-value < %.2f | Gene bodies dominate differential methylation", Q_THRESHOLD),
    x = "", y = "Number of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_region, file.path(OUTPUT_DIR, "03a_dmr_by_region"),
                        width = 10, height = 7)

# mC vs hmC comparison for gene bodies
cat("Creating mC vs hmC comparison...\n")

mc_hmc_compare <- data.frame(
  Type = c("5mC", "5hmC"),
  Total = c(nrow(mc_dmr), nrow(hmc_dmr)),
  Significant = c(sum(mc_dmr$significant), sum(hmc_dmr$significant)),
  Percentage = c(100 * sum(mc_dmr$significant) / nrow(mc_dmr),
                 100 * sum(hmc_dmr$significant) / nrow(hmc_dmr))
)

p_mc_hmc <- ggplot(mc_hmc_compare, aes(x = Type, y = Significant, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Significant, Percentage)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = COLORS$methylation) +
  labs(
    title = "Gene Body DMRs: 5mC vs 5hmC",
    subtitle = "5hmC shows more differential genes than 5mC",
    x = "", y = "Number of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Context comparison (CG vs CHG vs CHH)
cat("Creating methylation context comparison...\n")

context_data <- data.frame(
  Context = c("CG (CpG)", "CHG", "CHH"),
  Baseline_Methylation = c(72, 0.65, 0.87),  # From analysis_summary
  Significant_DMRs = c(sum(mc_dmr$significant), 0, 0),  # CHG/CHH have 0 significant
  Label = c("Primary Signal", "No Signal", "No Signal")
)

p_context_meth <- ggplot(context_data, aes(x = Context, y = Baseline_Methylation, fill = Context)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Baseline_Methylation)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("CG (CpG)" = "#E41A1C", "CHG" = "#999999", "CHH" = "#CCCCCC")) +
  labs(
    title = "Baseline Methylation by Context",
    subtitle = "CpG methylation is ~100x more abundant than non-CpG",
    x = "Methylation Context", y = "Methylation Level (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

p_context_dmr <- ggplot(context_data, aes(x = Context, y = Significant_DMRs, fill = Context)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Label), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("CG (CpG)" = "#E41A1C", "CHG" = "#999999", "CHH" = "#CCCCCC")) +
  labs(
    title = "Significant DMRs by Context",
    subtitle = "No significant changes in non-CpG methylation",
    x = "Methylation Context", y = "Significant Gene Body DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Combine into panel
p_region_combined <- (p_region | p_mc_hmc) / (p_context_meth | p_context_dmr) +
  plot_annotation(
    title = "DMR Statistics Summary",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_region_combined, file.path(OUTPUT_DIR, "03_dmr_region_statistics"),
                        width = 14, height = 12)

# -------------------------------------------------------------------------
# FIGURE 3b: Direction comparison - mC vs hmC asymmetry
# -------------------------------------------------------------------------
cat("Creating Figure 3b: Methylation direction comparison...\n")

# Calculate direction statistics for mC
mc_sig_stats <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::summarise(
    Hypermethylated = sum(direction == "Hypermethylated"),
    Hypomethylated = sum(direction == "Hypomethylated"),
    Total = n()
  ) %>%
 dplyr::mutate(
    Hyper_pct = 100 * Hypermethylated / Total,
    Hypo_pct = 100 * Hypomethylated / Total
  )

# Calculate direction statistics for hmC
hmc_sig_stats <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::summarise(
    Increased = sum(direction == "Hypermethylated"),
    Decreased = sum(direction == "Hypomethylated"),
    Total = n()
  ) %>%
  dplyr::mutate(
    Inc_pct = 100 * Increased / Total,
    Dec_pct = 100 * Decreased / Total
  )

# Create data frame for plotting
direction_data <- data.frame(
  Modification = c("5mC", "5mC", "5hmC", "5hmC"),
  Direction = c("Increased", "Decreased", "Increased", "Decreased"),
  Count = c(mc_sig_stats$Hypermethylated, mc_sig_stats$Hypomethylated,
            hmc_sig_stats$Increased, hmc_sig_stats$Decreased),
  Percentage = c(mc_sig_stats$Hyper_pct, mc_sig_stats$Hypo_pct,
                 hmc_sig_stats$Inc_pct, hmc_sig_stats$Dec_pct)
)

direction_data$Modification <- factor(direction_data$Modification, levels = c("5mC", "5hmC"))
direction_data$Direction <- factor(direction_data$Direction, levels = c("Increased", "Decreased"))

# Create the bar chart
p_direction <- ggplot(direction_data, aes(x = Modification, y = Percentage, fill = Direction)) +
 geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Increased" = "#D7191C", "Decreased" = "#2C7BB6"),
                    name = "Direction in\nBAP1-KO Mutant") +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    title = "Asymmetric Direction of Methylation Changes",
    subtitle = sprintf("5mC: %d significant genes | 5hmC: %d significant genes",
                       mc_sig_stats$Total, hmc_sig_stats$Total),
    x = "Methylation Type",
    y = "Percentage of Significant DMRs (%)"
  ) +
  theme_biomodal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16)
  ) +
  # Add annotation for key finding
  annotate("text", x = 1.5, y = 98,
           label = "mC increases while hmC decreases\n= Blocked TET-mediated demethylation",
           size = 4, fontface = "italic", color = "grey30")

save_multiformat_ggplot(p_direction, file.path(OUTPUT_DIR, "03b_direction_comparison"),
                        width = 10, height = 8)

cat("Section 3 complete.\n\n")

# =============================================================================
# SECTION 4: VOLCANO PLOTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 4: VOLCANO PLOTS\n")
cat("================================================================================\n\n")

# Function to create volcano plot
create_volcano <- function(df, title, subtitle, meth_type = "5mC") {
  # Identify genes to label
  df$label <- ""
  df$label[df$gene %in% KEY_GENES & df$significant] <- df$gene[df$gene %in% KEY_GENES & df$significant]

  # Color based on significance and direction
  df$color_group <- case_when(
    !df$significant ~ "Not Significant",
    df$mod_difference > 0 ~ "Hypermethylated",
    TRUE ~ "Hypomethylated"
  )

  # Cap -log10(q) for visualization
  df$neg_log10_q_capped <- pmin(df$neg_log10_q, 300)

  p <- ggplot(df, aes(x = mod_difference * 100, y = neg_log10_q_capped)) +
    geom_point(aes(color = color_group), alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    geom_text_repel(
      aes(label = label),
      size = 3.5,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = c("Hypermethylated" = "#D7191C",
                 "Hypomethylated" = "#2C7BB6",
                 "Not Significant" = "grey70"),
      name = "Direction"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste(meth_type, "Difference (Mutant - Control, %)"),
      y = expression(-log[10](q-value))
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  return(p)
}

cat("Creating 5mC volcano plot...\n")
p_volcano_mc <- create_volcano(
  mc_dmr,
  "Gene Body 5mC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f)", sum(mc_dmr$significant), Q_THRESHOLD),
  "5mC"
)

save_multiformat_ggplot(p_volcano_mc, file.path(OUTPUT_DIR, "04a_volcano_mc"),
                        width = 10, height = 8)

cat("Creating 5hmC volcano plot...\n")
p_volcano_hmc <- create_volcano(
  hmc_dmr,
  "Gene Body 5hmC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f)", sum(hmc_dmr$significant), Q_THRESHOLD),
  "5hmC"
)

save_multiformat_ggplot(p_volcano_hmc, file.path(OUTPUT_DIR, "04b_volcano_hmc"),
                        width = 10, height = 8)

# Combined volcano panel
cat("Creating combined volcano panel...\n")
p_volcano_combined <- p_volcano_mc + p_volcano_hmc +
  plot_annotation(
    title = "Differential Methylation: 5mC vs 5hmC",
    subtitle = "Key genes labeled: Syt1, Zbtb20, Trpm3, Cntnap2",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5))
  )

save_multiformat_ggplot(p_volcano_combined, file.path(OUTPUT_DIR, "04_volcano_plots"),
                        width = 16, height = 8)

cat("Section 4 complete.\n\n")

# =============================================================================
# SECTION 5: COORDINATED CHANGES ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 5: COORDINATED CHANGES ANALYSIS (KEY INSIGHT)\n")
cat("================================================================================\n\n")

# Merge mC and hmC data by gene
cat("Identifying genes with coordinated mC/hmC changes...\n")

mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_ctrl = mean_mod_group1, mc_mut = mean_mod_group2)

hmc_sig <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2)

coordinated <- inner_join(mc_sig, hmc_sig, by = "gene")

# Identify genes with opposite directions (mC up, hmC down)
coordinated <- coordinated %>%
  mutate(
    coordinated_pattern = (mc_diff > 0 & hmc_diff < 0),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  ) %>%
  arrange(desc(combined_effect))

cat(sprintf("  %d genes significant in both mC and hmC\n", nrow(coordinated)))
cat(sprintf("  %d genes show coordinated mC↑/hmC↓ pattern (%.1f%%)\n",
            sum(coordinated$coordinated_pattern),
            100 * sum(coordinated$coordinated_pattern) / nrow(coordinated)))

# Save coordinated changes table
write.table(coordinated, file.path(TABLES_DIR, "coordinated_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved coordinated_changes.tsv\n")

# Scatter plot: mC vs hmC change
cat("Creating mC vs hmC scatter plot...\n")

coordinated$label <- ""
coordinated$label[coordinated$gene %in% KEY_GENES] <- coordinated$gene[coordinated$gene %in% KEY_GENES]

p_scatter <- ggplot(coordinated, aes(x = mc_diff * 100, y = hmc_diff * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = coordinated_pattern), alpha = 0.6, size = 2) +
  geom_text_repel(
    aes(label = label),
    size = 4,
    max.overlaps = 15,
    box.padding = 0.5,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("TRUE" = "#8E44AD", "FALSE" = "grey60"),
    labels = c("TRUE" = "mC↑ / hmC↓", "FALSE" = "Other"),
    name = "Pattern"
  ) +
  labs(
    title = "Coordinated 5mC and 5hmC Changes",
    subtitle = sprintf("%d genes significant in both | %.0f%% show mC↑/hmC↓ pattern",
                       nrow(coordinated), 100 * mean(coordinated$coordinated_pattern)),
    x = "5mC Change (%)",
    y = "5hmC Change (%)"
  ) +
  theme_biomodal() +
  annotate("text", x = 15, y = -12, label = "COORDINATED\nmC↑ / hmC↓",
           color = "#8E44AD", fontface = "bold", size = 4) +
  annotate("rect", xmin = 0, xmax = 25, ymin = -20, ymax = 0,
           alpha = 0.1, fill = "#8E44AD")

save_multiformat_ggplot(p_scatter, file.path(OUTPUT_DIR, "05a_mc_hmc_scatter"),
                        width = 10, height = 9)

# Top coordinated genes bar chart
cat("Creating top coordinated genes bar chart...\n")

top_coordinated <- coordinated %>%
  filter(coordinated_pattern) %>%
  head(20) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

top_coord_long <- top_coordinated %>%
  dplyr::select(gene, mc_diff, hmc_diff) %>%
  tidyr::pivot_longer(cols = c(mc_diff, hmc_diff), names_to = "type", values_to = "change") %>%
  dplyr::mutate(
    type = ifelse(type == "mc_diff", "5mC", "5hmC"),
    change_pct = change * 100
  )

p_top_coord <- ggplot(top_coord_long, aes(x = gene, y = change_pct, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$methylation, name = "Modification") +
  coord_flip() +
  labs(
    title = "Top 20 Genes with Coordinated mC↑/hmC↓ Pattern",
    subtitle = "Consistent with impaired TET-mediated demethylation",
    x = "", y = "Change (Mutant - Control, %)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_top_coord, file.path(OUTPUT_DIR, "05b_top_coordinated_genes"),
                        width = 10, height = 10)

# Detailed locus plots for key genes
cat("Creating detailed locus plots...\n")

create_locus_plot <- function(gene_name, mc_data, hmc_data) {
  mc_gene <- mc_data %>% filter(gene == gene_name)
  hmc_gene <- hmc_data %>% filter(gene == gene_name)

  if (nrow(mc_gene) == 0 || nrow(hmc_gene) == 0) {
    return(NULL)
  }

  # Create data for plotting
  plot_data <- data.frame(
    Condition = rep(c("Control", "Mutant"), 2),
    Type = c("5mC", "5mC", "5hmC", "5hmC"),
    Methylation = c(mc_gene$mean_mod_group1 * 100, mc_gene$mean_mod_group2 * 100,
                    hmc_gene$mean_mod_group1 * 100, hmc_gene$mean_mod_group2 * 100)
  )

  # Calculate changes
  mc_change <- (mc_gene$mean_mod_group2 - mc_gene$mean_mod_group1) * 100
  hmc_change <- (hmc_gene$mean_mod_group2 - hmc_gene$mean_mod_group1) * 100

  p <- ggplot(plot_data, aes(x = Type, y = Methylation, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = COLORS$condition) +
    labs(
      title = sprintf("%s Gene Body Methylation", gene_name),
      subtitle = sprintf("mC: %+.1f%% | hmC: %+.1f%% | Coordinated pattern", mc_change, hmc_change),
      x = "", y = "Methylation Level (%)"
    ) +
    theme_biomodal() +
    # Add significance annotations
    annotate("text", x = 1, y = max(plot_data$Methylation) * 1.1,
             label = sprintf("q = %.1e", mc_gene$dmr_qvalue), size = 3) +
    annotate("text", x = 2, y = max(plot_data$Methylation) * 1.1,
             label = sprintf("q = %.1e", hmc_gene$dmr_qvalue), size = 3)

  return(p)
}

# Create Syt1 plot (highest effect)
p_syt1 <- create_locus_plot("Syt1", mc_dmr, hmc_dmr)
if (!is.null(p_syt1)) {
  save_multiformat_ggplot(p_syt1, file.path(OUTPUT_DIR, "05c_syt1_detail"),
                          width = 8, height = 6)
}

# Create Zbtb20 plot
p_zbtb20 <- create_locus_plot("Zbtb20", mc_dmr, hmc_dmr)
if (!is.null(p_zbtb20)) {
  save_multiformat_ggplot(p_zbtb20, file.path(OUTPUT_DIR, "05d_zbtb20_detail"),
                          width = 8, height = 6)
}

# Combined key genes panel
key_plots <- list()
for (gene in c("Syt1", "Zbtb20", "Trpm3", "Cntnap2")) {
  p <- create_locus_plot(gene, mc_dmr, hmc_dmr)
  if (!is.null(p)) {
    key_plots[[gene]] <- p
  }
}

if (length(key_plots) == 4) {
  p_key_genes <- (key_plots[[1]] | key_plots[[2]]) / (key_plots[[3]] | key_plots[[4]]) +
    plot_annotation(
      title = "Key Genes: Coordinated Methylation Changes",
      subtitle = "5mC increase + 5hmC decrease pattern in BAP1-KO",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                    plot.subtitle = element_text(hjust = 0.5))
    )

  save_multiformat_ggplot(p_key_genes, file.path(OUTPUT_DIR, "05_coordinated_changes"),
                          width = 14, height = 12)
}

cat("Section 5 complete.\n\n")

# =============================================================================
# SECTION 6: TOP DIFFERENTIALLY METHYLATED GENES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 6: TOP DIFFERENTIALLY METHYLATED GENES\n")
cat("================================================================================\n\n")

# Top 20 mC DMRs
cat("Creating top mC DMRs bar chart...\n")

top_mc <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::arrange(dmr_qvalue) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  head(20) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    direction = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
  )

p_top_mc <- ggplot(top_mc, aes(x = gene, y = mod_difference * 100, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  coord_flip() +
  labs(
    title = "Top 20 Gene Body 5mC DMRs",
    subtitle = "Ranked by q-value | BAP1-KO vs Control",
    x = "", y = "5mC Change (%)"
  ) +
  theme_biomodal()

# Top 20 hmC DMRs
cat("Creating top hmC DMRs bar chart...\n")

top_hmc <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::arrange(dmr_qvalue) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  head(20) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    direction = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
  )

p_top_hmc <- ggplot(top_hmc, aes(x = gene, y = mod_difference * 100, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  coord_flip() +
  labs(
    title = "Top 20 Gene Body 5hmC DMRs",
    subtitle = "Ranked by q-value | BAP1-KO vs Control",
    x = "", y = "5hmC Change (%)"
  ) +
  theme_biomodal()

p_top_combined <- p_top_mc | p_top_hmc
save_multiformat_ggplot(p_top_combined, file.path(OUTPUT_DIR, "06a_top_dmrs"),
                        width = 14, height = 10)

# Save top gene lists
write.table(top_mc, file.path(TABLES_DIR, "top_mc_dmrs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_hmc, file.path(TABLES_DIR, "top_hmc_dmrs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved top_mc_dmrs.tsv and top_hmc_dmrs.tsv\n")

# Venn diagram of overlap
cat("Creating overlap Venn diagram...\n")

mc_genes <- mc_dmr$gene[mc_dmr$significant]
hmc_genes <- hmc_dmr$gene[hmc_dmr$significant]

venn_list <- list(
  "5mC Significant" = mc_genes,
  "5hmC Significant" = hmc_genes
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#4DBBD5") +
  labs(
    title = "Overlap of Significant Genes",
    subtitle = sprintf("5mC: %d | 5hmC: %d | Both: %d",
                       length(mc_genes), length(hmc_genes), length(intersect(mc_genes, hmc_genes)))
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn, file.path(OUTPUT_DIR, "06b_venn_overlap"),
                        width = 8, height = 7)

p_top_all <- (p_top_mc | p_top_hmc) / p_venn +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Top Differentially Methylated Genes",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_top_all, file.path(OUTPUT_DIR, "06_top_genes"),
                        width = 14, height = 14)

cat("Section 6 complete.\n\n")

# =============================================================================
# SECTION 7: EFFECT SIZE DISTRIBUTIONS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 7: EFFECT SIZE DISTRIBUTIONS\n")
cat("================================================================================\n\n")

cat("Creating effect size distribution plots...\n")

# Histogram of mC differences
p_hist_mc <- ggplot(mc_dmr, aes(x = mod_difference * 100, fill = significant)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
                    name = "") +
  labs(
    title = "5mC Effect Size Distribution",
    subtitle = sprintf("Mean change: %+.2f%% (significant genes)",
                       mean(mc_dmr$mod_difference[mc_dmr$significant]) * 100),
    x = "5mC Change (Mutant - Control, %)", y = "Count"
  ) +
  theme_biomodal()

# Histogram of hmC differences
p_hist_hmc <- ggplot(hmc_dmr, aes(x = mod_difference * 100, fill = significant)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
                    name = "") +
  labs(
    title = "5hmC Effect Size Distribution",
    subtitle = sprintf("Mean change: %+.2f%% (significant genes)",
                       mean(hmc_dmr$mod_difference[hmc_dmr$significant]) * 100),
    x = "5hmC Change (Mutant - Control, %)", y = "Count"
  ) +
  theme_biomodal()

# Violin plot comparison
violin_data <- data.frame(
  Type = c(rep("5mC", sum(mc_dmr$significant)), rep("5hmC", sum(hmc_dmr$significant))),
  Change = c(mc_dmr$mod_difference[mc_dmr$significant] * 100,
             hmc_dmr$mod_difference[hmc_dmr$significant] * 100)
)

p_violin <- ggplot(violin_data, aes(x = Type, y = Change, fill = Type)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = COLORS$methylation) +
  labs(
    title = "Effect Size Comparison: 5mC vs 5hmC",
    subtitle = "5mC shifts positive (hypermethylation) | 5hmC shifts negative (hypomethylation)",
    x = "", y = "Change (Mutant - Control, %)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

p_effect_combined <- (p_hist_mc | p_hist_hmc) / p_violin +
  plot_annotation(
    title = "Effect Size Distributions",
    subtitle = "Significant genes show opposite directions: 5mC↑ / 5hmC↓",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5))
  )

save_multiformat_ggplot(p_effect_combined, file.path(OUTPUT_DIR, "07_effect_size_distributions"),
                        width = 14, height = 12)

cat("Section 7 complete.\n\n")

# =============================================================================
# SECTION 8: GO/KEGG ENRICHMENT ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 8: GO/KEGG ENRICHMENT ANALYSIS\n")
cat("================================================================================\n\n")

# Get hypermethylated genes (mC up in mutant)
hyper_genes <- mc_dmr %>%
  filter(significant & mod_difference > 0) %>%
  pull(gene) %>%
  unique()

cat(sprintf("Running enrichment on %d hypermethylated genes...\n", length(hyper_genes)))

# Convert gene symbols to Entrez IDs
tryCatch({
  gene_ids <- bitr(hyper_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

  if (nrow(gene_ids) > 10) {
    cat(sprintf("  Converted %d genes to Entrez IDs\n", nrow(gene_ids)))

    # GO Biological Process
    cat("Running GO Biological Process enrichment...\n")
    go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      p_go_bp <- dotplot(go_bp, showCategory = 15) +
        labs(title = "GO Biological Process Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_bp, file.path(OUTPUT_DIR, "08a_enrichment_go_bp"),
                              width = 10, height = 10)

      write.table(go_bp@result, file.path(TABLES_DIR, "enrichment_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO BP results\n")
    }

    # GO Cellular Component
    cat("Running GO Cellular Component enrichment...\n")
    go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      p_go_cc <- dotplot(go_cc, showCategory = 15) +
        labs(title = "GO Cellular Component Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_cc, file.path(OUTPUT_DIR, "08b_enrichment_go_cc"),
                              width = 10, height = 10)

      write.table(go_cc@result, file.path(TABLES_DIR, "enrichment_go_cc.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO CC results\n")
    }

    # GO Molecular Function
    cat("Running GO Molecular Function enrichment...\n")
    go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      p_go_mf <- dotplot(go_mf, showCategory = 15) +
        labs(title = "GO Molecular Function Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_mf, file.path(OUTPUT_DIR, "08c_enrichment_go_mf"),
                              width = 10, height = 10)

      write.table(go_mf@result, file.path(TABLES_DIR, "enrichment_go_mf.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO MF results\n")
    }

    # KEGG Pathway
    cat("Running KEGG pathway enrichment...\n")
    kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                       organism = "mmu",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1)

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      p_kegg <- dotplot(kegg, showCategory = 15) +
        labs(title = "KEGG Pathway Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_kegg, file.path(OUTPUT_DIR, "08d_enrichment_kegg"),
                              width = 10, height = 10)

      write.table(kegg@result, file.path(TABLES_DIR, "enrichment_kegg.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved KEGG results\n")
    }
  } else {
    cat("  Not enough genes converted for enrichment analysis\n")
  }
}, error = function(e) {
  cat(sprintf("  Enrichment analysis error: %s\n", e$message))
})

cat("Section 8 complete.\n\n")

# =============================================================================
# SECTION 9: SUMMARY STATISTICS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 9: SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

# Generate summary statistics
summary_text <- sprintf("
================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION ANALYSIS SUMMARY
================================================================================

Analysis Date: %s
Genome: GRCm38 (mm10)
Comparison: BAP1-KO Mutant vs Wildtype Control
Samples: 4 (2 Control, 2 Mutant)

--------------------------------------------------------------------------------
DIFFERENTIAL METHYLATION STATISTICS (CG Context, Gene Bodies)
--------------------------------------------------------------------------------

5mC (5-methylcytosine):
  Total genes tested: %d
  Significant (q < %.2f): %d (%.1f%%)
  Hypermethylated: %d
  Hypomethylated: %d
  Mean change (significant): %+.2f%%

5hmC (5-hydroxymethylcytosine):
  Total genes tested: %d
  Significant (q < %.2f): %d (%.1f%%)
  Increased: %d
  Decreased: %d
  Mean change (significant): %+.2f%%

--------------------------------------------------------------------------------
COORDINATED CHANGES
--------------------------------------------------------------------------------

Genes significant in both mC and hmC: %d
Genes with mC↑ / hmC↓ pattern: %d (%.1f%%)

This pattern indicates impaired TET-mediated active demethylation:
  5mC → (TET1/2/3) → 5hmC → 5fC → 5caC → C (unmethylated)

BAP1 loss → Increased H2AK119ub → Altered chromatin → TET access restricted

--------------------------------------------------------------------------------
TOP AFFECTED GENES (Coordinated Pattern)
--------------------------------------------------------------------------------

Gene      mC Change    hmC Change   mC q-value     hmC q-value
----      ---------    ----------   ----------     -----------
%s

--------------------------------------------------------------------------------
KEY FINDINGS
--------------------------------------------------------------------------------

1. PRIMARY AFFECTED REGIONS: Gene bodies (NOT promoters or CpG islands)
2. COORDINATED PATTERN: 5mC increase + 5hmC decrease at same loci
3. NEURONAL GENES: Top hits include synaptic/neuronal genes (Syt1, Cntnap2, Dlgap1)
4. NON-CpG: No significant changes in CHG/CHH contexts
5. MECHANISM: Consistent with blocked TET-mediated active demethylation

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------

Plots: %s
Tables: %s

Generated by: biomodal_visualizations.R
================================================================================
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(mc_dmr), Q_THRESHOLD, sum(mc_dmr$significant),
  100 * sum(mc_dmr$significant) / nrow(mc_dmr),
  sum(mc_dmr$significant & mc_dmr$mod_difference > 0),
  sum(mc_dmr$significant & mc_dmr$mod_difference < 0),
  mean(mc_dmr$mod_difference[mc_dmr$significant]) * 100,
  nrow(hmc_dmr), Q_THRESHOLD, sum(hmc_dmr$significant),
  100 * sum(hmc_dmr$significant) / nrow(hmc_dmr),
  sum(hmc_dmr$significant & hmc_dmr$mod_difference > 0),
  sum(hmc_dmr$significant & hmc_dmr$mod_difference < 0),
  mean(hmc_dmr$mod_difference[hmc_dmr$significant]) * 100,
  nrow(coordinated),
  sum(coordinated$coordinated_pattern),
  100 * mean(coordinated$coordinated_pattern),
  paste(sprintf("%-10s %+8.1f%%    %+9.1f%%   %.1e    %.1e",
                head(coordinated$gene[coordinated$coordinated_pattern], 7),
                head(coordinated$mc_diff[coordinated$coordinated_pattern], 7) * 100,
                head(coordinated$hmc_diff[coordinated$coordinated_pattern], 7) * 100,
                head(coordinated$mc_q[coordinated$coordinated_pattern], 7),
                head(coordinated$hmc_q[coordinated$coordinated_pattern], 7)),
        collapse = "\n"),
  OUTPUT_DIR,
  TABLES_DIR
)

# Write summary file
writeLines(summary_text, file.path(TABLES_DIR, "summary_statistics.txt"))
cat("Saved summary_statistics.txt\n")

# Print summary
cat(summary_text)

cat("Section 9 complete.\n\n")

# =============================================================================
# SECTION 10: ChIP-seq CHROMATIN STATE ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 10: ChIP-seq CHROMATIN STATE ANALYSIS\n")
cat("================================================================================\n\n")

cat("Loading ChIP-seq peak files...\n")

# Load all ChIP-seq peak files
chip_peaks <- list(
  ctcf = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  h3k27ac = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  h3k4me1 = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  h3k4me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)

# Check if all peaks loaded successfully
if (any(sapply(chip_peaks, is.null))) {
  cat("  WARNING: Some ChIP-seq peak files not found. Skipping Section 10.\n")
} else {
  cat("\nCreating GRanges from mC DMR data...\n")

  # Create GRanges from mC DMR data (significant genes only for analysis)
  mc_sig <- mc_dmr %>% dplyr::filter(significant)
  mc_gr <- dmr_to_granges(mc_sig)
  cat(sprintf("  %d significant mC DMRs\n", length(mc_gr)))

  # Compute TSS distances
  cat("\nComputing TSS distances...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss_gr <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d genes\n", length(genes)))

  nearest <- distanceToNearest(mc_gr, tss_gr)
  distance_to_tss <- rep(NA_real_, length(mc_gr))
  if (length(nearest) > 0) {
    distance_to_tss[queryHits(nearest)] <- mcols(nearest)$distance
  }

  cat(sprintf("  DMRs <= 2kb from TSS: %d (%.1f%%)\n",
              sum(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE),
              100 * mean(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE)))

  # Compute ChIP-seq overlaps
  cat("\nComputing ChIP-seq overlaps...\n")
  mc_overlaps <- compute_chip_overlaps(mc_gr, chip_peaks)

  cat("  Overlap summary:\n")
  cat(sprintf("    CTCF:     %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$CTCF_overlap), 100 * mean(mc_overlaps$CTCF_overlap)))
  cat(sprintf("    H3K27ac:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K27ac_overlap), 100 * mean(mc_overlaps$H3K27ac_overlap)))
  cat(sprintf("    H3K27me3: %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K27me3_overlap), 100 * mean(mc_overlaps$H3K27me3_overlap)))
  cat(sprintf("    H3K4me1:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K4me1_overlap), 100 * mean(mc_overlaps$H3K4me1_overlap)))
  cat(sprintf("    H3K4me3:  %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$H3K4me3_overlap), 100 * mean(mc_overlaps$H3K4me3_overlap)))
  cat(sprintf("    Bivalent: %d DMRs (%.1f%%)\n",
              sum(mc_overlaps$Bivalent_overlap), 100 * mean(mc_overlaps$Bivalent_overlap)))

  # Classify chromatin states
  cat("\nClassifying chromatin states (7 categories)...\n")
  mc_sig$chromatin_state <- classify_chromatin_state(mc_overlaps, distance_to_tss, TSS_THRESHOLD)

  # Add overlap columns to dataframe
  mc_sig$CTCF_overlap <- mc_overlaps$CTCF_overlap
  mc_sig$H3K27ac_overlap <- mc_overlaps$H3K27ac_overlap
  mc_sig$H3K27me3_overlap <- mc_overlaps$H3K27me3_overlap
  mc_sig$H3K4me1_overlap <- mc_overlaps$H3K4me1_overlap
  mc_sig$H3K4me3_overlap <- mc_overlaps$H3K4me3_overlap
  mc_sig$Bivalent_overlap <- mc_overlaps$Bivalent_overlap
  mc_sig$distance_to_tss <- distance_to_tss

  cat("\n  Chromatin state distribution (significant mC DMRs):\n")
  for (state in CHROMATIN_STATE_ORDER) {
    count <- sum(mc_sig$chromatin_state == state)
    pct <- 100 * count / nrow(mc_sig)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", state, count, pct))
  }

  # -----------------------------------------------------------------------
  # FIGURE 10a: Overall chromatin state distribution
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 10a: Chromatin state distribution...\n")

  state_summary <- mc_sig %>%
    group_by(chromatin_state, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction) %>%
    mutate(percentage = 100 * count / sum(count))

  # Bar chart
  p_10a_bar <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.3, size = 3) +
    facet_wrap(~direction, ncol = 2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.15), expand = c(0, 0)) +
    labs(
      title = "Chromatin State Distribution of Significant mC DMRs",
      subtitle = sprintf("Hypermethylated: n=%d | Hypomethylated: n=%d",
                         sum(mc_sig$direction == "Hypermethylated"),
                         sum(mc_sig$direction == "Hypomethylated")),
      x = "Chromatin State", y = "Percentage (%)"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Pie chart - overall
  state_overall <- mc_sig %>%
    group_by(chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  p_10a_pie <- ggplot(state_overall, aes(x = "", y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    labs(title = "Overall Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    ) +
    geom_text(aes(label = ifelse(percentage > 4, sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.5)

  p_10a <- p_10a_bar + p_10a_pie +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "Chromatin State Analysis of Differentially Methylated Genes",
      subtitle = "Based on ChIP-seq peak overlaps (Late timepoint)",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_10a, file.path(OUTPUT_DIR, "10a_chromatin_state_distribution"),
                          width = 14, height = 7)

  # -----------------------------------------------------------------------
  # FIGURE 10b: Chromatin state by methylation direction comparison
  # -----------------------------------------------------------------------
  cat("Creating Figure 10b: Chromatin state by direction...\n")

  # Stacked bar comparison
  p_10b_stacked <- ggplot(state_summary, aes(x = direction, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Chromatin State by Methylation Direction",
      x = "Direction", y = "Percentage (%)"
    ) +
    theme_biomodal()

  # Grouped bar comparison
  p_10b_grouped <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = COLORS$direction, name = "Direction") +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(
      title = "Comparison by Chromatin State",
      subtitle = "Do hypermethylated genes have specific chromatin signatures?",
      x = "Chromatin State", y = "Percentage (%)"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  p_10b <- p_10b_stacked / p_10b_grouped +
    plot_annotation(
      title = "Methylation Direction vs Chromatin State",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

  save_multiformat_ggplot(p_10b, file.path(OUTPUT_DIR, "10b_chromatin_by_methylation_direction"),
                          width = 12, height = 10)

  # -----------------------------------------------------------------------
  # FIGURE 10c: ChIP mark overlap heatmap
  # -----------------------------------------------------------------------
  cat("Creating Figure 10c: ChIP mark overlap heatmap...\n")

  # Calculate overlap percentages by direction
  mark_overlap <- mc_sig %>%
    group_by(direction) %>%
    summarise(
      CTCF = 100 * mean(CTCF_overlap),
      H3K27ac = 100 * mean(H3K27ac_overlap),
      H3K27me3 = 100 * mean(H3K27me3_overlap),
      H3K4me1 = 100 * mean(H3K4me1_overlap),
      H3K4me3 = 100 * mean(H3K4me3_overlap),
      Bivalent = 100 * mean(Bivalent_overlap),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -direction, names_to = "Mark", values_to = "Percentage")

  # Add overall row
  mark_overall <- mc_sig %>%
    summarise(
      CTCF = 100 * mean(CTCF_overlap),
      H3K27ac = 100 * mean(H3K27ac_overlap),
      H3K27me3 = 100 * mean(H3K27me3_overlap),
      H3K4me1 = 100 * mean(H3K4me1_overlap),
      H3K4me3 = 100 * mean(H3K4me3_overlap),
      Bivalent = 100 * mean(Bivalent_overlap)
    ) %>%
    mutate(direction = "All Significant") %>%
    pivot_longer(cols = -direction, names_to = "Mark", values_to = "Percentage")

  mark_overlap <- bind_rows(mark_overall, mark_overlap)
  mark_overlap$direction <- factor(mark_overlap$direction,
                                   levels = c("All Significant", "Hypermethylated", "Hypomethylated"))

  p_10c <- ggplot(mark_overlap, aes(x = Mark, y = direction, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 4) +
    scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                         midpoint = 50, name = "% Overlap", limits = c(0, 100)) +
    labs(
      title = "ChIP-seq Mark Overlap at Significant mC DMRs",
      subtitle = "% of DMRs overlapping each ChIP-seq mark",
      x = "ChIP-seq Mark", y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank()
    )

  save_multiformat_ggplot(p_10c, file.path(OUTPUT_DIR, "10c_chip_mark_overlap_heatmap"),
                          width = 10, height = 5)

  # -----------------------------------------------------------------------
  # FIGURE 10d: Chromatin state of coordinated genes
  # -----------------------------------------------------------------------
  cat("Creating Figure 10d: Coordinated genes chromatin state...\n")

  # Get chromatin states for coordinated genes (mC up / hmC down)
  coordinated_mc <- mc_sig %>%
    dplyr::filter(gene %in% coordinated$gene[coordinated$coordinated_pattern])

  if (nrow(coordinated_mc) > 0) {
    coord_state_summary <- coordinated_mc %>%
      group_by(chromatin_state) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(percentage = 100 * count / sum(count))

    p_10d_bar <- ggplot(coord_state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, count)), vjust = -0.1, size = 3) +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
      scale_y_continuous(limits = c(0, max(coord_state_summary$percentage) * 1.2), expand = c(0, 0)) +
      labs(
        title = "Chromatin State of Coordinated mC↑/hmC↓ Genes",
        subtitle = sprintf("n = %d genes with coordinated methylation changes", nrow(coordinated_mc)),
        x = "Chromatin State", y = "Percentage (%)"
      ) +
      theme_biomodal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    p_10d_pie <- ggplot(coord_state_summary, aes(x = "", y = percentage, fill = chromatin_state)) +
      geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
      labs(title = "Distribution") +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "right"
      ) +
      geom_text(aes(label = ifelse(percentage > 5, sprintf("%.0f%%", percentage), "")),
                position = position_stack(vjust = 0.5), size = 3)

    p_10d <- p_10d_bar + p_10d_pie +
      plot_layout(widths = c(2, 1)) +
      plot_annotation(
        title = "Chromatin State of Genes with Coordinated Methylation Changes",
        subtitle = "Genes showing mC increase + hmC decrease pattern",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
        )
      )

    save_multiformat_ggplot(p_10d, file.path(OUTPUT_DIR, "10d_coordinated_genes_chromatin"),
                            width = 13, height = 6)
  }

  # -----------------------------------------------------------------------
  # FIGURE 10e: Top genes with chromatin annotations
  # -----------------------------------------------------------------------
  cat("Creating Figure 10e: Top genes with chromatin annotations...\n")

  # Get top 20 coordinated genes with chromatin state
  top_coord_chrom <- coordinated %>%
    dplyr::filter(coordinated_pattern) %>%
    head(20) %>%
    left_join(mc_sig %>% dplyr::select(gene, chromatin_state), by = "gene")

  if (nrow(top_coord_chrom) > 0) {
    top_coord_chrom <- top_coord_chrom %>%
      mutate(gene = factor(gene, levels = rev(unique(gene))))

    p_10e <- ggplot(top_coord_chrom, aes(x = gene, y = combined_effect * 100, fill = chromatin_state)) +
      geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
      geom_text(aes(label = as.character(chromatin_state)),
                hjust = -0.1, size = 3, color = "grey30") +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
      coord_flip() +
      labs(
        title = "Top 20 Coordinated Genes with Chromatin State Annotation",
        subtitle = "Sorted by combined effect (|mC change| + |hmC change|)",
        x = "", y = "Combined Effect (%)"
      ) +
      theme_biomodal() +
      theme(legend.position = "bottom")

    save_multiformat_ggplot(p_10e, file.path(OUTPUT_DIR, "10e_top_genes_chromatin_annotation"),
                            width = 12, height = 10)
  }

  # -----------------------------------------------------------------------
  # Save chromatin state annotations
  # -----------------------------------------------------------------------
  cat("\nSaving chromatin state annotation tables...\n")

  # Full annotation table
  write.table(mc_sig, file.path(TABLES_DIR, "dmr_chromatin_state_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: dmr_chromatin_state_annotation.tsv\n")

  # Summary table
  chromatin_summary <- mc_sig %>%
    group_by(chromatin_state, direction) %>%
    summarise(
      count = n(),
      mean_mc_change = mean(mod_difference) * 100,
      median_mc_change = median(mod_difference) * 100,
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = direction,
                values_from = c(count, mean_mc_change, median_mc_change))

  write.table(chromatin_summary, file.path(TABLES_DIR, "chromatin_state_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: chromatin_state_summary.tsv\n")
}

cat("\nSection 10 complete.\n\n")

cat("\n")
cat("================================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Tables directory: %s\n", TABLES_DIR))
cat("\nGenerated files:\n")
system(sprintf("ls -la %s/*.pdf 2>/dev/null | head -20", OUTPUT_DIR))
cat("\n")
