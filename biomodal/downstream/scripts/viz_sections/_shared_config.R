# biomodal/downstream/scripts/viz_sections/_shared_config.R
# Shared configuration, packages, helper functions, and data loading
# Sourced by each individual section script

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base paths - run from downstream/ directory
BASE_DIR <- getwd()

# Data file paths
DATA_PATHS <- list(
  # Gene body DMRs (run-3, deep-seq)
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260221_190322/DMR_mc_control__mutant_20260221_190322.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260221_190322/DMR_hmc_control__mutant_20260221_190322.bed"),
  bioqc_json = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/BioQC_20260221_184206/biological_qc_report_4_samples_20260221_184206.json"),
  upstream_csv = file.path(BASE_DIR, "../upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv"),
  # Regional mC DMR files
  cpg_islands_mc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260221_185728/DMR_mc_control__mutant_20260221_185728.bed"),
  cpg_shores_mc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260221_190122/DMR_mc_control__mutant_20260221_190122.bed"),
  cpg_shelves_mc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260221_185914/DMR_mc_control__mutant_20260221_185914.bed"),
  promoters_mc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260221_190510/DMR_mc_control__mutant_20260221_190510.bed"),
  tss_mc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260221_190652/DMR_mc_control__mutant_20260221_190652.bed"),
  # Regional hmC DMR files
  cpg_islands_hmc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260221_185728/DMR_hmc_control__mutant_20260221_185728.bed"),
  cpg_shores_hmc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260221_190122/DMR_hmc_control__mutant_20260221_190122.bed"),
  cpg_shelves_hmc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260221_185914/DMR_hmc_control__mutant_20260221_185914.bed"),
  promoters_hmc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260221_190510/DMR_hmc_control__mutant_20260221_190510.bed"),
  tss_hmc = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260221_190652/DMR_hmc_control__mutant_20260221_190652.bed")
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

# MeCP2 differential binding files (from peaks/mecp2/)
MECP2_FILES <- list(
  annotated = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/mecp2/MeCP2_annotated.txt",
  up   = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/mecp2/MeCP2_up.bed",
  down = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/mecp2/MeCP2_down.bed"
)

# ATAC-seq differential and consensus peak files (from peaks/atac_seq/)
ATAC_FILES <- list(
  up   = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/atac_seq/ATAC_up.bed",
  down = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/atac_seq/ATAC_down.bed",
  consensus_ctrl = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/atac_seq/consensus_control.bed",
  consensus_mut  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/atac_seq/consensus_mutant.bed"
)

# H2AK119ub condition-specific peak files (from peaks/intersect/)
K119UB_FILES <- list(
  ctrl = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/intersect/P51_K119ub_ctrl_intersect.bed",
  mut  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/intersect/P51_K119ub_mut_intersect.bed"
)

# H3K27ac condition-specific peak files (from peaks/intersect/)
H3K27AC_FILES <- list(
  ctrl = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/intersect/P60_K27ac_ctrl_intersect.bed",
  mut  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/intersect/P60_K27ac_mut_intersect.bed"
)

# DiffBind quantitative differential binding results (4 marks)
DIFFBIND_FILES <- list(
  atac   = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/diffbind/ATAC_allATAC_diffbind_results_summit_appended_ap.txt",
  k27ac  = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt",
  k27me3 = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt",
  k119ub = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt"
)

# Hi-C loop annotation files (from mariner pipeline)
LOOP_FILES <- list(
  late = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"
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
  significant = c("Significant" = "#E41A1C", "Not Significant" = "grey70"),
  mecp2 = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3", "Not Significant" = "grey70"),
  atac = c("ATAC Up" = "#E6AB02", "ATAC Down" = "#66A61E", "Not Significant" = "grey70"),
  k119ub = c("K119ub Gained" = "#756BB1", "K119ub Lost" = "#74C476",
             "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k27ac = c("H3K27ac Gained" = "#FF7F00", "H3K27ac Lost" = "#1F78B4",
              "Shared" = "grey70", "Not Significant" = "grey70")
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
