# biomodal/downstream/scripts/viz_sections/_shared_config.R
# Shared configuration, packages, helper functions, and data loading
# Sourced by each individual section script

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base paths - run from downstream/ directory
BASE_DIR <- getwd()
REPO_ROOT <- normalizePath(file.path(BASE_DIR, "../.."))  # mariner_hi-c/

# Data file paths
DATA_PATHS <- list(
  # Gene body DMRs (run-5, 8 samples with sex covariate)
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_mc_control__mutant_20260402_191818.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_hmc_control__mutant_20260402_191818.bed"),
  bioqc_json = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/BioQC_20260402_185215/biological_qc_report_8_samples_20260402_185215.json"),
  upstream_csv = file.path(BASE_DIR, "../upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv"),
  # Regional mC DMR files
  cpg_islands_mc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260402_191006/DMR_mc_control__mutant_20260402_191006.bed"),
  cpg_shores_mc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260402_191531/DMR_mc_control__mutant_20260402_191531.bed"),
  cpg_shelves_mc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260402_191227/DMR_mc_control__mutant_20260402_191227.bed"),
  promoters_mc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260402_192045/DMR_mc_control__mutant_20260402_192045.bed"),
  tss_mc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260402_192302/DMR_mc_control__mutant_20260402_192302.bed"),
  # Regional hmC DMR files
  cpg_islands_hmc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260402_191006/DMR_hmc_control__mutant_20260402_191006.bed"),
  cpg_shores_hmc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260402_191531/DMR_hmc_control__mutant_20260402_191531.bed"),
  cpg_shelves_hmc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260402_191227/DMR_hmc_control__mutant_20260402_191227.bed"),
  promoters_hmc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260402_192045/DMR_hmc_control__mutant_20260402_192045.bed"),
  tss_hmc = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260402_192302/DMR_hmc_control__mutant_20260402_192302.bed")
)

# CG feature extraction paths (per-sample regional fractions)
EXTRACT_PATHS <- list(
  extract_dir = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260402_190326"),
  mc_regional_frac = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260402_190326/Extract_mc_regional-frac_20260402_190326.tsv.gz"),
  hmc_regional_frac = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260402_190326/Extract_hmc_regional-frac_20260402_190326.tsv.gz")
)

# CHG context DMR paths (gene body only — minimal signal expected)
CHG_DATA_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_222845/DMR_mc_control__mutant_20260402_222845.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_222845/DMR_hmc_control__mutant_20260402_222845.bed")
)

# CHH context DMR paths (gene body only — minimal signal expected)
CHH_DATA_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHH/Results/gencode.vM25.mouse.genes.annotation/DMR_20260403_162803/DMR_mc_control__mutant_20260403_162803.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHH/Results/gencode.vM25.mouse.genes.annotation/DMR_20260403_162803/DMR_hmc_control__mutant_20260403_162803.bed")
)

# CHG/CHH per-sample feature extraction paths (regional fractions at gene bodies)
CHG_EXTRACT_PATHS <- list(
  mc_regional_frac = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260402_211353/Extract_mc_regional-frac_20260402_211353.tsv.gz")
)
CHH_EXTRACT_PATHS <- list(
  mc_regional_frac = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CHH/Results/gencode.vM25.mouse.genes.annotation/Extract_20260403_151546/Extract_mc_regional-frac_20260403_151546.tsv.gz")
)

# Histone mark peak file paths (from peaks/beds/)
# Use Late timepoint peaks to match the adult BAP1-KO analysis
CHIP_PEAK_FILES <- list(
  ctcf     = file.path(REPO_ROOT, "peaks/CTCF.bed"),
  h3k27ac  = file.path(REPO_ROOT, "peaks/beds/H3K27acCerebellumLate2.bed"),
  h3k27me3 = file.path(REPO_ROOT, "peaks/beds/H3K27me3CerebellumLate1.bed"),
  h3k4me1  = file.path(REPO_ROOT, "peaks/beds/H3K4me1CerebellumLate1.bed"),
  h3k4me3  = file.path(REPO_ROOT, "peaks/beds/H3K4me3CerebellumLate2.bed"),
  bivalent = file.path(REPO_ROOT, "peaks/beds/Bivalent_Cerebellum_Late.bed")
)

# MeCP2 differential binding files (from peaks/mecp2/)
MECP2_FILES <- list(
  annotated = file.path(REPO_ROOT, "peaks/mecp2/MeCP2_annotated.txt"),
  up      = file.path(REPO_ROOT, "peaks/mecp2/MeCP2_up.bed"),
  down    = file.path(REPO_ROOT, "peaks/mecp2/MeCP2_down.bed"),
  ctrl_bw = "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Ctrl.bw",
  mut_bw  = "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Mut.bw"
)

# Merged methylation BigWigs (group-averaged, from sdsc/bigwigs/methylation/merged/)
METHYLATION_BIGWIGS <- list(
  cg_mc_ctrl  = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_ctrl.bw",
  cg_mc_mut   = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_mut.bw",
  cg_hmc_ctrl = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_ctrl.bw",
  cg_hmc_mut  = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_mut.bw"
)

# ATAC-seq differential and consensus peak files (from peaks/atac_seq/)
ATAC_FILES <- list(
  up   = file.path(REPO_ROOT, "peaks/atac_seq/ATAC_up.bed"),
  down = file.path(REPO_ROOT, "peaks/atac_seq/ATAC_down.bed"),
  consensus_ctrl = file.path(REPO_ROOT, "peaks/atac_seq/consensus_control.bed"),
  consensus_mut  = file.path(REPO_ROOT, "peaks/atac_seq/consensus_mutant.bed")
)

# H2AK119ub condition-specific peak files (from peaks/intersect/)
K119UB_FILES <- list(
  ctrl = file.path(REPO_ROOT, "peaks/intersect/P51_K119ub_ctrl_intersect.bed"),
  mut  = file.path(REPO_ROOT, "peaks/intersect/P51_K119ub_mut_intersect.bed")
)

# H3K27ac condition-specific peak files (from peaks/intersect/)
H3K27AC_FILES <- list(
  ctrl = file.path(REPO_ROOT, "peaks/intersect/P60_K27ac_ctrl_intersect.bed"),
  mut  = file.path(REPO_ROOT, "peaks/intersect/P60_K27ac_mut_intersect.bed")
)

# DiffBind quantitative differential binding results (4 marks)
DIFFBIND_FILES <- list(
  atac   = file.path(REPO_ROOT, "peaks/diffbind/ATAC_allATAC_diffbind_results_summit_appended_ap.txt"),
  k27ac  = file.path(REPO_ROOT, "peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt"),
  k27me3 = file.path(REPO_ROOT, "peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt"),
  k119ub = file.path(REPO_ROOT, "peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt")
)

# H3K36me2 DiffBind results (MACS2 narrow peaks, 12 samples: 6 ctrl + 6 mut)
H3K36ME2_FILES <- list(
  master = file.path(BASE_DIR, "peaks/03_2026_H3K36me2 files/260221_H3K36me2_macs2_narrow_diffbind_withsummits_results copy.txt"),
  up     = file.path(BASE_DIR, "peaks/03_2026_H3K36me2 files/diffbind  - using macs2 narrow/allreps_H3K36me2_up.bed"),
  down   = file.path(BASE_DIR, "peaks/03_2026_H3K36me2 files/diffbind  - using macs2 narrow/allreps_H3K36me2_down.bed")
)

# H3K36me3 DiffBind results (SEACR stringent peaks, 12 samples: 6 ctrl + 6 mut)
H3K36ME3_FILES <- list(
  master_annotated = file.path(BASE_DIR, "peaks/03_2026_H3K36me3 files/260221_seacr_stringent_summits_diffbind_results.txt"),
  sig_all          = file.path(BASE_DIR, "peaks/03_2026_H3K36me3 files/All Regions - DiffBind with summits/diffbind_H3K36me3_all_results.txt"),
  sig_pc           = file.path(BASE_DIR, "peaks/03_2026_H3K36me3 files/Protein Coding Regions ONLY DiffBind/diffbind_H3K36me3_PC_results.txt"),
  up               = file.path(BASE_DIR, "peaks/03_2026_H3K36me3 files/All Regions - DiffBind no summits/all regions H3K36me3 figures/allreps_H3K36me3_up.bed"),
  down             = file.path(BASE_DIR, "peaks/03_2026_H3K36me3 files/All Regions - DiffBind no summits/all regions H3K36me3 figures/allreps_H3K36me3_down.bed")
)

# Hi-C loop annotation files (from mariner pipeline)
LOOP_FILES <- list(
  late = file.path(REPO_ROOT, "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv")
)

# Output directory
OUTPUT_DIR <- file.path(BASE_DIR, "plots/visualizations")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Key genes to highlight
KEY_GENES <- c("Syt1", "Zbtb20", "Trpm3", "Epha3", "Mcu", "Cntnap2", "Lpp", "Dlgap1", "Arhgap26", "Cdh8")

# Statistical thresholds
Q_THRESHOLD <- 0.05

# Sample metadata (8 samples: batch 1 + batch 2 replicates)
SAMPLES <- data.frame(
  sample_id = c("evoC-Bap1-ctrl-F", "evoC-Bap1-ctrl-M", "evoC-Bap1-mut-F", "evoC-Bap1-mut-M",
                "evoC-Bap1-ctrl-F-B2", "evoC-Bap1-ctrl-M-B2", "evoC-Bap1-mut-F-B2", "evoC-Bap1-mut-M-B2"),
  condition = c("Control", "Control", "Mutant", "Mutant",
                "Control", "Control", "Mutant", "Mutant"),
  sex = c("Female", "Male", "Female", "Male",
          "Female", "Male", "Female", "Male"),
  batch = c(1, 1, 1, 1, 2, 2, 2, 2),
  short_name = c("ctrl-F", "ctrl-M", "mut-F", "mut-M",
                 "ctrl-F-B2", "ctrl-M-B2", "mut-F-B2", "mut-M-B2"),
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
              "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k36me2 = c("H3K36me2 Gained" = "#E6AB02", "H3K36me2 Lost" = "#66A61E",
               "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k36me3 = c("H3K36me3 Gained" = "#D95F02", "H3K36me3 Lost" = "#1B9E77",
               "Shared" = "grey70", "Not Significant" = "grey70"),
  h3k36_combined = c("H3K36me2" = "#E6AB02", "H3K36me3" = "#D95F02", "Both" = "#984EA3")
)

# Chromatin state classification (consistent with annotate_loops_extended.R)
CHROMATIN_STATE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                           "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Unmarked")

CHROMATIN_STATE_COLORS <- c(
  "Active_Promoter" = "#e41a1c",     # Red - active transcription
  "Repressed_Promoter" = "#756bb1",  # Purple - Polycomb-silenced promoter
  "Bivalent_Promoter" = "#984ea3",   # Magenta - developmental poised
  "Polycomb" = "#4daf4a",            # Green - distal repressive
  "Active_Enhancer" = "#377eb8",     # Blue - active enhancer
  "Poised_Enhancer" = "#ff7f00",     # Orange - primed enhancer
  "Unmarked" = "#999999"              # Gray - no histone mark overlap
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
  # Histone mark analysis packages (Section 10)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# Source multi-format output utility
util_path <- file.path(BASE_DIR, "scripts/utils/multi_format_output.R")
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

dedup_by_gene <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_min(dmr_qvalue, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

# =============================================================================
# HISTONE MARK HELPER FUNCTIONS (Section 10)
# =============================================================================

#' Load histone mark peaks from BED file as GRanges
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

#' Compute histone mark overlaps for DMR GRanges
#' @param dmr_gr GRanges with DMR coordinates
#' @param chip_peaks List of GRanges for each histone mark
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

#' Load DiffBind results with flexible column schema
#' Handles both summit-appended format (Summit_Chr/Summit_Start/Summit_End)
#' and raw DiffBind format (seqnames/start/end) as used by H3K36me2/me3 files.
#' @param filepath Path to DiffBind TSV
#' @param mark_name Display name for logging
#' @param fdr_threshold FDR threshold for significance counts
#' @return data.frame with standardized columns: Chr, Start, End, Fold, FDR, p.value
load_diffbind_flex <- function(filepath, mark_name = "Mark", fdr_threshold = 0.05) {
  stopifnot(file.exists(filepath))
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", fill = TRUE)

  # Detect schema and standardize coordinate columns
  if ("Summit_Chr" %in% colnames(df)) {
    df$Chr <- df$Summit_Chr
    df$Start <- df$Summit_Start
    df$End <- df$Summit_End
  } else if ("seqnames" %in% colnames(df)) {
    df$Chr <- df$seqnames
    df$Start <- df$start
    df$End <- df$end
  } else if ("Chr" %in% colnames(df)) {
    # Already standard
  } else {
    stop(sprintf("Unrecognized DiffBind column schema in %s. Expected Summit_Chr, seqnames, or Chr.", filepath))
  }

  # Normalize p-value column name (some files use p-value with hyphen)
  if ("p.value" %in% colnames(df)) {
    # already correct
  } else if ("p-value" %in% colnames(df)) {
    df$p.value <- df[["p-value"]]
  }

  # Validate required columns
  required <- c("Chr", "Start", "End", "Fold", "FDR")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns in %s: %s", filepath, paste(missing, collapse = ", ")))
  }

  n_up   <- sum(df$FDR < fdr_threshold & df$Fold > 0, na.rm = TRUE)
  n_down <- sum(df$FDR < fdr_threshold & df$Fold < 0, na.rm = TRUE)
  cat(sprintf("  %s: %d peaks (%d sig up, %d sig down at FDR<%.2f)\n",
              mark_name, nrow(df), n_up, n_down, fdr_threshold))
  df
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
  chromatin_state <- rep("Unmarked", n)

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

  # 7. Unmarked: Default (no histone marks)
  return(factor(chromatin_state, levels = CHROMATIN_STATE_ORDER))
}

# =============================================================================
# regioneReloaded BUG FIX: chooseHclustMet crashes with <=2 rows
# =============================================================================

#' Patch regioneReloaded::chooseHclustMet for small matrices
#'
#' chooseHclustMet computes cophenetic correlation to pick the best hclust method.
#' With <=2 rows, dist() returns a single value, cor() returns NA (undefined for
#' n=1), and which(NA == max(NA)) yields integer(0), crashing the [[ indexer.
#' This patch short-circuits to use the first specified method when nrow <= 2.
#'
#' Call this AFTER library(regioneReloaded) in permutation sections (34-36).
patch_chooseHclustMet <- function() {
  safe_fn <- function(GM, scale = TRUE, vecMet = NULL, distHC = "euclidean") {
    if (scale == TRUE) GM <- scale(GM)
    if (is.null(vecMet)) {
      vecMet <- c("complete", "average", "single", "ward.D2",
                   "median", "centroid", "mcquitty")
    }
    mat_dist <- stats::dist(x = GM, method = distHC)

    # With <=2 rows, cophenetic correlation is undefined (single distance value).
    # Just use the first specified method directly.
    if (nrow(GM) <= 2) {
      model <- stats::hclust(d = mat_dist, method = vecMet[1])
      methods::show(paste0("method for hclustering (<=2 elements, skipping cophenetic): ", vecMet[1]))
      return(model)
    }

    # Original logic for >2 rows
    resMetList <- lapply(seq_along(vecMet), FUN = function(i, mat_dist, vecMet) {
      stats::hclust(d = mat_dist, method = vecMet[[i]])
    }, mat_dist, vecMet)
    names(resMetList) <- vecMet
    resMetVec <- unlist(lapply(seq_along(resMetList), FUN = function(i, mat_dist, resMetList) {
      stats::cor(x = mat_dist, stats::cophenetic(resMetList[[i]]))
    }, mat_dist, resMetList))
    names(resMetVec) <- vecMet
    name_model <- vecMet[which(resMetVec == max(resMetVec))]
    if (length(name_model) > 1) name_model <- name_model[1]
    model <- resMetList[[name_model]]
    methods::show(paste0("method selected for hclustering: ", name_model))
    methods::show(resMetVec)
    return(model)
  }
  assignInNamespace("chooseHclustMet", safe_fn, ns = "regioneReloaded")
  cat("  Patched regioneReloaded::chooseHclustMet for small-matrix safety.\n")
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("================================================================================\n")
cat("LOADING DATA\n")
cat("================================================================================\n\n")

# Load gene body DMR data (deduplicate: one row per gene, lowest q-value wins)
cat("Loading gene body mC DMRs...\n")
mc_dmr_raw <- load_dmr_bed(DATA_PATHS$mc_dmr)
mc_dmr <- dedup_by_gene(mc_dmr_raw)
cat(sprintf("  Loaded %d rows, deduplicated to %d unique genes\n",
            nrow(mc_dmr_raw), nrow(mc_dmr)))

cat("Loading gene body hmC DMRs...\n")
hmc_dmr_raw <- load_dmr_bed(DATA_PATHS$hmc_dmr)
hmc_dmr <- dedup_by_gene(hmc_dmr_raw)
cat(sprintf("  Loaded %d rows, deduplicated to %d unique genes\n",
            nrow(hmc_dmr_raw), nrow(hmc_dmr)))

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
