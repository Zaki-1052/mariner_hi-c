# biomodal/downstream/scripts/preprocess_tet_ko_bigwig.R
# HPC preprocessing: Extract gene-body 5mC/5hmC signal from TET-KO bigWigs
#
# Computes per-gene 5mC and 5hmC from BS-seq and OxBS-seq bigWigs for both
# WT (2 replicates, averaged) and TET triple-KO (1 replicate) from GSE166423.
# Derives demethylation efficiency ratio = 5hmC/(5mC+5hmC) per gene.
#
# GSE166423 sample structure (P56 Purkinje neurons):
#   WT OxBS-Seq:  GSM5070963, GSM5070964 (2 replicates)
#   WT BS-Seq:    GSM5070969, GSM5070970 (2 replicates)
#   TKO OxBS-Seq: GSM5070971 (1 replicate, Pcp2-Cre conditional)
#   TKO BS-Seq:   GSM5070972 (1 replicate, Pcp2-Cre conditional)
#
# Chemistry:
#   5mC = OxBS signal (oxidative bisulfite converts 5hmC -> U, preserves 5mC)
#   5hmC = BS - OxBS (standard bisulfite preserves both 5mC and 5hmC as C)
#
# All 6 bigWig files are REQUIRED. WT replicates are averaged per gene.
#
# Usage (on Expanse HPC):
#   Rscript scripts/preprocess_tet_ko_bigwig.R \
#     --bs-wt-rep1 /path/to/BS_WT_rep1.bw \
#     --bs-wt-rep2 /path/to/BS_WT_rep2.bw \
#     --oxbs-wt-rep1 /path/to/OxBS_WT_rep1.bw \
#     --oxbs-wt-rep2 /path/to/OxBS_WT_rep2.bw \
#     --bs-ko /path/to/BS_KO.bw \
#     --oxbs-ko /path/to/OxBS_KO.bw \
#     --output data/tet_ko_gene_signal.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list(
    bs_wt_rep1   = NULL,
    bs_wt_rep2   = NULL,
    oxbs_wt_rep1 = NULL,
    oxbs_wt_rep2 = NULL,
    bs_ko        = NULL,
    oxbs_ko      = NULL,
    gene_bed     = "/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz",
    output       = "data/tet_ko_gene_signal.tsv"
  )

  usage_msg <- paste0(
    "Usage: Rscript preprocess_tet_ko_bigwig.R \\\n",
    "  --bs-wt-rep1 <bw> --bs-wt-rep2 <bw> \\\n",
    "  --oxbs-wt-rep1 <bw> --oxbs-wt-rep2 <bw> \\\n",
    "  --bs-ko <bw> --oxbs-ko <bw> \\\n",
    "  [--gene-bed <bed.gz>] [--output <tsv>]"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--bs-wt-rep1" && i < length(args)) {
      params$bs_wt_rep1 <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--bs-wt-rep2" && i < length(args)) {
      params$bs_wt_rep2 <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--oxbs-wt-rep1" && i < length(args)) {
      params$oxbs_wt_rep1 <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--oxbs-wt-rep2" && i < length(args)) {
      params$oxbs_wt_rep2 <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--bs-ko" && i < length(args)) {
      params$bs_ko <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--oxbs-ko" && i < length(args)) {
      params$oxbs_ko <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--gene-bed" && i < length(args)) {
      params$gene_bed <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      params$output <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\n%s", args[i], usage_msg))
    }
  }

  # All 6 bigWig files are required — no fallbacks
  if (is.null(params$bs_wt_rep1))   stop("--bs-wt-rep1 is required\n", usage_msg)
  if (is.null(params$bs_wt_rep2))   stop("--bs-wt-rep2 is required\n", usage_msg)
  if (is.null(params$oxbs_wt_rep1)) stop("--oxbs-wt-rep1 is required\n", usage_msg)
  if (is.null(params$oxbs_wt_rep2)) stop("--oxbs-wt-rep2 is required\n", usage_msg)
  if (is.null(params$bs_ko))        stop("--bs-ko is required\n", usage_msg)
  if (is.null(params$oxbs_ko))      stop("--oxbs-ko is required\n", usage_msg)

  params
}

params <- parse_args()

cat("================================================================================\n")
cat("TET-KO BIGWIG SIGNAL PREPROCESSING\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Inputs (6 bigWig files):\n")
cat("  BS-WT rep1:   ", params$bs_wt_rep1, "\n")
cat("  BS-WT rep2:   ", params$bs_wt_rep2, "\n")
cat("  OxBS-WT rep1: ", params$oxbs_wt_rep1, "\n")
cat("  OxBS-WT rep2: ", params$oxbs_wt_rep2, "\n")
cat("  BS-KO:        ", params$bs_ko, "\n")
cat("  OxBS-KO:      ", params$oxbs_ko, "\n")
cat("Gene BED:       ", params$gene_bed, "\n")
cat("Output:         ", params$output, "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating all 6 bigWig files exist...\n")
stopifnot("BS-WT rep1 not found"   = file.exists(params$bs_wt_rep1))
stopifnot("BS-WT rep2 not found"   = file.exists(params$bs_wt_rep2))
stopifnot("OxBS-WT rep1 not found" = file.exists(params$oxbs_wt_rep1))
stopifnot("OxBS-WT rep2 not found" = file.exists(params$oxbs_wt_rep2))
stopifnot("BS-KO not found"        = file.exists(params$bs_ko))
stopifnot("OxBS-KO not found"      = file.exists(params$oxbs_ko))
stopifnot("Gene BED not found"     = file.exists(params$gene_bed))
cat("  All 7 input files validated.\n\n")

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# LOAD GENE COORDINATES
# =============================================================================

cat("Loading gene coordinates from GENCODE vM25 BED...\n")

gene_bed <- read.table(gzfile(params$gene_bed), header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "#")
cat(sprintf("  Read %d rows from BED file\n", nrow(gene_bed)))

# BED format: chr, start (0-based), end, gene_name, ...
gene_gr <- GRanges(
  seqnames = gene_bed$V1,
  ranges = IRanges(start = gene_bed$V2 + 1, end = gene_bed$V3),
  gene = gene_bed$V4
)

# Standard chromosomes only, minimum gene length, deduplicate
std_chroms <- paste0("chr", c(1:19, "X"))
gene_gr <- keepSeqlevels(gene_gr, std_chroms, pruning.mode = "coarse")
gene_gr <- gene_gr[width(gene_gr) >= 500]

dups <- duplicated(mcols(gene_gr)$gene)
gene_gr <- gene_gr[!dups]
names(gene_gr) <- mcols(gene_gr)$gene

cat(sprintf("  After filtering: %d genes (standard chroms, >=500bp, unique)\n",
            length(gene_gr)))

# =============================================================================
# COMPUTE MEAN SIGNAL PER REGION
# =============================================================================

compute_mean_signal <- function(bw_path, regions, label) {
  cat(sprintf("  Extracting: %s (%d regions)...\n", label, length(regions)))

  bw <- import.bw(bw_path, which = regions)
  cov <- coverage(bw, weight = "score")

  ranges_by_chr <- split(ranges(regions), seqnames(regions))
  shared_chrs <- intersect(names(cov), names(ranges_by_chr))
  v <- Views(cov[shared_chrs], ranges_by_chr[shared_chrs])
  means <- viewMeans(v)

  result <- numeric(length(regions))
  chr_names <- as.character(seqnames(regions))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }

  cat(sprintf("    Range: [%.4f, %.4f], Median: %.4f, Non-zero: %d / %d\n",
              min(result), max(result), median(result),
              sum(result > 0), length(result)))
  result
}

cat("\n--- Extracting signal from all 6 bigWigs ---\n")

bs_wt_r1   <- compute_mean_signal(params$bs_wt_rep1,   gene_gr, "BS-WT rep1 (GSM5070969)")
bs_wt_r2   <- compute_mean_signal(params$bs_wt_rep2,   gene_gr, "BS-WT rep2 (GSM5070970)")
oxbs_wt_r1 <- compute_mean_signal(params$oxbs_wt_rep1, gene_gr, "OxBS-WT rep1 (GSM5070963)")
oxbs_wt_r2 <- compute_mean_signal(params$oxbs_wt_rep2, gene_gr, "OxBS-WT rep2 (GSM5070964)")
bs_ko      <- compute_mean_signal(params$bs_ko,        gene_gr, "BS-KO (GSM5070972)")
oxbs_ko    <- compute_mean_signal(params$oxbs_ko,      gene_gr, "OxBS-KO (GSM5070971)")

# Average WT replicates
bs_wt_signal   <- (bs_wt_r1 + bs_wt_r2) / 2
oxbs_wt_signal <- (oxbs_wt_r1 + oxbs_wt_r2) / 2

cat("\n  WT replicate correlations (sanity check):\n")
bs_wt_cor <- cor(bs_wt_r1, bs_wt_r2, method = "spearman")
oxbs_wt_cor <- cor(oxbs_wt_r1, oxbs_wt_r2, method = "spearman")
cat(sprintf("    BS-WT rep1 vs rep2:   Spearman rho = %.4f\n", bs_wt_cor))
cat(sprintf("    OxBS-WT rep1 vs rep2: Spearman rho = %.4f\n", oxbs_wt_cor))
stopifnot("BS-WT replicates poorly correlated (rho < 0.8)" = bs_wt_cor > 0.8)
stopifnot("OxBS-WT replicates poorly correlated (rho < 0.8)" = oxbs_wt_cor > 0.8)
cat("  WT replicates averaged.\n")

# =============================================================================
# COMPUTE 5mC, 5hmC, AND RATIO
# =============================================================================

cat("\n--- Computing 5mC, 5hmC, and demethylation ratio ---\n")

# 5mC = OxBS signal
mc_wt <- oxbs_wt_signal
mc_ko <- oxbs_ko

# 5hmC = BS - OxBS (floor at 0)
hmc_wt <- pmax(bs_wt_signal - oxbs_wt_signal, 0)
hmc_ko <- pmax(bs_ko - oxbs_ko, 0)

cat(sprintf("  WT 5mC: median=%.4f, range=[%.4f, %.4f]\n",
            median(mc_wt), min(mc_wt), max(mc_wt)))
cat(sprintf("  WT 5hmC: median=%.4f, range=[%.4f, %.4f]\n",
            median(hmc_wt), min(hmc_wt), max(hmc_wt)))
cat(sprintf("  KO 5mC: median=%.4f, range=[%.4f, %.4f]\n",
            median(mc_ko), min(mc_ko), max(mc_ko)))
cat(sprintf("  KO 5hmC: median=%.4f, range=[%.4f, %.4f]\n",
            median(hmc_ko), min(hmc_ko), max(hmc_ko)))
cat(sprintf("  5hmC negative values floored to 0: WT=%d, KO=%d\n",
            sum(bs_wt_signal - oxbs_wt_signal < 0),
            sum(bs_ko - oxbs_ko < 0)))

# Demethylation efficiency ratio = 5hmC / (5mC + 5hmC)
MIN_DENOM <- 0.01

total_wt <- mc_wt + hmc_wt
total_ko <- mc_ko + hmc_ko

ratio_wt <- hmc_wt / pmax(total_wt, MIN_DENOM)
ratio_ko <- hmc_ko / pmax(total_ko, MIN_DENOM)

# Set ratio to 0 where denominator < minimum threshold
ratio_wt[total_wt < MIN_DENOM] <- 0
ratio_ko[total_ko < MIN_DENOM] <- 0

delta_ratio_tet <- ratio_ko - ratio_wt

cat(sprintf("  WT ratio: median=%.4f, range=[%.4f, %.4f]\n",
            median(ratio_wt), min(ratio_wt), max(ratio_wt)))
cat(sprintf("  KO ratio: median=%.4f, range=[%.4f, %.4f]\n",
            median(ratio_ko), min(ratio_ko), max(ratio_ko)))
cat(sprintf("  Delta ratio: median=%.4f, range=[%.4f, %.4f]\n",
            median(delta_ratio_tet), min(delta_ratio_tet), max(delta_ratio_tet)))
cat(sprintf("  Genes with decreased ratio: %d / %d (%.1f%%)\n",
            sum(delta_ratio_tet < 0), length(delta_ratio_tet),
            100 * mean(delta_ratio_tet < 0)))

# =============================================================================
# BUILD OUTPUT TABLE
# =============================================================================

cat("\n--- Building output table ---\n")

gene_coords <- as.data.frame(gene_gr)

output_df <- data.frame(
  gene       = mcols(gene_gr)$gene,
  chr        = gene_coords$seqnames,
  start      = gene_coords$start,
  end        = gene_coords$end,
  mc_wt      = round(mc_wt, 6),
  mc_ko      = round(mc_ko, 6),
  hmc_wt     = round(hmc_wt, 6),
  hmc_ko     = round(hmc_ko, 6),
  ratio_wt   = round(ratio_wt, 6),
  ratio_ko   = round(ratio_ko, 6),
  delta_ratio_tet = round(delta_ratio_tet, 6),
  stringsAsFactors = FALSE
)

n_before <- nrow(output_df)
output_df <- output_df[!is.na(output_df$gene) & output_df$gene != "", ]
cat(sprintf("  Removed %d genes without names (%d remaining)\n",
            n_before - nrow(output_df), nrow(output_df)))

# =============================================================================
# WRITE OUTPUT
# =============================================================================

cat(sprintf("\nWriting output to: %s\n", params$output))
write.table(output_df, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d genes written\n", nrow(output_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Total genes:         %d\n", nrow(output_df)))
cat(sprintf("  WT replicates:       2 (averaged)\n"))
cat(sprintf("  KO replicates:       1\n"))
cat(sprintf("  BS-WT rep cor:       %.4f (Spearman)\n", bs_wt_cor))
cat(sprintf("  OxBS-WT rep cor:     %.4f (Spearman)\n", oxbs_wt_cor))
cat(sprintf("  Min denominator:     %.2f\n", MIN_DENOM))
cat(sprintf("  WT 5mC median:       %.4f\n", median(mc_wt)))
cat(sprintf("  WT 5hmC median:      %.4f\n", median(hmc_wt)))
cat(sprintf("  WT ratio median:     %.4f\n", median(ratio_wt)))
cat(sprintf("  KO 5mC median:       %.4f\n", median(mc_ko)))
cat(sprintf("  KO 5hmC median:      %.4f\n", median(hmc_ko)))
cat(sprintf("  KO ratio median:     %.4f\n", median(ratio_ko)))
cat(sprintf("  Delta ratio median:  %+.4f\n", median(delta_ratio_tet)))
cat(sprintf("  %% decreased:         %.1f%%\n", 100 * mean(delta_ratio_tet < 0)))

cat("\n  Sanity checks:\n")
if (median(ratio_wt) < 0.05 || median(ratio_wt) > 0.6) {
  cat(sprintf("  WARNING: WT ratio (%.4f) outside expected range [0.05, 0.6]\n",
              median(ratio_wt)))
}
if (median(hmc_wt) < 0) {
  cat("  WARNING: Negative median WT 5hmC — check BS/OxBS file assignment\n")
}
if (mean(delta_ratio_tet < 0) < 0.5) {
  cat("  WARNING: Fewer than 50%% of genes show decreased ratio — expected for TET-KO\n")
}
cat("  (Check manually: TET-KO should show near-complete 5hmC loss)\n")

cat("\nPreprocessing complete.\n")
cat("\nNext steps:\n")
cat(sprintf("  scp expanse:$(pwd)/%s biomodal/downstream/data/\n", params$output))
cat("  cd biomodal/downstream/\n")
cat("  Rscript scripts/viz_sections/section_26_tet_ko_comparison.R\n")
