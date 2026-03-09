# biomodal/downstream/scripts/preprocess_tet_ko_bigwig.R
# HPC preprocessing: Extract gene-body 5mC/5hmC signal from TET-KO bigWigs
#
# Computes per-gene 5mC and 5hmC from pre-processed bigWigs for both WT and
# TET triple-KO conditions from GSE166423. The bigWigs already contain separated
# methylCpG (5mC) and hydroxymethylCpG (5hmC) signal — no BS-OxBS subtraction needed.
# Derives demethylation efficiency ratio = 5hmC/(5mC+5hmC) per gene.
#
# GSE166423 bigWig files (P56 Purkinje neurons):
#   purkinje_adult_methyl_CpG_merged.bw    — WT 5mC (merged across 2 reps)
#   purkinje_adult_hydroxy_CpG_merged.bw   — WT 5hmC (merged across 2 reps)
#   3228_tet123_pcp2_CRE_POS_bs_rep2.methylCpG.bw       — TET-TKO 5mC (n=1)
#   3228_tet123_pcp2_CRE_POS_bs_rep2.hydroxymethylCpG.bw — TET-TKO 5hmC (n=1)
#
# All 4 bigWig files are REQUIRED.
#
# Usage (on Expanse HPC):
#   Rscript scripts/preprocess_tet_ko_bigwig.R \
#     --mc-wt /path/to/purkinje_adult_methyl_CpG_merged.bw \
#     --hmc-wt /path/to/purkinje_adult_hydroxy_CpG_merged.bw \
#     --mc-ko /path/to/3228_tet123_pcp2_CRE_POS_bs_rep2.methylCpG.bw \
#     --hmc-ko /path/to/3228_tet123_pcp2_CRE_POS_bs_rep2.hydroxymethylCpG.bw \
#     --output data/tet_ko_gene_signal.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list(
    mc_wt    = NULL,
    hmc_wt   = NULL,
    mc_ko    = NULL,
    hmc_ko   = NULL,
    gene_bed = "/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz",
    output   = "data/tet_ko_gene_signal.tsv"
  )

  usage_msg <- paste0(
    "Usage: Rscript preprocess_tet_ko_bigwig.R \\\n",
    "  --mc-wt <bw> --hmc-wt <bw> \\\n",
    "  --mc-ko <bw> --hmc-ko <bw> \\\n",
    "  [--gene-bed <bed.gz>] [--output <tsv>]"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--mc-wt" && i < length(args)) {
      params$mc_wt <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--hmc-wt" && i < length(args)) {
      params$hmc_wt <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--mc-ko" && i < length(args)) {
      params$mc_ko <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--hmc-ko" && i < length(args)) {
      params$hmc_ko <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--gene-bed" && i < length(args)) {
      params$gene_bed <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      params$output <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\n%s", args[i], usage_msg))
    }
  }

  # All 4 bigWig files are required
  if (is.null(params$mc_wt))  stop("--mc-wt is required\n", usage_msg)
  if (is.null(params$hmc_wt)) stop("--hmc-wt is required\n", usage_msg)
  if (is.null(params$mc_ko))  stop("--mc-ko is required\n", usage_msg)
  if (is.null(params$hmc_ko)) stop("--hmc-ko is required\n", usage_msg)

  params
}

params <- parse_args()

cat("================================================================================\n")
cat("TET-KO BIGWIG SIGNAL PREPROCESSING (GSE166423)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Inputs (4 pre-processed bigWig files):\n")
cat("  5mC WT:  ", params$mc_wt, "\n")
cat("  5hmC WT: ", params$hmc_wt, "\n")
cat("  5mC KO:  ", params$mc_ko, "\n")
cat("  5hmC KO: ", params$hmc_ko, "\n")
cat("Gene BED:  ", params$gene_bed, "\n")
cat("Output:    ", params$output, "\n\n")

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

cat("Validating all input files exist...\n")
stopifnot("5mC WT bigwig not found"  = file.exists(params$mc_wt))
stopifnot("5hmC WT bigwig not found" = file.exists(params$hmc_wt))
stopifnot("5mC KO bigwig not found"  = file.exists(params$mc_ko))
stopifnot("5hmC KO bigwig not found" = file.exists(params$hmc_ko))
stopifnot("Gene BED not found"       = file.exists(params$gene_bed))
cat("  All 5 input files validated.\n\n")

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

cat("\n--- Extracting signal from all 4 bigWigs ---\n")

mc_wt  <- compute_mean_signal(params$mc_wt,  gene_gr, "5mC WT (purkinje_adult_methyl_CpG_merged)")
hmc_wt <- compute_mean_signal(params$hmc_wt, gene_gr, "5hmC WT (purkinje_adult_hydroxy_CpG_merged)")
mc_ko  <- compute_mean_signal(params$mc_ko,  gene_gr, "5mC KO (tet123_methylCpG)")
hmc_ko <- compute_mean_signal(params$hmc_ko, gene_gr, "5hmC KO (tet123_hydroxymethylCpG)")

# =============================================================================
# COMPUTE DEMETHYLATION RATIO
# =============================================================================

cat("\n--- Computing demethylation efficiency ratio ---\n")

cat(sprintf("  WT 5mC: median=%.4f, range=[%.4f, %.4f]\n",
            median(mc_wt), min(mc_wt), max(mc_wt)))
cat(sprintf("  WT 5hmC: median=%.4f, range=[%.4f, %.4f]\n",
            median(hmc_wt), min(hmc_wt), max(hmc_wt)))
cat(sprintf("  KO 5mC: median=%.4f, range=[%.4f, %.4f]\n",
            median(mc_ko), min(mc_ko), max(mc_ko)))
cat(sprintf("  KO 5hmC: median=%.4f, range=[%.4f, %.4f]\n",
            median(hmc_ko), min(hmc_ko), max(hmc_ko)))

# ratio = 5hmC / (5mC + 5hmC), minimum denominator 0.01
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
cat(sprintf("  Delta ratio: median=%+.4f, range=[%+.4f, %+.4f]\n",
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
cat(sprintf("  WT source:           merged across 2 replicates (from GEO)\n"))
cat(sprintf("  KO source:           single replicate (Pcp2-Cre conditional TKO)\n"))
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
if (median(hmc_ko) > median(hmc_wt)) {
  cat("  WARNING: KO 5hmC > WT 5hmC — unexpected for TET-KO\n")
}
if (mean(delta_ratio_tet < 0) < 0.5) {
  cat("  WARNING: Fewer than 50%% of genes show decreased ratio — unexpected for TET-KO\n")
}

cat("\nPreprocessing complete.\n")
cat("\nNext steps:\n")
cat(sprintf("  scp expanse:$(pwd)/%s biomodal/downstream/data/\n", params$output))
cat("  cd biomodal/downstream/\n")
cat("  Rscript scripts/viz_sections/section_26_tet_ko_comparison.R\n")
