# scripts/preprocess_k119ub_anchor_signal.R
# HPC preprocessing: Extract K119ub continuous signal from merged bigwigs per loop anchor
#
# Computes mean K119ub signal over unique loop anchor regions for ctrl/mut,
# applies median-ratio normalization, and outputs a portable TSV.
#
# Mirrors biomodal/downstream/scripts/preprocess_k119ub_bigwig.R but for
# loop anchors instead of gene bodies/promoters.
#
# Usage (on Expanse HPC):
#   Rscript scripts/preprocess_k119ub_anchor_signal.R \
#     --ctrl /expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_ctrl_merged.bw \
#     --mut /expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_mut_merged.bw \
#     --loops 25042-late_outputs/merged_loops/characterized_loops.tsv \
#     --output data/k119ub_anchor_signal.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    ctrl   = "/expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_ctrl_merged.bw",
    mut    = "/expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_mut_merged.bw",
    loops  = "25042-late_outputs/merged_loops/characterized_loops.tsv",
    output = "data/k119ub_anchor_signal.tsv"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--ctrl" && i < length(args)) {
      defaults$ctrl <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--mut" && i < length(args)) {
      defaults$mut <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--loops" && i < length(args)) {
      defaults$loops <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      defaults$output <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\nUsage: Rscript preprocess_k119ub_anchor_signal.R --ctrl <bw> --mut <bw> --loops <tsv> --output <tsv>", args[i]))
    }
  }
  defaults
}

params <- parse_args()

cat("================================================================================\n")
cat("K119UB ANCHOR SIGNAL PREPROCESSING\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Ctrl bigwig:", params$ctrl, "\n")
cat("Mut bigwig: ", params$mut, "\n")
cat("Loops file: ", params$loops, "\n")
cat("Output:     ", params$output, "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

if (!file.exists(params$ctrl)) stop("Ctrl bigwig not found: ", params$ctrl)
if (!file.exists(params$mut))  stop("Mut bigwig not found: ", params$mut)
if (!file.exists(params$loops)) stop("Loops file not found: ", params$loops)

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# EXTRACT UNIQUE ANCHOR REGIONS
# =============================================================================

cat("Loading loop data and extracting unique anchors...\n")
loops <- read.table(params$loops, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d loops\n", nrow(loops)))

# Extract anchor1 and anchor2 coordinates
anchor1 <- data.frame(
  chr   = loops$anchor1_chr,
  start = loops$anchor1_start,
  end   = loops$anchor1_end,
  stringsAsFactors = FALSE
)
anchor2 <- data.frame(
  chr   = loops$anchor2_chr,
  start = loops$anchor2_start,
  end   = loops$anchor2_end,
  stringsAsFactors = FALSE
)

# Combine and deduplicate
all_anchors <- unique(rbind(anchor1, anchor2))
cat(sprintf("  Unique anchor regions: %d\n", nrow(all_anchors)))

# Create GRanges
anchor_gr <- GRanges(
  seqnames = all_anchors$chr,
  ranges = IRanges(start = all_anchors$start, end = all_anchors$end)
)

# Filter to standard chromosomes (chr1-19, chrX; chrY absent in loop data)
std_chroms <- intersect(paste0("chr", c(1:19, "X")), seqlevels(anchor_gr))
anchor_gr <- keepSeqlevels(anchor_gr, std_chroms, pruning.mode = "coarse")
cat(sprintf("  After filtering to standard chroms: %d anchors\n", length(anchor_gr)))

# Reduce overlapping anchors for signal computation
anchor_reduced <- reduce(anchor_gr)
cat(sprintf("  After reduce (merge overlapping): %d non-overlapping regions\n", length(anchor_reduced)))

# Keep original (non-reduced) for output mapping
anchor_ids <- paste0(seqnames(anchor_gr), ":", start(anchor_gr), "-", end(anchor_gr))

# =============================================================================
# COMPUTE MEAN SIGNAL PER ANCHOR REGION
# =============================================================================

compute_mean_signal <- function(bw_path, regions, label) {
  cat(sprintf("  Computing mean signal from %s (%d regions)...\n", label, length(regions)))

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

cat("\n--- Anchor signal extraction ---\n")
ctrl_raw <- compute_mean_signal(params$ctrl, anchor_gr, "ctrl anchors")
mut_raw  <- compute_mean_signal(params$mut,  anchor_gr, "mut anchors")

# =============================================================================
# BETWEEN-CONDITION NORMALIZATION (median ratio)
# =============================================================================

cat("\n--- Between-condition normalization ---\n")

ctrl_nonzero <- ctrl_raw[ctrl_raw > 0]
mut_nonzero  <- mut_raw[mut_raw > 0]

median_ctrl <- median(ctrl_nonzero)
median_mut  <- median(mut_nonzero)
scale_factor <- median_ctrl / median_mut

cat(sprintf("  Ctrl median (non-zero): %.6f\n", median_ctrl))
cat(sprintf("  Mut median (non-zero):  %.6f\n", median_mut))
cat(sprintf("  Scale factor (ctrl/mut): %.6f\n", scale_factor))

mut_norm <- mut_raw * scale_factor

# =============================================================================
# COMPUTE LOG2FC WITH ADAPTIVE PSEUDOCOUNT
# =============================================================================

cat("\n--- Computing log2 fold changes ---\n")

all_nonzero <- c(ctrl_raw[ctrl_raw > 0], mut_norm[mut_norm > 0])
pseudocount <- max(quantile(all_nonzero, 0.05), 0.01)
cat(sprintf("  Adaptive pseudocount: %.6f\n", pseudocount))

log2fc <- rep(NA_real_, length(ctrl_raw))
signal_class <- rep("no_signal", length(ctrl_raw))

both_above <- (ctrl_raw > pseudocount) & (mut_norm > pseudocount)
log2fc[both_above] <- log2((mut_norm[both_above] + pseudocount) / (ctrl_raw[both_above] + pseudocount))
signal_class[both_above] <- "quantifiable"

one_above <- xor(ctrl_raw > pseudocount, mut_norm > pseudocount)
log2fc[one_above] <- log2((mut_norm[one_above] + pseudocount) / (ctrl_raw[one_above] + pseudocount))
signal_class[one_above] <- "one_condition"

cat("  Signal classes:\n")
for (cls in c("quantifiable", "one_condition", "no_signal")) {
  n <- sum(signal_class == cls)
  cat(sprintf("    %-15s: %d (%.1f%%)\n", cls, n, 100 * n / length(signal_class)))
}

# =============================================================================
# BUILD OUTPUT TABLE
# =============================================================================

cat("\n--- Building output table ---\n")

output_df <- data.frame(
  anchor_id    = anchor_ids,
  chr          = as.character(seqnames(anchor_gr)),
  start        = start(anchor_gr),
  end          = end(anchor_gr),
  ctrl_signal  = round(ctrl_raw, 6),
  mut_signal   = round(mut_norm, 6),
  log2fc       = round(log2fc, 6),
  signal_class = signal_class,
  stringsAsFactors = FALSE
)

# =============================================================================
# WRITE OUTPUT
# =============================================================================

cat(sprintf("\nWriting output to: %s\n", params$output))
write.table(output_df, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d anchors written\n", nrow(output_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Input loops:           %d\n", nrow(loops)))
cat(sprintf("  Unique anchors:        %d\n", nrow(output_df)))
cat(sprintf("  Pseudocount:           %.6f\n", pseudocount))
cat(sprintf("  Normalization factor:  %.6f\n", scale_factor))

q <- log2fc[signal_class == "quantifiable"]
if (length(q) > 0) {
  cat(sprintf("\n  Log2FC (quantifiable, n=%d):\n", length(q)))
  cat(sprintf("    Median: %+.4f\n", median(q)))
  cat(sprintf("    Mean:   %+.4f\n", mean(q)))
  cat(sprintf("    SD:     %.4f\n", sd(q)))
  cat(sprintf("    Range:  [%+.4f, %+.4f]\n", min(q), max(q)))
  cat(sprintf("    Positive (mut > ctrl): %d (%.1f%%)\n",
              sum(q > 0), 100 * sum(q > 0) / length(q)))
}

cat("\nPreprocessing complete.\n")
