# abc/scripts/preprocess_k119ub_enhancer_signal.R
# HPC preprocessing: Extract per-replicate K119ub signal from 8 bigWigs at ABC enhancer regions
#
# Computes mean K119ub signal over unique enhancer regions (from ABC delta pairs)
# for 4 ctrl + 4 mut replicates, applies median-ratio normalization across all 8
# samples, and outputs a portable TSV for local correlation analysis.
#
# Mirrors scripts/preprocess_k119ub_anchor_signal.R but adapted for:
#   - Per-replicate bigWigs (8 files) instead of merged (2 files)
#   - ABC enhancer coordinates (chr, start, end) instead of loop anchors
#
# Usage (on Expanse HPC):
#   Rscript abc/scripts/preprocess_k119ub_enhancer_signal.R \
#     --abc-pairs /expanse/lustre/projects/csd940/zalibhai/abc/results/delta_abc_all_pairs.tsv \
#     --bw-dir /expanse/lustre/projects/csd940/zalibhai/loop-class/heatmaps/ \
#     --output abc/results/k119ub_enhancer_signal.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    abc_pairs = "/expanse/lustre/projects/csd940/zalibhai/abc/results/delta_abc_all_pairs.tsv",
    bw_dir    = "/expanse/lustre/projects/csd940/zalibhai/loop-class/heatmaps/",
    output    = "abc/results/k119ub_enhancer_signal.tsv"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--abc-pairs" && i < length(args)) {
      defaults$abc_pairs <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--bw-dir" && i < length(args)) {
      defaults$bw_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      defaults$output <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf(
        "Unknown argument: %s\nUsage: Rscript preprocess_k119ub_enhancer_signal.R --abc-pairs <tsv> --bw-dir <dir> --output <tsv>",
        args[i]
      ))
    }
  }
  defaults
}

params <- parse_args()

cat("================================================================================\n")
cat("K119UB ENHANCER SIGNAL PREPROCESSING (PER-REPLICATE)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("ABC pairs:  ", params$abc_pairs, "\n")
cat("BigWig dir: ", params$bw_dir, "\n")
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

stopifnot(file.exists(params$abc_pairs))
stopifnot(dir.exists(params$bw_dir))

# Define per-replicate bigWig files
bw_samples <- c(
  ctrl1 = "H2AK119ubCtrl1.bw",
  ctrl2 = "H2AK119ubCtrl2.bw",
  ctrl3 = "H2AK119ubCtrl3.bw",
  ctrl4 = "H2AK119ubCtrl4.bw",
  mut1  = "H2AK119ubMut1.bw",
  mut2  = "H2AK119ubMut2.bw",
  mut3  = "H2AK119ubMut3.bw",
  mut4  = "H2AK119ubMut4.bw"
)

bw_paths <- file.path(params$bw_dir, bw_samples)
names(bw_paths) <- names(bw_samples)

cat("Validating bigWig files:\n")
for (nm in names(bw_paths)) {
  if (!file.exists(bw_paths[nm])) {
    stop(sprintf("BigWig not found: %s (%s)", bw_paths[nm], nm))
  }
  cat(sprintf("  [OK] %s: %s\n", nm, bw_paths[nm]))
}
cat("\n")

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# EXTRACT UNIQUE ENHANCER REGIONS
# =============================================================================

cat("Loading ABC pairs and extracting unique enhancers...\n")
abc <- read.table(params$abc_pairs, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d E-G pairs\n", nrow(abc)))

# Unique enhancer coordinates
enhancers <- unique(abc[, c("chr", "start", "end")])
cat(sprintf("  Unique enhancer regions: %d\n", nrow(enhancers)))

# Create GRanges
enh_gr <- GRanges(
  seqnames = enhancers$chr,
  ranges = IRanges(start = enhancers$start, end = enhancers$end)
)

# Filter to standard chromosomes
std_chroms <- intersect(paste0("chr", c(1:19, "X")), seqlevels(enh_gr))
enh_gr <- keepSeqlevels(enh_gr, std_chroms, pruning.mode = "coarse")
cat(sprintf("  After filtering to standard chroms: %d enhancers\n", length(enh_gr)))

# Build enhancer IDs for output
enh_ids <- paste0(seqnames(enh_gr), ":", start(enh_gr), "-", end(enh_gr))

# =============================================================================
# COMPUTE MEAN SIGNAL PER ENHANCER REGION (PER SAMPLE)
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

cat("\n--- Per-replicate signal extraction ---\n")
signal_mat <- matrix(NA_real_, nrow = length(enh_gr), ncol = length(bw_paths))
colnames(signal_mat) <- names(bw_paths)

for (samp in names(bw_paths)) {
  signal_mat[, samp] <- compute_mean_signal(bw_paths[samp], enh_gr, samp)
}

# =============================================================================
# BETWEEN-SAMPLE NORMALIZATION (median-ratio across all 8 samples)
# =============================================================================

cat("\n--- Median-ratio normalization across 8 samples ---\n")

# Compute per-sample median of non-zero signals
sample_medians <- apply(signal_mat, 2, function(x) median(x[x > 0]))
cat("  Per-sample medians (non-zero):\n")
for (nm in names(sample_medians)) {
  cat(sprintf("    %-8s: %.6f\n", nm, sample_medians[nm]))
}

# Reference: geometric mean of medians
ref_median <- exp(mean(log(sample_medians)))
cat(sprintf("  Reference (geometric mean of medians): %.6f\n", ref_median))

# Scale factors
scale_factors <- ref_median / sample_medians
cat("  Scale factors:\n")
for (nm in names(scale_factors)) {
  cat(sprintf("    %-8s: %.6f\n", nm, scale_factors[nm]))
}

# Apply normalization
signal_norm <- signal_mat
for (samp in names(scale_factors)) {
  signal_norm[, samp] <- signal_mat[, samp] * scale_factors[samp]
}

# =============================================================================
# COMPUTE CONDITION MEANS AND LOG2FC
# =============================================================================

cat("\n--- Computing condition means and log2FC ---\n")

ctrl_cols <- c("ctrl1", "ctrl2", "ctrl3", "ctrl4")
mut_cols <- c("mut1", "mut2", "mut3", "mut4")

ctrl_mean <- rowMeans(signal_norm[, ctrl_cols])
mut_mean <- rowMeans(signal_norm[, mut_cols])

# Adaptive pseudocount
all_nonzero <- c(signal_norm[signal_norm > 0])
pseudocount <- max(quantile(all_nonzero, 0.05), 0.01)
cat(sprintf("  Adaptive pseudocount: %.6f\n", pseudocount))

# Log2FC and signal classification
log2fc <- rep(NA_real_, length(ctrl_mean))
signal_class <- rep("no_signal", length(ctrl_mean))

both_above <- (ctrl_mean > pseudocount) & (mut_mean > pseudocount)
log2fc[both_above] <- log2((mut_mean[both_above] + pseudocount) / (ctrl_mean[both_above] + pseudocount))
signal_class[both_above] <- "quantifiable"

one_above <- xor(ctrl_mean > pseudocount, mut_mean > pseudocount)
log2fc[one_above] <- log2((mut_mean[one_above] + pseudocount) / (ctrl_mean[one_above] + pseudocount))
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
  enhancer_id  = enh_ids,
  chr          = as.character(seqnames(enh_gr)),
  start        = start(enh_gr),
  end          = end(enh_gr),
  ctrl1        = round(signal_norm[, "ctrl1"], 6),
  ctrl2        = round(signal_norm[, "ctrl2"], 6),
  ctrl3        = round(signal_norm[, "ctrl3"], 6),
  ctrl4        = round(signal_norm[, "ctrl4"], 6),
  mut1         = round(signal_norm[, "mut1"], 6),
  mut2         = round(signal_norm[, "mut2"], 6),
  mut3         = round(signal_norm[, "mut3"], 6),
  mut4         = round(signal_norm[, "mut4"], 6),
  ctrl_mean    = round(ctrl_mean, 6),
  mut_mean     = round(mut_mean, 6),
  log2fc       = round(log2fc, 6),
  signal_class = signal_class,
  stringsAsFactors = FALSE
)

# =============================================================================
# WRITE OUTPUT
# =============================================================================

cat(sprintf("\nWriting output to: %s\n", params$output))
write.table(output_df, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d enhancers written\n", nrow(output_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Input E-G pairs:       %d\n", nrow(abc)))
cat(sprintf("  Unique enhancers:      %d\n", nrow(output_df)))
cat(sprintf("  BigWig samples:        %d (4 ctrl + 4 mut)\n", length(bw_paths)))
cat(sprintf("  Pseudocount:           %.6f\n", pseudocount))
cat(sprintf("  Reference median:      %.6f\n", ref_median))

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
