# biomodal/downstream/scripts/tad_methylation/add_mecp2_tad_signal.R
# Augment the existing TAD methylation signal TSV with MeCP2 ChIP-seq metrics.
#
# Runs LOCALLY (not on HPC) using merged MeCP2 BigWigs.
# Adds 13 columns: 10 signal metrics (same pattern as methylation contexts)
# plus 3 peak overlap counts from DiffBind MeCP2 peaks.
#
# Usage (from downstream/ directory):
#   Rscript scripts/tad_methylation/add_mecp2_tad_signal.R \
#     --input  data/tad_methylation_signal_late.tsv \
#     --output data/tad_methylation_signal_late.tsv \
#     --ctrl-bw /Users/zakiralibhai/sdsc/bigwigs/MeCP2Ctrl.bw \
#     --mut-bw  /Users/zakiralibhai/sdsc/bigwigs/MeCP2Mut.bw \
#     --peaks-up  ../../peaks/mecp2/MeCP2_up.bed \
#     --peaks-down ../../peaks/mecp2/MeCP2_down.bed

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input     = "data/tad_methylation_signal_late.tsv",
    output    = "data/tad_methylation_signal_late.tsv",
    ctrl_bw   = "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Ctrl.bw",
    mut_bw    = "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Mut.bw",
    peaks_up  = "../../peaks/mecp2/MeCP2_up.bed",
    peaks_down = "../../peaks/mecp2/MeCP2_down.bed"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" && i < length(args)) {
      defaults$input <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      defaults$output <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--ctrl-bw" && i < length(args)) {
      defaults$ctrl_bw <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--mut-bw" && i < length(args)) {
      defaults$mut_bw <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--peaks-up" && i < length(args)) {
      defaults$peaks_up <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--peaks-down" && i < length(args)) {
      defaults$peaks_down <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s", args[i]))
    }
  }
  defaults
}

params <- parse_args()

cat("================================================================================\n")
cat("ADD MeCP2 SIGNAL TO TAD METHYLATION TSV\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input:     ", params$input, "\n")
cat("Output:    ", params$output, "\n")
cat("Ctrl BW:   ", params$ctrl_bw, "\n")
cat("Mut BW:    ", params$mut_bw, "\n")
cat("Peaks Up:  ", params$peaks_up, "\n")
cat("Peaks Down:", params$peaks_down, "\n\n")

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

if (!file.exists(params$input))   stop("Input TSV not found: ", params$input)
if (!file.exists(params$ctrl_bw)) stop("Ctrl BigWig not found: ", params$ctrl_bw)
if (!file.exists(params$mut_bw))  stop("Mut BigWig not found: ", params$mut_bw)
if (!file.exists(params$peaks_up))   stop("Peaks-up BED not found: ", params$peaks_up)
if (!file.exists(params$peaks_down)) stop("Peaks-down BED not found: ", params$peaks_down)

# =============================================================================
# LOAD EXISTING TSV AND RECONSTRUCT DOMAIN GRanges
# =============================================================================

cat("Loading existing TAD signal TSV...\n")
tad <- read.table(params$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  %d domains x %d columns\n", nrow(tad), ncol(tad)))

if ("mecp2_ctrl_mean" %in% names(tad)) {
  cat("  WARNING: MeCP2 columns already exist. They will be overwritten.\n")
  drop_cols <- grep("^mecp2_", names(tad))
  tad <- tad[, -drop_cols]
  cat(sprintf("  Dropped %d existing MeCP2 columns. Now %d columns.\n",
              length(drop_cols), ncol(tad)))
}

domain_gr <- GRanges(
  seqnames = tad$chr,
  ranges = IRanges(start = tad$start, end = tad$end)
)
cat(sprintf("  Reconstructed %d domain GRanges\n\n", length(domain_gr)))

# =============================================================================
# SIGNAL EXTRACTION FUNCTIONS
# =============================================================================

compute_mean_signal <- function(bw_path, regions, label) {
  cat(sprintf("    %s (%d regions)...", label, length(regions)))

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

  cat(sprintf(" median=%.6f, non-zero=%d/%d\n",
              median(result, na.rm = TRUE),
              sum(result > 0, na.rm = TRUE), length(result)))
  result
}

compute_domain_cv <- function(bw_path, domain_gr, bin_size = 10000) {
  bw <- import.bw(bw_path, which = domain_gr)
  cov <- coverage(bw, weight = "score")

  cv_vals <- numeric(length(domain_gr))
  chr_names <- as.character(seqnames(domain_gr))

  for (i in seq_along(domain_gr)) {
    dom <- domain_gr[i]
    chr <- chr_names[i]
    dom_start <- start(dom)
    dom_end <- end(dom)

    tile_starts <- seq(dom_start, dom_end - 1, by = bin_size)
    tile_ends <- pmin(tile_starts + bin_size - 1, dom_end)
    tiles_ir <- IRanges(start = tile_starts, end = tile_ends)

    if (!(chr %in% names(cov))) {
      cv_vals[i] <- NA
      next
    }

    cov_chr <- cov[[chr]]
    max_pos <- length(cov_chr)
    tiles_ir <- tiles_ir[end(tiles_ir) <= max_pos]
    if (length(tiles_ir) < 2) {
      cv_vals[i] <- NA
      next
    }

    v <- Views(cov_chr, tiles_ir)
    bin_means <- as.numeric(viewMeans(v))
    domain_mean <- mean(bin_means, na.rm = TRUE)
    domain_sd <- sd(bin_means, na.rm = TRUE)

    cv_vals[i] <- if (!is.na(domain_mean) && domain_mean > 0) {
      domain_sd / domain_mean
    } else {
      NA
    }
  }
  cv_vals
}

compute_within_domain_var <- function(bw_path, domain_gr, bin_size = 10000) {
  bw <- import.bw(bw_path, which = domain_gr)
  cov <- coverage(bw, weight = "score")

  var_vals <- numeric(length(domain_gr))
  chr_names <- as.character(seqnames(domain_gr))

  for (i in seq_along(domain_gr)) {
    dom <- domain_gr[i]
    chr <- chr_names[i]
    dom_start <- start(dom)
    dom_end <- end(dom)

    tile_starts <- seq(dom_start, dom_end - 1, by = bin_size)
    tile_ends <- pmin(tile_starts + bin_size - 1, dom_end)
    tiles_ir <- IRanges(start = tile_starts, end = tile_ends)

    if (!(chr %in% names(cov))) {
      var_vals[i] <- NA
      next
    }

    cov_chr <- cov[[chr]]
    max_pos <- length(cov_chr)
    tiles_ir <- tiles_ir[end(tiles_ir) <= max_pos]
    if (length(tiles_ir) < 2) {
      var_vals[i] <- NA
      next
    }

    v <- Views(cov_chr, tiles_ir)
    bin_means <- as.numeric(viewMeans(v))
    var_vals[i] <- var(bin_means, na.rm = TRUE)
  }
  var_vals
}

compute_boundary_insulation <- function(bw_path, domain_gr, flank_kb = 50) {
  flank_size <- flank_kb * 1000

  left_inside <- GRanges(
    seqnames = seqnames(domain_gr),
    ranges = IRanges(
      start = start(domain_gr),
      end = pmin(start(domain_gr) + flank_size - 1, end(domain_gr))
    )
  )
  left_outside <- GRanges(
    seqnames = seqnames(domain_gr),
    ranges = IRanges(
      start = pmax(1, start(domain_gr) - flank_size),
      end = start(domain_gr) - 1
    )
  )
  right_inside <- GRanges(
    seqnames = seqnames(domain_gr),
    ranges = IRanges(
      start = pmax(start(domain_gr), end(domain_gr) - flank_size + 1),
      end = end(domain_gr)
    )
  )
  right_outside <- GRanges(
    seqnames = seqnames(domain_gr),
    ranges = IRanges(
      start = end(domain_gr) + 1,
      end = end(domain_gr) + flank_size
    )
  )

  li <- compute_mean_signal(bw_path, left_inside, "left_inside")
  lo <- compute_mean_signal(bw_path, left_outside, "left_outside")
  ri <- compute_mean_signal(bw_path, right_inside, "right_inside")
  ro <- compute_mean_signal(bw_path, right_outside, "right_outside")

  left_insulation <- abs(li - lo)
  right_insulation <- abs(ri - ro)
  mean_insulation <- (left_insulation + right_insulation) / 2

  list(left = left_insulation, right = right_insulation, mean = mean_insulation)
}

# =============================================================================
# COMPUTE MeCP2 SIGNAL METRICS
# =============================================================================

cat("================================================================================\n")
cat("MeCP2 SIGNAL EXTRACTION\n")
cat("================================================================================\n\n")

cat("--- Mean signal per domain ---\n")
cat("  Ctrl:\n")
mecp2_ctrl_mean <- compute_mean_signal(params$ctrl_bw, domain_gr, "MeCP2 ctrl")
cat("  Mut:\n")
mecp2_mut_mean <- compute_mean_signal(params$mut_bw, domain_gr, "MeCP2 mut")

all_nonzero <- c(mecp2_ctrl_mean[mecp2_ctrl_mean > 0],
                 mecp2_mut_mean[mecp2_mut_mean > 0])
pseudocount <- if (length(all_nonzero) > 0) {
  max(quantile(all_nonzero, 0.05), 1e-6)
} else {
  1e-6
}
mecp2_log2fc <- log2((mecp2_mut_mean + pseudocount) / (mecp2_ctrl_mean + pseudocount))
cat(sprintf("  Pseudocount: %.2e\n\n", pseudocount))

cat("--- Intra-TAD CV (10kb bins) ---\n")
cat("  Ctrl:\n")
mecp2_ctrl_cv <- compute_domain_cv(params$ctrl_bw, domain_gr)
cat(sprintf("    median CV=%.4f, non-NA=%d/%d\n",
            median(mecp2_ctrl_cv, na.rm = TRUE),
            sum(!is.na(mecp2_ctrl_cv)), length(mecp2_ctrl_cv)))
cat("  Mut:\n")
mecp2_mut_cv <- compute_domain_cv(params$mut_bw, domain_gr)
cat(sprintf("    median CV=%.4f, non-NA=%d/%d\n\n",
            median(mecp2_mut_cv, na.rm = TRUE),
            sum(!is.na(mecp2_mut_cv)), length(mecp2_mut_cv)))

cat("--- Within-domain variance (10kb bins) ---\n")
cat("  Ctrl:\n")
mecp2_ctrl_within_var <- compute_within_domain_var(params$ctrl_bw, domain_gr)
cat(sprintf("    median within-var=%.2e\n", median(mecp2_ctrl_within_var, na.rm = TRUE)))
cat("  Mut:\n")
mecp2_mut_within_var <- compute_within_domain_var(params$mut_bw, domain_gr)
cat(sprintf("    median within-var=%.2e\n\n", median(mecp2_mut_within_var, na.rm = TRUE)))

cat("--- Boundary insulation (50kb flanks) ---\n")
cat("  Ctrl:\n")
mecp2_ctrl_ins <- compute_boundary_insulation(params$ctrl_bw, domain_gr, flank_kb = 50)
cat("  Mut:\n")
mecp2_mut_ins <- compute_boundary_insulation(params$mut_bw, domain_gr, flank_kb = 50)

# =============================================================================
# COMPUTE MeCP2 PEAK OVERLAP COUNTS
# =============================================================================

cat("\n--- Peak overlap counts ---\n")

peaks_up <- read.table(params$peaks_up, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE)
peaks_up_gr <- GRanges(seqnames = peaks_up$V1,
                       ranges = IRanges(start = peaks_up$V2, end = peaks_up$V3))
cat(sprintf("  MeCP2-Up peaks loaded: %d\n", length(peaks_up_gr)))

peaks_down <- read.table(params$peaks_down, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
peaks_down_gr <- GRanges(seqnames = peaks_down$V1,
                         ranges = IRanges(start = peaks_down$V2, end = peaks_down$V3))
cat(sprintf("  MeCP2-Down peaks loaded: %d\n", length(peaks_down_gr)))

mecp2_up_count   <- countOverlaps(domain_gr, peaks_up_gr)
mecp2_down_count <- countOverlaps(domain_gr, peaks_down_gr)
mecp2_peak_count <- mecp2_up_count + mecp2_down_count

cat(sprintf("  TADs with any MeCP2 peak: %d / %d (%.1f%%)\n",
            sum(mecp2_peak_count > 0), length(mecp2_peak_count),
            100 * sum(mecp2_peak_count > 0) / length(mecp2_peak_count)))
cat(sprintf("  Total Up overlaps: %d, Down overlaps: %d\n\n",
            sum(mecp2_up_count), sum(mecp2_down_count)))

# =============================================================================
# ASSEMBLE AND WRITE OUTPUT
# =============================================================================

tad$mecp2_ctrl_mean      <- round(mecp2_ctrl_mean, 8)
tad$mecp2_mut_mean       <- round(mecp2_mut_mean, 8)
tad$mecp2_log2fc         <- round(mecp2_log2fc, 6)
tad$mecp2_pseudocount    <- pseudocount
tad$mecp2_ctrl_cv        <- round(mecp2_ctrl_cv, 6)
tad$mecp2_mut_cv         <- round(mecp2_mut_cv, 6)
tad$mecp2_ctrl_within_var <- mecp2_ctrl_within_var
tad$mecp2_mut_within_var  <- mecp2_mut_within_var
tad$mecp2_ctrl_insulation <- round(mecp2_ctrl_ins$mean, 8)
tad$mecp2_mut_insulation  <- round(mecp2_mut_ins$mean, 8)
tad$mecp2_peak_count     <- mecp2_peak_count
tad$mecp2_up_count       <- mecp2_up_count
tad$mecp2_down_count     <- mecp2_down_count

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Writing output to: %s\n", params$output))
write.table(tad, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d domains x %d columns written\n\n", nrow(tad), ncol(tad)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("================================================================================\n")
cat("MeCP2 AUGMENTATION SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Domains:           %d\n", nrow(tad)))
cat(sprintf("  New columns:       13 (10 signal + 3 peak counts)\n"))
cat(sprintf("  Total columns:     %d\n", ncol(tad)))
cat(sprintf("  Ctrl mean median:  %.6f\n", median(mecp2_ctrl_mean)))
cat(sprintf("  Mut mean median:   %.6f\n", median(mecp2_mut_mean)))
cat(sprintf("  Ctrl CV median:    %.4f\n", median(mecp2_ctrl_cv, na.rm = TRUE)))
cat(sprintf("  Mut CV median:     %.4f\n", median(mecp2_mut_cv, na.rm = TRUE)))

valid_ctrl <- !is.na(mecp2_ctrl_mean) & !is.na(mecp2_ctrl_within_var) &
              mecp2_ctrl_within_var > 0
if (sum(valid_ctrl) > 10) {
  vr <- var(mecp2_ctrl_mean[valid_ctrl]) / mean(mecp2_ctrl_within_var[valid_ctrl])
  cat(sprintf("  Ctrl variance ratio: %.4f (n=%d)\n", vr, sum(valid_ctrl)))
}

cat(sprintf("  Peak overlaps:     %d Up, %d Down across %d TADs\n",
            sum(mecp2_up_count), sum(mecp2_down_count),
            sum(mecp2_peak_count > 0)))
cat("\nDone.\n")
