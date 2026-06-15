# biomodal/downstream/scripts/tad_methylation/add_mecp2_tad_signal.R
# Augment the existing TAD methylation signal TSV with MeCP2 ChIP-seq metrics
# and compute bin-level residual variances for the CG-corrected proxy analysis.
#
# Designed for HPC (Expanse). Adds 13 MeCP2 signal columns to the main TSV
# and produces a separate bin-level variance TSV for section 55.
#
# Usage (from downstream/ directory on Expanse):
#   sbatch scripts/tad_methylation/add_mecp2_tad_signal.sb
#
# Or manually:
#   Rscript scripts/tad_methylation/add_mecp2_tad_signal.R \
#     --input  data/tad_methylation_signal_late.tsv \
#     --output data/tad_methylation_signal_late.tsv \
#     --ctrl-bw /expanse/.../bigwigs/MeCP2Ctrl.bw \
#     --mut-bw  /expanse/.../bigwigs/MeCP2Mut.bw \
#     --peaks-up  ../../peaks/mecp2/MeCP2_up.bed \
#     --peaks-down ../../peaks/mecp2/MeCP2_down.bed \
#     --cg-mc-ctrl /expanse/.../merged/CG_mc_ctrl.bw \
#     --cg-mc-mut  /expanse/.../merged/CG_mc_mut.bw \
#     --cg-hmc-ctrl /expanse/.../merged/CG_hmc_ctrl.bw \
#     --cg-hmc-mut  /expanse/.../merged/CG_hmc_mut.bw \
#     --binvar-output data/tad_mecp2_binlevel_variance.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

BW_BASE <- "/expanse/lustre/projects/csd940/zalibhai/bigwigs"
METH_BASE <- "/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/exports/bigwig/merged"

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input        = "data/tad_methylation_signal_late.tsv",
    output       = "data/tad_methylation_signal_late.tsv",
    ctrl_bw      = file.path(BW_BASE, "MeCP2Ctrl.bw"),
    mut_bw       = file.path(BW_BASE, "MeCP2Mut.bw"),
    peaks_up     = "../../peaks/mecp2/MeCP2_up.bed",
    peaks_down   = "../../peaks/mecp2/MeCP2_down.bed",
    cg_mc_ctrl   = file.path(METH_BASE, "CG_mc_ctrl.bw"),
    cg_mc_mut    = file.path(METH_BASE, "CG_mc_mut.bw"),
    cg_hmc_ctrl  = file.path(METH_BASE, "CG_hmc_ctrl.bw"),
    cg_hmc_mut   = file.path(METH_BASE, "CG_hmc_mut.bw"),
    binvar_output = "data/tad_mecp2_binlevel_variance.tsv"
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
    } else if (args[i] == "--cg-mc-ctrl" && i < length(args)) {
      defaults$cg_mc_ctrl <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--cg-mc-mut" && i < length(args)) {
      defaults$cg_mc_mut <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--cg-hmc-ctrl" && i < length(args)) {
      defaults$cg_hmc_ctrl <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--cg-hmc-mut" && i < length(args)) {
      defaults$cg_hmc_mut <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--binvar-output" && i < length(args)) {
      defaults$binvar_output <- args[i + 1]; i <- i + 2
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
cat("Input:         ", params$input, "\n")
cat("Output:        ", params$output, "\n")
cat("BinVar output: ", params$binvar_output, "\n")
cat("MeCP2 Ctrl BW: ", params$ctrl_bw, "\n")
cat("MeCP2 Mut BW:  ", params$mut_bw, "\n")
cat("CG mC Ctrl BW: ", params$cg_mc_ctrl, "\n")
cat("CG mC Mut BW:  ", params$cg_mc_mut, "\n")
cat("CG hmC Ctrl BW:", params$cg_hmc_ctrl, "\n")
cat("CG hmC Mut BW: ", params$cg_hmc_mut, "\n")
cat("Peaks Up:      ", params$peaks_up, "\n")
cat("Peaks Down:    ", params$peaks_down, "\n\n")

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

if (!file.exists(params$input))      stop("Input TSV not found: ", params$input)
if (!file.exists(params$ctrl_bw))    stop("MeCP2 Ctrl BigWig not found: ", params$ctrl_bw)
if (!file.exists(params$mut_bw))     stop("MeCP2 Mut BigWig not found: ", params$mut_bw)
if (!file.exists(params$peaks_up))   stop("Peaks-up BED not found: ", params$peaks_up)
if (!file.exists(params$peaks_down)) stop("Peaks-down BED not found: ", params$peaks_down)
if (!file.exists(params$cg_mc_ctrl))  stop("CG mC Ctrl BigWig not found: ", params$cg_mc_ctrl)
if (!file.exists(params$cg_mc_mut))   stop("CG mC Mut BigWig not found: ", params$cg_mc_mut)
if (!file.exists(params$cg_hmc_ctrl)) stop("CG hmC Ctrl BigWig not found: ", params$cg_hmc_ctrl)
if (!file.exists(params$cg_hmc_mut))  stop("CG hmC Mut BigWig not found: ", params$cg_hmc_mut)

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

cat(sprintf("Writing augmented TSV to: %s\n", params$output))
write.table(tad, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d domains x %d columns written\n\n", nrow(tad), ncol(tad)))

# =============================================================================
# BIN-LEVEL RESIDUAL VARIANCE FOR SECTION 55
# =============================================================================

cat("================================================================================\n")
cat("BIN-LEVEL RESIDUAL VARIANCE COMPUTATION\n")
cat("================================================================================\n\n")

cat("Fitting regression models on per-TAD means...\n")
valid_rows <- complete.cases(tad[, c("mecp2_ctrl_mean", "cg_mc_ctrl_mean",
                                      "cg_hmc_ctrl_mean")])
cat(sprintf("  Complete cases: %d / %d\n", sum(valid_rows), nrow(tad)))

model_a <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean, data = tad, subset = valid_rows)
model_b <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean + cg_hmc_ctrl_mean,
              data = tad, subset = valid_rows)

coefs_a <- coef(model_a)
coefs_b <- coef(model_b)
cat(sprintf("  Model A (MeCP2 ~ CG_mC): R²=%.4f, slope=%.4f\n",
            summary(model_a)$r.squared, coefs_a[2]))
cat(sprintf("  Model B (MeCP2 ~ CG_mC + CG_hmC): R²=%.4f\n\n",
            summary(model_b)$r.squared))

compute_bin_means <- function(bw_path, domain_gr, bin_size = 10000) {
  bw <- import.bw(bw_path, which = domain_gr)
  cov <- coverage(bw, weight = "score")
  chr_names <- as.character(seqnames(domain_gr))

  result <- vector("list", length(domain_gr))
  for (i in seq_along(domain_gr)) {
    dom <- domain_gr[i]
    chr <- chr_names[i]
    dom_start <- start(dom)
    dom_end <- end(dom)

    tile_starts <- seq(dom_start, dom_end - 1, by = bin_size)
    tile_ends <- pmin(tile_starts + bin_size - 1, dom_end)
    tiles_ir <- IRanges(start = tile_starts, end = tile_ends)

    if (!(chr %in% names(cov))) {
      result[[i]] <- rep(NA_real_, length(tiles_ir))
      next
    }

    cov_chr <- cov[[chr]]
    max_pos <- length(cov_chr)
    valid_mask <- end(tiles_ir) <= max_pos
    bin_means <- rep(NA_real_, length(tiles_ir))
    if (any(valid_mask)) {
      v <- Views(cov_chr, tiles_ir[valid_mask])
      bin_means[valid_mask] <- as.numeric(viewMeans(v))
    }
    result[[i]] <- bin_means
  }
  result
}

cat("Extracting bin-level signal from 6 BigWigs (10kb bins)...\n")
cat("  MeCP2 ctrl...\n")
mecp2_ctrl_bins <- compute_bin_means(params$ctrl_bw, domain_gr)
cat("  MeCP2 mut...\n")
mecp2_mut_bins <- compute_bin_means(params$mut_bw, domain_gr)
cat("  CG mC ctrl...\n")
cgmc_ctrl_bins <- compute_bin_means(params$cg_mc_ctrl, domain_gr)
cat("  CG mC mut...\n")
cgmc_mut_bins <- compute_bin_means(params$cg_mc_mut, domain_gr)
cat("  CG hmC ctrl...\n")
cghmc_ctrl_bins <- compute_bin_means(params$cg_hmc_ctrl, domain_gr)
cat("  CG hmC mut...\n")
cghmc_mut_bins <- compute_bin_means(params$cg_hmc_mut, domain_gr)

cat("\nComputing per-domain residual variances...\n")

within_var_raw_ctrl     <- numeric(nrow(tad))
within_var_raw_mut      <- numeric(nrow(tad))
within_var_resid_a_ctrl <- numeric(nrow(tad))
within_var_resid_a_mut  <- numeric(nrow(tad))
within_var_resid_b_ctrl <- numeric(nrow(tad))
within_var_resid_b_mut  <- numeric(nrow(tad))

for (i in seq_len(nrow(tad))) {
  m_ctrl <- mecp2_ctrl_bins[[i]]
  c_ctrl <- cgmc_ctrl_bins[[i]]
  h_ctrl <- cghmc_ctrl_bins[[i]]
  m_mut  <- mecp2_mut_bins[[i]]
  c_mut  <- cgmc_mut_bins[[i]]
  h_mut  <- cghmc_mut_bins[[i]]

  valid_ctrl <- !is.na(m_ctrl)
  valid_mut  <- !is.na(m_mut)

  within_var_raw_ctrl[i] <- if (sum(valid_ctrl) >= 2) {
    var(m_ctrl[valid_ctrl], na.rm = TRUE)
  } else NA_real_

  within_var_raw_mut[i] <- if (sum(valid_mut) >= 2) {
    var(m_mut[valid_mut], na.rm = TRUE)
  } else NA_real_

  valid_a_ctrl <- !is.na(m_ctrl) & !is.na(c_ctrl)
  if (sum(valid_a_ctrl) >= 2) {
    pred <- coefs_a[1] + coefs_a[2] * c_ctrl[valid_a_ctrl]
    within_var_resid_a_ctrl[i] <- var(m_ctrl[valid_a_ctrl] - pred)
  } else within_var_resid_a_ctrl[i] <- NA_real_

  valid_a_mut <- !is.na(m_mut) & !is.na(c_mut)
  if (sum(valid_a_mut) >= 2) {
    pred <- coefs_a[1] + coefs_a[2] * c_mut[valid_a_mut]
    within_var_resid_a_mut[i] <- var(m_mut[valid_a_mut] - pred)
  } else within_var_resid_a_mut[i] <- NA_real_

  valid_b_ctrl <- !is.na(m_ctrl) & !is.na(c_ctrl) & !is.na(h_ctrl)
  if (sum(valid_b_ctrl) >= 2) {
    pred <- coefs_b[1] + coefs_b[2] * c_ctrl[valid_b_ctrl] +
            coefs_b[3] * h_ctrl[valid_b_ctrl]
    within_var_resid_b_ctrl[i] <- var(m_ctrl[valid_b_ctrl] - pred)
  } else within_var_resid_b_ctrl[i] <- NA_real_

  valid_b_mut <- !is.na(m_mut) & !is.na(c_mut) & !is.na(h_mut)
  if (sum(valid_b_mut) >= 2) {
    pred <- coefs_b[1] + coefs_b[2] * c_mut[valid_b_mut] +
            coefs_b[3] * h_mut[valid_b_mut]
    within_var_resid_b_mut[i] <- var(m_mut[valid_b_mut] - pred)
  } else within_var_resid_b_mut[i] <- NA_real_
}

binvar_df <- data.frame(
  chr = tad$chr, start = tad$start, end = tad$end,
  within_var_raw_ctrl     = within_var_raw_ctrl,
  within_var_raw_mut      = within_var_raw_mut,
  within_var_resid_a_ctrl = within_var_resid_a_ctrl,
  within_var_resid_a_mut  = within_var_resid_a_mut,
  within_var_resid_b_ctrl = within_var_resid_b_ctrl,
  within_var_resid_b_mut  = within_var_resid_b_mut,
  model_a_intercept = coefs_a[1],
  model_a_slope_mc  = coefs_a[2],
  model_a_r_squared = summary(model_a)$r.squared,
  model_b_intercept  = coefs_b[1],
  model_b_slope_mc   = coefs_b[2],
  model_b_slope_hmc  = coefs_b[3],
  model_b_r_squared  = summary(model_b)$r.squared,
  stringsAsFactors = FALSE
)

binvar_dir <- dirname(params$binvar_output)
if (binvar_dir != "" && binvar_dir != ".") {
  dir.create(binvar_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("\nWriting bin-level variance TSV to: %s\n", params$binvar_output))
write.table(binvar_df, params$binvar_output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d domains x %d columns\n", nrow(binvar_df), ncol(binvar_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("MeCP2 AUGMENTATION SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Domains:           %d\n", nrow(tad)))
cat(sprintf("  New columns:       13 (10 signal + 3 peak counts)\n"))
cat(sprintf("  Total columns:     %d\n", ncol(tad)))
cat(sprintf("  Ctrl mean median:  %.6f\n", median(mecp2_ctrl_mean)))
cat(sprintf("  Mut mean median:   %.6f\n", median(mecp2_mut_mean)))
cat(sprintf("  Ctrl CV median:    %.4f\n", median(mecp2_ctrl_cv, na.rm = TRUE)))
cat(sprintf("  Mut CV median:     %.4f\n", median(mecp2_mut_cv, na.rm = TRUE)))

valid_ctrl <- !is.na(mecp2_ctrl_mean) & !is.na(within_var_raw_ctrl) &
              within_var_raw_ctrl > 0
if (sum(valid_ctrl) > 10) {
  vr <- var(mecp2_ctrl_mean[valid_ctrl]) / mean(within_var_raw_ctrl[valid_ctrl])
  cat(sprintf("  Ctrl variance ratio (raw): %.4f (n=%d)\n", vr, sum(valid_ctrl)))
}

valid_resid <- !is.na(tad$mecp2_ctrl_mean) & !is.na(within_var_resid_a_ctrl) &
               within_var_resid_a_ctrl > 0
if (sum(valid_resid) > 10) {
  resid_ctrl <- residuals(model_a)
  resid_means <- numeric(nrow(tad))
  resid_means[valid_rows] <- resid_ctrl
  vr_r <- var(resid_means[valid_resid]) / mean(within_var_resid_a_ctrl[valid_resid])
  cat(sprintf("  Ctrl variance ratio (resid A): %.4f (n=%d)\n", vr_r, sum(valid_resid)))
}

cat(sprintf("  Peak overlaps:     %d Up, %d Down across %d TADs\n",
            sum(mecp2_up_count), sum(mecp2_down_count),
            sum(mecp2_peak_count > 0)))

cat(sprintf("\nOutputs:\n  %s\n  %s\n", params$output, params$binvar_output))
cat("\nDone.\n")
