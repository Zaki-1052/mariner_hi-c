# biomodal/downstream/scripts/tad_methylation/compute_tad_methylation.R
# HPC preprocessing: Extract methylation signal over TAD domains and boundaries
#
# Reconstructs TAD domain intervals from TADCompare boundary coordinates,
# then computes per-sample mean methylation, intra-TAD CV, and boundary
# insulation scores for CG mC, CG hmC, CHG mC, and CHH mC.
#
# Usage (on Expanse HPC):
#   Rscript scripts/tad_methylation/compute_tad_methylation.R \
#     --tad-file ../../tads/results/late/final/tadcompare_final_filtered.tsv \
#     --timepoint late \
#     --bw-dir /expanse/lustre/projects/csd940/zalibhai/biomodal/modality/exports/bigwig \
#     --output data/tad_methylation_signal.tsv \
#     --min-domain-kb 100

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    tad_file = "../../tads/results/late/final/tadcompare_final_filtered.tsv",
    timepoint = "late",
    bw_dir = "/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/exports/bigwig",
    output = "data/tad_methylation_signal.tsv",
    min_domain_kb = 100
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--tad-file" && i < length(args)) {
      defaults$tad_file <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--timepoint" && i < length(args)) {
      defaults$timepoint <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--bw-dir" && i < length(args)) {
      defaults$bw_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      defaults$output <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--min-domain-kb" && i < length(args)) {
      defaults$min_domain_kb <- as.numeric(args[i + 1]); i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\nUsage: Rscript compute_tad_methylation.R --tad-file <tsv> --bw-dir <dir> --output <tsv>", args[i]))
    }
  }
  defaults
}

params <- parse_args()

cat("================================================================================\n")
cat("TAD METHYLATION SIGNAL PREPROCESSING\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("TAD file:  ", params$tad_file, "\n")
cat("Timepoint: ", params$timepoint, "\n")
cat("BigWig dir:", params$bw_dir, "\n")
cat("Output:    ", params$output, "\n")
cat("Min domain:", params$min_domain_kb, "kb\n\n")

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

if (!file.exists(params$tad_file)) stop("TAD file not found: ", params$tad_file)
if (!dir.exists(params$bw_dir)) stop("BigWig directory not found: ", params$bw_dir)

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# SAMPLE AND BIGWIG CONFIGURATION
# =============================================================================

CTRL_SAMPLES <- c("evoC-Bap1-ctrl-F", "evoC-Bap1-ctrl-M",
                   "evoC-Bap1-ctrl-F-B2", "evoC-Bap1-ctrl-M-B2")
MUT_SAMPLES  <- c("evoC-Bap1-mut-F", "evoC-Bap1-mut-M",
                   "evoC-Bap1-mut-F-B2", "evoC-Bap1-mut-M-B2")
ALL_SAMPLES  <- c(CTRL_SAMPLES, MUT_SAMPLES)

CONTEXTS <- list(
  cg_mc  = list(subdir = "mc",     suffix = ".genome.mc_bedgraph.bw"),
  cg_hmc = list(subdir = "hmc",    suffix = ".genome.hmc_bedgraph.bw"),
  chg_mc = list(subdir = "CHG_mc", suffix = ".genome.mc_bedgraph.bw"),
  chh_mc = list(subdir = "CHH_mc", suffix = ".genome.mc_bedgraph.bw")
)

build_bw_path <- function(bw_dir, context_info, sample_id) {
  file.path(bw_dir, context_info$subdir,
            paste0(sample_id, context_info$suffix))
}

cat("Checking BigWig availability:\n")
for (ctx_name in names(CONTEXTS)) {
  ctx <- CONTEXTS[[ctx_name]]
  found <- sum(sapply(ALL_SAMPLES, function(s)
    file.exists(build_bw_path(params$bw_dir, ctx, s))))
  cat(sprintf("  %s: %d / %d samples found\n", ctx_name, found, length(ALL_SAMPLES)))
}
cat("\n")

# =============================================================================
# RECONSTRUCT TAD DOMAINS FROM BOUNDARY COORDINATES
# =============================================================================

cat("Reconstructing TAD domains from boundaries...\n")

AUTOSOMES <- paste0("chr", 1:19)

reconstruct_tad_domains <- function(tad_file, min_kb, chromosomes = AUTOSOMES) {
  tad_raw <- read.table(tad_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tad_raw <- tad_raw[tad_raw$chr %in% chromosomes, ]
  cat(sprintf("  Raw boundaries: %d (autosomes)\n", nrow(tad_raw)))

  domain_list <- list()
  for (chr in chromosomes) {
    chr_bounds <- tad_raw[tad_raw$chr == chr, ]
    chr_bounds <- chr_bounds[order(chr_bounds$Boundary), ]
    if (nrow(chr_bounds) < 2) next

    for (i in seq_len(nrow(chr_bounds) - 1)) {
      b1 <- chr_bounds[i, ]
      b2 <- chr_bounds[i + 1, ]
      width_bp <- b2$Boundary - b1$Boundary
      if (width_bp < min_kb * 1000) next

      domain_list[[length(domain_list) + 1]] <- data.frame(
        chr = chr,
        start = b1$Boundary,
        end = b2$Boundary,
        width_kb = width_bp / 1000,
        boundary_left_pos = b1$Boundary,
        boundary_right_pos = b2$Boundary,
        boundary_left_type = b1$Type,
        boundary_right_type = b2$Type,
        boundary_left_differential = b1$Differential,
        boundary_right_differential = b2$Differential,
        boundary_left_gap_score = b1$Gap_Score,
        boundary_right_gap_score = b2$Gap_Score,
        boundary_left_tad_score1 = b1$TAD_Score1,
        boundary_left_tad_score2 = b1$TAD_Score2,
        boundary_right_tad_score1 = b2$TAD_Score1,
        boundary_right_tad_score2 = b2$TAD_Score2,
        boundary_left_enriched = b1$Enriched_In,
        boundary_right_enriched = b2$Enriched_In,
        stringsAsFactors = FALSE
      )
    }
  }

  domains <- do.call(rbind, domain_list)
  domains$is_differential <- (domains$boundary_left_differential == "Differential" |
                               domains$boundary_right_differential == "Differential")
  domains$mean_tad_score <- (domains$boundary_left_tad_score1 +
                              domains$boundary_right_tad_score1) / 2
  domains$mean_gap_score <- (domains$boundary_left_gap_score +
                              domains$boundary_right_gap_score) / 2
  domains
}

domains <- reconstruct_tad_domains(params$tad_file, params$min_domain_kb)
cat(sprintf("  Reconstructed domains (>= %d kb): %d\n", params$min_domain_kb, nrow(domains)))
cat(sprintf("  Differential: %d (%.1f%%)\n",
            sum(domains$is_differential),
            100 * sum(domains$is_differential) / nrow(domains)))
cat(sprintf("  Width range: %.0f - %.0f kb (median: %.0f kb)\n",
            min(domains$width_kb), max(domains$width_kb), median(domains$width_kb)))

domain_gr <- GRanges(
  seqnames = domains$chr,
  ranges = IRanges(start = domains$start, end = domains$end)
)

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
# EXTRACT SIGNAL FOR ALL CONTEXTS AND SAMPLES
# =============================================================================

cat("\n================================================================================\n")
cat("SIGNAL EXTRACTION\n")
cat("================================================================================\n\n")

result_df <- domains

for (ctx_name in names(CONTEXTS)) {
  ctx <- CONTEXTS[[ctx_name]]
  cat(sprintf("\n--- Context: %s ---\n", ctx_name))

  ctrl_means_mat <- matrix(NA, nrow = nrow(domains), ncol = length(CTRL_SAMPLES))
  mut_means_mat  <- matrix(NA, nrow = nrow(domains), ncol = length(MUT_SAMPLES))

  # Per-sample mean signal
  cat("  Per-sample domain means:\n")
  for (j in seq_along(CTRL_SAMPLES)) {
    bw_path <- build_bw_path(params$bw_dir, ctx, CTRL_SAMPLES[j])
    if (!file.exists(bw_path)) {
      cat(sprintf("    SKIP (not found): %s\n", CTRL_SAMPLES[j]))
      next
    }
    ctrl_means_mat[, j] <- compute_mean_signal(bw_path, domain_gr,
                                                paste0("ctrl:", CTRL_SAMPLES[j]))
  }
  for (j in seq_along(MUT_SAMPLES)) {
    bw_path <- build_bw_path(params$bw_dir, ctx, MUT_SAMPLES[j])
    if (!file.exists(bw_path)) {
      cat(sprintf("    SKIP (not found): %s\n", MUT_SAMPLES[j]))
      next
    }
    mut_means_mat[, j] <- compute_mean_signal(bw_path, domain_gr,
                                               paste0("mut:", MUT_SAMPLES[j]))
  }

  ctrl_mean <- rowMeans(ctrl_means_mat, na.rm = TRUE)
  mut_mean  <- rowMeans(mut_means_mat, na.rm = TRUE)

  # log2FC with pseudocount
  all_nonzero <- c(ctrl_mean[ctrl_mean > 0], mut_mean[mut_mean > 0])
  if (length(all_nonzero) > 0) {
    pseudocount <- max(quantile(all_nonzero, 0.05), 1e-6)
  } else {
    pseudocount <- 1e-6
  }
  log2fc <- log2((mut_mean + pseudocount) / (ctrl_mean + pseudocount))

  result_df[[paste0(ctx_name, "_ctrl_mean")]] <- round(ctrl_mean, 8)
  result_df[[paste0(ctx_name, "_mut_mean")]]  <- round(mut_mean, 8)
  result_df[[paste0(ctx_name, "_log2fc")]]     <- round(log2fc, 6)
  result_df[[paste0(ctx_name, "_pseudocount")]] <- pseudocount

  cat(sprintf("  %s: ctrl_median=%.6f, mut_median=%.6f, pseudocount=%.2e\n",
              ctx_name, median(ctrl_mean), median(mut_mean), pseudocount))

  # Intra-TAD CV (use first ctrl sample as representative — CV is condition-independent)
  cat("  Intra-TAD CV:\n")
  ctrl_cv_bw <- build_bw_path(params$bw_dir, ctx, CTRL_SAMPLES[1])
  mut_cv_bw  <- build_bw_path(params$bw_dir, ctx, MUT_SAMPLES[1])

  if (file.exists(ctrl_cv_bw)) {
    cat(sprintf("    Computing ctrl CV from %s\n", CTRL_SAMPLES[1]))
    ctrl_cv <- compute_domain_cv(ctrl_cv_bw, domain_gr)
    result_df[[paste0(ctx_name, "_ctrl_cv")]] <- round(ctrl_cv, 6)
    cat(sprintf("    ctrl CV: median=%.4f, non-NA=%d/%d\n",
                median(ctrl_cv, na.rm = TRUE),
                sum(!is.na(ctrl_cv)), length(ctrl_cv)))
  }

  if (file.exists(mut_cv_bw)) {
    cat(sprintf("    Computing mut CV from %s\n", MUT_SAMPLES[1]))
    mut_cv <- compute_domain_cv(mut_cv_bw, domain_gr)
    result_df[[paste0(ctx_name, "_mut_cv")]] <- round(mut_cv, 6)
    cat(sprintf("    mut CV: median=%.4f, non-NA=%d/%d\n",
                median(mut_cv, na.rm = TRUE),
                sum(!is.na(mut_cv)), length(mut_cv)))
  }

  # Within-domain variance (for inter/intra ratio computation downstream)
  cat("  Within-domain variance:\n")
  if (file.exists(ctrl_cv_bw)) {
    ctrl_within_var <- compute_within_domain_var(ctrl_cv_bw, domain_gr)
    result_df[[paste0(ctx_name, "_ctrl_within_var")]] <- ctrl_within_var
    cat(sprintf("    ctrl within-var: median=%.2e\n", median(ctrl_within_var, na.rm = TRUE)))
  }
  if (file.exists(mut_cv_bw)) {
    mut_within_var <- compute_within_domain_var(mut_cv_bw, domain_gr)
    result_df[[paste0(ctx_name, "_mut_within_var")]] <- mut_within_var
    cat(sprintf("    mut within-var: median=%.2e\n", median(mut_within_var, na.rm = TRUE)))
  }

  # Boundary insulation
  cat("  Boundary insulation (50kb flanks):\n")
  if (file.exists(ctrl_cv_bw)) {
    ctrl_ins <- compute_boundary_insulation(ctrl_cv_bw, domain_gr, flank_kb = 50)
    result_df[[paste0(ctx_name, "_ctrl_insulation")]] <- round(ctrl_ins$mean, 8)
  }
  if (file.exists(mut_cv_bw)) {
    mut_ins <- compute_boundary_insulation(mut_cv_bw, domain_gr, flank_kb = 50)
    result_df[[paste0(ctx_name, "_mut_insulation")]] <- round(mut_ins$mean, 8)
  }
}

# =============================================================================
# WRITE OUTPUT
# =============================================================================

cat(sprintf("\n\nWriting output to: %s\n", params$output))
write.table(result_df, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d domains x %d columns written\n", nrow(result_df), ncol(result_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Timepoint:      %s\n", params$timepoint))
cat(sprintf("  Total domains:  %d (>= %d kb)\n", nrow(result_df), params$min_domain_kb))
cat(sprintf("  Differential:   %d (%.1f%%)\n",
            sum(result_df$is_differential),
            100 * sum(result_df$is_differential) / nrow(result_df)))

cat("\n  Per-context domain mean summaries (median across domains):\n")
for (ctx_name in names(CONTEXTS)) {
  ctrl_col <- paste0(ctx_name, "_ctrl_mean")
  mut_col  <- paste0(ctx_name, "_mut_mean")
  cv_col   <- paste0(ctx_name, "_ctrl_cv")
  ins_col  <- paste0(ctx_name, "_ctrl_insulation")

  if (ctrl_col %in% names(result_df)) {
    cat(sprintf("    %s: ctrl=%.6f, mut=%.6f",
                ctx_name,
                median(result_df[[ctrl_col]], na.rm = TRUE),
                median(result_df[[mut_col]], na.rm = TRUE)))
    if (cv_col %in% names(result_df)) {
      cat(sprintf(", CV=%.4f", median(result_df[[cv_col]], na.rm = TRUE)))
    }
    if (ins_col %in% names(result_df)) {
      cat(sprintf(", insulation=%.6f", median(result_df[[ins_col]], na.rm = TRUE)))
    }
    cat("\n")
  }
}

# Inter/intra variance ratio preview
cat("\n  Inter/Intra TAD variance ratio (higher = more TAD-organized):\n")
for (ctx_name in names(CONTEXTS)) {
  mean_col <- paste0(ctx_name, "_ctrl_mean")
  var_col  <- paste0(ctx_name, "_ctrl_within_var")
  if (mean_col %in% names(result_df) && var_col %in% names(result_df)) {
    domain_means <- result_df[[mean_col]]
    within_vars  <- result_df[[var_col]]
    valid <- !is.na(domain_means) & !is.na(within_vars) & within_vars > 0
    if (sum(valid) > 10) {
      inter_var <- var(domain_means[valid])
      intra_var <- mean(within_vars[valid])
      ratio <- inter_var / intra_var
      cat(sprintf("    %s: %.4f (inter=%.2e, intra=%.2e, n=%d)\n",
                  ctx_name, ratio, inter_var, intra_var, sum(valid)))
    }
  }
}

cat("\nPreprocessing complete.\n")
