# loops/scripts/ctcf_motif_analysis.R
#
# 2. Loop Anchor CTCF Motif Analysis
#
# Refines loop anchor annotation by centering on CTCF motifs within anchors.
# Current 10 kb anchors are coarse; this script examines motif-level signal.
# Sub-analyses:
#   2a. Re-center chromatin state annotation on individual CTCF motifs
#   2b. CTCF motif filter for loop type classification
#   2c. "Other-end" anchor type for CTCF-anchored loops
#
# Usage:
#   Rscript loops/scripts/ctcf_motif_analysis.R --timepoint late
#   Rscript loops/scripts/ctcf_motif_analysis.R --timepoint early
#   Rscript loops/scripts/ctcf_motif_analysis.R --timepoint both
#
# Output:
#   loops/output/ctcf_motif_analysis/{timepoint}/

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(RColorBrewer)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))

MOTIF_WINDOW_BP <- 500L
TSS_THRESHOLD <- 2000L
VALID_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

PEAK_FILES <- list(
  early = list(
    h3k27ac  = file.path(BASE_DIR, "peaks/beds/H3K27acCerebellumEarly2.bed"),
    h3k27me3 = file.path(BASE_DIR, "peaks/beds/H3K27me3CerebellumEarly1.bed"),
    h3k4me1  = file.path(BASE_DIR, "peaks/beds/H3K4me1CerebellumEarly1.bed"),
    h3k4me3  = file.path(BASE_DIR, "peaks/beds/H3K4me3CerebellumEarly2.bed"),
    bivalent = file.path(BASE_DIR, "peaks/beds/Bivalent_Cerebellum_Early.bed")
  ),
  late = list(
    h3k27ac  = file.path(BASE_DIR, "peaks/beds/H3K27acCerebellumLate2.bed"),
    h3k27me3 = file.path(BASE_DIR, "peaks/beds/H3K27me3CerebellumLate1.bed"),
    h3k4me1  = file.path(BASE_DIR, "peaks/beds/H3K4me1CerebellumLate1.bed"),
    h3k4me3  = file.path(BASE_DIR, "peaks/beds/H3K4me3CerebellumLate2.bed"),
    bivalent = file.path(BASE_DIR, "peaks/beds/Bivalent_Cerebellum_Late.bed")
  )
)

INPUT_FILES <- list(
  early = file.path(BASE_DIR, "peaks/loop_annotation_extended/early/extended_characterized_loops.tsv"),
  late  = file.path(BASE_DIR, "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv")
)

CTCF_MOTIF_FILE <- file.path(BASE_DIR, "peaks/ctcf_motifs_mm10.bed")
CTCF_CHIP_FILE  <- file.path(BASE_DIR, "peaks/CTCF.bed")

ANCHOR_TYPE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                       "Polycomb", "Active_Enhancer", "Poised_Enhancer", "CTCF_Site", "Other")

ANCHOR_COLORS <- c(
  "Active_Promoter"    = "#e41a1c",
  "Repressed_Promoter" = "#756bb1",
  "Bivalent_Promoter"  = "#984ea3",
  "Polycomb"           = "#4daf4a",
  "Active_Enhancer"    = "#377eb8",
  "Poised_Enhancer"    = "#ff7f00",
  "CTCF_Site"          = "#a65628",
  "Other"              = "#999999"
)

CTCF_COUNT_COLORS <- c("0" = "#bdbdbd", "1" = "#6baed6", "2" = "#08519c")

# =============================================================================
# 2. ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  timepoint <- "both"
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--help") {
      cat("Usage: Rscript loops/scripts/ctcf_motif_analysis.R [--timepoint early|late|both]\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }
  return(list(timepoint = timepoint))
}

parsed <- parse_args(args)

# =============================================================================
# 3. CORE FUNCTIONS
# =============================================================================

load_loops <- function(timepoint) {
  input_file <- INPUT_FILES[[timepoint]]
  if (!file.exists(input_file)) stop(sprintf("Input file not found: %s", input_file))

  loops_df <- readr::read_tsv(input_file, show_col_types = FALSE)

  required_cols <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
                     "direction", "anchor1_CTCF_motif_overlap", "anchor2_CTCF_motif_overlap",
                     "anchor1_type_extended", "anchor2_type_extended", "loop_type_extended")
  missing <- setdiff(required_cols, colnames(loops_df))
  if (length(missing) > 0) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))

  n_up <- sum(loops_df$direction == "up_in_mutant")
  n_down <- sum(loops_df$direction == "down_in_mutant")
  cat(sprintf("  Loaded %d loops (%d up, %d down) for %s\n", nrow(loops_df), n_up, n_down, timepoint))
  return(loops_df)
}

load_chip_peaks_local <- function(bed_path, mark_name) {
  if (!file.exists(bed_path)) stop(sprintf("%s BED not found: %s", mark_name, bed_path))

  bed <- read.table(bed_path, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
  bed <- bed[bed[, 1] %in% VALID_CHROMS, ]

  gr <- GRanges(seqnames = bed[, 1],
                ranges = IRanges(start = as.integer(bed[, 2]), end = as.integer(bed[, 3])))
  cat(sprintf("  Loaded %d %s peaks\n", length(gr), mark_name))
  return(gr)
}

load_ctcf_motifs <- function() {
  if (!file.exists(CTCF_MOTIF_FILE)) stop(sprintf("CTCF motif file not found: %s", CTCF_MOTIF_FILE))

  bed <- read.table(CTCF_MOTIF_FILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  bed <- bed[bed[, 1] %in% VALID_CHROMS, ]

  gr <- GRanges(seqnames = bed[, 1],
                ranges = IRanges(start = as.integer(bed[, 2]), end = as.integer(bed[, 3])),
                motif_name = bed[, 4],
                score = as.integer(bed[, 5]),
                strand = bed[, 6])
  cat(sprintf("  Loaded %d CTCF motifs\n", length(gr)))
  return(gr)
}

anchors_to_granges <- function(loops_df, anchor_num) {
  chr_col   <- paste0("chr", anchor_num)
  start_col <- paste0("start", anchor_num)
  end_col   <- paste0("end", anchor_num)

  GRanges(seqnames = loops_df[[chr_col]],
          ranges = IRanges(start = loops_df[[start_col]], end = loops_df[[end_col]]),
          loop_id    = loops_df$loop_id,
          direction  = loops_df$direction,
          anchor_num = anchor_num)
}

find_motifs_in_anchors <- function(anchor_gr, motif_gr) {
  hits <- findOverlaps(anchor_gr, motif_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(tibble(loop_id = character(), anchor_num = integer(), direction = character(),
                  anchor_chr = character(), anchor_start = integer(), anchor_end = integer(),
                  motif_chr = character(), motif_start = integer(), motif_end = integer(),
                  motif_score = integer(), motif_strand = character()))
  }

  a_idx <- queryHits(hits)
  m_idx <- subjectHits(hits)

  tibble(
    loop_id      = mcols(anchor_gr)$loop_id[a_idx],
    anchor_num   = mcols(anchor_gr)$anchor_num[a_idx],
    direction    = mcols(anchor_gr)$direction[a_idx],
    anchor_chr   = as.character(seqnames(anchor_gr))[a_idx],
    anchor_start = start(anchor_gr)[a_idx],
    anchor_end   = end(anchor_gr)[a_idx],
    motif_chr    = as.character(seqnames(motif_gr))[m_idx],
    motif_start  = start(motif_gr)[m_idx],
    motif_end    = end(motif_gr)[m_idx],
    motif_score  = mcols(motif_gr)$score[m_idx],
    motif_strand = as.character(strand(motif_gr))[m_idx]
  )
}

build_motif_windows <- function(motif_df, window_bp = MOTIF_WINDOW_BP) {
  motif_center <- as.integer((motif_df$motif_start + motif_df$motif_end) / 2)
  window_start <- pmax(1L, motif_center - window_bp)
  window_end   <- motif_center + window_bp

  GRanges(seqnames    = motif_df$motif_chr,
          ranges      = IRanges(start = window_start, end = window_end),
          loop_id     = motif_df$loop_id,
          anchor_num  = motif_df$anchor_num,
          direction   = motif_df$direction,
          motif_start = motif_df$motif_start,
          motif_end   = motif_df$motif_end)
}

annotate_chip_overlaps_local <- function(window_gr, k27ac_gr, k27me3_gr,
                                         k4me1_gr, k4me3_gr, bivalent_gr,
                                         ctcf_chip_gr) {
  data.frame(
    H3K27ac_overlap           = countOverlaps(window_gr, k27ac_gr, ignore.strand = TRUE) > 0,
    H3K27me3_overlap          = countOverlaps(window_gr, k27me3_gr, ignore.strand = TRUE) > 0,
    H3K4me1_overlap           = countOverlaps(window_gr, k4me1_gr, ignore.strand = TRUE) > 0,
    H3K4me3_overlap           = countOverlaps(window_gr, k4me3_gr, ignore.strand = TRUE) > 0,
    Bivalent_Promoter_overlap = countOverlaps(window_gr, bivalent_gr, ignore.strand = TRUE) > 0,
    CTCF_chip_overlap         = countOverlaps(window_gr, ctcf_chip_gr, ignore.strand = TRUE) > 0
  )
}

compute_tss_distances <- function(window_gr, tss_gr) {
  d <- distanceToNearest(window_gr, tss_gr, ignore.strand = TRUE)
  dist_vec <- rep(NA_real_, length(window_gr))
  if (length(d) > 0) {
    dist_vec[queryHits(d)] <- mcols(d)$distance
  }
  return(dist_vec)
}

classify_anchor_type_local <- function(h3k27ac_overlap, h3k27me3_overlap,
                                       h3k4me1_overlap, h3k4me3_overlap,
                                       bivalent_overlap, ctcf_motif_overlap,
                                       distance_to_tss,
                                       tss_threshold = TSS_THRESHOLD) {
  n <- length(h3k27ac_overlap)
  anchor_type <- rep("Other", n)

  is_ap <- h3k4me3_overlap & !h3k27me3_overlap &
           !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_ap] <- "Active_Promoter"

  is_rp <- !is_ap & h3k27me3_overlap & !h3k27ac_overlap &
           !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_rp] <- "Repressed_Promoter"

  is_bv <- !is_ap & !is_rp & bivalent_overlap
  anchor_type[is_bv] <- "Bivalent_Promoter"

  is_pc <- !is_ap & !is_rp & !is_bv & h3k27me3_overlap &
           (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_pc] <- "Polycomb"

  is_ae <- !is_ap & !is_rp & !is_bv & !is_pc & h3k27ac_overlap &
           (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_ae] <- "Active_Enhancer"

  is_pe <- !is_ap & !is_rp & !is_bv & !is_pc & !is_ae &
           h3k4me1_overlap & !h3k27ac_overlap & !h3k27me3_overlap &
           (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_pe] <- "Poised_Enhancer"

  is_ctcf <- !is_ap & !is_rp & !is_bv & !is_pc & !is_ae & !is_pe & ctcf_motif_overlap
  anchor_type[is_ctcf] <- "CTCF_Site"

  return(anchor_type)
}

classify_loop_type_local <- function(anchor1_type, anchor2_type) {
  n <- length(anchor1_type)
  loop_type <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    t1 <- anchor1_type[i]
    t2 <- anchor2_type[i]
    if (t1 == t2) {
      loop_type[i] <- sprintf("%s-%s", t1, t2)
    } else {
      pos1 <- match(t1, ANCHOR_TYPE_ORDER)
      pos2 <- match(t2, ANCHOR_TYPE_ORDER)
      if (pos1 < pos2) {
        loop_type[i] <- sprintf("%s-%s", t1, t2)
      } else {
        loop_type[i] <- sprintf("%s-%s", t2, t1)
      }
    }
  }
  return(loop_type)
}

add_ctcf_anchor_count <- function(loops_df) {
  loops_df %>%
    dplyr::mutate(ctcf_anchor_count = as.integer(anchor1_CTCF_motif_overlap) +
                                      as.integer(anchor2_CTCF_motif_overlap))
}

extract_other_end_types <- function(loops_df) {
  loops_1ctcf <- loops_df %>%
    dplyr::filter(ctcf_anchor_count == 1) %>%
    dplyr::mutate(
      ctcf_anchor    = ifelse(anchor1_CTCF_motif_overlap, 1L, 2L),
      ctcf_end_type  = ifelse(ctcf_anchor == 1L, anchor1_type_extended, anchor2_type_extended),
      other_end_type = ifelse(ctcf_anchor == 1L, anchor2_type_extended, anchor1_type_extended)
    ) %>%
    dplyr::select(loop_id, direction, ctcf_anchor_count, ctcf_anchor,
                  ctcf_end_type, other_end_type, loop_type_extended)

  loops_2ctcf <- loops_df %>%
    dplyr::filter(ctcf_anchor_count == 2) %>%
    dplyr::select(loop_id, direction, ctcf_anchor_count,
                  anchor1_type_extended, anchor2_type_extended, loop_type_extended)

  list(one_ctcf = loops_1ctcf, two_ctcf = loops_2ctcf)
}

reannotate_at_window <- function(motif_df, window_bp, k27ac_gr, k27me3_gr,
                                  k4me1_gr, k4me3_gr, bivalent_gr,
                                  ctcf_chip_gr, tss_gr) {
  windows_gr <- build_motif_windows(motif_df, window_bp)
  chip <- annotate_chip_overlaps_local(windows_gr, k27ac_gr, k27me3_gr,
                                       k4me1_gr, k4me3_gr, bivalent_gr, ctcf_chip_gr)
  tss_dist <- compute_tss_distances(windows_gr, tss_gr)
  types <- classify_anchor_type_local(
    chip$H3K27ac_overlap, chip$H3K27me3_overlap,
    chip$H3K4me1_overlap, chip$H3K4me3_overlap,
    chip$Bivalent_Promoter_overlap,
    rep(TRUE, nrow(chip)),
    tss_dist
  )

  tibble(
    loop_id          = mcols(windows_gr)$loop_id,
    anchor_num       = mcols(windows_gr)$anchor_num,
    direction        = mcols(windows_gr)$direction,
    motif_start      = mcols(windows_gr)$motif_start,
    motif_end        = mcols(windows_gr)$motif_end,
    window_chr       = as.character(seqnames(windows_gr)),
    window_start     = start(windows_gr),
    window_end       = end(windows_gr),
    H3K27ac_overlap  = chip$H3K27ac_overlap,
    H3K27me3_overlap = chip$H3K27me3_overlap,
    H3K4me1_overlap  = chip$H3K4me1_overlap,
    H3K4me3_overlap  = chip$H3K4me3_overlap,
    Bivalent_overlap = chip$Bivalent_Promoter_overlap,
    CTCF_chip_overlap = chip$CTCF_chip_overlap,
    distance_to_tss  = tss_dist,
    type_motif_window = types
  )
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

plot_2a_motif_count_per_anchor <- function(anchor_summary, output_dir) {
  plot_df <- anchor_summary %>%
    dplyr::mutate(
      n_motifs_bin = factor(dplyr::case_when(
        n_motifs == 0 ~ "0",
        n_motifs == 1 ~ "1",
        n_motifs == 2 ~ "2",
        TRUE          ~ "3+"
      ), levels = c("0", "1", "2", "3+")),
      type_10kb = factor(type_10kb, levels = ANCHOR_TYPE_ORDER)
    )

  p <- ggplot(plot_df, aes(x = type_10kb, fill = n_motifs_bin)) +
    geom_bar(position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_brewer(palette = "Blues", name = "CTCF motifs\nper anchor") +
    labs(x = "Anchor Type (10 kb)", y = "Number of Anchors",
         title = "CTCF Motif Count per Anchor by Chromatin State") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_multiformat_ggplot(p, file.path(output_dir, "2a_motif_count_per_anchor"),
                          width = 10, height = 7)
}

plot_2a_state_concordance <- function(anchor_summary, output_dir) {
  motif_anchors <- anchor_summary %>%
    dplyr::filter(n_motifs > 0, !is.na(type_motif_consensus))

  if (nrow(motif_anchors) < 3) {
    cat("  Skipping concordance plot: too few motif-containing anchors\n")
    return(invisible(NULL))
  }

  concordance_long <- dplyr::bind_rows(
    motif_anchors %>%
      dplyr::transmute(type = factor(type_10kb, levels = ANCHOR_TYPE_ORDER),
                       source = "10 kb Anchor"),
    motif_anchors %>%
      dplyr::transmute(type = factor(type_motif_consensus, levels = ANCHOR_TYPE_ORDER),
                       source = "Motif Window")
  ) %>%
    dplyr::mutate(source = factor(source, levels = c("10 kb Anchor", "Motif Window")))

  n_concord <- sum(motif_anchors$is_concordant, na.rm = TRUE)
  n_total <- nrow(motif_anchors)
  pct <- round(100 * n_concord / n_total, 1)

  p <- ggplot(concordance_long, aes(x = source, fill = type)) +
    geom_bar(position = "fill", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = ANCHOR_COLORS, name = "Anchor Type") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "", y = "Proportion",
         title = sprintf("Chromatin State: 10 kb vs Motif Window (%d bp)",
                         MOTIF_WINDOW_BP * 2),
         subtitle = sprintf("Concordance: %d / %d anchors (%.1f%%)",
                            n_concord, n_total, pct)) +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, "2a_state_concordance"),
                          width = 8, height = 7)
}

plot_2a_no_motif_anchors <- function(anchor_summary, output_dir) {
  no_motif <- anchor_summary %>%
    dplyr::filter(n_motifs == 0) %>%
    dplyr::count(type_10kb) %>%
    dplyr::mutate(type_10kb = factor(type_10kb, levels = ANCHOR_TYPE_ORDER))

  if (nrow(no_motif) == 0) {
    cat("  Skipping no-motif plot: all anchors have CTCF motifs\n")
    return(invisible(NULL))
  }

  p <- ggplot(no_motif, aes(x = type_10kb, y = n, fill = type_10kb)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
    labs(x = "Anchor Type (10 kb)", y = "Count",
         title = "Anchors Without Any CTCF Motif",
         subtitle = sprintf("%d / %d total anchors", sum(no_motif$n),
                            nrow(anchor_summary))) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_multiformat_ggplot(p, file.path(output_dir, "2a_no_motif_anchor_types"),
                          width = 8, height = 6)
}

plot_2b_ctcf_filter_comparison <- function(dist_all, dist_filtered, output_dir) {
  plot_df <- dplyr::bind_rows(dist_all, dist_filtered)

  top_types <- plot_df %>%
    dplyr::group_by(loop_type_extended) %>%
    dplyr::summarize(total = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::slice_head(n = 12) %>%
    dplyr::pull(loop_type_extended)

  plot_df <- plot_df %>%
    dplyr::mutate(
      loop_type_plot = ifelse(loop_type_extended %in% top_types,
                              loop_type_extended, "Other combinations"),
      filter_label = factor(filter_label, levels = c("All loops", "≥1 CTCF anchor"))
    )

  type_colors <- c(setNames(brewer.pal(12, "Set3"), top_types), "Other combinations" = "#cccccc")

  p <- ggplot(plot_df, aes(x = filter_label, y = n, fill = loop_type_plot)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = type_colors, name = "Loop Type") +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~direction, labeller = labeller(
      direction = c("up_in_mutant" = "Up in Mutant", "down_in_mutant" = "Down in Mutant"))) +
    labs(x = "", y = "Proportion",
         title = "Loop Type Distribution: All vs CTCF-Filtered") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  save_multiformat_ggplot(p, file.path(output_dir, "2b_ctcf_filter_comparison"),
                          width = 12, height = 7)
}

plot_2b_ctcf_count_distribution <- function(loops_with_count, output_dir) {
  count_summary <- loops_with_count %>%
    dplyr::count(direction, ctcf_anchor_count) %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(
      pct = n / sum(n),
      label = sprintf("%d\n(%.0f%%)", n, 100 * pct),
      ctcf_anchor_count = factor(ctcf_anchor_count)
    ) %>%
    dplyr::ungroup()

  p <- ggplot(count_summary, aes(x = direction, y = n, fill = ctcf_anchor_count)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
    geom_text(aes(label = label), position = position_fill(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = CTCF_COUNT_COLORS,
                      name = "CTCF motif\nanchors",
                      labels = c("0" = "Neither", "1" = "One", "2" = "Both")) +
    scale_x_discrete(labels = c("up_in_mutant" = "Up", "down_in_mutant" = "Down")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Loop Direction", y = "Proportion",
         title = "Loops by Number of CTCF-Motif-Containing Anchors") +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, "2b_ctcf_count_distribution"),
                          width = 7, height = 6)
}

plot_2c_other_end_types <- function(one_ctcf_df, output_dir) {
  if (nrow(one_ctcf_df) < 3) {
    cat("  Skipping other-end plot: too few 1-CTCF loops\n")
    return(invisible(NULL))
  }

  plot_df <- one_ctcf_df %>%
    dplyr::mutate(other_end_type = factor(other_end_type, levels = rev(ANCHOR_TYPE_ORDER)))

  p <- ggplot(plot_df, aes(x = direction, fill = other_end_type)) +
    geom_bar(position = "fill", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = ANCHOR_COLORS, name = "Other-End\nAnchor Type") +
    scale_x_discrete(labels = c("up_in_mutant" = "Up", "down_in_mutant" = "Down")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Loop Direction", y = "Proportion",
         title = "Other-End Anchor Type for 1-CTCF Loops",
         subtitle = sprintf("n = %d loops with exactly one CTCF-motif anchor",
                            nrow(one_ctcf_df))) +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, "2c_other_end_types"),
                          width = 8, height = 7)
}

plot_2c_both_ctcf_heatmap <- function(two_ctcf_df, output_dir) {
  if (nrow(two_ctcf_df) < 3) {
    cat("  Skipping both-CTCF heatmap: too few 2-CTCF loops\n")
    return(invisible(NULL))
  }

  pair_counts <- two_ctcf_df %>%
    dplyr::count(anchor1_type_extended, anchor2_type_extended) %>%
    tidyr::complete(
      anchor1_type_extended = ANCHOR_TYPE_ORDER,
      anchor2_type_extended = ANCHOR_TYPE_ORDER,
      fill = list(n = 0L)
    ) %>%
    dplyr::mutate(
      anchor1_type_extended = factor(anchor1_type_extended, levels = ANCHOR_TYPE_ORDER),
      anchor2_type_extended = factor(anchor2_type_extended, levels = rev(ANCHOR_TYPE_ORDER))
    )

  p <- ggplot(pair_counts, aes(x = anchor1_type_extended, y = anchor2_type_extended, fill = n)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = 4) +
    scale_fill_gradient(low = "white", high = "#b2182b", name = "Count") +
    labs(x = "Anchor 1 Type", y = "Anchor 2 Type",
         title = "Anchor Type Pairs for 2-CTCF Loops",
         subtitle = sprintf("n = %d loops with both anchors containing a CTCF motif",
                            nrow(two_ctcf_df))) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_multiformat_ggplot(p, file.path(output_dir, "2c_both_ctcf_heatmap"),
                          width = 9, height = 8)
}

# =============================================================================
# 5. MAIN ANALYSIS
# =============================================================================

run_analysis <- function(timepoint) {
  cat(sprintf("\n========== CTCF Motif Analysis: %s ==========\n\n", timepoint))

  output_dir <- file.path(BASE_DIR, "loops/output/ctcf_motif_analysis", timepoint)
  tables_dir <- file.path(output_dir, "tables")
  stats_dir  <- file.path(output_dir, "statistics")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir,  recursive = TRUE, showWarnings = FALSE)

  # --- Step 1: Load annotated loops ---
  cat("Step 1: Loading annotated loops...\n")
  loops_df <- load_loops(timepoint)

  # --- Step 2: Load genomic data ---
  cat("\nStep 2: Loading ChIP-seq peaks and CTCF motifs...\n")
  peaks <- PEAK_FILES[[timepoint]]
  k27ac_gr    <- load_chip_peaks_local(peaks$h3k27ac, "H3K27ac")
  k27me3_gr   <- load_chip_peaks_local(peaks$h3k27me3, "H3K27me3")
  k4me1_gr    <- load_chip_peaks_local(peaks$h3k4me1, "H3K4me1")
  k4me3_gr    <- load_chip_peaks_local(peaks$h3k4me3, "H3K4me3")
  bivalent_gr <- load_chip_peaks_local(peaks$bivalent, "Bivalent")
  ctcf_chip_gr <- load_chip_peaks_local(CTCF_CHIP_FILE, "CTCF_ChIP")
  motif_gr    <- load_ctcf_motifs()

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  all_genes <- genes(txdb)
  tss_gr <- resize(all_genes, width = 1, fix = "start")
  cat(sprintf("  TSS reference: %d gene starts\n", length(tss_gr)))

  # --- Step 3: Find motifs in anchors (2a) ---
  cat("\nStep 3: Finding CTCF motifs within loop anchors...\n")
  anchor1_gr <- anchors_to_granges(loops_df, 1)
  anchor2_gr <- anchors_to_granges(loops_df, 2)
  motifs_a1 <- find_motifs_in_anchors(anchor1_gr, motif_gr)
  motifs_a2 <- find_motifs_in_anchors(anchor2_gr, motif_gr)
  all_motifs <- dplyr::bind_rows(motifs_a1, motifs_a2)

  n_anchors_total <- nrow(loops_df) * 2
  n_anchors_with_motif <- all_motifs %>%
    dplyr::distinct(loop_id, anchor_num) %>% nrow()
  cat(sprintf("  Motif-anchor pairs: %d\n", nrow(all_motifs)))
  cat(sprintf("  Anchors with motif: %d / %d (%.1f%%)\n",
              n_anchors_with_motif, n_anchors_total,
              100 * n_anchors_with_motif / n_anchors_total))

  # --- Step 4: Re-annotate motif-centered windows (2a) ---
  cat(sprintf("\nStep 4: Re-annotating motif-centered windows (+/-%d bp)...\n", MOTIF_WINDOW_BP))

  window_annotation_df <- NULL
  if (nrow(all_motifs) > 0) {
    window_annotation_df <- reannotate_at_window(
      all_motifs, MOTIF_WINDOW_BP, k27ac_gr, k27me3_gr,
      k4me1_gr, k4me3_gr, bivalent_gr, ctcf_chip_gr, tss_gr)
    cat(sprintf("  Re-annotated %d motif windows\n", nrow(window_annotation_df)))
  } else {
    cat("  No motifs found in anchors; skipping re-annotation\n")
  }

  # --- Step 5: Per-anchor consensus (2a) ---
  cat("\nStep 5: Computing per-anchor motif consensus...\n")

  anchor_types_10kb <- dplyr::bind_rows(
    loops_df %>% dplyr::transmute(loop_id, anchor_num = 1L, type_10kb = anchor1_type_extended,
                                  direction),
    loops_df %>% dplyr::transmute(loop_id, anchor_num = 2L, type_10kb = anchor2_type_extended,
                                  direction)
  )

  if (!is.null(window_annotation_df) && nrow(window_annotation_df) > 0) {
    per_anchor_motif <- window_annotation_df %>%
      dplyr::group_by(loop_id, anchor_num) %>%
      dplyr::summarize(
        n_motifs           = n(),
        type_motif_consensus = ANCHOR_TYPE_ORDER[min(match(type_motif_window, ANCHOR_TYPE_ORDER))],
        n_chip_bound       = sum(CTCF_chip_overlap),
        .groups = "drop"
      )

    anchor_summary <- anchor_types_10kb %>%
      dplyr::left_join(per_anchor_motif, by = c("loop_id", "anchor_num")) %>%
      dplyr::mutate(
        n_motifs = dplyr::coalesce(as.integer(n_motifs), 0L),
        n_chip_bound = dplyr::coalesce(as.integer(n_chip_bound), 0L),
        is_concordant = (type_10kb == type_motif_consensus)
      )
  } else {
    anchor_summary <- anchor_types_10kb %>%
      dplyr::mutate(n_motifs = 0L, type_motif_consensus = NA_character_,
                    n_chip_bound = 0L, is_concordant = NA)
  }

  n_with <- sum(anchor_summary$n_motifs > 0)
  n_concordant <- sum(anchor_summary$is_concordant, na.rm = TRUE)
  if (n_with > 0) {
    cat(sprintf("  Concordance (10kb == motif window): %d / %d (%.1f%%)\n",
                n_concordant, n_with, 100 * n_concordant / n_with))
  }

  # --- Step 6: Window size sensitivity (2a) ---
  cat("\nStep 6: Window size sensitivity analysis...\n")

  sensitivity_results <- tibble(window_bp = integer(), n_windows = integer(),
                                concordance_rate = numeric())

  if (nrow(all_motifs) > 0) {
    for (wb in c(250L, 500L, 1000L)) {
      sens_windows <- reannotate_at_window(
        all_motifs, wb, k27ac_gr, k27me3_gr,
        k4me1_gr, k4me3_gr, bivalent_gr, ctcf_chip_gr, tss_gr)

      sens_consensus <- sens_windows %>%
        dplyr::group_by(loop_id, anchor_num) %>%
        dplyr::summarize(
          type_motif = ANCHOR_TYPE_ORDER[min(match(type_motif_window, ANCHOR_TYPE_ORDER))],
          .groups = "drop"
        ) %>%
        dplyr::left_join(anchor_types_10kb %>% dplyr::select(loop_id, anchor_num, type_10kb),
                         by = c("loop_id", "anchor_num"))

      conc_rate <- mean(sens_consensus$type_motif == sens_consensus$type_10kb, na.rm = TRUE)
      sensitivity_results <- dplyr::bind_rows(sensitivity_results,
        tibble(window_bp = as.integer(wb * 2), n_windows = nrow(sens_windows),
               concordance_rate = conc_rate))
      cat(sprintf("  +/-%d bp (total %d bp): concordance = %.1f%%\n",
                  wb, wb * 2, 100 * conc_rate))
    }
  }

  # --- Step 7: CTCF motif filter (2b) ---
  cat("\nStep 7: CTCF motif filter analysis...\n")

  loops_counted <- add_ctcf_anchor_count(loops_df)
  loops_filtered <- loops_counted %>% dplyr::filter(ctcf_anchor_count >= 1)

  cat(sprintf("  All loops: %d | >=1 CTCF anchor: %d | Reduction: %.1f%%\n",
              nrow(loops_counted), nrow(loops_filtered),
              100 * (1 - nrow(loops_filtered) / nrow(loops_counted))))

  compute_dist <- function(df, label) {
    df %>%
      dplyr::count(direction, loop_type_extended) %>%
      dplyr::group_by(direction) %>%
      dplyr::mutate(pct = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(filter_label = label)
  }

  dist_all      <- compute_dist(loops_counted, "All loops")
  dist_filtered <- compute_dist(loops_filtered, ">=1 CTCF anchor")

  fisher_2b_results <- tibble(loop_type = character(), odds_ratio = numeric(),
                              ci_low = numeric(), ci_high = numeric(),
                              p_value = numeric())

  all_types <- unique(c(dist_all$loop_type_extended, dist_filtered$loop_type_extended))
  for (lt in all_types) {
    n_filt_type  <- sum(loops_filtered$loop_type_extended == lt)
    n_filt_other <- nrow(loops_filtered) - n_filt_type
    n_excl_type  <- sum(loops_counted$ctcf_anchor_count == 0 & loops_counted$loop_type_extended == lt)
    n_excl_other <- sum(loops_counted$ctcf_anchor_count == 0) - n_excl_type

    if (n_filt_type + n_excl_type == 0) next

    ft <- fisher.test(matrix(c(n_filt_type, n_excl_type, n_filt_other, n_excl_other), nrow = 2))
    fisher_2b_results <- dplyr::bind_rows(fisher_2b_results, tibble(
      loop_type = lt, odds_ratio = ft$estimate,
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
      p_value = ft$p.value))
  }

  # --- Step 8: Other-end analysis (2c) ---
  cat("\nStep 8: Other-end anchor type analysis...\n")

  other_end_data <- extract_other_end_types(loops_counted)
  one_ctcf_df <- other_end_data$one_ctcf
  two_ctcf_df <- other_end_data$two_ctcf

  cat(sprintf("  1-CTCF loops: %d | 2-CTCF loops: %d\n",
              nrow(one_ctcf_df), nrow(two_ctcf_df)))

  fisher_2c_results <- tibble(other_end_type = character(), odds_ratio = numeric(),
                              ci_low = numeric(), ci_high = numeric(),
                              p_value = numeric())

  if (nrow(one_ctcf_df) >= 5) {
    for (ot in unique(one_ctcf_df$other_end_type)) {
      n_up_type   <- sum(one_ctcf_df$direction == "up_in_mutant" & one_ctcf_df$other_end_type == ot)
      n_up_other  <- sum(one_ctcf_df$direction == "up_in_mutant") - n_up_type
      n_dn_type   <- sum(one_ctcf_df$direction == "down_in_mutant" & one_ctcf_df$other_end_type == ot)
      n_dn_other  <- sum(one_ctcf_df$direction == "down_in_mutant") - n_dn_type

      ft <- fisher.test(matrix(c(n_up_type, n_dn_type, n_up_other, n_dn_other), nrow = 2))
      fisher_2c_results <- dplyr::bind_rows(fisher_2c_results, tibble(
        other_end_type = ot, odds_ratio = ft$estimate,
        ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
        p_value = ft$p.value))
    }
  }

  # --- Step 9: Write statistics ---
  cat("\nStep 9: Writing statistics...\n")

  stats_lines <- c(
    sprintf("=== CTCF Motif Analysis: %s ===", timepoint),
    sprintf("Date: %s", Sys.Date()),
    "",
    "--- 2a: Re-Centering on CTCF Motifs ---",
    sprintf("Total loops: %d", nrow(loops_df)),
    sprintf("Total anchors: %d", n_anchors_total),
    sprintf("Anchors with >=1 CTCF motif: %d (%.1f%%)",
            n_anchors_with_motif, 100 * n_anchors_with_motif / n_anchors_total),
    sprintf("Anchors with 0 CTCF motifs: %d (%.1f%%)",
            n_anchors_total - n_anchors_with_motif,
            100 * (1 - n_anchors_with_motif / n_anchors_total)),
    ""
  )

  motif_count_tab <- anchor_summary %>%
    dplyr::mutate(n_bin = dplyr::case_when(n_motifs == 0 ~ "0", n_motifs == 1 ~ "1",
                                           n_motifs == 2 ~ "2", TRUE ~ "3+")) %>%
    dplyr::count(n_bin) %>%
    dplyr::mutate(pct = sprintf("%.1f%%", 100 * n / sum(n)))
  stats_lines <- c(stats_lines, "Motif count distribution:")
  for (i in seq_len(nrow(motif_count_tab))) {
    stats_lines <- c(stats_lines, sprintf("  %s motifs: %d (%s)",
                                          motif_count_tab$n_bin[i],
                                          motif_count_tab$n[i],
                                          motif_count_tab$pct[i]))
  }

  if (n_with > 0) {
    stats_lines <- c(stats_lines, "",
      sprintf("Concordance (10kb type == motif-window consensus): %d / %d (%.1f%%)",
              n_concordant, n_with, 100 * n_concordant / n_with))

    transition_tab <- anchor_summary %>%
      dplyr::filter(n_motifs > 0, !is.na(type_motif_consensus)) %>%
      dplyr::count(type_10kb, type_motif_consensus) %>%
      dplyr::arrange(type_10kb, dplyr::desc(n))
    stats_lines <- c(stats_lines, "", "Type transition matrix (10kb -> motif window):")
    for (i in seq_len(nrow(transition_tab))) {
      stats_lines <- c(stats_lines, sprintf("  %s -> %s: %d",
                                            transition_tab$type_10kb[i],
                                            transition_tab$type_motif_consensus[i],
                                            transition_tab$n[i]))
    }
  }

  if (nrow(sensitivity_results) > 0) {
    stats_lines <- c(stats_lines, "", "Window size sensitivity:")
    for (i in seq_len(nrow(sensitivity_results))) {
      stats_lines <- c(stats_lines, sprintf("  %d bp window: concordance = %.1f%% (%d windows)",
                                            sensitivity_results$window_bp[i],
                                            100 * sensitivity_results$concordance_rate[i],
                                            sensitivity_results$n_windows[i]))
    }
  }

  n0 <- sum(loops_counted$ctcf_anchor_count == 0)
  n1 <- sum(loops_counted$ctcf_anchor_count == 1)
  n2 <- sum(loops_counted$ctcf_anchor_count == 2)
  stats_lines <- c(stats_lines, "",
    "--- 2b: CTCF Motif Filter ---",
    sprintf("All loops: %d", nrow(loops_counted)),
    sprintf("  0-CTCF anchor loops: %d (%.1f%%) | up=%d down=%d", n0, 100 * n0 / nrow(loops_counted),
            sum(loops_counted$ctcf_anchor_count == 0 & loops_counted$direction == "up_in_mutant"),
            sum(loops_counted$ctcf_anchor_count == 0 & loops_counted$direction == "down_in_mutant")),
    sprintf("  1-CTCF anchor loops: %d (%.1f%%) | up=%d down=%d", n1, 100 * n1 / nrow(loops_counted),
            sum(loops_counted$ctcf_anchor_count == 1 & loops_counted$direction == "up_in_mutant"),
            sum(loops_counted$ctcf_anchor_count == 1 & loops_counted$direction == "down_in_mutant")),
    sprintf("  2-CTCF anchor loops: %d (%.1f%%) | up=%d down=%d", n2, 100 * n2 / nrow(loops_counted),
            sum(loops_counted$ctcf_anchor_count == 2 & loops_counted$direction == "up_in_mutant"),
            sum(loops_counted$ctcf_anchor_count == 2 & loops_counted$direction == "down_in_mutant")),
    sprintf(">=1 CTCF filtered: %d (%.1f%% retained)", nrow(loops_filtered),
            100 * nrow(loops_filtered) / nrow(loops_counted))
  )

  if (nrow(fisher_2b_results) > 0) {
    stats_lines <- c(stats_lines, "", "Fisher tests (filtered vs excluded, per loop type):")
    sig_2b <- fisher_2b_results %>% dplyr::filter(p_value < 0.05) %>% dplyr::arrange(p_value)
    if (nrow(sig_2b) > 0) {
      for (i in seq_len(nrow(sig_2b))) {
        stats_lines <- c(stats_lines, sprintf("  %s: OR=%.2f [%.2f-%.2f] p=%.4g",
                                              sig_2b$loop_type[i], sig_2b$odds_ratio[i],
                                              sig_2b$ci_low[i], sig_2b$ci_high[i],
                                              sig_2b$p_value[i]))
      }
    } else {
      stats_lines <- c(stats_lines, "  No significant enrichments (p < 0.05)")
    }
  }

  stats_lines <- c(stats_lines, "",
    "--- 2c: Other-End Anchor Types ---",
    sprintf("1-CTCF loops: %d", nrow(one_ctcf_df)),
    sprintf("2-CTCF loops: %d", nrow(two_ctcf_df))
  )

  if (nrow(one_ctcf_df) > 0) {
    oe_tab <- one_ctcf_df %>% dplyr::count(other_end_type) %>% dplyr::arrange(dplyr::desc(n))
    stats_lines <- c(stats_lines, "", "Other-end type distribution (1-CTCF loops):")
    for (i in seq_len(nrow(oe_tab))) {
      stats_lines <- c(stats_lines, sprintf("  %s: %d (%.1f%%)",
                                            oe_tab$other_end_type[i], oe_tab$n[i],
                                            100 * oe_tab$n[i] / nrow(one_ctcf_df)))
    }
  }

  if (nrow(fisher_2c_results) > 0) {
    stats_lines <- c(stats_lines, "", "Fisher tests (other-end type, up vs down):")
    for (i in seq_len(nrow(fisher_2c_results))) {
      stats_lines <- c(stats_lines, sprintf("  %s: OR=%.2f [%.2f-%.2f] p=%.4g",
                                            fisher_2c_results$other_end_type[i],
                                            fisher_2c_results$odds_ratio[i],
                                            fisher_2c_results$ci_low[i],
                                            fisher_2c_results$ci_high[i],
                                            fisher_2c_results$p_value[i]))
    }
  }

  if (nrow(two_ctcf_df) > 0) {
    pair_tab <- two_ctcf_df %>%
      dplyr::count(anchor1_type_extended, anchor2_type_extended) %>%
      dplyr::arrange(dplyr::desc(n))
    stats_lines <- c(stats_lines, "", "2-CTCF loop anchor type pairs:")
    for (i in seq_len(min(nrow(pair_tab), 15))) {
      stats_lines <- c(stats_lines, sprintf("  %s - %s: %d",
                                            pair_tab$anchor1_type_extended[i],
                                            pair_tab$anchor2_type_extended[i],
                                            pair_tab$n[i]))
    }
  }

  stats_file <- file.path(stats_dir, "ctcf_motif_statistics.txt")
  writeLines(stats_lines, stats_file)
  cat(sprintf("  Statistics written to: %s\n", stats_file))

  # --- Step 10: Write tables ---
  cat("\nStep 10: Writing tables...\n")

  if (!is.null(window_annotation_df) && nrow(window_annotation_df) > 0) {
    readr::write_tsv(window_annotation_df,
                     file.path(tables_dir, "motif_window_annotations.tsv"))
  }

  readr::write_tsv(anchor_summary, file.path(tables_dir, "anchor_motif_summary.tsv"))
  readr::write_tsv(loops_counted, file.path(tables_dir, "loops_ctcf_count.tsv"))
  readr::write_tsv(dplyr::bind_rows(dist_all, dist_filtered),
                   file.path(tables_dir, "loop_type_filtered_comparison.tsv"))

  if (nrow(one_ctcf_df) > 0) {
    readr::write_tsv(one_ctcf_df, file.path(tables_dir, "other_end_types.tsv"))
  }
  if (nrow(fisher_2b_results) > 0) {
    readr::write_tsv(fisher_2b_results, file.path(tables_dir, "fisher_2b_results.tsv"))
  }
  if (nrow(fisher_2c_results) > 0) {
    readr::write_tsv(fisher_2c_results, file.path(tables_dir, "fisher_2c_results.tsv"))
  }

  cat(sprintf("  Tables written to: %s\n", tables_dir))

  # --- Step 11: Generate plots ---
  cat("\nStep 11: Generating plots...\n")

  plot_2a_motif_count_per_anchor(anchor_summary, output_dir)
  plot_2a_state_concordance(anchor_summary, output_dir)
  plot_2a_no_motif_anchors(anchor_summary, output_dir)
  plot_2b_ctcf_filter_comparison(dist_all, dist_filtered, output_dir)
  plot_2b_ctcf_count_distribution(loops_counted, output_dir)
  plot_2c_other_end_types(one_ctcf_df, output_dir)
  plot_2c_both_ctcf_heatmap(two_ctcf_df, output_dir)

  cat(sprintf("\n  All plots saved to: %s\n", output_dir))
  cat(sprintf("\n========== %s complete ==========\n", timepoint))
}

# =============================================================================
# 6. DISPATCH
# =============================================================================

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  run_analysis(tp)
}

cat("\nDone.\n")
