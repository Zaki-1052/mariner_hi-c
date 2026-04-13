# scripts/cpg_ctcf_analysis.R
# CpG Island / GC Content Analysis at CTCF-Anchored Loop Anchors
# Author: Zakir Alibhai
# Date: 2026-04-12
#
# Purpose:
#   Test whether CpG islands and GC content differ at CTCF anchors of
#   lost vs gained chromatin loops in BAP1-KO. Connects to Brad Bernstein's
#   IDH mutant glioma work where repressive chromatin disrupts CTCF binding
#   at CpG-rich regions. Jesse Dixon suggested this analysis (04/10/26).
#
# Biological hypothesis:
#   Lost CTCF loops are longer than gained (810kb vs 325kb). If Polycomb
#   repression targets CpG island-associated CTCF sites, lost anchors may
#   be enriched for CpG islands relative to gained anchors.
#
# Usage:
#   Rscript scripts/cpg_ctcf_analysis.R
#
# Output:
#   output/cpg_ctcf_analysis/
#     ├── tables/
#     ├── plots/
#     └── cpg_ctcf_analysis_summary.txt

# ==============================================================================
# SECTION 1: LIBRARIES & CONFIGURATION
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("scripts/utils/multi_format_output.R")

cat("=", rep("=", 79) |> paste(collapse = ""), "\n")
cat("CpG Island / GC Content at CTCF-Anchored Loop Anchors\n")
cat("=", rep("=", 79) |> paste(collapse = ""), "\n\n")

OUTPUT_DIR <- "output/cpg_ctcf_analysis"
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")

INPUT_FILE <- "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"
CPG_ISLAND_FILE <- "biomodal/downstream/modality/mm10/gencode.vM25.mouse.cpg_islands.annotation.bed"
CTCF_PEAK_FILE <- "peaks/CTCF.bed"

COLORS <- list(
  lost = "#d73027",
  gained = "#4575b4",
  lost_light = "#f4a582",
  gained_light = "#92c5de",
  ns = "grey60"
)

DISTANCE_BINS <- c(0, 200e3, 500e3, 1e6, Inf)
DISTANCE_LABELS <- c("<200kb", "200-500kb", "500kb-1Mb", ">1Mb")

STANDARD_CHROMS <- paste0("chr", c(1:19, "X"))

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 2: HELPER FUNCTIONS
# ==============================================================================

#' Load CpG island BED with chr-prefix correction
load_cpg_islands <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$Chromosome <- paste0("chr", df$Chromosome)
  gr <- GRanges(seqnames = df$Chromosome,
                ranges = IRanges(start = df$Start, end = df$End))
  gr <- gr[seqnames(gr) %in% STANDARD_CHROMS]
  gr <- sort(gr)
  cat(sprintf("  Loaded %d CpG islands (%d on standard chromosomes)\n",
              nrow(df), length(gr)))
  gr
}

#' Fisher's exact test for CpG island overlap enrichment
run_fisher_cpg <- function(test_cpg, test_total, ref_cpg, ref_total,
                           test_label, ref_label) {
  test_no_cpg <- test_total - test_cpg
  ref_no_cpg <- ref_total - ref_cpg

  mat <- matrix(c(test_cpg, ref_cpg,
                  test_no_cpg, ref_no_cpg),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("CpG_island", "No_CpG"),
                               c(test_label, ref_label)))

  fisher_result <- tryCatch({
    fisher.test(mat)
  }, error = function(e) {
    list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  })

  list(
    test_label = test_label,
    ref_label = ref_label,
    test_n = test_total,
    test_cpg = test_cpg,
    test_pct = 100 * test_cpg / test_total,
    ref_n = ref_total,
    ref_cpg = ref_cpg,
    ref_pct = 100 * ref_cpg / ref_total,
    odds_ratio = as.numeric(fisher_result$estimate),
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    p_value = fisher_result$p.value
  )
}

#' Compute GC content for a set of genomic ranges
compute_gc_content <- function(anchor_gr, genome) {
  seqs <- getSeq(genome, anchor_gr)
  gc_frac <- letterFrequency(seqs, letters = "GC", as.prob = TRUE)[, 1]
  gc_frac
}

#' Extract unique CTCF+ anchors from loops data frame
#' Returns a data frame with anchor coordinates, direction, and loop distance
extract_ctcf_anchors <- function(loops_df, ctcf_col_prefix = "CTCF_overlap") {
  a1_col <- paste0("anchor1_", ctcf_col_prefix)
  a2_col <- paste0("anchor2_", ctcf_col_prefix)

  # Anchor 1 entries where CTCF+
  a1 <- loops_df %>%
    filter(.data[[a1_col]] == TRUE) %>%
    transmute(chr = chr1, start = start1, end = end1,
              direction = direction, loop_distance = loop_distance,
              loop_id = loop_id, anchor_num = 1L)

  # Anchor 2 entries where CTCF+
  a2 <- loops_df %>%
    filter(.data[[a2_col]] == TRUE) %>%
    transmute(chr = chr2, start = start2, end = end2,
              direction = direction, loop_distance = loop_distance,
              loop_id = loop_id, anchor_num = 2L)

  anchors <- bind_rows(a1, a2)

  # Deduplicate: same coordinates + direction → keep one, take max loop_distance
  anchors <- anchors %>%
    group_by(chr, start, end, direction) %>%
    summarise(
      loop_distance = max(loop_distance),
      n_loops = n(),
      .groups = "drop"
    )

  anchors
}

# ==============================================================================
# SECTION 3: DATA LOADING & CTCF IDENTIFICATION
# ==============================================================================

cat("\n[1] Loading input data...\n")

loops <- read.table(INPUT_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d loops from extended annotation\n", nrow(loops)))
cat(sprintf("  Direction: %d lost (down_in_mutant), %d gained (up_in_mutant)\n",
            sum(loops$direction == "down_in_mutant"),
            sum(loops$direction == "up_in_mutant")))

cpg_islands <- load_cpg_islands(CPG_ISLAND_FILE)

ctcf_peaks <- read.table(CTCF_PEAK_FILE, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
ctcf_peaks_gr <- GRanges(seqnames = ctcf_peaks$V1,
                         ranges = IRanges(start = ctcf_peaks$V2, end = ctcf_peaks$V3))
ctcf_peaks_gr <- ctcf_peaks_gr[seqnames(ctcf_peaks_gr) %in% STANDARD_CHROMS]
cat(sprintf("  Loaded %d CTCF ChIP-seq peaks (for GC baseline)\n", length(ctcf_peaks_gr)))

genome <- BSgenome.Mmusculus.UCSC.mm10

# Define loop subsets
cat("\n[2] Identifying CTCF-anchored loops...\n")

# CTCF-CTCF loops (both anchors CTCF ChIP+)
ctcf_ctcf <- loops %>%
  filter(anchor1_CTCF_overlap == TRUE & anchor2_CTCF_overlap == TRUE)
ctcf_ctcf_lost <- ctcf_ctcf %>% filter(direction == "down_in_mutant")
ctcf_ctcf_gained <- ctcf_ctcf %>% filter(direction == "up_in_mutant")

# Any-CTCF loops (at least one anchor CTCF ChIP+)
any_ctcf <- loops %>%
  filter(anchor1_CTCF_overlap == TRUE | anchor2_CTCF_overlap == TRUE)
any_ctcf_lost <- any_ctcf %>% filter(direction == "down_in_mutant")
any_ctcf_gained <- any_ctcf %>% filter(direction == "up_in_mutant")

cat(sprintf("  CTCF-CTCF loops: %d lost, %d gained\n",
            nrow(ctcf_ctcf_lost), nrow(ctcf_ctcf_gained)))
cat(sprintf("  Any-CTCF loops:  %d lost, %d gained\n",
            nrow(any_ctcf_lost), nrow(any_ctcf_gained)))

# Extract unique CTCF+ anchors (ChIP-seq based)
anchors_chip <- extract_ctcf_anchors(loops, "CTCF_overlap")
anchors_chip_lost <- anchors_chip %>% filter(direction == "down_in_mutant")
anchors_chip_gained <- anchors_chip %>% filter(direction == "up_in_mutant")

cat(sprintf("  Unique CTCF anchors (ChIP): %d lost, %d gained\n",
            nrow(anchors_chip_lost), nrow(anchors_chip_gained)))

# Extract unique CTCF+ anchors (motif-based, sensitivity analysis)
anchors_motif <- extract_ctcf_anchors(loops, "CTCF_motif_overlap")
anchors_motif_lost <- anchors_motif %>% filter(direction == "down_in_mutant")
anchors_motif_gained <- anchors_motif %>% filter(direction == "up_in_mutant")

cat(sprintf("  Unique CTCF anchors (motif): %d lost, %d gained\n",
            nrow(anchors_motif_lost), nrow(anchors_motif_gained)))

# ==============================================================================
# SECTION 4: CpG ISLAND OVERLAP ANALYSIS
# ==============================================================================

cat("\n[3] Computing CpG island overlaps...\n")

# Helper: compute CpG overlap for a set of loops (any anchor in CpG)
compute_loop_cpg_overlap <- function(loop_df, cpg_gr) {
  a1_gr <- GRanges(seqnames = loop_df$chr1,
                   ranges = IRanges(start = loop_df$start1, end = loop_df$end1))
  a2_gr <- GRanges(seqnames = loop_df$chr2,
                   ranges = IRanges(start = loop_df$start2, end = loop_df$end2))
  a1_hit <- countOverlaps(a1_gr, cpg_gr) > 0
  a2_hit <- countOverlaps(a2_gr, cpg_gr) > 0
  any_hit <- a1_hit | a2_hit
  list(any = any_hit, a1 = a1_hit, a2 = a2_hit)
}

# Helper: compute CpG overlap for anchor data frame
compute_anchor_cpg_overlap <- function(anchor_df, cpg_gr) {
  gr <- GRanges(seqnames = anchor_df$chr,
                ranges = IRanges(start = anchor_df$start, end = anchor_df$end))
  countOverlaps(gr, cpg_gr) > 0
}

# Loop-level CpG overlaps
cpg_ctcf_ctcf_lost <- compute_loop_cpg_overlap(ctcf_ctcf_lost, cpg_islands)
cpg_ctcf_ctcf_gained <- compute_loop_cpg_overlap(ctcf_ctcf_gained, cpg_islands)
cpg_any_ctcf_lost <- compute_loop_cpg_overlap(any_ctcf_lost, cpg_islands)
cpg_any_ctcf_gained <- compute_loop_cpg_overlap(any_ctcf_gained, cpg_islands)

# Anchor-level CpG overlaps (ChIP-seq)
anchors_chip_lost$cpg_overlap <- compute_anchor_cpg_overlap(anchors_chip_lost, cpg_islands)
anchors_chip_gained$cpg_overlap <- compute_anchor_cpg_overlap(anchors_chip_gained, cpg_islands)

# Anchor-level CpG overlaps (motif)
anchors_motif_lost$cpg_overlap <- compute_anchor_cpg_overlap(anchors_motif_lost, cpg_islands)
anchors_motif_gained$cpg_overlap <- compute_anchor_cpg_overlap(anchors_motif_gained, cpg_islands)

cat("  CpG overlap rates (any anchor in CpG island):\n")
cat(sprintf("    CTCF-CTCF lost:  %.1f%% (%d/%d)\n",
            100 * sum(cpg_ctcf_ctcf_lost$any) / nrow(ctcf_ctcf_lost),
            sum(cpg_ctcf_ctcf_lost$any), nrow(ctcf_ctcf_lost)))
cat(sprintf("    CTCF-CTCF gained: %.1f%% (%d/%d)\n",
            100 * sum(cpg_ctcf_ctcf_gained$any) / nrow(ctcf_ctcf_gained),
            sum(cpg_ctcf_ctcf_gained$any), nrow(ctcf_ctcf_gained)))
cat(sprintf("    Any-CTCF lost:   %.1f%% (%d/%d)\n",
            100 * sum(cpg_any_ctcf_lost$any) / nrow(any_ctcf_lost),
            sum(cpg_any_ctcf_lost$any), nrow(any_ctcf_lost)))
cat(sprintf("    Any-CTCF gained: %.1f%% (%d/%d)\n",
            100 * sum(cpg_any_ctcf_gained$any) / nrow(any_ctcf_gained),
            sum(cpg_any_ctcf_gained$any), nrow(any_ctcf_gained)))

# --- Fisher's exact tests ---
cat("\n[4] Running Fisher's exact tests...\n")

fisher_results <- list()

# Test 1: CTCF-CTCF loop-level (ChIP)
fisher_results[[1]] <- run_fisher_cpg(
  sum(cpg_ctcf_ctcf_lost$any), nrow(ctcf_ctcf_lost),
  sum(cpg_ctcf_ctcf_gained$any), nrow(ctcf_ctcf_gained),
  "CTCF-CTCF_lost", "CTCF-CTCF_gained")

# Test 2: Any-CTCF loop-level (ChIP)
fisher_results[[2]] <- run_fisher_cpg(
  sum(cpg_any_ctcf_lost$any), nrow(any_ctcf_lost),
  sum(cpg_any_ctcf_gained$any), nrow(any_ctcf_gained),
  "Any-CTCF_lost", "Any-CTCF_gained")

# Test 3: Anchor-level (ChIP)
fisher_results[[3]] <- run_fisher_cpg(
  sum(anchors_chip_lost$cpg_overlap), nrow(anchors_chip_lost),
  sum(anchors_chip_gained$cpg_overlap), nrow(anchors_chip_gained),
  "CTCF_anchor_lost", "CTCF_anchor_gained")

# Tests 4-6: Repeat with motif-based definition (sensitivity)
# Need motif-based loop subsets
ctcf_motif_ctcf <- loops %>%
  filter(anchor1_CTCF_motif_overlap == TRUE & anchor2_CTCF_motif_overlap == TRUE)
ctcf_motif_lost <- ctcf_motif_ctcf %>% filter(direction == "down_in_mutant")
ctcf_motif_gained <- ctcf_motif_ctcf %>% filter(direction == "up_in_mutant")

cpg_motif_lost <- compute_loop_cpg_overlap(ctcf_motif_lost, cpg_islands)
cpg_motif_gained <- compute_loop_cpg_overlap(ctcf_motif_gained, cpg_islands)

fisher_results[[4]] <- run_fisher_cpg(
  sum(cpg_motif_lost$any), nrow(ctcf_motif_lost),
  sum(cpg_motif_gained$any), nrow(ctcf_motif_gained),
  "Motif-CTCF_lost", "Motif-CTCF_gained")

any_motif <- loops %>%
  filter(anchor1_CTCF_motif_overlap == TRUE | anchor2_CTCF_motif_overlap == TRUE)
any_motif_lost <- any_motif %>% filter(direction == "down_in_mutant")
any_motif_gained <- any_motif %>% filter(direction == "up_in_mutant")
cpg_any_motif_lost <- compute_loop_cpg_overlap(any_motif_lost, cpg_islands)
cpg_any_motif_gained <- compute_loop_cpg_overlap(any_motif_gained, cpg_islands)

fisher_results[[5]] <- run_fisher_cpg(
  sum(cpg_any_motif_lost$any), nrow(any_motif_lost),
  sum(cpg_any_motif_gained$any), nrow(any_motif_gained),
  "Any-motif_lost", "Any-motif_gained")

fisher_results[[6]] <- run_fisher_cpg(
  sum(anchors_motif_lost$cpg_overlap), nrow(anchors_motif_lost),
  sum(anchors_motif_gained$cpg_overlap), nrow(anchors_motif_gained),
  "Motif_anchor_lost", "Motif_anchor_gained")

fisher_df <- bind_rows(lapply(fisher_results, as_tibble))
fisher_df$fdr <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$significant <- fisher_df$fdr < 0.05
fisher_df$enrichment <- case_when(
  fisher_df$odds_ratio > 1 & fisher_df$significant ~ "enriched",
  fisher_df$odds_ratio < 1 & fisher_df$significant ~ "depleted",
  TRUE ~ "ns"
)
fisher_df$ctcf_definition <- ifelse(grepl("motif|Motif", fisher_df$test_label),
                                    "motif", "ChIP-seq")

cat("  Primary Fisher's test results:\n")
for (i in 1:nrow(fisher_df)) {
  row <- fisher_df[i, ]
  cat(sprintf("    %s vs %s: OR=%.2f [%.2f-%.2f], p=%.2e, FDR=%.2e %s\n",
              row$test_label, row$ref_label,
              row$odds_ratio, row$ci_lower, row$ci_upper,
              row$p_value, row$fdr,
              ifelse(row$significant, paste0("(", row$enrichment, ")"), "(ns)")))
}

# --- Distance-stratified Fisher's tests ---
cat("\n[5] Distance-stratified CpG enrichment (controlling for loop length)...\n")

ctcf_ctcf_lost_dist <- ctcf_ctcf_lost %>%
  mutate(cpg_any = cpg_ctcf_ctcf_lost$any,
         distance_bin = cut(loop_distance, breaks = DISTANCE_BINS,
                            labels = DISTANCE_LABELS, include.lowest = TRUE))

ctcf_ctcf_gained_dist <- ctcf_ctcf_gained %>%
  mutate(cpg_any = cpg_ctcf_ctcf_gained$any,
         distance_bin = cut(loop_distance, breaks = DISTANCE_BINS,
                            labels = DISTANCE_LABELS, include.lowest = TRUE))

strat_fisher_results <- list()
strat_idx <- 1

for (bin in DISTANCE_LABELS) {
  lost_bin <- ctcf_ctcf_lost_dist %>% filter(distance_bin == bin)
  gained_bin <- ctcf_ctcf_gained_dist %>% filter(distance_bin == bin)

  if (nrow(lost_bin) >= 5 & nrow(gained_bin) >= 5) {
    result <- run_fisher_cpg(
      sum(lost_bin$cpg_any), nrow(lost_bin),
      sum(gained_bin$cpg_any), nrow(gained_bin),
      paste0("lost_", bin), paste0("gained_", bin))
    result$distance_bin <- bin
    strat_fisher_results[[strat_idx]] <- as_tibble(result)
    strat_idx <- strat_idx + 1

    cat(sprintf("    %s: lost %d/%d (%.1f%%) vs gained %d/%d (%.1f%%), OR=%.2f, p=%.2e\n",
                bin,
                sum(lost_bin$cpg_any), nrow(lost_bin),
                100 * sum(lost_bin$cpg_any) / nrow(lost_bin),
                sum(gained_bin$cpg_any), nrow(gained_bin),
                100 * sum(gained_bin$cpg_any) / nrow(gained_bin),
                result$odds_ratio, result$p_value))
  } else {
    cat(sprintf("    %s: skipped (lost n=%d, gained n=%d — too few)\n",
                bin, nrow(lost_bin), nrow(gained_bin)))
  }
}

if (length(strat_fisher_results) > 0) {
  strat_fisher_df <- bind_rows(strat_fisher_results)
  strat_fisher_df$fdr <- p.adjust(strat_fisher_df$p_value, method = "BH")
  strat_fisher_df$significant <- strat_fisher_df$fdr < 0.05
} else {
  strat_fisher_df <- tibble()
}

# ==============================================================================
# SECTION 5: GC CONTENT ANALYSIS
# ==============================================================================

cat("\n[6] Computing GC content at CTCF anchors...\n")

# GC at lost CTCF anchors (ChIP)
lost_gr <- GRanges(seqnames = anchors_chip_lost$chr,
                   ranges = IRanges(start = anchors_chip_lost$start,
                                    end = anchors_chip_lost$end))
anchors_chip_lost$gc_content <- compute_gc_content(lost_gr, genome)

# GC at gained CTCF anchors (ChIP)
gained_gr <- GRanges(seqnames = anchors_chip_gained$chr,
                     ranges = IRanges(start = anchors_chip_gained$start,
                                      end = anchors_chip_gained$end))
anchors_chip_gained$gc_content <- compute_gc_content(gained_gr, genome)

# GC at genome-wide CTCF peaks (baseline)
ctcf_baseline_gc <- compute_gc_content(ctcf_peaks_gr, genome)

cat(sprintf("  GC content (median [IQR]):\n"))
cat(sprintf("    Lost CTCF anchors:   %.3f [%.3f-%.3f] (n=%d)\n",
            median(anchors_chip_lost$gc_content),
            quantile(anchors_chip_lost$gc_content, 0.25),
            quantile(anchors_chip_lost$gc_content, 0.75),
            nrow(anchors_chip_lost)))
cat(sprintf("    Gained CTCF anchors: %.3f [%.3f-%.3f] (n=%d)\n",
            median(anchors_chip_gained$gc_content),
            quantile(anchors_chip_gained$gc_content, 0.25),
            quantile(anchors_chip_gained$gc_content, 0.75),
            nrow(anchors_chip_gained)))
cat(sprintf("    All CTCF peaks:      %.3f [%.3f-%.3f] (n=%d)\n",
            median(ctcf_baseline_gc),
            quantile(ctcf_baseline_gc, 0.25),
            quantile(ctcf_baseline_gc, 0.75),
            length(ctcf_baseline_gc)))

# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(anchors_chip_lost$gc_content,
                             anchors_chip_gained$gc_content,
                             alternative = "two.sided")
n_total <- nrow(anchors_chip_lost) + nrow(anchors_chip_gained)
z_stat <- qnorm(wilcox_result$p.value / 2, lower.tail = FALSE)
effect_size_r <- z_stat / sqrt(n_total)

cat(sprintf("\n  Wilcoxon test (lost vs gained): W=%.0f, p=%.2e, r=%.3f\n",
            wilcox_result$statistic, wilcox_result$p.value, effect_size_r))

# Distance-stratified Wilcoxon tests
cat("\n[7] Distance-stratified GC content comparison...\n")

anchors_chip_lost$distance_bin <- cut(anchors_chip_lost$loop_distance,
                                      breaks = DISTANCE_BINS,
                                      labels = DISTANCE_LABELS,
                                      include.lowest = TRUE)
anchors_chip_gained$distance_bin <- cut(anchors_chip_gained$loop_distance,
                                        breaks = DISTANCE_BINS,
                                        labels = DISTANCE_LABELS,
                                        include.lowest = TRUE)

strat_wilcox_results <- list()
strat_wilcox_idx <- 1

for (bin in DISTANCE_LABELS) {
  lost_gc <- anchors_chip_lost %>% filter(distance_bin == bin)
  gained_gc <- anchors_chip_gained %>% filter(distance_bin == bin)

  if (nrow(lost_gc) >= 5 & nrow(gained_gc) >= 5) {
    wt <- wilcox.test(lost_gc$gc_content, gained_gc$gc_content,
                      alternative = "two.sided")
    n_bin <- nrow(lost_gc) + nrow(gained_gc)
    z_bin <- qnorm(wt$p.value / 2, lower.tail = FALSE)
    r_bin <- z_bin / sqrt(n_bin)

    strat_wilcox_results[[strat_wilcox_idx]] <- tibble(
      distance_bin = bin,
      lost_n = nrow(lost_gc),
      lost_median_gc = median(lost_gc$gc_content),
      gained_n = nrow(gained_gc),
      gained_median_gc = median(gained_gc$gc_content),
      W = wt$statistic,
      p_value = wt$p.value,
      effect_size_r = r_bin
    )
    strat_wilcox_idx <- strat_wilcox_idx + 1

    cat(sprintf("    %s: lost %.3f (n=%d) vs gained %.3f (n=%d), p=%.2e\n",
                bin,
                median(lost_gc$gc_content), nrow(lost_gc),
                median(gained_gc$gc_content), nrow(gained_gc),
                wt$p.value))
  } else {
    cat(sprintf("    %s: skipped (lost n=%d, gained n=%d)\n",
                bin, nrow(lost_gc), nrow(gained_gc)))
  }
}

if (length(strat_wilcox_results) > 0) {
  strat_wilcox_df <- bind_rows(strat_wilcox_results)
  strat_wilcox_df$fdr <- p.adjust(strat_wilcox_df$p_value, method = "BH")
} else {
  strat_wilcox_df <- tibble()
}

# ==============================================================================
# SECTION 6: VISUALIZATIONS
# ==============================================================================

cat("\n[8] Generating plots...\n")

# --- Plot 1: CpG overlap barplot ---
cat("  01_cpg_overlap_barplot\n")

overlap_summary <- tibble(
  category = rep(c("CTCF-CTCF", "Any-CTCF", "Anchor (ChIP)"), each = 2),
  direction = rep(c("Lost", "Gained"), 3),
  cpg_count = c(
    sum(cpg_ctcf_ctcf_lost$any), sum(cpg_ctcf_ctcf_gained$any),
    sum(cpg_any_ctcf_lost$any), sum(cpg_any_ctcf_gained$any),
    sum(anchors_chip_lost$cpg_overlap), sum(anchors_chip_gained$cpg_overlap)
  ),
  total = c(
    nrow(ctcf_ctcf_lost), nrow(ctcf_ctcf_gained),
    nrow(any_ctcf_lost), nrow(any_ctcf_gained),
    nrow(anchors_chip_lost), nrow(anchors_chip_gained)
  )
) %>%
  mutate(pct = 100 * cpg_count / total,
         label = sprintf("%d/%d", cpg_count, total),
         direction = factor(direction, levels = c("Lost", "Gained")),
         category = factor(category, levels = c("CTCF-CTCF", "Any-CTCF", "Anchor (ChIP)")))

p1 <- ggplot(overlap_summary, aes(x = category, y = pct, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                    name = "Loop Direction") +
  labs(x = NULL,
       y = "% with CpG island overlap",
       title = "CpG Island Overlap at CTCF Loop Anchors",
       subtitle = "Lost (down in mutant) vs Gained (up in mutant)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top") +
  coord_cartesian(ylim = c(0, max(overlap_summary$pct) * 1.15))

save_multiformat_ggplot(p1, file.path(PLOT_DIR, "01_cpg_overlap_barplot"),
                        width = 8, height = 6)

# --- Plot 2: Fisher's test forest plot ---
cat("  02_cpg_enrichment_forest\n")

fisher_plot_df <- fisher_df %>%
  mutate(
    comparison = paste(test_label, "vs", ref_label),
    log2_or = log2(odds_ratio),
    log2_ci_lower = log2(ci_lower),
    log2_ci_upper = log2(ci_upper),
    neg_log10_fdr = -log10(fdr),
    comparison = factor(comparison, levels = rev(comparison))
  )

p2 <- ggplot(fisher_plot_df, aes(x = log2_or, y = comparison)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(xmin = log2_ci_lower, xmax = log2_ci_upper),
                width = 0.25, color = "grey30", orientation = "y") +
  geom_point(aes(size = neg_log10_fdr,
                 color = ifelse(significant, enrichment, "ns")),
             shape = 16) +
  scale_color_manual(values = c("enriched" = COLORS$lost,
                                "depleted" = COLORS$gained,
                                "ns" = COLORS$ns),
                     name = "Direction") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
  labs(x = "log2(Odds Ratio)",
       y = NULL,
       title = "CpG Island Enrichment at CTCF Anchors",
       subtitle = "Lost vs Gained loops, Fisher's exact test") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

save_multiformat_ggplot(p2, file.path(PLOT_DIR, "02_cpg_enrichment_forest"),
                        width = 10, height = 5)

# --- Plot 3: GC content box + violin ---
cat("  03_gc_content_boxplot\n")

gc_df <- bind_rows(
  anchors_chip_lost %>% mutate(group = "Lost"),
  anchors_chip_gained %>% mutate(group = "Gained")
) %>%
  mutate(group = factor(group, levels = c("Lost", "Gained")))

baseline_gc_median <- median(ctcf_baseline_gc)

p_label <- ifelse(wilcox_result$p.value < 2.2e-16, "p < 2.2e-16",
                  sprintf("p = %.2e", wilcox_result$p.value))

p3 <- ggplot(gc_df, aes(x = group, y = gc_content, fill = group)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = baseline_gc_median, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = 2.4, y = baseline_gc_median + 0.005,
           label = sprintf("All CTCF peaks: %.3f", baseline_gc_median),
           size = 3, color = "grey40", hjust = 1) +
  scale_fill_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                    guide = "none") +
  labs(x = "CTCF Anchor Direction",
       y = "GC Content",
       title = "GC Content at CTCF Loop Anchors",
       subtitle = sprintf("Wilcoxon rank-sum: %s, r = %.3f", p_label, effect_size_r)) +
  theme_bw(base_size = 12)

save_multiformat_ggplot(p3, file.path(PLOT_DIR, "03_gc_content_boxplot"),
                        width = 6, height = 6)

# --- Plot 4: GC content density ---
cat("  04_gc_content_density\n")

p4 <- ggplot(gc_df, aes(x = gc_content, fill = group, color = group)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = median(anchors_chip_lost$gc_content),
             color = COLORS$lost, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = median(anchors_chip_gained$gc_content),
             color = COLORS$gained, linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                    name = "Direction") +
  scale_color_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                     name = "Direction") +
  labs(x = "GC Content",
       y = "Density",
       title = "GC Content Distribution at CTCF Anchors",
       subtitle = sprintf("Lost median=%.3f, Gained median=%.3f",
                          median(anchors_chip_lost$gc_content),
                          median(anchors_chip_gained$gc_content))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

save_multiformat_ggplot(p4, file.path(PLOT_DIR, "04_gc_content_density"),
                        width = 8, height = 6)

# --- Plot 5: CpG overlap by distance ---
cat("  05_cpg_by_distance\n")

# Build per-distance-bin summary for CTCF-CTCF loops
cpg_by_dist <- bind_rows(
  ctcf_ctcf_lost_dist %>%
    group_by(distance_bin) %>%
    summarise(cpg_pct = 100 * mean(cpg_any), n = n(), .groups = "drop") %>%
    mutate(direction = "Lost"),
  ctcf_ctcf_gained_dist %>%
    group_by(distance_bin) %>%
    summarise(cpg_pct = 100 * mean(cpg_any), n = n(), .groups = "drop") %>%
    mutate(direction = "Gained")
) %>%
  filter(!is.na(distance_bin)) %>%
  mutate(direction = factor(direction, levels = c("Lost", "Gained")))

p5 <- ggplot(cpg_by_dist, aes(x = distance_bin, y = cpg_pct,
                               color = direction, group = direction)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = n), shape = 16) +
  geom_text(aes(label = sprintf("n=%d", n)),
            vjust = -1.2, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                     name = "Direction") +
  scale_size_continuous(name = "Loop count", range = c(2, 5)) +
  labs(x = "Loop Distance",
       y = "% with CpG Island at Anchor",
       title = "CpG Island Overlap by Loop Distance (CTCF-CTCF Loops)",
       subtitle = "Controls for loop length confound") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, max(cpg_by_dist$cpg_pct, na.rm = TRUE) * 1.2))

save_multiformat_ggplot(p5, file.path(PLOT_DIR, "05_cpg_by_distance"),
                        width = 8, height = 6)

# --- Plot 6: GC vs distance scatter ---
cat("  06_gc_vs_distance_scatter\n")

gc_dist_df <- gc_df %>%
  filter(loop_distance > 0) %>%
  mutate(log10_distance = log10(loop_distance))

p6 <- ggplot(gc_dist_df, aes(x = log10_distance, y = gc_content, color = group)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2) +
  scale_color_manual(values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
                     name = "Direction") +
  scale_x_continuous(
    breaks = log10(c(50e3, 100e3, 500e3, 1e6, 5e6)),
    labels = c("50kb", "100kb", "500kb", "1Mb", "5Mb")
  ) +
  labs(x = "Loop Distance",
       y = "GC Content at CTCF Anchor",
       title = "GC Content vs Loop Distance at CTCF Anchors",
       subtitle = "Loess smoothed trends by direction") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

save_multiformat_ggplot(p6, file.path(PLOT_DIR, "06_gc_vs_distance_scatter"),
                        width = 8, height = 6)

# ==============================================================================
# SECTION 7: SAVE TABLES & SUMMARY REPORT
# ==============================================================================

cat("\n[9] Saving output tables...\n")

# CpG overlap summary
write.table(overlap_summary, file.path(TABLE_DIR, "cpg_overlap_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Fisher test results
write.table(fisher_df, file.path(TABLE_DIR, "fisher_test_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Distance-stratified Fisher results
if (nrow(strat_fisher_df) > 0) {
  write.table(strat_fisher_df, file.path(TABLE_DIR, "fisher_test_distance_stratified.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# GC content summary
gc_summary <- tibble(
  group = c("Lost CTCF anchors", "Gained CTCF anchors", "All CTCF peaks (baseline)"),
  n = c(nrow(anchors_chip_lost), nrow(anchors_chip_gained), length(ctcf_baseline_gc)),
  median_gc = c(median(anchors_chip_lost$gc_content),
                median(anchors_chip_gained$gc_content),
                median(ctcf_baseline_gc)),
  q25_gc = c(quantile(anchors_chip_lost$gc_content, 0.25),
             quantile(anchors_chip_gained$gc_content, 0.25),
             quantile(ctcf_baseline_gc, 0.25)),
  q75_gc = c(quantile(anchors_chip_lost$gc_content, 0.75),
             quantile(anchors_chip_gained$gc_content, 0.75),
             quantile(ctcf_baseline_gc, 0.75)),
  mean_gc = c(mean(anchors_chip_lost$gc_content),
              mean(anchors_chip_gained$gc_content),
              mean(ctcf_baseline_gc))
)
write.table(gc_summary, file.path(TABLE_DIR, "gc_content_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Wilcoxon results
wilcox_summary <- tibble(
  comparison = "Lost_vs_Gained_CTCF_anchors",
  W = wilcox_result$statistic,
  p_value = wilcox_result$p.value,
  effect_size_r = effect_size_r,
  lost_median = median(anchors_chip_lost$gc_content),
  gained_median = median(anchors_chip_gained$gc_content),
  lost_n = nrow(anchors_chip_lost),
  gained_n = nrow(anchors_chip_gained)
)
if (nrow(strat_wilcox_df) > 0) {
  wilcox_all <- bind_rows(
    wilcox_summary %>% mutate(distance_bin = "overall"),
    strat_wilcox_df %>% transmute(
      comparison = paste0("Lost_vs_Gained_", distance_bin),
      distance_bin = distance_bin,
      W = W, p_value = p_value, effect_size_r = effect_size_r,
      lost_median = lost_median_gc, gained_median = gained_median_gc,
      lost_n = lost_n, gained_n = gained_n
    )
  )
} else {
  wilcox_all <- wilcox_summary %>% mutate(distance_bin = "overall")
}
write.table(wilcox_all, file.path(TABLE_DIR, "wilcoxon_test_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Full annotated anchor table
anchors_annotated <- bind_rows(
  anchors_chip_lost %>% mutate(group = "lost"),
  anchors_chip_gained %>% mutate(group = "gained")
)
write.table(anchors_annotated, file.path(TABLE_DIR, "ctcf_anchors_annotated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  Saved all tables to ", TABLE_DIR, "\n")

# --- Summary report ---
cat("\n[10] Writing summary report...\n")

report_lines <- c(
  "CpG ISLAND / GC CONTENT ANALYSIS AT CTCF-ANCHORED LOOP ANCHORS",
  "================================================================",
  sprintf("Date: %s", Sys.time()),
  "",
  "MOTIVATION:",
  "-----------",
  "  Jesse Dixon (04/10/26) connected lost CTCF-CTCF loop length asymmetry",
  "  to Brad Bernstein's IDH mutant glioma work (repressive chromatin disrupts",
  "  CTCF binding at CpG-rich regions). This analysis tests whether CpG islands",
  "  and GC content differ at CTCF anchors of lost vs gained loops.",
  "",
  "INPUT DATA:",
  "-----------",
  sprintf("  Loops: %s (%d loops)", INPUT_FILE, nrow(loops)),
  sprintf("  CpG islands: %s (%d islands on standard chroms)", CPG_ISLAND_FILE, length(cpg_islands)),
  sprintf("  CTCF peaks: %s (%d peaks)", CTCF_PEAK_FILE, length(ctcf_peaks_gr)),
  "",
  "LOOP SUBSETS:",
  "-------------",
  sprintf("  CTCF-CTCF (both anchors ChIP+): %d lost, %d gained", nrow(ctcf_ctcf_lost), nrow(ctcf_ctcf_gained)),
  sprintf("  Any-CTCF (>=1 anchor ChIP+):    %d lost, %d gained", nrow(any_ctcf_lost), nrow(any_ctcf_gained)),
  sprintf("  Unique anchors (ChIP):          %d lost, %d gained", nrow(anchors_chip_lost), nrow(anchors_chip_gained)),
  "",
  "CpG ISLAND OVERLAP RATES:",
  "-------------------------",
  sprintf("  CTCF-CTCF lost:   %.1f%% (%d/%d)",
          100 * sum(cpg_ctcf_ctcf_lost$any) / nrow(ctcf_ctcf_lost),
          sum(cpg_ctcf_ctcf_lost$any), nrow(ctcf_ctcf_lost)),
  sprintf("  CTCF-CTCF gained: %.1f%% (%d/%d)",
          100 * sum(cpg_ctcf_ctcf_gained$any) / nrow(ctcf_ctcf_gained),
          sum(cpg_ctcf_ctcf_gained$any), nrow(ctcf_ctcf_gained)),
  sprintf("  Any-CTCF lost:    %.1f%% (%d/%d)",
          100 * sum(cpg_any_ctcf_lost$any) / nrow(any_ctcf_lost),
          sum(cpg_any_ctcf_lost$any), nrow(any_ctcf_lost)),
  sprintf("  Any-CTCF gained:  %.1f%% (%d/%d)",
          100 * sum(cpg_any_ctcf_gained$any) / nrow(any_ctcf_gained),
          sum(cpg_any_ctcf_gained$any), nrow(any_ctcf_gained)),
  ""
)

# Fisher's tests
report_lines <- c(report_lines,
  "FISHER'S EXACT TESTS (CpG island enrichment, lost vs gained):",
  "--------------------------------------------------------------"
)
for (i in 1:nrow(fisher_df)) {
  row <- fisher_df[i, ]
  report_lines <- c(report_lines,
    sprintf("  %s vs %s [%s]: OR=%.2f [%.2f-%.2f], p=%.2e, FDR=%.2e %s",
            row$test_label, row$ref_label, row$ctcf_definition,
            row$odds_ratio, row$ci_lower, row$ci_upper,
            row$p_value, row$fdr,
            ifelse(row$significant, paste0("** ", row$enrichment, " **"), "(ns)")))
}

# Distance-stratified Fisher
if (nrow(strat_fisher_df) > 0) {
  report_lines <- c(report_lines, "",
    "DISTANCE-STRATIFIED FISHER'S TESTS (CTCF-CTCF loops):",
    "------------------------------------------------------"
  )
  for (i in 1:nrow(strat_fisher_df)) {
    row <- strat_fisher_df[i, ]
    report_lines <- c(report_lines,
      sprintf("  %s: OR=%.2f, p=%.2e, FDR=%.2e, lost %.1f%% vs gained %.1f%%",
              row$distance_bin, row$odds_ratio, row$p_value, row$fdr,
              row$test_pct, row$ref_pct))
  }
}

# GC content
report_lines <- c(report_lines, "",
  "GC CONTENT AT CTCF ANCHORS:",
  "---------------------------",
  sprintf("  Lost CTCF anchors:   median=%.3f [IQR %.3f-%.3f] (n=%d)",
          median(anchors_chip_lost$gc_content),
          quantile(anchors_chip_lost$gc_content, 0.25),
          quantile(anchors_chip_lost$gc_content, 0.75),
          nrow(anchors_chip_lost)),
  sprintf("  Gained CTCF anchors: median=%.3f [IQR %.3f-%.3f] (n=%d)",
          median(anchors_chip_gained$gc_content),
          quantile(anchors_chip_gained$gc_content, 0.25),
          quantile(anchors_chip_gained$gc_content, 0.75),
          nrow(anchors_chip_gained)),
  sprintf("  All CTCF peaks:      median=%.3f [IQR %.3f-%.3f] (n=%d, baseline)",
          median(ctcf_baseline_gc),
          quantile(ctcf_baseline_gc, 0.25),
          quantile(ctcf_baseline_gc, 0.75),
          length(ctcf_baseline_gc)),
  "",
  "WILCOXON RANK-SUM TEST (GC content, lost vs gained):",
  "----------------------------------------------------",
  sprintf("  Overall: W=%.0f, p=%.2e, effect size r=%.3f",
          wilcox_result$statistic, wilcox_result$p.value, effect_size_r)
)

# Distance-stratified Wilcoxon
if (nrow(strat_wilcox_df) > 0) {
  report_lines <- c(report_lines, "",
    "  Distance-stratified:")
  for (i in 1:nrow(strat_wilcox_df)) {
    row <- strat_wilcox_df[i, ]
    report_lines <- c(report_lines,
      sprintf("    %s: lost %.3f vs gained %.3f, p=%.2e, FDR=%.2e",
              row$distance_bin, row$lost_median_gc, row$gained_median_gc,
              row$p_value, row$fdr))
  }
}

report_lines <- c(report_lines, "",
  "INTERPRETATION:",
  "---------------",
  "  - OR > 1 means CpG islands are more common at LOST loop anchors",
  "  - OR < 1 means CpG islands are more common at GAINED loop anchors",
  "  - If distance-stratified tests lose significance, the CpG/GC effect",
  "    is confounded with loop distance (lost loops are ~2.5x longer)",
  "  - If significant within distance bins, the effect is independent of length",
  "  - Bernstein IDH prediction: lost anchors enriched for CpG islands (OR > 1)",
  "",
  "CAVEATS:",
  "--------",
  "  - CTCF ChIP-seq peaks are from reference (not differential binding)",
  "  - CpG islands defined by UCSC/GENCODE annotation (not methylation state)",
  "  - Anchor width is 5-25kb depending on resolution; CpG islands are ~1kb",
  "  - ChIP-seq based CTCF definition used (not motifs); 81% of motifs are unbound",
  "",
  "OUTPUT FILES:",
  "-------------",
  "  tables/cpg_overlap_summary.tsv",
  "  tables/fisher_test_results.tsv",
  "  tables/fisher_test_distance_stratified.tsv",
  "  tables/gc_content_summary.tsv",
  "  tables/wilcoxon_test_results.tsv",
  "  tables/ctcf_anchors_annotated.tsv",
  "  plots/01_cpg_overlap_barplot/{pdf,svg,jpg}",
  "  plots/02_cpg_enrichment_forest/{pdf,svg,jpg}",
  "  plots/03_gc_content_boxplot/{pdf,svg,jpg}",
  "  plots/04_gc_content_density/{pdf,svg,jpg}",
  "  plots/05_cpg_by_distance/{pdf,svg,jpg}",
  "  plots/06_gc_vs_distance_scatter/{pdf,svg,jpg}"
)

writeLines(report_lines, file.path(OUTPUT_DIR, "cpg_ctcf_analysis_summary.txt"))
cat("  Saved summary report\n")

cat("\n", "=", rep("=", 50) |> paste(collapse = ""), "\n")
cat("  Analysis complete. Results in:", OUTPUT_DIR, "\n")
cat("=", rep("=", 50) |> paste(collapse = ""), "\n")
