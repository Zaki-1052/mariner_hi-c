# abc/scripts/step14_loop_type_at_enhancer_subsets.R
# Loop type decomposition at enhancer phenotypic subsets.
#
# Uses ALL loops (not just differential) from the merged results, annotates
# them with 7-category chromatin states using ChIP-seq peaks, then breaks
# down loop types present at each enhancer class (Activity_Lost, K119ub_Only,
# Activity_Gain, Stable). Focuses on the "partner anchor" — what the enhancer
# connects to.
#
# Key question: if enhancer-TSS (promoter) changes are minimal at K119ub_Only
# sites, what other loop types are present and contributing to the logFC?
#
# Inputs (all relative to abc/ working directory):
#   enhancer_subsets/enhancer_classes_*.tsv    -- 4 enhancer class files
#   ../outputs/250402-late_outputs/merged_loops/merged_all_results.tsv
#                                              -- ALL loops (39K, not just
#                                                 differential)
#   ../peaks/beds/*.bed                        -- ChIP-seq peaks for annotation
#   reference/mm10_tss.bed                     -- TSS positions with gene names
#   results/delta_abc_all_pairs.tsv            -- ABC enhancer-gene pairs
#
# Outputs:
#   all_characterized_loops.tsv                -- cached annotated all-loops file
#   results/loop_type_at_subsets/              -- TSVs + plots (PDF + SVG + JPG)
#
# Usage:
#   cd abc && Rscript scripts/step14_loop_type_at_enhancer_subsets.R
#   cd abc && Rscript scripts/step14_loop_type_at_enhancer_subsets.R --force-annotate

# =============================================================================
# CONFIGURATION
# =============================================================================

ENHANCER_FILES <- c(
  Activity_Lost = "enhancer_subsets/enhancer_classes_activity_lost.tsv",
  K119ub_Only   = "enhancer_subsets/enhancer_classes_k119ub_only.tsv",
  Activity_Gain = "enhancer_subsets/enhancer_classes_activity_gain.tsv",
  Stable        = "enhancer_subsets/enhancer_classes_stable.tsv"
)

ALL_LOOPS_FILE  <- "../outputs/250402-late_outputs/merged_loops/merged_all_results.tsv"
ANNOTATED_CACHE <- "all_characterized_loops.tsv"

PEAK_FILES <- list(
  h3k27ac  = "../peaks/beds/H3K27acCerebellumLate2.bed",
  h3k27me3 = "../peaks/beds/H3K27me3CerebellumLate1.bed",
  h3k4me1  = "../peaks/beds/H3K4me1CerebellumLate1.bed",
  h3k4me3  = "../peaks/beds/H3K4me3CerebellumLate2.bed",
  bivalent = "../peaks/beds/Bivalent_Cerebellum_Late.bed"
)

TSS_FILE       <- "reference/mm10_tss.bed"
ABC_PAIRS_FILE <- "results/delta_abc_all_pairs.tsv"
TSS_THRESHOLD  <- 2000

OUTPUT_DIR <- "results/loop_type_at_subsets"

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

PARTNER_COLORS <- c(
  Active_Promoter    = "#E41A1C",
  Repressed_Promoter = "#377EB8",
  Bivalent_Promoter  = "#FF7F00",
  Polycomb           = "#984EA3",
  Active_Enhancer    = "#4DAF4A",
  Poised_Enhancer    = "#A6D854",
  Other              = "#999999"
)
PARTNER_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                   "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")

# Used by classify_loop_type_extended for canonical ordering
ANCHOR_TYPE_ORDER <- PARTNER_ORDER

FUNCTIONAL_GROUPS <- c(
  Active_Promoter    = "Promoter",
  Repressed_Promoter = "Promoter",
  Bivalent_Promoter  = "Promoter",
  Polycomb           = "Polycomb",
  Active_Enhancer    = "Enhancer",
  Poised_Enhancer    = "Enhancer",
  Other              = "Structural"
)
FUNC_COLORS <- c(
  Promoter   = "#E41A1C",
  Polycomb   = "#984EA3",
  Enhancer   = "#4DAF4A",
  Structural = "#999999"
)
FUNC_ORDER <- c("Promoter", "Polycomb", "Enhancer", "Structural")

# Promoter types for enhancer-TSS contact analysis
PROMOTER_TYPES <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter")

cat("================================================================================\n")
cat("STEP 14: LOOP TYPE DECOMPOSITION AT ENHANCER SUBSETS (ALL LOOPS)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(scales)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
for (nm in names(ENHANCER_FILES)) {
  stopifnot(file.exists(ENHANCER_FILES[[nm]]))
  cat(sprintf("  [OK] %s: %s\n", nm, ENHANCER_FILES[[nm]]))
}
stopifnot(file.exists(ALL_LOOPS_FILE))
cat(sprintf("  [OK] All loops: %s\n", ALL_LOOPS_FILE))

for (nm in names(PEAK_FILES)) {
  stopifnot(file.exists(PEAK_FILES[[nm]]))
  cat(sprintf("  [OK] Peak %s: %s\n", nm, PEAK_FILES[[nm]]))
}

stopifnot(file.exists(TSS_FILE))
cat(sprintf("  [OK] TSS reference: %s\n", TSS_FILE))

stopifnot(file.exists(ABC_PAIRS_FILE))
cat(sprintf("  [OK] ABC pairs: %s\n", ABC_PAIRS_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# ANNOTATION FUNCTIONS (from scripts/annotate_loops_extended.R)
# =============================================================================

load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    stop(sprintf("%s bed file not found: %s", mark_name, bed_path))
  }
  cat(sprintf("  Loading %s peaks from: %s\n", mark_name, bed_path))

  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  gr
}

annotate_chip_overlaps_extended <- function(anchor_gr, k27ac_gr, k27me3_gr,
                                            k4me1_gr, k4me3_gr, bivalent_gr) {
  data.frame(
    H3K27ac_overlap = countOverlaps(anchor_gr, k27ac_gr) > 0,
    H3K27me3_overlap = countOverlaps(anchor_gr, k27me3_gr) > 0,
    H3K4me1_overlap = countOverlaps(anchor_gr, k4me1_gr) > 0,
    H3K4me3_overlap = countOverlaps(anchor_gr, k4me3_gr) > 0,
    Bivalent_Promoter_overlap = countOverlaps(anchor_gr, bivalent_gr) > 0
  )
}

classify_anchor_type_extended <- function(h3k27ac_overlap, h3k27me3_overlap,
                                          h3k4me1_overlap, h3k4me3_overlap,
                                          bivalent_overlap, distance_to_tss,
                                          tss_threshold = 2000) {
  n <- length(h3k27ac_overlap)
  anchor_type <- rep("Other", n)

  # 1. Active_Promoter: H3K4me3+ AND NOT H3K27me3 AND <=2kb from TSS
  is_active_promoter <- h3k4me3_overlap & !h3k27me3_overlap &
                        !is.na(distance_to_tss) &
                        distance_to_tss <= tss_threshold
  anchor_type[is_active_promoter] <- "Active_Promoter"

  # 2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND <=2kb from TSS
  is_repressed_promoter <- !is_active_promoter &
                           h3k27me3_overlap & !h3k27ac_overlap &
                           !is.na(distance_to_tss) &
                           distance_to_tss <= tss_threshold
  anchor_type[is_repressed_promoter] <- "Repressed_Promoter"

  # 3. Bivalent_Promoter: K4me3+K27me3 overlap (not already classified)
  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent_overlap
  anchor_type[is_bivalent] <- "Bivalent_Promoter"

  # 4. Polycomb: H3K27me3+ AND >2kb from TSS (distal repressive)
  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
                 h3k27me3_overlap &
                 (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_polycomb] <- "Polycomb"

  # 5. Active_Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb &
                        h3k27ac_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  # 6. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb & !is_active_enhancer &
                        h3k4me1_overlap & !h3k27ac_overlap & !h3k27me3_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

  # 7. Other: default (no ChIP-seq marks / structural elements)
  return(anchor_type)
}

classify_loop_type_extended <- function(anchor1_type, anchor2_type) {
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

# =============================================================================
# LOAD & ANNOTATE ALL LOOPS
# =============================================================================

cat("Loading and annotating all loops...\n")

force_annotate <- "--force-annotate" %in% commandArgs(trailingOnly = TRUE)

expected_annotation_cols <- c("anchor1_type", "anchor2_type", "loop_type",
                              "anchor1_distance_to_tss", "anchor2_distance_to_tss",
                              "anchor1_nearest_gene", "anchor2_nearest_gene")

if (!force_annotate && file.exists(ANNOTATED_CACHE)) {
  cat(sprintf("  Loading cached annotated loops from %s...\n", ANNOTATED_CACHE))
  all_loops <- read.table(ANNOTATED_CACHE, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
  missing <- setdiff(expected_annotation_cols, colnames(all_loops))
  if (length(missing) > 0) {
    stop(sprintf("Cache file %s is missing annotation columns: %s. Re-run with --force-annotate.",
                 ANNOTATED_CACHE, paste(missing, collapse = ", ")))
  }
  cat(sprintf("  Loaded %d annotated loops from cache.\n", nrow(all_loops)))
} else {
  if (force_annotate) cat("  --force-annotate: re-annotating de novo.\n")
  cat(sprintf("  Annotating all loops from %s...\n", ALL_LOOPS_FILE))

  # 1. Load raw loops
  all_loops <- read.table(ALL_LOOPS_FILE, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
  stopifnot(nrow(all_loops) > 0)
  cat(sprintf("  Loaded %d loops\n", nrow(all_loops)))

  # 2. Load peak files
  k27ac_gr   <- load_chip_peaks(PEAK_FILES$h3k27ac, "H3K27ac")
  k27me3_gr  <- load_chip_peaks(PEAK_FILES$h3k27me3, "H3K27me3")
  k4me1_gr   <- load_chip_peaks(PEAK_FILES$h3k4me1, "H3K4me1")
  k4me3_gr   <- load_chip_peaks(PEAK_FILES$h3k4me3, "H3K4me3")
  bivalent_gr <- load_chip_peaks(PEAK_FILES$bivalent, "Bivalent")

  # 3. Load TSS reference for distance computation AND gene identification
  tss_df <- read.table(TSS_FILE, sep = "\t", header = TRUE,
                       comment.char = "", stringsAsFactors = FALSE)
  # Header is #chr — R reads as X.chr
  chr_col <- grep("chr", colnames(tss_df), value = TRUE)[1]
  tss_gr <- GRanges(
    seqnames = tss_df[[chr_col]],
    ranges = IRanges(start = tss_df$start + 250, width = 1),
    gene_name = tss_df$name,
    gene_type = tss_df$gene_type,
    ensembl_id = tss_df$Ensembl_ID
  )
  cat(sprintf("  Loaded %d TSS positions from %s\n", length(tss_gr), TSS_FILE))

  # 4. Build anchor GRanges
  anchor1_gr <- GRanges(seqnames = all_loops$chr1,
                        ranges = IRanges(start = all_loops$start1,
                                         end = all_loops$end1))
  anchor2_gr <- GRanges(seqnames = all_loops$chr2,
                        ranges = IRanges(start = all_loops$start2,
                                         end = all_loops$end2))

  # 5. Compute TSS distances and nearest gene
  cat("  Computing TSS distances...\n")
  nearest1 <- distanceToNearest(anchor1_gr, tss_gr)
  a1_dist <- rep(NA_real_, length(anchor1_gr))
  a1_gene <- rep(NA_character_, length(anchor1_gr))
  a1_dist[queryHits(nearest1)] <- mcols(nearest1)$distance
  a1_gene[queryHits(nearest1)] <- tss_gr$gene_name[subjectHits(nearest1)]

  nearest2 <- distanceToNearest(anchor2_gr, tss_gr)
  a2_dist <- rep(NA_real_, length(anchor2_gr))
  a2_gene <- rep(NA_character_, length(anchor2_gr))
  a2_dist[queryHits(nearest2)] <- mcols(nearest2)$distance
  a2_gene[queryHits(nearest2)] <- tss_gr$gene_name[subjectHits(nearest2)]

  # 6. Annotate ChIP-seq overlaps
  cat("  Annotating ChIP-seq overlaps...\n")
  a1_chip <- annotate_chip_overlaps_extended(anchor1_gr, k27ac_gr, k27me3_gr,
                                              k4me1_gr, k4me3_gr, bivalent_gr)
  a2_chip <- annotate_chip_overlaps_extended(anchor2_gr, k27ac_gr, k27me3_gr,
                                              k4me1_gr, k4me3_gr, bivalent_gr)

  # 7. Classify anchor types (7-category system)
  cat("  Classifying anchor types...\n")
  anchor1_type <- classify_anchor_type_extended(
    a1_chip$H3K27ac_overlap, a1_chip$H3K27me3_overlap,
    a1_chip$H3K4me1_overlap, a1_chip$H3K4me3_overlap,
    a1_chip$Bivalent_Promoter_overlap, a1_dist, TSS_THRESHOLD)

  anchor2_type <- classify_anchor_type_extended(
    a2_chip$H3K27ac_overlap, a2_chip$H3K27me3_overlap,
    a2_chip$H3K4me1_overlap, a2_chip$H3K4me3_overlap,
    a2_chip$Bivalent_Promoter_overlap, a2_dist, TSS_THRESHOLD)

  # 8. Classify loop types
  loop_type <- classify_loop_type_extended(anchor1_type, anchor2_type)

  # 9. Add annotation columns
  all_loops$anchor1_type <- anchor1_type
  all_loops$anchor2_type <- anchor2_type
  all_loops$loop_type    <- loop_type
  all_loops$anchor1_distance_to_tss <- a1_dist
  all_loops$anchor2_distance_to_tss <- a2_dist
  all_loops$anchor1_nearest_gene <- a1_gene
  all_loops$anchor2_nearest_gene <- a2_gene

  # 10. Save cache
  write.table(all_loops, ANNOTATED_CACHE, sep = "\t", quote = FALSE,
              row.names = FALSE)
  cat(sprintf("  Saved annotated cache: %s (%d loops, %d columns)\n",
              ANNOTATED_CACHE, nrow(all_loops), ncol(all_loops)))
}

# Add convenience column
all_loops$is_differential <- all_loops$direction != "unchanged"

# Print annotation summary
cat(sprintf("\n  Total loops: %d\n", nrow(all_loops)))
cat(sprintf("    Differential (up/down): %d\n", sum(all_loops$is_differential)))
cat(sprintf("    Unchanged:             %d\n", sum(!all_loops$is_differential)))

cat("\n  Anchor1 type distribution:\n")
for (tp in PARTNER_ORDER) {
  n <- sum(all_loops$anchor1_type == tp)
  cat(sprintf("    %-22s %5d (%.1f%%)\n", tp, n, 100 * n / nrow(all_loops)))
}

cat("\n  Anchor2 type distribution:\n")
for (tp in PARTNER_ORDER) {
  n <- sum(all_loops$anchor2_type == tp)
  cat(sprintf("    %-22s %5d (%.1f%%)\n", tp, n, 100 * n / nrow(all_loops)))
}

# Validate annotation
valid_types <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                 "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")
actual_types <- unique(c(all_loops$anchor1_type, all_loops$anchor2_type))
unexpected <- setdiff(actual_types, valid_types)
if (length(unexpected) > 0) {
  stop(sprintf("Unexpected anchor types: %s", paste(unexpected, collapse = ", ")))
}

cat(sprintf("\n  Unique loop types: %d\n", length(unique(all_loops$loop_type))))
cat("\n")


# =============================================================================
# LOAD ENHANCER DATA
# =============================================================================

cat("Loading enhancer data...\n")

enh_list <- lapply(names(ENHANCER_FILES), function(cls) {
  df <- read.table(ENHANCER_FILES[[cls]], sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
  df$enhancer_class <- cls
  df
})
enh_all <- do.call(rbind, enh_list)
enh_all$enhancer_class <- factor(enh_all$enhancer_class, levels = CLASS_ORDER)
cat(sprintf("  Enhancers: %d total across %d classes\n", nrow(enh_all),
            length(CLASS_ORDER)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_all$enhancer_class == cls)))
}
cat("\n")


# #############################################################################
# PART A: ENHANCER-LOOP OVERLAP WITH PARTNER ANCHOR CLASSIFICATION
# #############################################################################

cat("=== PART A: Enhancer-Loop Overlap & Partner Classification ===\n\n")

# Build GRanges
enh_gr <- GRanges(
  seqnames = enh_all$chr,
  ranges = IRanges(start = enh_all$start, end = enh_all$end),
  enhancer_class = enh_all$enhancer_class,
  enh_idx = seq_len(nrow(enh_all))
)

loop_anchor1_gr <- GRanges(
  seqnames = all_loops$chr1,
  ranges = IRanges(start = all_loops$start1, end = all_loops$end1),
  loop_idx = seq_len(nrow(all_loops))
)

loop_anchor2_gr <- GRanges(
  seqnames = all_loops$chr2,
  ranges = IRanges(start = all_loops$start2, end = all_loops$end2),
  loop_idx = seq_len(nrow(all_loops))
)

cat("GRanges built.\n")

# Find overlaps: enhancer <-> loop anchors
hits_a1 <- findOverlaps(enh_gr, loop_anchor1_gr, type = "any")
hits_a2 <- findOverlaps(enh_gr, loop_anchor2_gr, type = "any")

cat(sprintf("  Enhancer <-> anchor1 overlaps: %d\n", length(hits_a1)))
cat(sprintf("  Enhancer <-> anchor2 overlaps: %d\n", length(hits_a2)))

# Build enhancer-loop pair table
# When enhancer overlaps anchor1 -> partner is anchor2 (and vice versa)
pairs_a1 <- data.frame(
  enh_idx       = queryHits(hits_a1),
  loop_idx      = subjectHits(hits_a1),
  enh_on_anchor = 1L,
  partner_type  = all_loops$anchor2_type[subjectHits(hits_a1)],
  own_type      = all_loops$anchor1_type[subjectHits(hits_a1)],
  stringsAsFactors = FALSE
)

pairs_a2 <- data.frame(
  enh_idx       = queryHits(hits_a2),
  loop_idx      = subjectHits(hits_a2),
  enh_on_anchor = 2L,
  partner_type  = all_loops$anchor1_type[subjectHits(hits_a2)],
  own_type      = all_loops$anchor2_type[subjectHits(hits_a2)],
  stringsAsFactors = FALSE
)

enh_loop_pairs <- rbind(pairs_a1, pairs_a2)

# Deduplicate: if enhancer overlaps BOTH anchors of same loop, keep first
enh_loop_pairs <- enh_loop_pairs[
  !duplicated(enh_loop_pairs[, c("enh_idx", "loop_idx")]), ]

cat(sprintf("  Unique enhancer-loop pairs (after dedup): %d\n",
            nrow(enh_loop_pairs)))

# Attach metadata
enh_loop_pairs$enhancer_class <- enh_all$enhancer_class[enh_loop_pairs$enh_idx]
enh_loop_pairs$loop_type      <- all_loops$loop_type[enh_loop_pairs$loop_idx]
enh_loop_pairs$logFC           <- all_loops$logFC[enh_loop_pairs$loop_idx]
enh_loop_pairs$FDR             <- all_loops$FDR[enh_loop_pairs$loop_idx]
enh_loop_pairs$direction       <- all_loops$direction[enh_loop_pairs$loop_idx]
enh_loop_pairs$is_differential <- all_loops$is_differential[enh_loop_pairs$loop_idx]

# Factor levels
enh_loop_pairs$enhancer_class <- factor(enh_loop_pairs$enhancer_class,
                                         levels = CLASS_ORDER)
enh_loop_pairs$partner_type <- factor(enh_loop_pairs$partner_type,
                                       levels = PARTNER_ORDER)

# Functional grouping
enh_loop_pairs$partner_group <- factor(
  FUNCTIONAL_GROUPS[as.character(enh_loop_pairs$partner_type)],
  levels = FUNC_ORDER
)

cat(sprintf("  Enhancers with >=1 loop: %d / %d (%.1f%%)\n",
            length(unique(enh_loop_pairs$enh_idx)), nrow(enh_all),
            100 * length(unique(enh_loop_pairs$enh_idx)) / nrow(enh_all)))

# Per-class summary
cat("\n  Per-class enhancer-loop pair counts:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls, ]
  n_enh <- length(unique(sub$enh_idx))
  n_total_enh <- sum(enh_all$enhancer_class == cls)
  n_diff <- sum(sub$is_differential)
  cat(sprintf("    %s: %d pairs (%d differential), %d unique enhancers (%.1f%% of %d)\n",
              cls, nrow(sub), n_diff, n_enh, 100 * n_enh / n_total_enh, n_total_enh))
}

# Sanity check: overall logFC by class
cat("\n  Sanity check - overall median logFC per class:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls, ]
  cat(sprintf("    %s: median logFC = %.4f (n=%d)\n",
              cls, median(sub$logFC), nrow(sub)))
}

# Save detailed pairs table
detailed_pairs <- enh_loop_pairs
detailed_pairs$enh_chr   <- enh_all$chr[detailed_pairs$enh_idx]
detailed_pairs$enh_start <- enh_all$start[detailed_pairs$enh_idx]
detailed_pairs$enh_end   <- enh_all$end[detailed_pairs$enh_idx]
detailed_pairs$loop_chr1   <- all_loops$chr1[detailed_pairs$loop_idx]
detailed_pairs$loop_start1 <- all_loops$start1[detailed_pairs$loop_idx]
detailed_pairs$loop_end1   <- all_loops$end1[detailed_pairs$loop_idx]
detailed_pairs$loop_chr2   <- all_loops$chr2[detailed_pairs$loop_idx]
detailed_pairs$loop_start2 <- all_loops$start2[detailed_pairs$loop_idx]
detailed_pairs$loop_end2   <- all_loops$end2[detailed_pairs$loop_idx]

write.table(detailed_pairs,
            file.path(OUTPUT_DIR, "detailed_enhancer_loop_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n",
            file.path(OUTPUT_DIR, "detailed_enhancer_loop_pairs.tsv")))


# #############################################################################
# PART B: LOOP TYPE COMPOSITION BY ENHANCER CLASS
# #############################################################################

cat("\n=== PART B: Loop Type Composition ===\n\n")

# --- Composition table (7-cat partner types) ---
comp_table <- as.data.frame(table(
  enhancer_class = enh_loop_pairs$enhancer_class,
  partner_type = enh_loop_pairs$partner_type
))
comp_table$proportion <- ave(comp_table$Freq, comp_table$enhancer_class,
                              FUN = function(x) x / sum(x))

# Print summary
cat("  Partner type composition by enhancer class:\n")
for (cls in CLASS_ORDER) {
  sub <- comp_table[comp_table$enhancer_class == cls, ]
  sub <- sub[order(-sub$proportion), ]
  n_total <- sum(sub$Freq)
  cat(sprintf("\n  %s (n=%d loops):\n", cls, n_total))
  for (i in seq_len(nrow(sub))) {
    if (sub$Freq[i] > 0) {
      cat(sprintf("    %-22s %4d (%5.1f%%)\n",
                  as.character(sub$partner_type[i]),
                  sub$Freq[i], 100 * sub$proportion[i]))
    }
  }
}

write.table(comp_table, file.path(OUTPUT_DIR, "partner_type_composition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Chi-squared test ---
contingency <- xtabs(Freq ~ enhancer_class + partner_type, data = comp_table)
chi_test <- chisq.test(contingency, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n  Chi-squared (partner type ~ enhancer class): X2=%.1f, %s\n",
            chi_test$statistic, fmt_p(chi_test$p.value)))

# --- Plot 01: Stacked bar — partner type by enhancer class (proportions) ---
cat("\nGenerating Part B plots...\n")

p01 <- ggplot(comp_table[comp_table$Freq > 0, ],
              aes(x = enhancer_class, y = Freq, fill = partner_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = PARTNER_COLORS) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Enhancer class", y = "Proportion of loops",
    fill = "Partner anchor type",
    title = "Loop partner composition by enhancer class (all loops)",
    subtitle = sprintf("7-category classification | %d total loops | Chi-sq %s",
                       nrow(all_loops), fmt_p(chi_test$p.value))
  ) +
  theme_pub +
  guides(fill = guide_legend(nrow = 2))
save_plot(p01, "01_partner_type_stacked_bar", w = 8, h = 6)

# --- Plot 02: Simplified functional grouping stacked bar ---
func_comp <- as.data.frame(table(
  enhancer_class = enh_loop_pairs$enhancer_class,
  partner_group = enh_loop_pairs$partner_group
))
func_comp$proportion <- ave(func_comp$Freq, func_comp$enhancer_class,
                             FUN = function(x) x / sum(x))

n_labels <- aggregate(Freq ~ enhancer_class, data = func_comp, FUN = sum)
colnames(n_labels) <- c("enhancer_class", "n_total")

p02 <- ggplot(func_comp[func_comp$Freq > 0, ],
              aes(x = enhancer_class, y = Freq, fill = partner_group)) +
  geom_col(position = "fill") +
  geom_text(data = n_labels, aes(x = enhancer_class, y = 1.03,
            label = paste0("n=", format(n_total, big.mark = ","))),
            inherit.aes = FALSE, size = 3.2) +
  scale_fill_manual(values = FUNC_COLORS) +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Enhancer class", y = "Proportion of loops",
    fill = "Partner category",
    title = "Functional loop partner categories by enhancer class (all loops)",
    subtitle = paste("Promoter (Active/Repressed/Bivalent) | Polycomb |",
                     "Enhancer (Active/Poised) | Structural (Other)")
  ) +
  theme_pub
save_plot(p02, "02_functional_group_stacked_bar", w = 8, h = 6)

# --- Plot 03: Heatmap — proportions (partner_type x enhancer_class) ---
heat_data <- comp_table[, c("enhancer_class", "partner_type",
                             "proportion", "Freq")]
heat_data$label <- ifelse(
  heat_data$Freq > 0,
  sprintf("%d\n(%.0f%%)", heat_data$Freq, 100 * heat_data$proportion),
  ""
)

p03 <- ggplot(heat_data, aes(x = enhancer_class, y = partner_type,
                              fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, NA),
                      labels = percent_format()) +
  scale_y_discrete(limits = rev(PARTNER_ORDER)) +
  labs(
    x = "Enhancer class", y = "Partner anchor type",
    fill = "Proportion",
    title = "Partner type proportions by enhancer class (all loops)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p03, "03_partner_type_heatmap", w = 8, h = 7)


# #############################################################################
# PART C: logFC BY LOOP TYPE x ENHANCER CLASS
# #############################################################################

cat("\n=== PART C: logFC by Loop Type x Enhancer Class ===\n\n")

# --- Build summary statistics table ---
logfc_summary <- do.call(rbind, lapply(CLASS_ORDER, function(cls) {
  do.call(rbind, lapply(PARTNER_ORDER, function(pt) {
    sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls &
                            enh_loop_pairs$partner_type == pt, ]
    n <- nrow(sub)
    if (n >= 3) {
      wt <- tryCatch(wilcox.test(sub$logFC, mu = 0),
                     error = function(e) list(p.value = NA))
      data.frame(
        enhancer_class = cls,
        partner_type = pt,
        n = n,
        median_logFC = median(sub$logFC),
        mean_logFC = mean(sub$logFC),
        sd_logFC = sd(sub$logFC),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        enhancer_class = cls,
        partner_type = pt,
        n = n,
        median_logFC = ifelse(n > 0, median(sub$logFC), NA_real_),
        mean_logFC = ifelse(n > 0, mean(sub$logFC), NA_real_),
        sd_logFC = NA_real_,
        wilcox_p = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }))
}))

# Print formatted table
cat("  logFC by partner type x enhancer class (n >= 5):\n\n")
cat(sprintf("  %-15s %-22s %5s %10s %12s\n",
            "Class", "Partner Type", "n", "med logFC", "Wilcoxon p"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in seq_len(nrow(logfc_summary))) {
  row <- logfc_summary[i, ]
  if (row$n >= 5) {
    cat(sprintf("  %-15s %-22s %5d %10.4f %12s\n",
                row$enhancer_class, row$partner_type, row$n,
                row$median_logFC, fmt_p(row$wilcox_p)))
  }
}

write.table(logfc_summary,
            file.path(OUTPUT_DIR, "loop_logfc_by_type_and_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plot 04: Faceted boxplot — logFC by partner type, faceted by class ---
cat("\nGenerating Part C plots...\n")

# Filter to partner types with sufficient data
pt_counts <- table(enh_loop_pairs$partner_type)
viable_pts <- names(pt_counts[pt_counts >= 5])
plot_data <- enh_loop_pairs[enh_loop_pairs$partner_type %in% viable_pts, ]
plot_data$partner_type <- factor(plot_data$partner_type, levels = PARTNER_ORDER)

p04 <- ggplot(plot_data, aes(x = partner_type, y = logFC, fill = partner_type)) +
  geom_boxplot(outlier.size = 0.3, notch = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ enhancer_class, scales = "free_x") +
  scale_fill_manual(values = PARTNER_COLORS) +
  labs(
    x = "Partner anchor type", y = "Loop logFC (KO/WT)",
    title = "Loop strength change by partner type and enhancer class (all loops)",
    subtitle = "Dashed line = no change"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
save_plot(p04, "04_logfc_by_partner_faceted", w = 12, h = 8)

# --- Plot 05: Heatmap — median logFC (partner_type x enhancer_class) ---
heat_logfc <- logfc_summary[logfc_summary$n >= 5, ]
heat_logfc$enhancer_class <- factor(heat_logfc$enhancer_class, levels = CLASS_ORDER)
heat_logfc$partner_type <- factor(heat_logfc$partner_type, levels = PARTNER_ORDER)

# Significance markers
heat_logfc$sig <- ""
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 0.05] <- "*"
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 0.001] <- "**"
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 1e-10] <- "***"

heat_logfc$label <- sprintf("%.3f\n(n=%d)%s",
                             heat_logfc$median_logFC,
                             heat_logfc$n,
                             heat_logfc$sig)

p05 <- ggplot(heat_logfc, aes(x = enhancer_class, y = partner_type,
                               fill = median_logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.2, 0.2), oob = squish,
                       name = "Median\nlogFC") +
  scale_y_discrete(limits = rev(PARTNER_ORDER)) +
  labs(
    x = "Enhancer class", y = "Partner anchor type",
    title = "Median loop logFC by partner type x enhancer class (all loops)",
    subtitle = "Wilcoxon vs 0: * p<0.05, ** p<0.001, *** p<1e-10"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p05, "05_logfc_heatmap", w = 8, h = 7)

# --- Plot 06: Contribution analysis ---
# For each class: contribution_i = (n_i / N_total) * median_logFC_i
contrib <- logfc_summary[logfc_summary$n >= 3, ]
contrib <- do.call(rbind, lapply(CLASS_ORDER, function(cls) {
  sub <- contrib[contrib$enhancer_class == cls, ]
  n_total <- sum(sub$n)
  sub$weight <- sub$n / n_total
  sub$contribution <- sub$weight * sub$median_logFC
  sub
}))

cat("\n  Contribution analysis (weighted logFC by partner type):\n")
for (cls in CLASS_ORDER) {
  sub <- contrib[contrib$enhancer_class == cls, ]
  sub <- sub[order(-abs(sub$contribution)), ]
  total_contrib <- sum(sub$contribution, na.rm = TRUE)
  cat(sprintf("\n  %s (total weighted logFC = %.4f):\n", cls, total_contrib))
  for (i in seq_len(min(nrow(sub), 7))) {
    row <- sub[i, ]
    if (!is.na(row$contribution)) {
      cat(sprintf("    %-22s weight=%.3f x median=%.4f = %.5f\n",
                  row$partner_type, row$weight, row$median_logFC,
                  row$contribution))
    }
  }
}

write.table(contrib, file.path(OUTPUT_DIR, "contribution_analysis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

contrib_plot <- contrib[!is.na(contrib$contribution), ]
contrib_plot$enhancer_class <- factor(contrib_plot$enhancer_class,
                                       levels = CLASS_ORDER)
contrib_plot$partner_type <- factor(contrib_plot$partner_type,
                                     levels = PARTNER_ORDER)

p06 <- ggplot(contrib_plot, aes(x = partner_type, y = contribution,
                                 fill = partner_type)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ enhancer_class) +
  scale_fill_manual(values = PARTNER_COLORS) +
  labs(
    x = "Partner anchor type",
    y = "Weighted logFC contribution\n(weight x median logFC)",
    title = "Loop type contribution to overall logFC by enhancer class (all loops)",
    subtitle = "Weight = fraction of loops with that partner type"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
save_plot(p06, "06_contribution_analysis", w = 12, h = 8)


# #############################################################################
# PART D: KEY COMPARISONS AND STATISTICS
# #############################################################################

cat("\n=== PART D: Key Comparisons ===\n\n")

# --- K119ub_Only vs Stable: Wilcoxon by partner type ---
cat("  K119ub_Only vs Stable logFC (by partner type):\n")
ks_comparisons <- list()
for (pt in PARTNER_ORDER) {
  k_vals <- enh_loop_pairs$logFC[enh_loop_pairs$enhancer_class == "K119ub_Only" &
                                   enh_loop_pairs$partner_type == pt]
  s_vals <- enh_loop_pairs$logFC[enh_loop_pairs$enhancer_class == "Stable" &
                                   enh_loop_pairs$partner_type == pt]
  if (length(k_vals) >= 5 && length(s_vals) >= 5) {
    wt <- wilcox.test(k_vals, s_vals)
    cat(sprintf("    %-22s K119ub(n=%d, med=%.4f) vs Stable(n=%d, med=%.4f): %s\n",
                pt, length(k_vals), median(k_vals),
                length(s_vals), median(s_vals), fmt_p(wt$p.value)))
    ks_comparisons[[pt]] <- list(
      partner_type = pt,
      k_n = length(k_vals), k_median = median(k_vals),
      s_n = length(s_vals), s_median = median(s_vals),
      p = wt$p.value
    )
  } else {
    cat(sprintf("    %-22s insufficient data (K=%d, S=%d)\n",
                pt, length(k_vals), length(s_vals)))
  }
}

# --- K119ub_Only vs Activity_Lost: Fisher's exact for partner type enrichment ---
cat("\n  K119ub_Only vs Activity_Lost: partner type enrichment (Fisher's exact):\n")
k_all <- enh_loop_pairs[enh_loop_pairs$enhancer_class == "K119ub_Only", ]
al_all <- enh_loop_pairs[enh_loop_pairs$enhancer_class == "Activity_Lost", ]

fisher_results <- list()
for (pt in PARTNER_ORDER) {
  k_yes <- sum(k_all$partner_type == pt)
  k_no  <- nrow(k_all) - k_yes
  al_yes <- sum(al_all$partner_type == pt)
  al_no  <- nrow(al_all) - al_yes

  mat <- matrix(c(k_yes, k_no, al_yes, al_no), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  cat(sprintf("    %-22s K119ub=%d/%d (%.1f%%) vs ActLost=%d/%d (%.1f%%) OR=%.2f %s\n",
              pt, k_yes, nrow(k_all), 100 * k_yes / nrow(k_all),
              al_yes, nrow(al_all), 100 * al_yes / nrow(al_all),
              ft$estimate, fmt_p(ft$p.value)))
  fisher_results[[pt]] <- list(
    partner_type = pt,
    k_count = k_yes, k_frac = k_yes / nrow(k_all),
    al_count = al_yes, al_frac = al_yes / nrow(al_all),
    OR = ft$estimate, p = ft$p.value
  )
}

# --- Functional group comparison ---
cat("\n  Functional group proportions by enhancer class:\n")
func_table <- table(enh_loop_pairs$enhancer_class, enh_loop_pairs$partner_group)
for (cls in CLASS_ORDER) {
  total <- sum(func_table[cls, ])
  props <- 100 * func_table[cls, ] / total
  cat(sprintf("    %s: Promoter=%.1f%%, Polycomb=%.1f%%, Enhancer=%.1f%%, Structural=%.1f%%\n",
              cls, props["Promoter"], props["Polycomb"],
              props["Enhancer"], props["Structural"]))
}

chi_func <- chisq.test(func_table, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n  Chi-squared (functional group ~ enhancer class): X2=%.1f, %s\n",
            chi_func$statistic, fmt_p(chi_func$p.value)))


# #############################################################################
# PART F: ENHANCER-TSS CONTACT ANALYSIS
# #############################################################################

cat("\n=== PART F: Enhancer-TSS Contact Analysis ===\n\n")

# Filter to enhancer->promoter contacts (partner is a promoter type)
enh_tss_pairs <- enh_loop_pairs[enh_loop_pairs$partner_type %in% PROMOTER_TYPES, ]
cat(sprintf("  Enhancer-promoter loop contacts: %d (of %d total pairs)\n",
            nrow(enh_tss_pairs), nrow(enh_loop_pairs)))

# Identify gene at promoter (partner) anchor
# When enhancer is on anchor1, partner/promoter is anchor2 -> use anchor2 gene
# When enhancer is on anchor2, partner/promoter is anchor1 -> use anchor1 gene
enh_tss_pairs$promoter_gene <- ifelse(
  enh_tss_pairs$enh_on_anchor == 1,
  all_loops$anchor2_nearest_gene[enh_tss_pairs$loop_idx],
  all_loops$anchor1_nearest_gene[enh_tss_pairs$loop_idx]
)

enh_tss_pairs$promoter_distance <- ifelse(
  enh_tss_pairs$enh_on_anchor == 1,
  all_loops$anchor2_distance_to_tss[enh_tss_pairs$loop_idx],
  all_loops$anchor1_distance_to_tss[enh_tss_pairs$loop_idx]
)

# Count per enhancer class
cat("\n  Enhancer-TSS contacts per enhancer class:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_tss_pairs[enh_tss_pairs$enhancer_class == cls, ]
  n_total_enh <- sum(enh_all$enhancer_class == cls)
  n_enh_with_contact <- length(unique(sub$enh_idx))
  n_genes <- length(unique(na.omit(sub$promoter_gene)))
  n_diff <- sum(sub$is_differential)
  cat(sprintf("    %s: %d E-P contacts (%d differential), %d unique enhancers (%.1f%%), %d unique genes\n",
              cls, nrow(sub), n_diff, n_enh_with_contact,
              100 * n_enh_with_contact / n_total_enh, n_genes))
}

# Summary by enhancer class x promoter type
etss_summary <- as.data.frame(table(
  enhancer_class = enh_tss_pairs$enhancer_class,
  promoter_type = enh_tss_pairs$partner_type
))

# Add per-enhancer rates
for (cls in CLASS_ORDER) {
  n_enh <- sum(enh_all$enhancer_class == cls)
  idx <- etss_summary$enhancer_class == cls
  etss_summary$n_enhancers[idx] <- n_enh
  etss_summary$contacts_per_enhancer[idx] <- etss_summary$Freq[idx] / n_enh
}

write.table(etss_summary,
            file.path(OUTPUT_DIR, "enhancer_tss_contacts_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n",
            file.path(OUTPUT_DIR, "enhancer_tss_contacts_summary.tsv")))

# --- Plot 07: Bar chart of E-P contact counts per enhancer class ---
cat("\nGenerating Part F plots...\n")

etss_class_total <- aggregate(Freq ~ enhancer_class, data = etss_summary, FUN = sum)
colnames(etss_class_total) <- c("enhancer_class", "n_contacts")
etss_class_total$n_enhancers <- sapply(etss_class_total$enhancer_class, function(cls) {
  sum(enh_all$enhancer_class == cls)
})
etss_class_total$contacts_per_enhancer <- etss_class_total$n_contacts /
                                           etss_class_total$n_enhancers

p07 <- ggplot(etss_summary[etss_summary$Freq > 0, ],
              aes(x = enhancer_class, y = Freq, fill = promoter_type)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  geom_text(data = etss_class_total,
            aes(x = enhancer_class, y = n_contacts,
                label = paste0("n=", format(n_contacts, big.mark = ","))),
            inherit.aes = FALSE, vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = PARTNER_COLORS[PROMOTER_TYPES]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Enhancer class", y = "Number of E-P loop contacts",
    fill = "Promoter type",
    title = "Enhancer-TSS loop contacts by enhancer class",
    subtitle = "All loops (not just differential) | Partner anchor = promoter type"
  ) +
  theme_pub
save_plot(p07, "07_enhancer_tss_contact_counts", w = 8, h = 6)

# --- Plot 08: Per-enhancer E-P contact rate ---
p08 <- ggplot(etss_class_total,
              aes(x = enhancer_class, y = contacts_per_enhancer,
                  fill = enhancer_class)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", contacts_per_enhancer)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = CLASS_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Enhancer class",
    y = "E-P contacts per enhancer",
    title = "Per-enhancer rate of promoter loop contacts",
    subtitle = "All loops | contacts = loops where partner is Active/Repressed/Bivalent Promoter"
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p08, "08_enhancer_tss_contact_rate", w = 7, h = 6)


# #############################################################################
# PART G: GENE-LEVEL TSV + ABC CROSS-REFERENCE
# #############################################################################

cat("\n=== PART G: Gene-Level Output & ABC Cross-Reference ===\n\n")

# Load ABC pairs
cat("  Loading ABC pairs...\n")
abc_pairs <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d\n", nrow(abc_pairs)))
cat(sprintf("  Unique ABC genes: %d\n", length(unique(abc_pairs$TargetGene))))

# Build gene-level table from enh_tss_pairs
gene_table <- enh_tss_pairs[!is.na(enh_tss_pairs$promoter_gene), ]

# Attach enhancer coordinates for ABC matching
gene_table$enh_chr   <- enh_all$chr[gene_table$enh_idx]
gene_table$enh_start <- enh_all$start[gene_table$enh_idx]
gene_table$enh_end   <- enh_all$end[gene_table$enh_idx]

cat(sprintf("  E-P contacts with gene assignment: %d\n", nrow(gene_table)))
cat(sprintf("  Unique genes in E-P contacts: %d\n",
            length(unique(gene_table$promoter_gene))))

# Cross-reference with ABC pairs using coordinate overlap + gene name match
cat("  Cross-referencing with ABC pairs...\n")
gene_table$in_abc <- FALSE
gene_table$abc_score_wt <- NA_real_
gene_table$abc_score_ko <- NA_real_
gene_table$delta_abc <- NA_real_

# Build lookup of ABC pairs by gene name for efficient matching
abc_by_gene <- split(seq_len(nrow(abc_pairs)), abc_pairs$TargetGene)

# Build GRanges for gene_table enhancers
gt_enh_gr <- GRanges(
  seqnames = gene_table$enh_chr,
  ranges = IRanges(start = gene_table$enh_start, end = gene_table$enh_end)
)

# Build GRanges for all ABC enhancers
abc_enh_gr <- GRanges(
  seqnames = abc_pairs$chr,
  ranges = IRanges(start = abc_pairs$start, end = abc_pairs$end)
)

# Match: for each unique gene, find overlapping enhancers in both tables
unique_genes <- unique(gene_table$promoter_gene)
n_matched <- 0L

for (gene in unique_genes) {
  gt_rows <- which(gene_table$promoter_gene == gene)
  abc_rows <- abc_by_gene[[gene]]
  if (is.null(abc_rows) || length(abc_rows) == 0) next

  # Find overlaps between this gene's enhancers in gene_table and ABC
  hits <- findOverlaps(gt_enh_gr[gt_rows], abc_enh_gr[abc_rows], type = "any")
  if (length(hits) == 0) next

  # Assign ABC scores (take first match per gene_table row)
  for (h_idx in seq_len(length(hits))) {
    q_local <- queryHits(hits)[h_idx]
    s_local <- subjectHits(hits)[h_idx]
    q_global <- gt_rows[q_local]
    s_global <- abc_rows[s_local]

    if (!gene_table$in_abc[q_global]) {
      gene_table$in_abc[q_global] <- TRUE
      gene_table$abc_score_wt[q_global] <- abc_pairs$ABC.Score_WT[s_global]
      gene_table$abc_score_ko[q_global] <- abc_pairs$ABC.Score_KO[s_global]
      gene_table$delta_abc[q_global] <- abc_pairs$delta_ABC[s_global]
      n_matched <- n_matched + 1L
    }
  }
}

cat(sprintf("  ABC-matched E-P contacts: %d / %d (%.1f%%)\n",
            n_matched, nrow(gene_table), 100 * n_matched / nrow(gene_table)))

# Per-class ABC match summary
cat("\n  ABC match rate by enhancer class:\n")
for (cls in CLASS_ORDER) {
  sub <- gene_table[gene_table$enhancer_class == cls, ]
  n_in <- sum(sub$in_abc)
  cat(sprintf("    %s: %d / %d matched (%.1f%%)\n",
              cls, n_in, nrow(sub),
              ifelse(nrow(sub) > 0, 100 * n_in / nrow(sub), 0)))
  if (n_in > 0) {
    cat(sprintf("      mean delta_ABC = %.5f, median delta_ABC = %.5f\n",
                mean(sub$delta_abc, na.rm = TRUE),
                median(sub$delta_abc, na.rm = TRUE)))
  }
}

# Save detailed E-P contacts with ABC annotation
write.table(gene_table,
            file.path(OUTPUT_DIR, "enhancer_tss_contacts_with_abc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n",
            file.path(OUTPUT_DIR, "enhancer_tss_contacts_with_abc.tsv")))

# Aggregate to gene level
gene_level <- do.call(rbind, lapply(
  split(gene_table, gene_table$promoter_gene), function(g) {
    data.frame(
      gene = g$promoter_gene[1],
      n_enhancer_contacts = nrow(g),
      n_classes = length(unique(as.character(g$enhancer_class))),
      classes = paste(sort(unique(as.character(g$enhancer_class))), collapse = ";"),
      n_differential = sum(g$is_differential),
      median_logFC = median(g$logFC),
      mean_logFC = mean(g$logFC),
      n_in_abc = sum(g$in_abc),
      mean_delta_abc = ifelse(any(g$in_abc),
                               mean(g$delta_abc, na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  }
))
gene_level <- gene_level[order(-gene_level$n_enhancer_contacts), ]
rownames(gene_level) <- NULL

cat(sprintf("  Gene-level summary: %d genes\n", nrow(gene_level)))
cat(sprintf("    Genes with ABC match: %d (%.1f%%)\n",
            sum(gene_level$n_in_abc > 0),
            100 * sum(gene_level$n_in_abc > 0) / nrow(gene_level)))

write.table(gene_level,
            file.path(OUTPUT_DIR, "enhancer_tss_gene_level.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n",
            file.path(OUTPUT_DIR, "enhancer_tss_gene_level.tsv")))

# --- Plot 09: Loop logFC vs delta_ABC for matched E-P contacts ---
cat("\nGenerating Part G plots...\n")

matched_data <- gene_table[gene_table$in_abc, ]
if (nrow(matched_data) >= 10) {
  cor_test <- cor.test(matched_data$logFC, matched_data$delta_abc,
                       method = "spearman")
  p09 <- ggplot(matched_data, aes(x = logFC, y = delta_abc,
                                    color = enhancer_class)) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = CLASS_COLORS) +
    labs(
      x = "Loop logFC (KO/WT)", y = "delta ABC score (KO - WT)",
      color = "Enhancer class",
      title = "Loop strength vs ABC score change at E-P contacts",
      subtitle = sprintf("n=%d matched | Spearman rho=%.3f, %s",
                          nrow(matched_data), cor_test$estimate,
                          fmt_p(cor_test$p.value))
    ) +
    theme_pub
  save_plot(p09, "09_logfc_vs_delta_abc", w = 8, h = 7)
} else {
  cat("  Skipping plot 09: <10 matched E-P contacts for correlation.\n")
}

# --- Plot 10: ABC match rate per enhancer class ---
abc_rate <- data.frame(
  enhancer_class = CLASS_ORDER,
  stringsAsFactors = FALSE
)
abc_rate$n_contacts <- sapply(abc_rate$enhancer_class, function(cls) {
  sum(gene_table$enhancer_class == cls)
})
abc_rate$n_matched <- sapply(abc_rate$enhancer_class, function(cls) {
  sum(gene_table$enhancer_class == cls & gene_table$in_abc)
})
abc_rate$match_rate <- ifelse(abc_rate$n_contacts > 0,
                               abc_rate$n_matched / abc_rate$n_contacts, 0)
abc_rate$enhancer_class <- factor(abc_rate$enhancer_class, levels = CLASS_ORDER)

p10 <- ggplot(abc_rate, aes(x = enhancer_class, y = match_rate,
                              fill = enhancer_class)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d/%d)",
                                100 * match_rate, n_matched, n_contacts)),
            vjust = -0.2, size = 3) +
  scale_fill_manual(values = CLASS_COLORS) +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Enhancer class",
    y = "ABC match rate",
    title = "Fraction of E-P loop contacts with ABC model support",
    subtitle = "Matched by enhancer coordinate overlap + gene name"
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p10, "10_abc_match_rate", w = 7, h = 6)


# #############################################################################
# PART E: SUMMARY
# #############################################################################

cat("\n=== PART E: Summary ===\n\n")

summary_file <- file.path(OUTPUT_DIR, "summary.txt")
sink(summary_file)

cat("================================================================================\n")
cat("STEP 14: LOOP TYPE DECOMPOSITION AT ENHANCER SUBSETS (ALL LOOPS)\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUTS:\n")
cat(sprintf("  All loops: %d (from %s)\n", nrow(all_loops), ALL_LOOPS_FILE))
cat(sprintf("    Differential (up/down): %d\n", sum(all_loops$is_differential)))
cat(sprintf("    Unchanged:             %d\n", sum(!all_loops$is_differential)))
cat(sprintf("  Enhancers: %d total\n", nrow(enh_all)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_all$enhancer_class == cls)))
}

cat(sprintf("\nANNOTATION:\n"))
cat(sprintf("  Peak files (Late timepoint):\n"))
for (nm in names(PEAK_FILES)) {
  cat(sprintf("    %s: %s\n", nm, PEAK_FILES[[nm]]))
}
cat(sprintf("  TSS reference: %s\n", TSS_FILE))
cat(sprintf("  TSS threshold: %d bp\n", TSS_THRESHOLD))

cat(sprintf("\nTOTAL ENHANCER-LOOP PAIRS: %d\n", nrow(enh_loop_pairs)))

cat("\nPARTNER TYPE COMPOSITION:\n")
for (cls in CLASS_ORDER) {
  sub <- comp_table[comp_table$enhancer_class == cls, ]
  sub <- sub[order(-sub$proportion), ]
  n_total <- sum(sub$Freq)
  cat(sprintf("\n  %s (n=%d loops):\n", cls, n_total))
  for (i in seq_len(nrow(sub))) {
    if (sub$Freq[i] > 0) {
      cat(sprintf("    %-22s %4d (%5.1f%%)\n",
                  as.character(sub$partner_type[i]),
                  sub$Freq[i], 100 * sub$proportion[i]))
    }
  }
}

cat("\nFUNCTIONAL GROUP PROPORTIONS:\n")
for (cls in CLASS_ORDER) {
  total <- sum(func_table[cls, ])
  props <- 100 * func_table[cls, ] / total
  cat(sprintf("  %s: Promoter=%.1f%%, Polycomb=%.1f%%, Enhancer=%.1f%%, Structural=%.1f%%\n",
              cls, props["Promoter"], props["Polycomb"],
              props["Enhancer"], props["Structural"]))
}

cat(sprintf("\nChi-squared (partner type ~ class): %s\n", fmt_p(chi_test$p.value)))
cat(sprintf("Chi-squared (functional group ~ class): %s\n", fmt_p(chi_func$p.value)))

cat("\nLOGFC BY PARTNER TYPE (K119ub_Only, n >= 5):\n")
k_summary <- logfc_summary[logfc_summary$enhancer_class == "K119ub_Only" &
                              logfc_summary$n >= 5, ]
k_summary <- k_summary[order(-abs(k_summary$median_logFC)), ]
for (i in seq_len(nrow(k_summary))) {
  row <- k_summary[i, ]
  cat(sprintf("  %-22s n=%d, median logFC=%.4f, %s\n",
              row$partner_type, row$n, row$median_logFC, fmt_p(row$wilcox_p)))
}

cat("\nCONTRIBUTION ANALYSIS (K119ub_Only):\n")
k_contrib <- contrib[contrib$enhancer_class == "K119ub_Only", ]
k_contrib <- k_contrib[order(-abs(k_contrib$contribution)), ]
total_k <- sum(k_contrib$contribution, na.rm = TRUE)
cat(sprintf("  Total weighted logFC: %.4f\n", total_k))
for (i in seq_len(nrow(k_contrib))) {
  row <- k_contrib[i, ]
  if (!is.na(row$contribution)) {
    cat(sprintf("  %-22s weight=%.3f x med=%.4f -> contrib=%.5f (%.0f%%)\n",
                row$partner_type, row$weight, row$median_logFC,
                row$contribution,
                ifelse(total_k != 0, 100 * row$contribution / total_k, 0)))
  }
}

cat("\nK119ub_Only vs Activity_Lost PARTNER ENRICHMENT:\n")
for (pt in names(fisher_results)) {
  r <- fisher_results[[pt]]
  cat(sprintf("  %-22s K119ub=%.1f%% vs ActLost=%.1f%%, OR=%.2f, %s\n",
              r$partner_type, 100 * r$k_frac, 100 * r$al_frac,
              r$OR, fmt_p(r$p)))
}

cat("\nENHANCER-TSS CONTACTS:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_tss_pairs[enh_tss_pairs$enhancer_class == cls, ]
  n_total_enh <- sum(enh_all$enhancer_class == cls)
  n_enh_with <- length(unique(sub$enh_idx))
  n_genes <- length(unique(na.omit(sub$promoter_gene)))
  cat(sprintf("  %s: %d contacts, %d enhancers (%.1f%%), %d genes\n",
              cls, nrow(sub), n_enh_with,
              100 * n_enh_with / n_total_enh, n_genes))
}

cat(sprintf("\nGENE-LEVEL SUMMARY: %d genes with E-P contacts\n", nrow(gene_level)))
cat(sprintf("  Genes with ABC match: %d (%.1f%%)\n",
            sum(gene_level$n_in_abc > 0),
            100 * sum(gene_level$n_in_abc > 0) / nrow(gene_level)))

cat("\nABC MATCH RATE BY CLASS:\n")
for (cls in CLASS_ORDER) {
  sub <- gene_table[gene_table$enhancer_class == cls, ]
  n_in <- sum(sub$in_abc)
  cat(sprintf("  %s: %d / %d (%.1f%%)\n",
              cls, n_in, nrow(sub),
              ifelse(nrow(sub) > 0, 100 * n_in / nrow(sub), 0)))
}

sink()
cat(sprintf("  Summary written to: %s\n", summary_file))

cat("\n================================================================================\n")
cat("STEP 14 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("Plots: 10 panels (PDF + SVG + JPG)\n")
cat("Tables: 7 TSVs + 1 summary\n")
cat(sprintf("Cached annotated loops: %s\n", ANNOTATED_CACHE))
cat("================================================================================\n")
