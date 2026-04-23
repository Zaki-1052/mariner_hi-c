#!/usr/bin/env Rscript
# scripts/stripe_loop_crossref_extended.R
#
# T1.4 - Stripe x loop cross-reference across all 7 anchor categories.
# Extension of scripts/ctcf_stripe_crossref.R, which only tested CTCF loops.
# Uses the anchor1_type / anchor2_type columns from characterized_loops.tsv
# (7 categories matching annotate_loops_extended.R taxonomy).
#
# Usage:
#   Rscript scripts/stripe_loop_crossref_extended.R --timepoint late --resolution 5000
#
# Inputs:
#   - Loops (late):  outputs/250402-late_outputs/merged_loops/characterized_loops.tsv
#     Loops (early): outputs/250831-early_outputs/merged_loops/characterized_loops.tsv
#     Columns: loop_id, chr1/2, start1/2, end1/2, logFC, FDR, direction,
#              anchor1_type, anchor2_type, loop_type, ...
#   - Stripes: stripes/stripenn/outputs/{tp}/visualizations/{tp}_annotated_stripes.tsv
#
# Outputs (outputs/stripe_integration/{label}/loop_crossref/):
#   - loop_stripe_overlap_by_type.tsv     Per-anchor-type contingency table
#   - permutation_pvalues.tsv             Circular-shuffle null enrichment p-values
#   - enrichment_heatmap                  Odds ratio heatmap (anchor type x stripe direction)
#   - loop_type_by_stripe_direction_bar   Bar plot
#   - summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

source("scripts/utils/multi_format_output.R")

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default)
  if (i == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i + 1]
}
TP_LABEL   <- parse_arg("--timepoint", "late")
RESOLUTION <- as.integer(parse_arg("--resolution", "5000"))
ANCHOR_TOL <- as.integer(parse_arg("--tolerance", "50000"))
N_PERM     <- as.integer(parse_arg("--n-permutations", "1000"))
TP_MAP <- list(late = "250402", early = "250831")
if (!TP_LABEL %in% names(TP_MAP)) stop("--timepoint must be 'late' or 'early'")
TP_ID <- TP_MAP[[TP_LABEL]]

STRIPE_FILE <- sprintf("stripes/stripenn/outputs/%s/visualizations/%s_annotated_stripes.tsv",
                      TP_ID, TP_ID)
LOOPS_FILE <- if (TP_LABEL == "late") {
  "outputs/250402-late_outputs/merged_loops/characterized_loops.tsv"
} else {
  "outputs/250831-early_outputs/merged_loops/characterized_loops.tsv"
}
OUTPUT_DIR <- sprintf("outputs/stripe_integration/%s/loop_crossref", TP_LABEL)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("T1.4 - Stripe x Loop Cross-Reference (all 7 anchor types)\n")
cat("==============================================================\n")
cat(sprintf("Timepoint    : %s (%s)\n", TP_LABEL, TP_ID))
cat(sprintf("Tolerance    : %d bp\n", ANCHOR_TOL))
cat(sprintf("Permutations : %d\n", N_PERM))
cat(sprintf("Stripes      : %s\n", STRIPE_FILE))
cat(sprintf("Loops        : %s\n", LOOPS_FILE))
cat(sprintf("Output       : %s\n\n", OUTPUT_DIR))

stopifnot(file.exists(STRIPE_FILE), file.exists(LOOPS_FILE))

# -------- LOAD --------
stripes <- read.table(STRIPE_FILE, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      na.strings = c("NA", ""))
loops <- read.table(LOOPS_FILE, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", comment.char = "",
                    na.strings = c("NA", ""))
cat(sprintf("Stripes: %d; Loops: %d\n", nrow(stripes), nrow(loops)))

# -------- BUILD GRanges --------
# Stripe anchor (narrow; the "foot" of the stripe)
stripe_anchor_gr <- GRanges(seqnames = stripes$chr,
                            ranges = IRanges(start = as.integer(stripes$pos1),
                                             end   = as.integer(stripes$pos2)))
mcols(stripe_anchor_gr) <- stripes[, c("stripe_id", "direction",
                                       "direction_confidence",
                                       "significant_FDR05",
                                       "anchor_type")]
names(stripe_anchor_gr) <- stripes$stripe_id

# Loop anchor 1 and anchor 2 as separate GRanges, both tagged with loop_id and loop direction
make_loop_anchor_gr <- function(loops, which_anchor = 1) {
  if (which_anchor == 1) {
    chr <- loops$anchor1_chr; s <- loops$anchor1_start; e <- loops$anchor1_end
    typ <- loops$anchor1_type
  } else {
    chr <- loops$anchor2_chr; s <- loops$anchor2_start; e <- loops$anchor2_end
    typ <- loops$anchor2_type
  }
  # fall back to chr1/start1 if anchorN_* not populated
  if (any(is.na(chr))) {
    if (which_anchor == 1) {
      chr <- ifelse(is.na(chr), loops$chr1, chr)
      s   <- ifelse(is.na(s),   loops$start1, s)
      e   <- ifelse(is.na(e),   loops$end1,   e)
    } else {
      chr <- ifelse(is.na(chr), loops$chr2, chr)
      s   <- ifelse(is.na(s),   loops$start2, s)
      e   <- ifelse(is.na(e),   loops$end2,   e)
    }
  }
  keep <- !is.na(chr) & !is.na(s) & !is.na(e)
  gr <- GRanges(seqnames = chr[keep],
                ranges = IRanges(start = as.integer(s[keep]),
                                 end   = as.integer(e[keep])))
  mcols(gr) <- data.frame(
    loop_id     = loops$loop_id[keep],
    loop_fdr    = as.numeric(loops$FDR)[keep],
    loop_logfc  = as.numeric(loops$logFC)[keep],
    loop_direction = loops$direction[keep],
    anchor_num  = which_anchor,
    anchor_type = typ[keep],
    stringsAsFactors = FALSE
  )
  gr
}
loop_a1 <- make_loop_anchor_gr(loops, 1)
loop_a2 <- make_loop_anchor_gr(loops, 2)
all_loop_anchors <- c(loop_a1, loop_a2)
cat(sprintf("Loop anchors: %d (a1=%d + a2=%d)\n",
            length(all_loop_anchors), length(loop_a1), length(loop_a2)))

# Categorical axes
STRIPE_DIRS <- c("gained", "lost", "unchanged")
# Anchor types come from the loops (7 categories): Active_Promoter, Active_Enhancer,
# Poised_Enhancer, Repressed_Promoter, Bivalent_Promoter, Polycomb, Other
ANCHOR_TYPES <- sort(unique(na.omit(mcols(all_loop_anchors)$anchor_type)))
cat("Loop anchor types present:\n")
for (t in ANCHOR_TYPES) {
  cat(sprintf("  %-22s %6d\n", t, sum(mcols(all_loop_anchors)$anchor_type == t, na.rm = TRUE)))
}

# Restrict stripes to analyzable direction subset (drop strengthened/weakened
# because populations are tiny per the late summary: 3 / 4 at FDR<0.05)
stripe_keep <- stripes$direction %in% STRIPE_DIRS
stripe_anchor_gr <- stripe_anchor_gr[stripe_keep]
cat(sprintf("Stripes used (direction in {%s}): %d\n",
            paste(STRIPE_DIRS, collapse = ","), length(stripe_anchor_gr)))

# -------- OBSERVED OVERLAPS --------
# Per (anchor_type, stripe_direction): count unique loops whose anchor falls within
# ANCHOR_TOL of a stripe anchor of the given direction.
observed_counts <- function(loop_gr, stripe_gr, stripe_direction_vec) {
  # loop_gr: GRanges with mcols$anchor_type
  out <- data.frame(
    anchor_type = character(),
    stripe_direction = character(),
    n_overlapping_loops = integer(),
    n_total_loops_of_type = integer(),
    stringsAsFactors = FALSE
  )
  for (at in ANCHOR_TYPES) {
    la <- loop_gr[!is.na(mcols(loop_gr)$anchor_type) & mcols(loop_gr)$anchor_type == at]
    if (length(la) == 0) next
    n_total_loops <- length(unique(mcols(la)$loop_id))
    for (sd in STRIPE_DIRS) {
      sg <- stripe_gr[stripe_direction_vec == sd]
      if (length(sg) == 0) {
        out <- rbind(out, data.frame(anchor_type = at, stripe_direction = sd,
                                     n_overlapping_loops = 0L,
                                     n_total_loops_of_type = n_total_loops,
                                     stringsAsFactors = FALSE))
        next
      }
      ov <- findOverlaps(la, sg, maxgap = ANCHOR_TOL)
      n_ov <- length(unique(mcols(la)$loop_id[queryHits(ov)]))
      out <- rbind(out, data.frame(anchor_type = at, stripe_direction = sd,
                                   n_overlapping_loops = n_ov,
                                   n_total_loops_of_type = n_total_loops,
                                   stringsAsFactors = FALSE))
    }
  }
  out
}
stripe_dir_vec <- stripes$direction[stripe_keep]
obs <- observed_counts(all_loop_anchors, stripe_anchor_gr, stripe_dir_vec)
obs$observed_fraction <- obs$n_overlapping_loops / obs$n_total_loops_of_type

# -------- PERMUTATION TEST (circular shuffle within chromosome) --------
# Shuffle stripe positions within their own chromosome; preserve widths and
# per-chromosome counts.
set.seed(42)

# Chromosome "size" estimates - max coordinate seen across stripes+loops per chr
all_coords <- data.frame(
  chr = c(as.character(seqnames(stripe_anchor_gr)),
          as.character(seqnames(all_loop_anchors))),
  endv = c(end(stripe_anchor_gr), end(all_loop_anchors)),
  stringsAsFactors = FALSE
)
chr_max <- tapply(all_coords$endv, all_coords$chr, max)

do_permutation <- function(n_perm) {
  s_chr   <- as.character(seqnames(stripe_anchor_gr))
  s_start <- start(stripe_anchor_gr)
  s_width <- width(stripe_anchor_gr)
  res_arr <- array(NA_real_, dim = c(length(ANCHOR_TYPES), length(STRIPE_DIRS), n_perm),
                   dimnames = list(ANCHOR_TYPES, STRIPE_DIRS, NULL))

  # Pre-subset loop anchors by type once (outside permutation loop).
  loop_by_type <- lapply(ANCHOR_TYPES, function(at) {
    all_loop_anchors[!is.na(mcols(all_loop_anchors)$anchor_type) &
                     mcols(all_loop_anchors)$anchor_type == at]
  })
  names(loop_by_type) <- ANCHOR_TYPES
  loop_ids_by_type <- lapply(loop_by_type, function(g) mcols(g)$loop_id)

  # Pre-compute per-chromosome index sets for faster shuffling.
  chr_uniq <- unique(s_chr)
  chr_idx <- lapply(chr_uniq, function(ch) which(s_chr == ch))
  names(chr_idx) <- chr_uniq
  chr_lens <- vapply(chr_uniq, function(ch) {
    v <- as.integer(chr_max[ch])
    if (is.na(v) || v <= 0) max(s_start[chr_idx[[ch]]]) + 1L else v
  }, integer(1))
  names(chr_lens) <- chr_uniq

  # Pre-compute stripe-direction masks once.
  dir_masks <- lapply(STRIPE_DIRS, function(sd) stripe_dir_vec == sd)
  names(dir_masks) <- STRIPE_DIRS

  for (i in seq_len(n_perm)) {
    shifted_start <- integer(length(s_start))
    for (ch in chr_uniq) {
      idx <- chr_idx[[ch]]
      chr_len <- chr_lens[[ch]]
      shift <- sample.int(chr_len, size = 1L) - 1L
      shifted_start[idx] <- ((s_start[idx] - 1L + shift) %% chr_len) + 1L
    }
    shifted_end <- shifted_start + s_width - 1L
    shuf_gr <- GRanges(seqnames = s_chr,
                       ranges = IRanges(start = shifted_start, end = shifted_end))

    for (ai in seq_along(ANCHOR_TYPES)) {
      la <- loop_by_type[[ai]]
      lids <- loop_ids_by_type[[ai]]
      if (length(la) == 0) next
      for (di in seq_along(STRIPE_DIRS)) {
        sg <- shuf_gr[dir_masks[[di]]]
        if (length(sg) == 0) { res_arr[ai, di, i] <- 0; next }
        ov <- findOverlaps(la, sg, maxgap = ANCHOR_TOL)
        res_arr[ai, di, i] <- length(unique(lids[queryHits(ov)]))
      }
    }
    if (i %% max(1L, n_perm %/% 10L) == 0L) {
      cat(sprintf("  permutation %d / %d\n", i, n_perm))
    }
  }
  res_arr
}

cat(sprintf("\nRunning %d circular-shuffle permutations...\n", N_PERM))
null_arr <- do_permutation(N_PERM)

# Derive one-tailed (enrichment) and two-tailed empirical p-values per cell
perm_rows <- list()
for (at in ANCHOR_TYPES) {
  for (sd in STRIPE_DIRS) {
    obs_val <- obs$n_overlapping_loops[obs$anchor_type == at & obs$stripe_direction == sd]
    null_vec <- as.numeric(null_arr[at, sd, ])
    null_vec <- null_vec[!is.na(null_vec)]
    if (length(null_vec) == 0) next
    null_mean <- mean(null_vec); null_sd <- sd(null_vec)
    p_enrich  <- (sum(null_vec >= obs_val) + 1L) / (length(null_vec) + 1L)
    p_deplete <- (sum(null_vec <= obs_val) + 1L) / (length(null_vec) + 1L)
    z         <- ifelse(null_sd > 0, (obs_val - null_mean) / null_sd, NA_real_)
    fold_enrich <- ifelse(null_mean > 0, obs_val / null_mean, NA_real_)
    perm_rows[[length(perm_rows) + 1L]] <- data.frame(
      anchor_type = at, stripe_direction = sd,
      observed = obs_val, null_mean = null_mean, null_sd = null_sd,
      z_score = z, fold_enrichment = fold_enrich,
      p_enrichment = p_enrich, p_depletion = p_deplete,
      stringsAsFactors = FALSE
    )
  }
}
perm_df <- do.call(rbind, perm_rows)
perm_df$p_enrichment_BH <- p.adjust(perm_df$p_enrichment, method = "BH")
perm_df$p_depletion_BH  <- p.adjust(perm_df$p_depletion,  method = "BH")

# Merge observed + permutation into one table
out_tbl <- merge(obs, perm_df[, c("anchor_type", "stripe_direction",
                                  "null_mean", "null_sd", "z_score", "fold_enrichment",
                                  "p_enrichment", "p_enrichment_BH",
                                  "p_depletion", "p_depletion_BH")],
                 by = c("anchor_type", "stripe_direction"))
write.table(out_tbl,
            file.path(OUTPUT_DIR, "loop_stripe_overlap_by_type.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(perm_df,
            file.path(OUTPUT_DIR, "permutation_pvalues.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------- ENRICHMENT HEATMAP --------
heat_df <- out_tbl %>%
  mutate(label = sprintf("%.2f\n(p=%.2g)", fold_enrichment, p_enrichment_BH),
         star = ifelse(p_enrichment_BH < 0.05, "*", ""))
p_heat <- ggplot(heat_df,
                 aes(x = stripe_direction, y = anchor_type,
                     fill = log2(pmax(fold_enrichment, 1e-3)))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = paste0(sprintf("%.2fx", fold_enrichment), star)),
            size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "log2(fold\nenrichment)") +
  labs(title = sprintf("Stripe-anchor overlap enrichment vs circular-shuffle null (%s)", TP_LABEL),
       subtitle = sprintf("n=%d permutations; * = BH FDR<0.05", N_PERM),
       x = "Stripe direction", y = "Loop anchor type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))
save_multiformat_ggplot(p_heat,
                        file.path(OUTPUT_DIR, "enrichment_heatmap"),
                        width = 8, height = 7, use_subfolders = TRUE)

# -------- LOOP TYPE x STRIPE DIRECTION BAR --------
# Distribution of loop anchor types among loops whose anchor overlaps a stripe
# of each direction.
bar_rows <- list()
for (sd in STRIPE_DIRS) {
  sg <- stripe_anchor_gr[stripe_dir_vec == sd]
  if (length(sg) == 0) next
  ov <- findOverlaps(all_loop_anchors, sg, maxgap = ANCHOR_TOL)
  hit_types <- mcols(all_loop_anchors)$anchor_type[queryHits(ov)]
  if (length(hit_types) == 0) next
  tab <- as.data.frame(table(hit_types), stringsAsFactors = FALSE)
  colnames(tab) <- c("anchor_type", "n")
  tab$stripe_direction <- sd
  bar_rows[[length(bar_rows) + 1L]] <- tab
}
if (length(bar_rows) > 0) {
  bar_df <- do.call(rbind, bar_rows) %>%
    group_by(stripe_direction) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  p_bar <- ggplot(bar_df,
                  aes(x = stripe_direction, y = pct, fill = anchor_type)) +
    geom_col() +
    labs(title = sprintf("Loop anchor types at stripe anchors by direction (%s)", TP_LABEL),
         x = "Stripe direction", y = "Percent of overlapping loop anchors",
         fill = "Loop anchor type") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_multiformat_ggplot(p_bar,
                          file.path(OUTPUT_DIR, "loop_type_by_stripe_direction_bar"),
                          width = 9, height = 6, use_subfolders = TRUE)
}

# -------- SUMMARY --------
summary_lines <- c(
  sprintf("Stripe x loop crossref (%s, res=%d)", TP_LABEL, RESOLUTION),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Tolerance: %d bp; Permutations: %d", ANCHOR_TOL, N_PERM),
  sprintf("Stripes analysed (direction in gained/lost/unchanged): %d",
          length(stripe_anchor_gr)),
  sprintf("Loop anchors: %d", length(all_loop_anchors)),
  ""
)
sig_rows <- perm_df %>%
  filter(p_enrichment_BH < 0.05) %>%
  arrange(p_enrichment_BH)
summary_lines <- c(summary_lines,
                   sprintf("Cells with BH-adj enrichment p<0.05: %d", nrow(sig_rows)),
                   capture.output(print(sig_rows, row.names = FALSE)))
writeLines(summary_lines, file.path(OUTPUT_DIR, "summary.txt"))
cat("\nDone.\n")
