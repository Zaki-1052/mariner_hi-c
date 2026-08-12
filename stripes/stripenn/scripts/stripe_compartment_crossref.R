#!/usr/bin/env Rscript
# stripes/stripenn/scripts/stripe_compartment_crossref.R
#
# T1.2 - Stripe x Compartment (PC1) cross-reference.
#
# Tests whether differential stripes (gained/lost/strengthened/weakened)
# associate with A/B compartment state at their anchor and body regions.
# Structural analog of scripts/loop_compartment_crossref.R, adapted for
# 1D stripe geometry (anchor pos1-pos2, body pos3-pos4).
#
# Usage:
#   Rscript scripts/stripe_compartment_crossref.R --timepoint late --resolution 5000
#
# Inputs:
#   - Stripes: stripes/stripenn/outputs/{tp}/visualizations/{tp}_annotated_stripes.tsv
#   - Compartments: tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv  (late)
#                   tads/tad-pc-analysis/output/compartment_analysis_early/...                       (early)
#     Schema: Region_ID, Chr, Start, End, ctrl_avg_PC1, mut_avg_PC1, Difference, adj_pvalue,
#             direction ('A_to_B_in_Mutant' | 'B_to_A_in_Mutant' | 'No_Change'),
#             direction_label, significant (in the _significant_standard subset)
#
# Outputs (outputs/stripe_integration/{label}/compartment_crossref/):
#   - stripe_compartment_assignment.tsv   Per-stripe: anchor/body mean PC1 (ctrl+mut), A/B calls
#   - enrichment_tests.tsv                Fisher + Wilcoxon test results by direction
#   - anchor_compartment_by_direction     Stacked-bar plot (pdf/svg/jpg)
#   - pc1_by_direction                    Violin of anchor mean PC1 by direction
#   - switched_bin_overlap                Counts of stripes on A<->B significant bins
#   - summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

if (dir.exists("/expanse/lustre/projects/csd940/zalibhai")) {
  CODE_DIR  <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
  DATA_DIR  <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
  REPO_ROOT <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c"
} else {
  CODE_DIR  <- normalizePath(file.path(getwd()))
  DATA_DIR  <- CODE_DIR
  REPO_ROOT <- normalizePath(file.path(CODE_DIR, "..", ".."))
}

source(file.path(REPO_ROOT, "scripts/utils/multi_format_output.R"))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default)
  if (i == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i + 1]
}
TP_LABEL   <- parse_arg("--timepoint", "late")
RESOLUTION <- as.integer(parse_arg("--resolution", "5000"))
TP_MAP <- list(late = "250402", early = "250831")
if (!TP_LABEL %in% names(TP_MAP)) stop("--timepoint must be 'late' or 'early'")
TP_ID <- TP_MAP[[TP_LABEL]]

STRIPE_FILE <- file.path(CODE_DIR, sprintf("outputs/%s/visualizations/%s_annotated_stripes.tsv",
                      TP_ID, TP_ID))

PC_BASE <- if (TP_LABEL == "late") {
  file.path(REPO_ROOT, "tads/tad-pc-analysis/output/compartment_analysis")
} else {
  file.path(REPO_ROOT, "tads/tad-pc-analysis/output/compartment_analysis_early")
}
PC_ALL_FILE <- file.path(PC_BASE, "compartment_all_annotated.tsv")
PC_SIG_FILE <- file.path(PC_BASE, "compartment_significant_standard.tsv")

OUTPUT_DIR <- file.path(REPO_ROOT, sprintf("outputs/stripe_integration/%s/compartment_crossref", TP_LABEL))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("T1.2 - Stripe x Compartment Cross-Reference\n")
cat("==============================================================\n")
cat(sprintf("Timepoint : %s (%s)\n", TP_LABEL, TP_ID))
cat(sprintf("Resolution: %d bp\n", RESOLUTION))
cat(sprintf("Stripes   : %s\n", STRIPE_FILE))
cat(sprintf("PC all    : %s\n", PC_ALL_FILE))
cat(sprintf("PC sig    : %s\n", PC_SIG_FILE))
cat(sprintf("Output    : %s\n\n", OUTPUT_DIR))

stopifnot(file.exists(STRIPE_FILE), file.exists(PC_ALL_FILE))

# -------- LOAD --------
stripes <- read.table(STRIPE_FILE, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      na.strings = c("NA", ""))
pc <- read.table(PC_ALL_FILE, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "", comment.char = "",
                 na.strings = c("NA", ""))
cat(sprintf("Stripes: %d rows\n", nrow(stripes)))
cat(sprintf("PC bins: %d rows\n", nrow(pc)))

# Coerce PC1 columns
pc$ctrl_avg_PC1 <- suppressWarnings(as.numeric(pc$ctrl_avg_PC1))
pc$mut_avg_PC1  <- suppressWarnings(as.numeric(pc$mut_avg_PC1))
pc$Difference   <- suppressWarnings(as.numeric(pc$Difference))
pc$adj_pvalue   <- suppressWarnings(as.numeric(pc$adj_pvalue))

# Drop bins with no signal
pc <- pc[!is.na(pc$ctrl_avg_PC1) & !is.na(pc$mut_avg_PC1), , drop = FALSE]
cat(sprintf("PC bins with signal: %d\n", nrow(pc)))

# -------- GRANGES --------
pc_gr <- GRanges(seqnames = pc$Chr,
                 ranges = IRanges(start = as.integer(pc$Start) + 1L,
                                  end   = as.integer(pc$End)))
mcols(pc_gr) <- pc[, c("Region_ID", "ctrl_avg_PC1", "mut_avg_PC1",
                       "Difference", "adj_pvalue", "direction")]

anchor_gr <- GRanges(seqnames = stripes$chr,
                     ranges = IRanges(start = as.integer(stripes$pos1),
                                      end   = as.integer(stripes$pos2)))
body_gr   <- GRanges(seqnames = stripes$chr,
                     ranges = IRanges(start = as.integer(stripes$pos3),
                                      end   = as.integer(stripes$pos4)))

# -------- AGGREGATE PC1 OVER ANCHORS AND BODIES --------
# For each stripe, mean PC1 across overlapping bins (query=stripe, subject=pc).
mean_pc_over <- function(query_gr) {
  ov <- findOverlaps(query_gr, pc_gr)
  ctrl_vals <- rep(NA_real_, length(query_gr))
  mut_vals  <- rep(NA_real_, length(query_gr))
  n_over    <- rep(0L, length(query_gr))
  if (length(ov) > 0) {
    qh <- queryHits(ov); sh <- subjectHits(ov)
    ctrl_bin <- mcols(pc_gr)$ctrl_avg_PC1[sh]
    mut_bin  <- mcols(pc_gr)$mut_avg_PC1[sh]
    agg_ctrl <- tapply(ctrl_bin, qh, mean, na.rm = TRUE)
    agg_mut  <- tapply(mut_bin,  qh, mean, na.rm = TRUE)
    agg_n    <- tapply(qh, qh, length)
    idx <- as.integer(names(agg_ctrl))
    ctrl_vals[idx] <- as.numeric(agg_ctrl)
    mut_vals[idx]  <- as.numeric(agg_mut)
    n_over[idx]    <- as.integer(agg_n)
  }
  data.frame(ctrl = ctrl_vals, mut = mut_vals, n_bins = n_over)
}

cat("\nAggregating PC1 over anchors...\n")
anchor_pc <- mean_pc_over(anchor_gr)
cat("Aggregating PC1 over stripe bodies...\n")
body_pc   <- mean_pc_over(body_gr)

# -------- CLASSIFY --------
# PC1 > 0 is canonically A; PC1 < 0 is B. Use ctrl average for baseline call.
classify_compartment <- function(pc1_ctrl) {
  cls <- rep(NA_character_, length(pc1_ctrl))
  cls[!is.na(pc1_ctrl) & pc1_ctrl >  0] <- "A"
  cls[!is.na(pc1_ctrl) & pc1_ctrl <  0] <- "B"
  cls[!is.na(pc1_ctrl) & pc1_ctrl == 0] <- "boundary"
  cls
}
anchor_class <- classify_compartment(anchor_pc$ctrl)
body_class   <- classify_compartment(body_pc$ctrl)

# -------- SWITCHED-BIN OVERLAP --------
switched_hit <- rep(FALSE, length(anchor_gr))
if (file.exists(PC_SIG_FILE)) {
  pc_sig <- read.table(PC_SIG_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "", comment.char = "",
                       na.strings = c("NA", ""))
  if (nrow(pc_sig) > 0) {
    pc_sig_gr <- GRanges(seqnames = pc_sig$Chr,
                         ranges = IRanges(start = as.integer(pc_sig$Start) + 1L,
                                          end   = as.integer(pc_sig$End)))
    switched_hit <- overlapsAny(anchor_gr, pc_sig_gr)
  }
  cat(sprintf("Significant-PC1 bins: %d; anchors on switched bins: %d\n",
              nrow(pc_sig), sum(switched_hit)))
}

# -------- ASSEMBLE OUTPUT TABLE --------
assign_tbl <- data.frame(
  stripe_id           = stripes$stripe_id,
  chr                 = stripes$chr,
  pos1                = stripes$pos1,
  pos2                = stripes$pos2,
  pos3                = stripes$pos3,
  pos4                = stripes$pos4,
  direction           = stripes$direction,
  direction_confidence = stripes$direction_confidence,
  significant_FDR05   = stripes$significant_FDR05,
  stripe_logFC        = stripes$logFC,
  stripe_FDR          = stripes$FDR,
  anchor_pc1_ctrl     = anchor_pc$ctrl,
  anchor_pc1_mut      = anchor_pc$mut,
  anchor_pc1_delta    = anchor_pc$mut - anchor_pc$ctrl,
  anchor_n_bins       = anchor_pc$n_bins,
  anchor_compartment  = anchor_class,
  body_pc1_ctrl       = body_pc$ctrl,
  body_pc1_mut        = body_pc$mut,
  body_pc1_delta      = body_pc$mut - body_pc$ctrl,
  body_n_bins         = body_pc$n_bins,
  body_compartment    = body_class,
  anchor_on_switched_bin = switched_hit,
  stringsAsFactors    = FALSE
)
write.table(assign_tbl,
            file.path(OUTPUT_DIR, "stripe_compartment_assignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# -------- STATISTICAL TESTS --------
tests <- list()

# 2x2 Fisher: anchor A vs anchor B  x  gained vs lost (high-conf only)
hc <- assign_tbl %>%
  filter(direction %in% c("gained", "lost"),
         direction_confidence == "high",
         significant_FDR05 == TRUE,
         !is.na(anchor_compartment),
         anchor_compartment %in% c("A", "B"))
count_cell <- function(df, dir, comp) {
  sum(df$direction == dir & df$anchor_compartment == comp, na.rm = TRUE)
}
if (nrow(hc) >= 10) {
  n_gA <- count_cell(hc, "gained", "A"); n_gB <- count_cell(hc, "gained", "B")
  n_lA <- count_cell(hc, "lost",   "A"); n_lB <- count_cell(hc, "lost",   "B")
  tab <- matrix(c(n_gA, n_lA, n_gB, n_lB), nrow = 2,
                dimnames = list(c("gained", "lost"), c("A", "B")))
  if (sum(tab) >= 10 && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
    ft <- fisher.test(tab)
    tests$fisher_anchor_AB_by_direction <- data.frame(
      test = "Fisher: anchor A vs B x gained vs lost (high-conf sig)",
      n_gained_A = n_gA, n_gained_B = n_gB,
      n_lost_A   = n_lA, n_lost_B   = n_lB,
      odds_ratio = as.numeric(ft$estimate),
      p_value    = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}

# Wilcoxon: anchor ctrl PC1 by direction (pairwise gained vs lost)
wlx <- function(df, col) {
  g <- df[df$direction == "gained", col]
  l <- df[df$direction == "lost",   col]
  g <- g[!is.na(g)]; l <- l[!is.na(l)]
  if (length(g) < 5 || length(l) < 5) return(NULL)
  w <- suppressWarnings(wilcox.test(g, l))
  data.frame(
    test       = sprintf("Wilcoxon: %s (gained vs lost, high-conf sig)", col),
    n_gained   = length(g), n_lost = length(l),
    median_gained = median(g), median_lost = median(l),
    W          = as.numeric(w$statistic),
    p_value    = w$p.value,
    stringsAsFactors = FALSE
  )
}
wres <- list()
for (col in c("anchor_pc1_ctrl", "anchor_pc1_delta",
              "body_pc1_ctrl", "body_pc1_delta")) {
  r <- wlx(hc, col); if (!is.null(r)) wres[[col]] <- r
}
if (length(wres) > 0) tests$wilcoxon <- do.call(rbind, wres)

# Enrichment of switched-bin anchors among gained/lost high-conf
if (nrow(hc) >= 10) {
  count_sw <- function(df, dir, on_sw) {
    sum(df$direction == dir & df$anchor_on_switched_bin == on_sw, na.rm = TRUE)
  }
  n_g_sw <- count_sw(hc, "gained", TRUE);  n_g_st <- count_sw(hc, "gained", FALSE)
  n_l_sw <- count_sw(hc, "lost",   TRUE);  n_l_st <- count_sw(hc, "lost",   FALSE)
  tab2 <- matrix(c(n_g_sw, n_l_sw, n_g_st, n_l_st), nrow = 2,
                 dimnames = list(c("gained", "lost"), c("switched", "stable")))
  if (sum(tab2) >= 10 && all(rowSums(tab2) > 0) && all(colSums(tab2) > 0)) {
    ft2 <- fisher.test(tab2)
    tests$fisher_switched_bin <- data.frame(
      test = "Fisher: switched-PC1 bin x gained vs lost (high-conf sig)",
      n_gained_switched = n_g_sw, n_gained_stable = n_g_st,
      n_lost_switched   = n_l_sw, n_lost_stable   = n_l_st,
      odds_ratio = as.numeric(ft2$estimate),
      p_value    = ft2$p.value,
      stringsAsFactors = FALSE
    )
  }
}

tests_long <- do.call(rbind, lapply(tests, function(x) {
  detail_cols <- setdiff(colnames(x), c("test", "p_value"))
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    vals <- vapply(detail_cols, function(k) format(x[[k]][i], digits = 4), character(1))
    data.frame(test = x$test[i],
               details = paste(detail_cols, vals, sep = "=", collapse = "; "),
               p_value = x$p_value[i],
               stringsAsFactors = FALSE)
  }))
}))
if (!is.null(tests_long)) {
  tests_long$p_adj_BH <- p.adjust(tests_long$p_value, method = "BH")
  write.table(tests_long,
              file.path(OUTPUT_DIR, "enrichment_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# -------- PLOTS --------
# Stacked bar: anchor compartment x direction
plot_df <- assign_tbl %>%
  filter(!is.na(direction), !is.na(anchor_compartment),
         anchor_compartment %in% c("A", "B"),
         direction %in% c("gained", "lost", "unchanged", "strengthened", "weakened")) %>%
  count(direction, anchor_compartment) %>%
  group_by(direction) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p1 <- ggplot(plot_df,
             aes(x = direction, y = pct, fill = anchor_compartment)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d", n)),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(A = "#E74C3C", B = "#3498DB")) +
  labs(title = sprintf("Anchor compartment by stripe direction (%s)", TP_LABEL),
       y = "Percent of stripes", x = "Stripe direction",
       fill = "Anchor\ncompartment") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_multiformat_ggplot(p1,
                        file.path(OUTPUT_DIR, "anchor_compartment_by_direction"),
                        width = 8, height = 6, use_subfolders = TRUE)

# Violin: anchor ctrl PC1 by direction
vdf <- assign_tbl %>%
  filter(!is.na(anchor_pc1_ctrl),
         direction %in% c("gained", "lost", "unchanged"))
if (nrow(vdf) >= 10) {
  p2 <- ggplot(vdf, aes(x = direction, y = anchor_pc1_ctrl, fill = direction)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c(gained = "#D73027", lost = "#4575B4",
                                 unchanged = "#999999")) +
    labs(title = sprintf("Anchor mean PC1 (ctrl) by direction (%s)", TP_LABEL),
         x = "Stripe direction", y = "Anchor mean PC1 (ctrl)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  save_multiformat_ggplot(p2,
                          file.path(OUTPUT_DIR, "pc1_by_direction"),
                          width = 7, height = 6, use_subfolders = TRUE)
}

# Switched-bin bar
sw_df <- assign_tbl %>%
  filter(direction %in% c("gained", "lost", "unchanged")) %>%
  count(direction, anchor_on_switched_bin) %>%
  group_by(direction) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
if (nrow(sw_df) > 0) {
  p3 <- ggplot(sw_df,
               aes(x = direction, y = pct,
                   fill = factor(anchor_on_switched_bin,
                                 levels = c(TRUE, FALSE),
                                 labels = c("switched", "stable")))) +
    geom_col() +
    scale_fill_manual(values = c(switched = "#9B59B6", stable = "#BDC3C7")) +
    labs(title = sprintf("Anchors on significant-PC1-switched bins (%s)", TP_LABEL),
         x = "Stripe direction", y = "Percent",
         fill = "Bin type") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_multiformat_ggplot(p3,
                          file.path(OUTPUT_DIR, "switched_bin_overlap"),
                          width = 8, height = 6, use_subfolders = TRUE)
}

# -------- SUMMARY --------
summary_lines <- c(
  sprintf("Stripe x compartment crossref (%s, res=%d)", TP_LABEL, RESOLUTION),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Stripes:      %d", nrow(assign_tbl)),
  sprintf("PC bins used: %d", nrow(pc)),
  sprintf("Anchors with PC1 signal: %d", sum(!is.na(assign_tbl$anchor_pc1_ctrl))),
  sprintf("Bodies with PC1 signal:  %d", sum(!is.na(assign_tbl$body_pc1_ctrl))),
  sprintf("Anchors on switched bins: %d", sum(assign_tbl$anchor_on_switched_bin, na.rm = TRUE)),
  ""
)
if (!is.null(tests_long)) {
  summary_lines <- c(summary_lines, "Test results (BH-adjusted):", "")
  summary_lines <- c(summary_lines,
                     capture.output(print(tests_long[, c("test", "p_value", "p_adj_BH")],
                                          row.names = FALSE)))
}
writeLines(summary_lines, file.path(OUTPUT_DIR, "summary.txt"))
cat("\nDone.\n")
